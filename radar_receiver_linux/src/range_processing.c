#include "processor_internal.h"

#include <math.h>
#include <string.h>

#define MF_LANES   24           /* beams per fast-convolution chunk (3 AVX registers) */

/* Fast convolution of the matched filter: FFT length and the FFT of the taps.
 * It is used when its estimated time is shorter than that of the direct
 * convolution (which is split over the threads beam by beam; the FFTs of
 * MF_LANES beams run in one thread; a single beam is split into the 8
 * interleaved sub-FFTs of a four-step FFT instead, which fill the vectors).
 * Costs per operation measured on an AVX2 Xeon: 0.3 ns per complex MAC
 * (direct), 0.6 ns per FFT unit (N log2 N + N per beam, or per sub-pulse
 * with the four-step FFT). */
int mf_fft_init(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int method, int nThreads)
{
    int N = 1;
    while (N < pc->ddcCount) N *= 2;
    int log2N = 0;
    while ((1 << log2N) < N) log2N++;
    const int nChunks = (cfg->nBeams + MF_LANES - 1) / MF_LANES, four = cfg->nBeams == 1 && N >= 64;
    const double direct = 0.3 * pc->nTaps * pc->nRange * cfg->nBeams / team_size(nThreads, cfg->nBeams);
    const double fast = 0.6 * ((double)N * log2N + N) * (four ? 1.0 : (double)nChunks * MF_LANES / team_size(nThreads, nChunks));
    s->mfFft = method == MF_FFT || (method == MF_AUTO && fast < direct);
    if (!s->mfFft) return 0;
    s->mfN = N;
    s->mfChunks = nChunks;
    s->mfFour = four;
    if (fft_plan_init(&s->mfPlan, N) != 0) return -1;
    s->mfHre = proc_calloc(N, sizeof(float));
    s->mfHim = proc_calloc(N, sizeof(float));
    for (int j = 0; j < pc->nTaps; j++) {                           /* bit-reversed input      */
        s->mfHre[s->mfPlan.bitrev[j]] = pc->tapsRe[j] / N;
        s->mfHim[s->mfPlan.bitrev[j]] = pc->tapsIm[j] / N;
    }
    fft_rows(&s->mfPlan, s->mfHre, s->mfHim, 1, 1);
    if (four) {
        const int Q = N / 8;
        if (fft_plan_init(&s->mfPlanQ, Q) != 0 || fft_plan_init(&s->mfPlan8, 8) != 0) return -1;
        s->mfTwRe = proc_calloc(N, sizeof(float));
        s->mfTwIm = proc_calloc(N, sizeof(float));
        for (int n1 = 0; n1 < 8; n1++)
            for (int k2 = 0; k2 < Q; k2++) {
                const double a = -2.0 * PROCESSOR_PI * n1 * k2 / N;
                s->mfTwRe[n1 * Q + k2] = (float)cos(a);
                s->mfTwIm[n1 * Q + k2] = (float)sin(a);
            }
        s->mfWork = proc_calloc(4 * (size_t)N, sizeof(float));
        return 0;
    }
    s->mfWork = proc_calloc((size_t)s->mfChunks * 4 * N * MF_LANES, sizeof(float));
    return 0;
}

/* The same valid convolution by FFT: the circular convolution of length
 * mfN >= ddcCount equals the linear one at the valid outputs r + nTaps - 1.
 * Inverse FFT = conj(FFT(conj(.))) (1/mfN is in the taps). The beams are
 * processed in chunks of MF_LANES that are the fast index of the FFT rows. */
static void matched_filter_fft(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int slot, int nThreads)
{
    const int N = s->mfN, n = pc->ddcCount, L = pc->nTaps, nR = s->nR, nB = cfg->nBeams;
    const int *rev = s->mfPlan.bitrev;
    (void)nThreads;                                         /* without OpenMP */
    #pragma omp parallel for schedule(static) num_threads(team_size(nThreads, s->mfChunks)) if(s->mfChunks > 1)
    for (int c = 0; c < s->mfChunks; c++) {
        const int b0 = c * MF_LANES, nb = (nB - b0 < MF_LANES) ? nB - b0 : MF_LANES;
        float *Xr = s->mfWork + (size_t)c * 4 * N * MF_LANES;
        float *Xi = Xr + (size_t)N * MF_LANES, *Yr = Xi + (size_t)N * MF_LANES, *Yi = Yr + (size_t)N * MF_LANES;
        for (int m = 0; m < N; m++) {                           /* beams -> bit-reversed rows */
            float *dr = Xr + (size_t)rev[m] * MF_LANES, *di = Xi + (size_t)rev[m] * MF_LANES;
            for (int j = 0; j < MF_LANES; j++) {
                const int ok = m < n && j < nb;
                dr[j] = ok ? s->beamRe[(size_t)(b0 + j) * n + m] : 0.0f;
                di[j] = ok ? s->beamIm[(size_t)(b0 + j) * n + m] : 0.0f;
            }
        }
        fft_rows(&s->mfPlan, Xr, Xi, MF_LANES, MF_LANES);
        for (int k = 0; k < N; k++) {                           /* conj(X * H), bit-reversed  */
            const float hr = s->mfHre[k], hi = s->mfHim[k];
            const float *xr = Xr + (size_t)k * MF_LANES, *xi = Xi + (size_t)k * MF_LANES;
            float *yr = Yr + (size_t)rev[k] * MF_LANES, *yi = Yi + (size_t)rev[k] * MF_LANES;
            for (int j = 0; j < MF_LANES; j++) {
                yr[j] = xr[j] * hr - xi[j] * hi;
                yi[j] = -(xr[j] * hi + xi[j] * hr);
            }
        }
        fft_rows(&s->mfPlan, Yr, Yi, MF_LANES, MF_LANES);
        for (int j = 0; j < nb; j++) {                          /* valid outputs, conjugated  */
            float *outR = s->cpiRe + (size_t)slot * s->nCells + (size_t)(b0 + j) * nR;
            float *outI = s->cpiIm + (size_t)slot * s->nCells + (size_t)(b0 + j) * nR;
            const float *yr = Yr + (size_t)(L - 1) * MF_LANES + j, *yi = Yi + (size_t)(L - 1) * MF_LANES + j;
            for (int r = 0; r < nR; r++) {
                outR[r] = yr[(size_t)r * MF_LANES];
                outI[r] = -yi[(size_t)r * MF_LANES];
            }
        }
    }
}

/* Four-step FFT of length N = 8Q of one signal: A holds x[n1 + 8 n2] in row
 * bitrev_Q(n2), lane n1 ([Q][8]; overwritten). The 8 interleaved FFTs of
 * length Q run as the lanes of one vector, they are multiplied by the
 * twiddles exp(-2 pi j n1 k2 / N), and the FFTs of length 8 over the lanes
 * give B[k1 Q + k2] = X[k2 + Q k1]: X in natural order. */
static void fft_four_step(const PulseState *s, float *Ar, float *Ai, float *restrict Br, float *restrict Bi)
{
    const int Q = s->mfN / 8;
    const int *rev8 = s->mfPlan8.bitrev;
    fft_rows(&s->mfPlanQ, Ar, Ai, 8, 8);
    for (int k0 = 0; k0 < Q; k0 += 8) {                    /* transpose 8 x 8 blocks, twiddles */
        float tr[8][8], ti[8][8];
        for (int k = 0; k < 8; k++)
            for (int n1 = 0; n1 < 8; n1++) {
                tr[n1][k] = Ar[(k0 + k) * 8 + n1];
                ti[n1][k] = Ai[(k0 + k) * 8 + n1];
            }
        for (int n1 = 0; n1 < 8; n1++) {
            float *restrict br = Br + (size_t)rev8[n1] * Q + k0, *restrict bi = Bi + (size_t)rev8[n1] * Q + k0;
            const float *restrict wr = s->mfTwRe + (size_t)n1 * Q + k0, *restrict wi = s->mfTwIm + (size_t)n1 * Q + k0;
            for (int k = 0; k < 8; k++) {
                br[k] = tr[n1][k] * wr[k] - ti[n1][k] * wi[k];
                bi[k] = tr[n1][k] * wi[k] + ti[n1][k] * wr[k];
            }
        }
    }
    fft_rows(&s->mfPlan8, Br, Bi, Q, Q);
}

/* Fast convolution of the single beam (no DBF) with four-step FFTs, as
 * matched_filter_fft: X = FFT(x), conj(FFT(conj(X H))) at the valid outputs. */
static void matched_filter_four(const PulseConfig *pc, PulseState *s, int slot)
{
    const int N = s->mfN, Q = N / 8, n = pc->ddcCount, L = pc->nTaps, nR = s->nR;
    const int *revQ = s->mfPlanQ.bitrev;
    float *Ar = s->mfWork, *Ai = Ar + N, *Br = Ai + N, *Bi = Br + N;
    for (int n2 = 0; n2 < Q; n2++) {                        /* zero padded beam into the rows */
        float *dr = Ar + (size_t)revQ[n2] * 8, *di = Ai + (size_t)revQ[n2] * 8;
        const int m0 = 8 * n2, k = n - m0 < 0 ? 0 : n - m0 > 8 ? 8 : n - m0;
        if (k > 0) {
            memcpy(dr, s->beamRe + m0, sizeof(float) * k);
            memcpy(di, s->beamIm + m0, sizeof(float) * k);
        }
        memset(dr + k, 0, sizeof(float) * (8 - k));
        memset(di + k, 0, sizeof(float) * (8 - k));
    }
    fft_four_step(s, Ar, Ai, Br, Bi);
    for (int n2 = 0; n2 < Q; n2++) {                        /* conj(X * H) into the rows      */
        float *dr = Ar + (size_t)revQ[n2] * 8, *di = Ai + (size_t)revQ[n2] * 8;
        for (int n1 = 0; n1 < 8; n1++) {
            const int k = 8 * n2 + n1;
            const float xr = Br[k], xi = Bi[k], hr = s->mfHre[k], hi = s->mfHim[k];
            dr[n1] = xr * hr - xi * hi;
            di[n1] = -(xr * hi + xi * hr);
        }
    }
    fft_four_step(s, Ar, Ai, Br, Bi);
    float *outR = s->cpiRe + (size_t)slot * s->nCells, *outI = s->cpiIm + (size_t)slot * s->nCells;
    for (int r = 0; r < nR; r++) {                          /* valid outputs, conjugated      */
        outR[r] = Br[r + L - 1];
        outI[r] = -Bi[r + L - 1];
    }
}

/* Valid convolution of every beam with the matched filter, written into the
 * CPI buffer slot: y[r] = sum_j h[j] * x[r + nTaps - 1 - j]. */
void matched_filter(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int slot, int nThreads)
{
    if (s->mfFour) {
        matched_filter_four(pc, s, slot);
        return;
    }
    if (s->mfFft) {
        matched_filter_fft(cfg, pc, s, slot, nThreads);
        return;
    }
    const int n = pc->ddcCount, L = pc->nTaps, nR = s->nR;
    #pragma omp parallel for schedule(static) num_threads(team_size(nThreads, cfg->nBeams)) if(cfg->nBeams > 1)
    for (int b = 0; b < cfg->nBeams; b++) {
        const float *restrict xr = s->beamRe + (size_t)b * n, *restrict xi = s->beamIm + (size_t)b * n;
        float *restrict yr = s->cpiRe + (size_t)slot * s->nCells + (size_t)b * nR;
        float *restrict yi = s->cpiIm + (size_t)slot * s->nCells + (size_t)b * nR;
        memset(yr, 0, sizeof(float) * nR);
        memset(yi, 0, sizeof(float) * nR);
        for (int j = 0; j < L; j++) {
            const float hr = pc->tapsRe[j], hi = pc->tapsIm[j];
            const float *restrict ar = xr + (L - 1 - j), *restrict ai = xi + (L - 1 - j);
            for (int r = 0; r < nR; r++) {
                yr[r] += hr * ar[r] - hi * ai[r];
                yi[r] += hr * ai[r] + hi * ar[r];
            }
        }
    }
}

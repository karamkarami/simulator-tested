/*
 * processor.c - receiver processing chain of one PRI, for the short and the
 * long sub-pulse separately (same algorithms as the MATLAB simulator):
 *
 *   DDC output (int16 I/Q) -> DBF -> matched filter -> CPI buffer -> MTI
 *   -> Doppler (window, FFT, non-coherent integration of rows, noise-gain
 *   equalisation) -> 2-D CA-CFAR -> pruning -> detection reports -> peak
 *   interpolation (Doppler, elevation) -> plot extraction.
 *
 * Speed: long matched filters run as fast convolution (FFT, the valid
 * convolution within float rounding; a single beam as a four-step FFT whose
 * 8 sub-FFTs fill a vector); the MTI output of every pulse is computed once
 * and kept in a ring buffer stored by cell blocks; the Doppler FFTs of all
 * beam/range cells run together (cells are the fast index, so the loops
 * vectorise) in blocks that fit the L1 cache (radix-2^2); the CFAR works
 * beam by beam, or by range chunks with one beam (range sums from partial
 * sums in a ring of a few Doppler rows, running Doppler sums) so its data
 * stays in the cache, and it marks only the cells that are also local peaks
 * (the first pruning step); the noise floor median is selected among the few
 * values of a bracket around the previous median; the plot extractor finds
 * the groups near a detection in a hash grid; DBF, matched filter, Doppler,
 * CFAR and the noise floors are split over threads with OpenMP. Work buffers
 * are allocated once.
 *
 * Array layouts (float, real and imaginary parts in separate arrays):
 *   cell             = beam * nRange + rangeCell
 *   CPI buffer       = [slot][cell]
 *   MTI buffer       = [cellBlock][slot][CELL_BLOCK]   (cell = cellBlock * CELL_BLOCK + c)
 *   power map P      = [dopplerBin][cell]   (Doppler after fftshift)
 */
#include "radar_rx.h"
#include "extractor.h"
#include "tracker.h"
#include "cfar.h"
#include "fft.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(__AVX2__)
#include <immintrin.h>
#endif

#define PI 3.14159265358979323846
#define CELL_BLOCK 32           /* cells per Doppler FFT block (32 KB: fits the L1 cache) */
#define MF_LANES   24           /* beams per fast-convolution chunk (3 AVX registers) */
#define MEDIAN_SAMPLE 2048      /* noise floor: sample that brackets the median */

#if defined(__clang__)
#  define IVDEP _Pragma("clang loop vectorize(assume_safety)")
#elif defined(__GNUC__)
#  define IVDEP _Pragma("GCC ivdep")
#elif defined(_MSC_VER)
#  define IVDEP __pragma(loop(ivdep))
#else
#  define IVDEP
#endif

typedef struct {
    int nR, nCells, nSlots, nMtiOut;
    float *segRe, *segIm;       /* DDC samples of this sub-pulse [ch][n]          */
    float *beamRe, *beamIm;     /* beamformer output [beam][n]                    */
    float *cpiRe, *cpiIm;       /* matched filter outputs of the newest nMtiTaps
                                   pulses (all the MTI needs) [mfSlot][cell]      */
    float *mtiRe, *mtiIm;       /* MTI outputs [cell block][slot][CELL_BLOCK]     */
    int    last, count, nPushed;            /* CPI buffer state                   */
    int    mfLast;                          /* newest slot of cpiRe/cpiIm         */
    int    sinceUpdate;                     /* pulses since the last performed update */
    int    mtiLast, mtiCount;               /* MTI ring state                     */
    uint16_t *staggerIdx;       /* stagger position of every buffered pulse        */
    double   *az, *time;
    uint32_t newestPulse;
    float *P;                   /* power map [bin][cell]                          */
    unsigned char *det;         /* CFAR detections that are local peaks [bin][cell] */
    double *cntInv;             /* CFAR: 1 / training cells of every range cell (more at the edges) */
    double *thrScale;           /* CFAR: alpha * cntInv of every range cell       */
    float  *noise;              /* CFAR noise estimate of the detected cells [bin][cell]
                                   (cfg->snrNoise == SNR_NOISE_CFAR, else NULL)   */
    float  *lastMedian;         /* noise floor: median of every beam's map at its last
                                   computation (bracket of the next one; 0 = none) */
    float *pkM, *pkA;           /* dense sidelobe blanking: peak map, its maxima [bin][cell] */
    float *pkE, *pkD, *pkT;     /* maxima over the beams [bin][r], over Doppler [beam][r], temp */
    Detection *reports;         /* pruned detections of the latest CPI           */
    int nReports, capReports;
    Extractor ext;
    /* matched filter by fast convolution (used when cheaper than direct) */
    int    mfFft, mfN, mfChunks;
    int    mfFour;              /* one beam: FFT of length mfN as 8 FFTs of mfN/8 (four-step) */
    FftPlan mfPlan, mfPlanQ, mfPlan8;   /* lengths mfN; mfN/8 and 8 (four-step)  */
    float *mfHre, *mfHim;       /* FFT of the zero-padded taps, divided by mfN    */
    float *mfTwRe, *mfTwIm;     /* four-step twiddles exp(-2 pi j n1 k2 / mfN) [8][mfN/8] */
    float *mfWork;              /* per chunk: Xre, Xim, Yre, Yim [mfN][MF_LANES]
                                   (four-step: Are, Aim, Bre, Bim [mfN])           */
} PulseState;

typedef struct {                /* work buffers of one thread                      */
    float  *xr, *xi;            /* Doppler: one cell block [nfft][CELL_BLOCK]     */
    CfarScratch cfar;           /* CFAR (cfar.c)                                  */
    float  *work;               /* noise floor: median candidates [nfft*maxNR + 8] */
} Scratch;

struct Processor {
    const RadarConfig *cfg;
    FftPlan fft;
    PulseState ps[N_PULSES];
    int nThreads;
    Scratch *scratch;           /* [nThreads] */
    int slowTime;               /* 0: skip the Doppler/CFAR updates (overload)    */
    double alpha;               /* CFAR threshold factor in use                   */
    Tracker tracker;            /* EKF tracker (cfg->trackEnable)                 */
    Detection *newPlots;        /* plots of the latest extractor update           */
    int nNewPlots, capNewPlots;
    Detection *trackPlots;      /* plots of this PRI for the tracker              */
    int nTrackPlots, capTrackPlots;
    uint32_t lastPulse;
    double lastTime;            /* time of the newest PRI                          */
    int havePulse;
    ProcessorTiming timing;
    /* display taps (processor_set_profile, processor_set_rdmap) */
    int profileKind, profileBeam;
    int rdEnabled, rdBeam;
    double rdPeriod_s, rdLast;
    int nBins;                  /* range bins of the taps (dR from 0 to maxRange)  */
    double dR;
    float *profile;             /* [nBins]                                         */
    int *profileCnt;            /* [nBins]                                         */
    float *rd;                  /* [nfft][nBins]                                   */
    double *rdVel;              /* [nfft]                                          */
    Detection *rdDets;
    int nRdDets, capRdDets, rdUpdated;
    float *rdWork;
};

/* ------------------------------------------------------------------------ */

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n ? n : 1, size);
    if (!p) abort();
    return p;
}

static void *xrealloc(void *p, size_t size)
{
    p = realloc(p, size ? size : 1);
    if (!p) abort();
    return p;
}

static void reset_buffers(PulseState *s)
{
    s->last = -1;
    s->mfLast = -1;
    s->sinceUpdate = 0;
    s->count = 0;
    s->nPushed = 0;
    s->mtiLast = -1;
    s->mtiCount = 0;
}

static int thread_num(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* Threads for a parallel loop of nWork items: no more than there are items
 * (and no parallel region for one item: if(nWork > 1)).
 * Every thread of a team has to be woken and reach the closing barrier even
 * without work; when the other cores are busy (display, UDP input) that wait
 * costs more than the work. With one beam (no DBF) the loops over beams run
 * in the calling thread. */
static int team_size(int nThreads, int nWork)
{
    return nWork < nThreads ? (nWork > 1 ? nWork : 1) : nThreads;
}

/* Noise floor median: lane numbers of the set bits of every 8-bit mask, in
 * order, and their count (left-packing of 8 values with AVX2). The same
 * values are written by every processor_create. */
static uint8_t g_packIdx[256][8], g_packCnt[256];

static void init_pack_table(void)
{
    for (int m = 0; m < 256; m++) {
        int k = 0;
        for (int lane = 0; lane < 8; lane++)
            if (m >> lane & 1) g_packIdx[m][k++] = (uint8_t)lane;
        g_packCnt[m] = (uint8_t)k;
    }
}

/* Fast convolution of the matched filter: FFT length and the FFT of the taps.
 * It is used when its estimated time is shorter than that of the direct
 * convolution (which is split over the threads beam by beam; the FFTs of
 * MF_LANES beams run in one thread; a single beam is split into the 8
 * interleaved sub-FFTs of a four-step FFT instead, which fill the vectors).
 * Costs per operation measured on an AVX2 Xeon: 0.3 ns per complex MAC
 * (direct), 0.6 ns per FFT unit (N log2 N + N per beam, or per sub-pulse
 * with the four-step FFT). */
static int mf_fft_init(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int method, int nThreads)
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
    s->mfHre = xcalloc(N, sizeof(float));
    s->mfHim = xcalloc(N, sizeof(float));
    for (int j = 0; j < pc->nTaps; j++) {                           /* bit-reversed input      */
        s->mfHre[s->mfPlan.bitrev[j]] = pc->tapsRe[j] / N;
        s->mfHim[s->mfPlan.bitrev[j]] = pc->tapsIm[j] / N;
    }
    fft_rows(&s->mfPlan, s->mfHre, s->mfHim, 1, 1);
    if (four) {
        const int Q = N / 8;
        if (fft_plan_init(&s->mfPlanQ, Q) != 0 || fft_plan_init(&s->mfPlan8, 8) != 0) return -1;
        s->mfTwRe = xcalloc(N, sizeof(float));
        s->mfTwIm = xcalloc(N, sizeof(float));
        for (int n1 = 0; n1 < 8; n1++)
            for (int k2 = 0; k2 < Q; k2++) {
                const double a = -2.0 * PI * n1 * k2 / N;
                s->mfTwRe[n1 * Q + k2] = (float)cos(a);
                s->mfTwIm[n1 * Q + k2] = (float)sin(a);
            }
        s->mfWork = xcalloc(4 * (size_t)N, sizeof(float));
        return 0;
    }
    s->mfWork = xcalloc((size_t)s->mfChunks * 4 * N * MF_LANES, sizeof(float));
    return 0;
}

Processor *processor_create(const RadarConfig *cfg, int nThreads)
{
#ifdef _OPENMP
    /* default: OMP_NUM_THREADS if set, else all cores but one (left for the
     * UDP input and the display) */
    if (nThreads <= 0 && !getenv("OMP_NUM_THREADS"))
        nThreads = omp_get_num_procs() > 1 ? omp_get_num_procs() - 1 : 1;
    if (nThreads > 0) omp_set_num_threads(nThreads);
    const int nThr = omp_get_max_threads();
#else
    (void)nThreads;
    const int nThr = 1;
#endif
    Processor *p = xcalloc(1, sizeof(Processor));
    p->cfg = cfg;
    if (fft_plan_init(&p->fft, cfg->nfft) != 0) { free(p); return NULL; }
    int maxNR = 1;
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        PulseState *s = &p->ps[i];
        if (!pc->enable) {                         /* sub-pulse switched off: no buffers */
            extractor_init(&s->ext, cfg->nfft, cfg->minHits, cfg->maxMiss);
            s->ext.nBeams = cfg->nBeams;
            reset_buffers(s);
            continue;
        }
        s->nR = pc->nRange;
        if (s->nR > maxNR) maxNR = s->nR;
        s->nCells = cfg->nBeams * s->nR;
        s->nSlots = cfg->cpiLength;
        s->nMtiOut = cfg->cpiLength - cfg->nMtiTaps + 1;
        size_t seg = (size_t)cfg->nCh * pc->ddcCount, beams = (size_t)cfg->nBeams * pc->ddcCount;
        size_t cpi = (size_t)cfg->nMtiTaps * s->nCells;
        size_t mti = (size_t)(s->nCells + CELL_BLOCK - 1) / CELL_BLOCK * s->nMtiOut * CELL_BLOCK;
        size_t map = (size_t)cfg->nfft * s->nCells;
        s->segRe = xcalloc(seg, sizeof(float));   s->segIm = xcalloc(seg, sizeof(float));
        s->beamRe = xcalloc(beams, sizeof(float)); s->beamIm = xcalloc(beams, sizeof(float));
        s->cpiRe = xcalloc(cpi, sizeof(float));   s->cpiIm = xcalloc(cpi, sizeof(float));
        s->mtiRe = xcalloc(mti, sizeof(float));   s->mtiIm = xcalloc(mti, sizeof(float));
        s->staggerIdx = xcalloc(s->nSlots, sizeof(uint16_t));
        s->az = xcalloc(s->nSlots, sizeof(double));
        s->time = xcalloc(s->nSlots, sizeof(double));
        s->P = xcalloc(map, sizeof(float));
        s->det = xcalloc(map, 1);
        s->cntInv = xcalloc(s->nR, sizeof(double));
        s->thrScale = xcalloc(s->nR, sizeof(double));
        if (cfg->snrNoise == SNR_NOISE_CFAR) s->noise = xcalloc(map, sizeof(float));
        s->lastMedian = xcalloc(cfg->nBeams, sizeof(float));
        {
            const int gR = cfg->guardRange, oR = gR + cfg->trainRange;
            const int gD = cfg->guardDoppler, oD = gD + cfg->trainDoppler;
            for (int r = 0; r < s->nR; r++) {          /* guard cells cut at the range edges */
                const int g0 = r - gR > 0 ? r - gR : 0, g1 = r + gR < s->nR - 1 ? r + gR : s->nR - 1;
                const int cnt = (2 * oR + 1) * (2 * oD + 1) - (g1 - g0 + 1) * (2 * gD + 1);
                s->cntInv[r] = 1.0 / cnt;
            }
        }
        if (mf_fft_init(cfg, pc, s, cfg->mfMethod, nThr) != 0) { processor_destroy(p); return NULL; }
        extractor_init(&s->ext, cfg->nfft, cfg->minHits, cfg->maxMiss);
        s->ext.nBeams = cfg->nBeams;
        reset_buffers(s);
    }
    p->nThreads = nThr;
    p->slowTime = 1;
    p->alpha = cfg->alpha;
    if (cfg->trackEnable) tracker_init(&p->tracker, cfg);
    init_pack_table();
    p->scratch = xcalloc(nThr, sizeof(Scratch));
    for (int t = 0; t < nThr; t++) {
        Scratch *w = &p->scratch[t];
        w->xr = xcalloc((size_t)cfg->nfft * CELL_BLOCK, sizeof(float));
        w->xi = xcalloc((size_t)cfg->nfft * CELL_BLOCK, sizeof(float));
        cfar_scratch_init(&w->cfar, cfg, maxNR);
        w->work = xcalloc((size_t)cfg->nfft * maxNR + 8, sizeof(float));
    }
    return p;
}

void processor_destroy(Processor *p)
{
    if (!p) return;
    free(p->profile); free(p->profileCnt); free(p->rd); free(p->rdVel); free(p->rdDets); free(p->rdWork);
    for (int i = 0; i < N_PULSES; i++) {
        PulseState *s = &p->ps[i];
        free(s->segRe); free(s->segIm); free(s->beamRe); free(s->beamIm);
        free(s->cpiRe); free(s->cpiIm); free(s->mtiRe); free(s->mtiIm);
        free(s->staggerIdx); free(s->az); free(s->time);
        free(s->P); free(s->det); free(s->cntInv); free(s->thrScale); free(s->noise); free(s->lastMedian);
        free(s->pkM); free(s->pkA); free(s->pkE); free(s->pkD); free(s->pkT);
        free(s->reports);
        if (s->mfFft) fft_plan_free(&s->mfPlan);
        if (s->mfFour) { fft_plan_free(&s->mfPlanQ); fft_plan_free(&s->mfPlan8); }
        free(s->mfHre); free(s->mfHim); free(s->mfTwRe); free(s->mfTwIm); free(s->mfWork);
        extractor_free(&s->ext);
    }
    for (int t = 0; p->scratch && t < p->nThreads; t++) {
        Scratch *w = &p->scratch[t];
        free(w->xr); free(w->xi); free(w->work);
        cfar_scratch_free(&w->cfar);
    }
    free(p->scratch);
    fft_plan_free(&p->fft);
    tracker_free(&p->tracker);
    free(p->newPlots);
    free(p->trackPlots);
    free(p);
}

void processor_set_slow_time(Processor *p, int enabled)
{
    p->slowTime = enabled;
}

void processor_set_cfar_alpha(Processor *p, double alpha)
{
    if (alpha > 0.0 && isfinite(alpha)) p->alpha = alpha;
}

double processor_cfar_alpha(const Processor *p)
{
    return p->alpha;
}

void processor_timing(const Processor *p, ProcessorTiming *t)
{
    *t = p->timing;
}

int processor_threads(const Processor *p)
{
    return p->nThreads;
}

/* ------------------------------------------------------------------------ */
/* Fast time: DBF and matched filter                                        */
/* ------------------------------------------------------------------------ */

/* DDC samples of the sub-pulse, as float [ch][n]; the samples outside its
 * receive interval (transmission, other sub-pulse, end of the record) are zero. */
static void load_segment(const RadarConfig *cfg, const PulseConfig *pc, const PriData *pri, PulseState *s)
{
    const int n = pc->ddcCount;
    const int k0 = pc->rxFirst > pc->ddcFirst ? pc->rxFirst - pc->ddcFirst : 0;      /* first inside   */
    const int rxEnd = pc->rxFirst + pc->rxCount - pc->ddcFirst;
    const int k1 = rxEnd < n ? rxEnd : n;                                            /* first after it */
    for (int ch = 0; ch < cfg->nCh; ch++) {
        const int16_t *src = pri->iq + 2 * ((size_t)ch * pri->nSamples + pc->ddcFirst + k0);   /* sample k0 */
        float *re = s->segRe + (size_t)ch * n, *im = s->segIm + (size_t)ch * n;
        memset(re, 0, sizeof(float) * k0);
        memset(im, 0, sizeof(float) * k0);
        for (int k = k0; k < k1; k++) {
            re[k] = src[2 * (k - k0)];
            im[k] = src[2 * (k - k0) + 1];
        }
        memset(re + k1, 0, sizeof(float) * (n - k1));
        memset(im + k1, 0, sizeof(float) * (n - k1));
    }
}

/* beam b = sum over channels of W[b][ch] * x[ch]  (or the single channel).
 * Three beams are updated per pass over a channel (fewer loads). */
static void beamform(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int nThreads)
{
    const int n = pc->ddcCount, nB = cfg->nBeams, nCh = cfg->nCh;
    (void)nThreads;                                         /* without OpenMP */
    if (!cfg->dbfEnabled) {
        memcpy(s->beamRe, s->segRe, sizeof(float) * n);
        memcpy(s->beamIm, s->segIm, sizeof(float) * n);
        return;
    }
    const int nGroups = (nB + 2) / 3;
    #pragma omp parallel for schedule(static) num_threads(team_size(nThreads, nGroups)) if(nGroups > 1)
    for (int g = 0; g < nGroups; g++) {
        const int b0 = 3 * g, nb = nB - b0 < 3 ? nB - b0 : 3;
        float *y0r = s->beamRe + (size_t)b0 * n, *y0i = s->beamIm + (size_t)b0 * n;
        memset(y0r, 0, sizeof(float) * nb * n);
        memset(y0i, 0, sizeof(float) * nb * n);
        if (nb < 3) {                                       /* last 1 or 2 beams */
            for (int b = b0; b < nB; b++) {
                float *restrict yr = s->beamRe + (size_t)b * n, *restrict yi = s->beamIm + (size_t)b * n;
                for (int ch = 0; ch < nCh; ch++) {
                    const float wr = cfg->wRe[b * nCh + ch], wi = cfg->wIm[b * nCh + ch];
                    const float *restrict xr = s->segRe + (size_t)ch * n, *restrict xi = s->segIm + (size_t)ch * n;
                    for (int k = 0; k < n; k++) {
                        yr[k] += wr * xr[k] - wi * xi[k];
                        yi[k] += wr * xi[k] + wi * xr[k];
                    }
                }
            }
            continue;
        }
        float *restrict a0r = y0r, *restrict a0i = y0i, *restrict a1r = y0r + n, *restrict a1i = y0i + n;
        float *restrict a2r = y0r + 2 * (size_t)n, *restrict a2i = y0i + 2 * (size_t)n;
        for (int ch = 0; ch < nCh; ch++) {
            const float w0r = cfg->wRe[b0 * nCh + ch], w0i = cfg->wIm[b0 * nCh + ch];
            const float w1r = cfg->wRe[(b0 + 1) * nCh + ch], w1i = cfg->wIm[(b0 + 1) * nCh + ch];
            const float w2r = cfg->wRe[(b0 + 2) * nCh + ch], w2i = cfg->wIm[(b0 + 2) * nCh + ch];
            const float *restrict xr = s->segRe + (size_t)ch * n, *restrict xi = s->segIm + (size_t)ch * n;
            IVDEP
            for (int k = 0; k < n; k++) {
                const float r = xr[k], i = xi[k];
                a0r[k] += w0r * r - w0i * i;  a0i[k] += w0r * i + w0i * r;
                a1r[k] += w1r * r - w1i * i;  a1i[k] += w1r * i + w1i * r;
                a2r[k] += w2r * r - w2i * i;  a2i[k] += w2r * i + w2i * r;
            }
        }
    }
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
static void matched_filter(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int slot, int nThreads)
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

/* ------------------------------------------------------------------------ */
/* Slow time: MTI and Doppler                                               */
/* ------------------------------------------------------------------------ */

/* Row of the MTI ring: the cells c0 ... c0+CELL_BLOCK-1 of cell block blk in
 * slot slot ([block][slot][CELL_BLOCK]: the Doppler processing of a block
 * reads contiguous memory). */
static float *mti_row(const PulseState *s, float *ring, int blk, int slot)
{
    return ring + ((size_t)blk * s->nMtiOut + slot) * CELL_BLOCK;
}

/* MTI output of the newest L pulses: sum_l c[stagger of the oldest][l] * x[oldest + l]. */
static void mti_newest(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s)
{
    const int L = cfg->nMtiTaps, N = s->nSlots;
    const int first = (s->last - (L - 1) + N) % N;               /* oldest pulse of the output */
    const double *c = cfg->mtiCoeffs + (size_t)(s->staggerIdx[first] - 1) * L;
    const int firstMf = (s->mfLast + 1) % L;                     /* its matched filter output  */
    const int nBlocks = (s->nCells + CELL_BLOCK - 1) / CELL_BLOCK;
    s->mtiLast = (s->mtiLast + 1) % s->nMtiOut;
    for (int blk = 0; blk < nBlocks; blk++) {
        const int c0 = blk * CELL_BLOCK, nc = s->nCells - c0 < CELL_BLOCK ? s->nCells - c0 : CELL_BLOCK;
        float *restrict yr = mti_row(s, s->mtiRe, blk, s->mtiLast);
        float *restrict yi = mti_row(s, s->mtiIm, blk, s->mtiLast);
        memset(yr, 0, sizeof(float) * CELL_BLOCK);
        memset(yi, 0, sizeof(float) * CELL_BLOCK);
        for (int l = 0; l < L; l++) {
            const float w = (float)c[l];
            const int slot = (firstMf + l) % L;
            const float *restrict xr = s->cpiRe + (size_t)slot * s->nCells + c0;
            const float *restrict xi = s->cpiIm + (size_t)slot * s->nCells + c0;
            for (int k = 0; k < nc; k++) {
                yr[k] += w * xr[k];
                yi[k] += w * xi[k];
            }
        }
        if (pc->mtiCells < s->nR) {              /* beyond cfg.mti.maxRange_m: the newest pulse, unfiltered */
            const int slot = (firstMf + L - 1) % L;
            const float *xr = s->cpiRe + (size_t)slot * s->nCells + c0;
            const float *xi = s->cpiIm + (size_t)slot * s->nCells + c0;
            for (int k = 0; k < nc; k++)
                if ((c0 + k) % s->nR >= pc->mtiCells) { yr[k] = xr[k]; yi[k] = xi[k]; }
        }
    }
    if (s->mtiCount < s->nMtiOut) s->mtiCount++;
}

/* Range-Doppler power map P[bin][cell] from the newest nRows*rowLen MTI outputs:
 * rows (fixed PRF: one row of consecutive outputs; staggered: de-interleaved by
 * stagger position), Hamming window, FFT, |X|^2 averaged over the rows, fftshift
 * and division by the noise gain of every bin. */
static void doppler(const Processor *p, const PulseConfig *pc, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const FftPlan *plan = &p->fft;
    const int nfft = cfg->nfft, nRows = cfg->nRows, rowLen = cfg->rowLen;
    const int nUsed = nRows * rowLen;
    const int nBlocks = (s->nCells + CELL_BLOCK - 1) / CELL_BLOCK;
    const float scale = 1.0f / nRows;

    #pragma omp parallel num_threads(team_size(p->nThreads, nBlocks)) if(nBlocks > 1)
    {
        float *xr = p->scratch[thread_num()].xr;
        float *xi = p->scratch[thread_num()].xi;
        #pragma omp for schedule(static)
        for (int blk = 0; blk < nBlocks; blk++) {
            const int c0 = blk * CELL_BLOCK;
            const int nc = (s->nCells - c0 < CELL_BLOCK) ? s->nCells - c0 : CELL_BLOCK;
            if (nRows > 1)
                for (int d = 0; d < nfft; d++)
                    memset(s->P + (size_t)d * s->nCells + c0, 0, sizeof(float) * nc);
            for (int row = 0; row < nRows; row++) {
                /* load the row in bit-reversed order, windowed, zero padded */
                for (int j = 0; j < nfft; j++) {
                    float *dr = xr + (size_t)plan->bitrev[j] * CELL_BLOCK;
                    float *di = xi + (size_t)plan->bitrev[j] * CELL_BLOCK;
                    if (j >= rowLen) {
                        memset(dr, 0, sizeof(float) * nc);
                        memset(di, 0, sizeof(float) * nc);
                        continue;
                    }
                    const int q = j * nRows + row;                       /* chronological index */
                    const int slot = (s->mtiLast - (nUsed - 1) + q + 2 * s->nMtiOut) % s->nMtiOut;
                    const float w = cfg->dopWindow[j];
                    const float *restrict sr = mti_row(s, s->mtiRe, blk, slot);
                    const float *restrict si = mti_row(s, s->mtiIm, blk, slot);
                    for (int c = 0; c < nc; c++) {
                        dr[c] = w * sr[c];
                        di[c] = w * si[c];
                    }
                }
                fft_rows(plan, xr, xi, CELL_BLOCK, nc);
                for (int d = 0; d < nfft; d++) {                       /* fftshift + power */
                    const int sh = (d + nfft / 2) % nfft;
                    float *restrict pd = s->P + (size_t)sh * s->nCells + c0;
                    const float *restrict ar = xr + (size_t)d * CELL_BLOCK, *restrict ai = xi + (size_t)d * CELL_BLOCK;
                    if (nRows == 1) {                                  /* one row: final value */
                        const float g = 1.0f / cfg->noiseGain[sh];
                        for (int c = 0; c < nc; c++) pd[c] = (ar[c] * ar[c] + ai[c] * ai[c]) * g;
                    } else {
                        for (int c = 0; c < nc; c++) pd[c] += ar[c] * ar[c] + ai[c] * ai[c];
                    }
                }
            }
            for (int d = 0; d < nfft && nRows > 1; d++) {          /* mean over rows, noise gain */
                const float g = scale / cfg->noiseGain[d];
                float *restrict pd = s->P + (size_t)d * s->nCells + c0;
                for (int c = 0; c < nc; c++) pd[c] *= g;
            }
            if (pc->mtiCells < s->nR) {           /* unfiltered cells: their own noise gain */
                for (int d = 0; d < nfft; d++) {
                    const float g = cfg->noiseGain[d] / cfg->noiseGainBypass[d];
                    float *restrict pd = s->P + (size_t)d * s->nCells + c0;
                    for (int c = 0; c < nc; c++)
                        if ((c0 + c) % s->nR >= pc->mtiCells) pd[c] *= g;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* CFAR (cfar.c)                                                            */
/* ------------------------------------------------------------------------ */

/* The work items are the beams, or, with fewer beams than threads (one
 * element), chunks of tested range cells of at least CFAR_MIN_CHUNK cells:
 * each chunk sums the range cells its windows need, so its sums and
 * decisions are those of the whole beam. */
#define CFAR_MIN_CHUNK 64

static void cfar(const Processor *p, PulseState *s, const PulseConfig *pc)
{
    const RadarConfig *cfg = p->cfg;
    const int nB = cfg->nBeams;
    const int f = pc->testFirst, l = pc->testFirst + pc->testCount - 1, nTest = l - f + 1;

    memset(s->det, 0, (size_t)cfg->nfft * s->nCells);
    if (l < f) return;
    for (int r = f; r <= l; r++) s->thrScale[r] = p->alpha * s->cntInv[r];
    int notchLastCell = -1;
    for (int r = 0; r < pc->nRange && pc->range_m[r] <= cfg->notchMaxRange_m; r++) notchLastCell = r;
    const CfarMap map = { s->P, s->det, s->noise, s->cntInv, s->thrScale, s->nR, s->nCells, nB, cfg->nfft };
    int nChunks = 1;
    if (nB < p->nThreads) {
        const int most = nTest / CFAR_MIN_CHUNK > 1 ? nTest / CFAR_MIN_CHUNK : 1;
        nChunks = (p->nThreads + nB - 1) / nB;
        if (nChunks > most) nChunks = most;
    }
    const int nItems = nB * nChunks;
    #pragma omp parallel num_threads(team_size(p->nThreads, nItems)) if(nItems > 1)
    {
        Scratch *w = &p->scratch[thread_num()];
        #pragma omp for schedule(dynamic)
        for (int it = 0; it < nItems; it++) {
            const int c = it % nChunks;
            cfar_chunk(cfg, &map, p->alpha, it / nChunks, f + (int)((long long)nTest * c / nChunks),
                       f + (int)((long long)nTest * (c + 1) / nChunks) - 1, notchLastCell, &w->cfar);
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Pruning and detection reports                                            */
/* ------------------------------------------------------------------------ */

typedef struct { int b, r, d; float pw; } Cell;

/* MATLAB's find() order: Doppler, then range, then beam. */
static int cell_order(const void *a, const void *b)
{
    const Cell *x = a, *y = b;
    if (x->d != y->d) return x->d - y->d;
    if (x->r != y->r) return x->r - y->r;
    return x->b - y->b;
}

typedef struct { double pw; int i; } Ranked;
static int ranked_desc(const void *a, const void *b)
{
    const double x = ((const Ranked *)a)->pw, y = ((const Ranked *)b)->pw;
    return (x < y) - (x > y);
}

static float P_at(const PulseState *s, int d, int b, int r)
{
    return s->P[(size_t)d * s->nCells + (size_t)b * s->nR + r];
}

/* k-th smallest value (0-based) of v[0..n-1]; v is reordered. */
static float select_kth(float *v, int n, int k)
{
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        float pivot = v[(lo + hi) / 2];
        int i = lo, j = hi;
        while (i <= j) {
            while (v[i] < pivot) i++;
            while (v[j] > pivot) j--;
            if (i <= j) { float t = v[i]; v[i] = v[j]; v[j] = t; i++; j--; }
        }
        if (k <= j) hi = j; else if (k >= i) lo = i; else return v[k];
    }
    return v[k];
}

/* Appends the values of x[0 ... n-1] inside [lo, hi] to out + *cnt (advancing
 * *cnt) and adds the number of values below lo to *below. out needs room for
 * 8 more floats than it gets values (AVX2 stores whole vectors: the values
 * are left-packed with a permutation). */
static void bracket_pass(const float *restrict x, int n, float lo, float hi, float *restrict out, int *cnt,
                         int *below)
{
    int c = *cnt, b = 0, i = 0;
#if defined(__AVX2__)
    const __m256 vlo = _mm256_set1_ps(lo), vhi = _mm256_set1_ps(hi);
    __m256i vb = _mm256_setzero_si256();
    for (; i + 8 <= n; i += 8) {
        const __m256 v = _mm256_loadu_ps(x + i);
        const __m256 lt = _mm256_cmp_ps(v, vlo, _CMP_LT_OQ);
        const int m = _mm256_movemask_ps(_mm256_andnot_ps(lt, _mm256_cmp_ps(v, vhi, _CMP_LE_OQ)));
        vb = _mm256_sub_epi32(vb, _mm256_castps_si256(lt));            /* below: lanes of -1 */
        const __m256i idx = _mm256_cvtepu8_epi32(_mm_loadl_epi64((const __m128i *)g_packIdx[m]));
        _mm256_storeu_ps(out + c, _mm256_permutevar8x32_ps(v, idx));
        c += g_packCnt[m];
    }
    int lanes[8];
    _mm256_storeu_si256((__m256i *)lanes, vb);
    for (int k = 0; k < 8; k++) b += lanes[k];
#endif
    for (; i < n; i++) {
        const float v = x[i];
        b += v < lo;
        out[c] = v;
        c += (v >= lo) & (v <= hi);
    }
    *cnt = c;
    *below += b;
}

/* Exact median of the nRows x rowLen values (row d at base + d*stride), as
 * MATLAB's median (mean of the two middle values when their number is even).
 * The middle values are selected among the few values of a bracket: one pass
 * counts the values below it and copies those inside it. The bracket is
 * +-1.6 % around the median of the previous update (*guess, 0 = none; the
 * map changes little from one update to the next: about 1 % of the values),
 * else, or when it misses the middle ranks, the range of +-3.5 standard
 * deviations of the middle ranks in a strided sample of MEDIAN_SAMPLE values
 * (about 8 % of the values). If that misses them too (an unusual map), all
 * values are searched: the result never depends on the bracket. *guess gets
 * the new median. work: nRows*rowLen + 8 floats. */
static double map_median(const float *base, int nRows, int rowLen, size_t stride, float *work, float *guess)
{
    const int n = nRows * rowLen, kHi = n / 2, kLo = (n - 1) / 2;    /* ranks of the middle values */
    int below = 0, cnt = 0, found = 0;
    if (*guess > 0.0f) {
        const float lo = *guess * (1.0f - 1.0f / 64), hi = *guess * (1.0f + 1.0f / 64);
        for (int d = 0; d < nRows; d++) bracket_pass(base + (size_t)d * stride, rowLen, lo, hi, work, &cnt, &below);
        found = below <= kLo && kHi < below + cnt;
    }
    if (!found) {
        const int S = n < MEDIAN_SAMPLE ? n : MEDIAN_SAMPLE;
        for (int i = 0; i < S; i++) {
            const long long q = (long long)i * n / S;
            work[i] = base[(size_t)(q / rowLen) * stride + (size_t)(q % rowLen)];
        }
        const int delta = (int)(1.75 * sqrt((double)S)) + 2;
        const int sLo = (int)((long long)kLo * S / n) - delta, sHi = (int)((long long)kHi * S / n) + 1 + delta;
        const float lo = sLo <= 0 ? -FLT_MAX : select_kth(work, S, sLo);
        const float hi = sHi >= S - 1 ? FLT_MAX : select_kth(work, S, sHi);
        below = cnt = 0;
        for (int d = 0; d < nRows; d++) bracket_pass(base + (size_t)d * stride, rowLen, lo, hi, work, &cnt, &below);
        found = below <= kLo && kHi < below + cnt;
    }
    if (!found) {
        cnt = below = 0;
        for (int d = 0; d < nRows; d++, cnt += rowLen)
            memcpy(work + cnt, base + (size_t)d * stride, sizeof(float) * rowLen);
    }
    const float vHi = select_kth(work, cnt, kHi - below);   /* work[0 ... kHi-below-1] <= vHi */
    *guess = vHi;
    if (kLo == kHi) return vHi;
    float vLo = work[0];
    for (int i = 1; i < kHi - below; i++) vLo = work[i] > vLo ? work[i] : vLo;
    return 0.5 * ((double)vLo + vHi);
}

/* Noise floor of a beam: median of its range-Doppler map / (median/mean of
 * the noise). */
static double beam_noise_floor(const RadarConfig *cfg, PulseState *s, int b, Scratch *w)
{
    return map_median(s->P + (size_t)b * s->nR, cfg->nfft, s->nR, (size_t)s->nCells, w->work, &s->lastMedian[b]) /
           cfg->medianToMean;
}

static void add_report(PulseState *s, const Detection *d)
{
    if (s->nReports == s->capReports) {
        s->capReports = s->capReports ? 2 * s->capReports : 64;
        s->reports = xrealloc(s->reports, sizeof(Detection) * s->capReports);
    }
    s->reports[s->nReports++] = *d;
}

/* y[j] = max of x over j-L ... j+L (zero outside the row; x >= 0): van Herk /
 * Gil-Werman, blocks of the window length with prefix and suffix maxima, so
 * the cost does not depend on L. buf: 2 * (n + 4L + 2) floats. */
static void range_max_row(const float *x, float *y, int n, int L, float *buf)
{
    if (L <= 0) { memcpy(y, x, sizeof(float) * n); return; }
    const int w = 2 * L + 1, nb = (n + 2 * L + w - 1) / w, len = nb * w;
    float *g = buf, *h = buf + len;
    for (int k = 0; k < nb; k++) {
        const int q0 = k * w;
        float m = 0.0f;
        for (int q = q0; q < q0 + w; q++) {
            const float v = (q >= L && q < L + n) ? x[q - L] : 0.0f;
            if (v > m) m = v;
            g[q] = m;
        }
        m = 0.0f;
        for (int q = q0 + w - 1; q >= q0; q--) {
            const float v = (q >= L && q < L + n) ? x[q - L] : 0.0f;
            if (v > m) m = v;
            h[q] = m;
        }
    }
    for (int j = 0; j < n; j++) y[j] = h[j] > g[j + w - 1] ? h[j] : g[j + w - 1];
}

/* Sidelobe blanking with maxima over the whole map (as prune_peaks.m with many
 * peaks): a peak is blanked when the largest peak power of one of its three
 * neighbourhoods is more than the margin above it. The cost depends on the
 * map size, not on the number of peaks. keep[] gets 0 for the blanked peaks. */
static void blank_dense(const Processor *p, PulseState *s, const Cell *cells, int n, int pulseCells,
                        double ratio, unsigned char *keep)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams, nR = s->nR;
    const size_t nCells = (size_t)s->nCells, map = (size_t)nfft * nCells;
    const int rows = nfft > nB ? nfft : nB;
    if (!s->pkM) {
        s->pkM = xcalloc(map, sizeof(float));
        s->pkA = xcalloc(map, sizeof(float));
        s->pkE = xcalloc((size_t)nfft * nR, sizeof(float));
        s->pkD = xcalloc((size_t)nB * nR, sizeof(float));
        s->pkT = xcalloc((size_t)rows * nR, sizeof(float));
    }
    float *M = s->pkM, *A = s->pkA, *E = s->pkE, *D = s->pkD, *T = s->pkT;
    memset(M, 0, sizeof(float) * map);
    for (int i = 0; i < n; i++) M[(size_t)cells[i].d * nCells + (size_t)cells[i].b * nR + cells[i].r] = cells[i].pw;
    const size_t bufLen = 2 * ((size_t)nR + 4 * (size_t)pulseCells + 2) + 8;

    /* range sidelobes: beam +-1, Doppler +-1 (circular), range +-pulseCells */
    #pragma omp parallel num_threads(team_size(p->nThreads, nfft)) if(nfft > 1)
    {
        float *row = xcalloc(nR, sizeof(float)), *buf = xcalloc(bufLen, sizeof(float));
        #pragma omp for schedule(static)
        for (int d = 0; d < nfft; d++) {
            for (int b = 0; b < nB; b++) {
                for (int r = 0; r < nR; r++) row[r] = 0.0f;
                for (int dd = -1; dd <= 1; dd++) {
                    const float *Md = M + (size_t)((d + dd + nfft) % nfft) * nCells;
                    for (int bb = b - 1; bb <= b + 1; bb++) {
                        if (bb < 0 || bb >= nB) continue;
                        const float *x = Md + (size_t)bb * nR;
                        for (int r = 0; r < nR; r++) if (x[r] > row[r]) row[r] = x[r];
                    }
                }
                range_max_row(row, A + (size_t)d * nCells + (size_t)b * nR, nR, pulseCells, buf);
            }
        }
        free(row); free(buf);
    }
    float *buf = xcalloc(bufLen, sizeof(float));
    /* elevation sidelobes: any beam, Doppler +-1, range +-1 */
    for (int d = 0; d < nfft; d++) {
        float *t = T + (size_t)d * nR;
        for (int r = 0; r < nR; r++) t[r] = 0.0f;
        for (int b = 0; b < nB; b++) {
            const float *x = M + (size_t)d * nCells + (size_t)b * nR;
            for (int r = 0; r < nR; r++) if (x[r] > t[r]) t[r] = x[r];
        }
    }
    for (int d = 0; d < nfft; d++) {
        float *e = E + (size_t)d * nR;
        const float *t0 = T + (size_t)((d + nfft - 1) % nfft) * nR, *t1 = T + (size_t)d * nR,
                    *t2 = T + (size_t)((d + 1) % nfft) * nR;
        for (int r = 0; r < nR; r++) {
            float m = t0[r] > t1[r] ? t0[r] : t1[r];
            e[r] = m > t2[r] ? m : t2[r];
        }
    }
    for (int d = 0; d < nfft; d++) {                                /* range +-1 */
        memcpy(T + (size_t)d * nR, E + (size_t)d * nR, sizeof(float) * nR);
        range_max_row(T + (size_t)d * nR, E + (size_t)d * nR, nR, 1, buf);
    }
    /* Doppler sidelobes: any Doppler bin, beam +-1, range +-1 */
    for (int b = 0; b < nB; b++) {
        float *t = T + (size_t)b * nR;
        for (int r = 0; r < nR; r++) t[r] = 0.0f;
        for (int d = 0; d < nfft; d++) {
            const float *x = M + (size_t)d * nCells + (size_t)b * nR;
            for (int r = 0; r < nR; r++) if (x[r] > t[r]) t[r] = x[r];
        }
    }
    for (int b = 0; b < nB; b++) {
        float *dst = D + (size_t)b * nR;
        for (int r = 0; r < nR; r++) {
            float m = T[(size_t)b * nR + r];
            if (b > 0 && T[(size_t)(b - 1) * nR + r] > m) m = T[(size_t)(b - 1) * nR + r];
            if (b < nB - 1 && T[(size_t)(b + 1) * nR + r] > m) m = T[(size_t)(b + 1) * nR + r];
            dst[r] = m;
        }
    }
    for (int b = 0; b < nB; b++) {                                  /* range +-1 */
        memcpy(T + (size_t)b * nR, D + (size_t)b * nR, sizeof(float) * nR);
        range_max_row(T + (size_t)b * nR, D + (size_t)b * nR, nR, 1, buf);
    }
    free(buf);
    for (int i = 0; i < n; i++) {
        const Cell *c = &cells[i];
        const double x = ratio * c->pw;
        keep[i] = !((double)A[(size_t)c->d * nCells + (size_t)c->b * nR + c->r] > x ||
                    (double)E[(size_t)c->d * nR + c->r] > x || (double)D[(size_t)c->b * nR + c->r] > x);
    }
}

/* Keeps the CFAR detections that are not sidelobes of a much stronger peak
 * (cfar() has kept only the local peaks: 3x3x3 beam x range x Doppler) and
 * converts them to reports. */
static void prune_and_report(const Processor *p, int pulse, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams, nR = s->nR;
    const PulseConfig *pc = &cfg->pulse[pulse];
    Cell *cells = NULL;
    int n = 0, cap = 0;

    /* 1) the local peaks above the threshold, scanned in memory order and (with
     *    several beams) sorted into the order of MATLAB's find() on [beam x
     *    range x Doppler]; with one beam the memory order is that order */
    for (int d = 0; d < nfft; d++) {
        for (int b = 0; b < nB; b++) {
            const unsigned char *detRow = s->det + (size_t)d * s->nCells + (size_t)b * nR, *hit;
            for (int r = 0; r < nR && (hit = memchr(detRow + r, 1, (size_t)(nR - r))) != NULL; r++) {
                r = (int)(hit - detRow);
                if (n == cap) {
                    cap = cap ? 2 * cap : 64;
                    cells = xrealloc(cells, sizeof(Cell) * cap);
                }
                cells[n++] = (Cell){ b, r, d, P_at(s, d, b, r) };
            }
        }
    }
    if (n > 1 && nB > 1) qsort(cells, n, sizeof(Cell), cell_order);

    /* 2) sidelobe blanking (range, elevation and Doppler sidelobes). Only the
     *    peaks more than the margin stronger can blank a peak: the "strong"
     *    peaks, more than the margin above the weakest one, sorted by power,
     *    and of them the first lo[i] (binary search) for peak i. With a low
     *    threshold the many weak peaks are not sorted. When these
     *    comparisons would cost more than the passes over the map (a low CFAR
     *    threshold: many weak peaks, many strong ones), the map maxima are
     *    used instead (blank_dense, same result), so the time is bounded by
     *    the map size and hardly depends on the number of detections. */
    const double ratio = pow(10.0, cfg->sidelobeMargin_dB / 10.0);
    const int pulseCells = pc->nTaps;
    unsigned char *keep = xcalloc(n > 0 ? n : 1, 1);
    Ranked *rank = xcalloc(n > 0 ? n : 1, sizeof(Ranked));
    int *loK = xcalloc(n > 0 ? n : 1, sizeof(int));
    double minPw = INFINITY;
    for (int i = 0; i < n; i++) minPw = cells[i].pw < minPw ? cells[i].pw : minPw;
    int nStrong = 0;                                        /* only these can blank a peak */
    for (int i = 0; i < n; i++)
        if (cells[i].pw > ratio * minPw) rank[nStrong++] = (Ranked){ cells[i].pw, i };
    if (nStrong > 1) qsort(rank, nStrong, sizeof(Ranked), ranked_desc);
    double work = 0.0;
    for (int i = 0; i < n; i++) {
        const double x = ratio * cells[i].pw;
        int lo = 0, hi = nStrong;                            /* rank[0..lo) > x, rank[hi..) <= x */
        while (lo < hi) {
            const int mid = (lo + hi) / 2;
            if (rank[mid].pw > x) lo = mid + 1; else hi = mid;
        }
        loK[i] = lo;
        work += lo;
    }
    if (work > 20.0 * nfft * s->nCells) {                  /* a dense pass: ~20 operations per map cell */
        blank_dense(p, s, cells, n, pulseCells, ratio, keep);
    } else {
        for (int i = 0; i < n; i++) {
            keep[i] = 1;
            for (int k = 0; k < loK[i] && keep[i]; k++) {
                const int j = rank[k].i;
                int dDop = abs(cells[i].d - cells[j].d);
                if (nfft - dDop < dDop) dDop = nfft - dDop;
                const int sameD = dDop <= 1, sameB = abs(cells[i].b - cells[j].b) <= 1;
                const int dR = abs(cells[i].r - cells[j].r), sameR = dR <= 1;
                if ((sameB && sameD && dR <= pulseCells) || (sameR && sameD) || (sameR && sameB)) keep[i] = 0;
            }
        }
    }
    free(loK);

    free(rank);

    /* 3) reports: range, velocity, elevation, CPI centre azimuth/time, SNR */
    s->nReports = 0;
    double sinSum = 0, cosSum = 0, tSum = 0;
    for (int k = 0; k < s->nSlots; k++) {
        sinSum += sin(s->az[k] * PI / 180);
        cosSum += cos(s->az[k] * PI / 180);
        tSum += s->time[k];
    }
    double azCenter = atan2(sinSum / s->nSlots, cosSum / s->nSlots) * 180 / PI;
    azCenter = fmod(azCenter + 360.0, 360.0);
    /* noise floors of the beams with detections, in parallel (not needed with
     * the CFAR noise of every cell) */
    double *floors = xcalloc(nB, sizeof(double));
    int *beams = xcalloc(nB, sizeof(int)), nNeed = 0;
    for (int i = 0; i < n && !s->noise; i++)
        if (keep[i] && floors[cells[i].b] == 0) { floors[cells[i].b] = -1; beams[nNeed++] = cells[i].b; }
    #pragma omp parallel for schedule(dynamic) num_threads(team_size(p->nThreads, nNeed)) if(nNeed > 1)
    for (int j = 0; j < nNeed; j++)
        floors[beams[j]] = beam_noise_floor(cfg, s, beams[j], &p->scratch[thread_num()]);
    free(beams);
    for (int i = 0; i < n; i++) {
        if (!keep[i]) continue;
        const Cell *c = &cells[i];
        Detection d;
        memset(&d, 0, sizeof d);
        d.time_s = tSum / s->nSlots;
        d.az_deg = azCenter;
        d.range_m = pc->range_m[c->r];
        d.vr_mps = -cfg->dopFreq_Hz[c->d] * pc->lambda_m / 2;
        d.el_deg = cfg->dbfEnabled ? cfg->tilt_deg + asin(pc->beamU[c->b]) * 180 / PI : NAN;   /* beam centre (on its carrier) */
        const double noise = s->noise ? s->noise[(size_t)c->d * s->nCells + (size_t)c->b * nR + c->r] : floors[c->b];
        double ratioSnr = c->pw / noise - 1;
        d.snr_dB = 10 * log10(ratioSnr > DBL_EPSILON ? ratioSnr : DBL_EPSILON);
        d.power = c->pw;
        d.beam = c->b;  d.rIdx = c->r;  d.dIdx = c->d;
        d.pulse = pulse;
        d.hits = 0;
        d.pulseNumber = s->newestPulse;
        add_report(s, &d);
    }
    free(floors);
    free(keep);
    free(cells);
}

/* Vertex offset of the parabola through the powers [dB] p0, p1 (peak), p2;
 * 0 when it has no maximum (or a power is zero: -Inf dB). */
static int parabola_vertex(double p0, double p1, double p2, double *delta)
{
    const double l0 = 10 * log10(p0), l1 = 10 * log10(p1), l2 = 10 * log10(p2);
    const double den = l0 - 2 * l1 + l2;
    *delta = 0.5 * (l0 - l2) / den;
    return den < 0 && isfinite(*delta);
}

/* Peak interpolation of the reports of the latest CPI update (as
 * peak_interpolation.m): the radial velocity from the parabola through the
 * Doppler bins d-1, d, d+1 (circular), the elevation from the parabola
 * through the beams b-1, b, b+1 (fractional beam -> u by linear interpolation
 * of beamU); otherwise the centres of the bin and of the beam stay. */
static void peak_interpolation(const Processor *p, int pulse, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams;
    const int doDop = cfg->dopplerInterp && nfft > 2;
    const int doEl = cfg->elevationInterp && cfg->dbfEnabled && nB > 2;
    if (!doDop && !doEl) return;
    const double df = cfg->dopFreq_Hz[1] - cfg->dopFreq_Hz[0];      /* Doppler bin width [Hz] */
    const double lambda = cfg->pulse[pulse].lambda_m;
    const double *beamU = cfg->pulse[pulse].beamU ? cfg->pulse[pulse].beamU : cfg->beamU;
    for (int i = 0; i < s->nReports; i++) {
        Detection *d = &s->reports[i];
        const int b = d->beam, r = d->rIdx, k = d->dIdx;
        double delta;
        if (doDop && parabola_vertex(P_at(s, (k + nfft - 1) % nfft, b, r), P_at(s, k, b, r),
                                     P_at(s, (k + 1) % nfft, b, r), &delta))
            d->vr_mps = -(cfg->dopFreq_Hz[k] + delta * df) * lambda / 2;
        if (doEl && b > 0 && b < nB - 1 &&
            parabola_vertex(P_at(s, k, b - 1, r), P_at(s, k, b, r), P_at(s, k, b + 1, r), &delta)) {
            const double pos = b + delta;                              /* fractional beam index */
            const int i0 = (int)floor(pos);
            const double f = pos - i0;
            const double u = (i0 + 1 < nB) ? beamU[i0] * (1 - f) + beamU[i0 + 1] * f : beamU[i0];
            d->el_deg = cfg->tilt_deg + asin(u) * 180 / PI;
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Display taps: profile of a beam after DBF / pulse compression, range-     */
/* Doppler map of a beam                                                     */
/* ------------------------------------------------------------------------ */

static void taps_alloc(Processor *p)
{
    if (p->profile) return;
    const RadarConfig *cfg = p->cfg;
    p->dR = 299792458.0 / (2.0 * cfg->fsOut_Hz);
    p->nBins = (int)ceil(cfg->maxRange_m / p->dR);
    if (p->nBins < 1) p->nBins = 1;
    p->profile = xcalloc((size_t)p->nBins, sizeof(float));
    p->profileCnt = xcalloc((size_t)p->nBins, sizeof(int));
    p->rd = xcalloc((size_t)p->nBins * cfg->nfft, sizeof(float));
    p->rdVel = xcalloc((size_t)cfg->nfft, sizeof(double));
    p->rdWork = xcalloc((size_t)p->nBins * cfg->nfft + 1, sizeof(float));
    const int vp = cfg->pulse[PULSE_LONG].enable ? PULSE_LONG : PULSE_SHORT;   /* velocity on its carrier */
    for (int d = 0; d < cfg->nfft; d++) p->rdVel[d] = -cfg->dopFreq_Hz[d] * cfg->pulse[vp].lambda_m / 2;
}

void processor_set_profile(Processor *p, int kind, int beam)
{
    if (beam < 0 || beam >= p->cfg->nBeams) beam = 0;
    p->profileKind = kind == PROFILE_DBF || kind == PROFILE_MF ? kind : PROFILE_NONE;
    p->profileBeam = beam;
    if (p->profileKind) taps_alloc(p);
}

void processor_set_rdmap(Processor *p, int enabled, int beam, double minPeriod_s)
{
    if (beam < 0 || beam >= p->cfg->nBeams) beam = 0;
    if (enabled && (!p->rdEnabled || beam != p->rdBeam)) {   /* a new map: no old cells */
        taps_alloc(p);
        memset(p->rd, 0, sizeof(float) * (size_t)p->nBins * p->cfg->nfft);
        p->nRdDets = 0;
        p->rdLast = -1e300;
    }
    p->rdEnabled = enabled != 0;
    p->rdBeam = beam;
    p->rdPeriod_s = minPeriod_s;
}

/* Range bin of a range (-1 outside the taps). */
static int tap_bin(const Processor *p, double r)
{
    const int b = (int)floor(r / p->dR);
    return b >= 0 && b < p->nBins ? b : -1;
}

/* Adds |x|^2 of the selected beam after DBF (pulse compression input inside the
 * receive interval) or after the pulse compression to the profile. */
static void profile_add(Processor *p, int pulse, const PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const PulseConfig *pc = &cfg->pulse[pulse];
    const double handover = cfg->handover_m;
    const int b = p->profileBeam;
    if (p->profileKind == PROFILE_DBF) {
        const int n = pc->ddcCount;
        const float *xr = s->beamRe + (size_t)b * n, *xi = s->beamIm + (size_t)b * n;
        for (int k = 0; k < pc->rxCount; k++) {
            const int j = pc->rxFirst + k - pc->ddcFirst;        /* segment index of rx sample k */
            const double r = pc->range_m[0] + k * p->dR;
            if (pulse == PULSE_SHORT && !(r < handover)) break;
            const int bin = tap_bin(p, r);
            if (bin < 0 || j < 0 || j >= n) continue;
            p->profile[bin] += xr[j] * xr[j] + xi[j] * xi[j];
            p->profileCnt[bin]++;
        }
    } else {
        const float *yr = s->cpiRe + (size_t)s->mfLast * s->nCells + (size_t)b * s->nR;
        const float *yi = s->cpiIm + (size_t)s->mfLast * s->nCells + (size_t)b * s->nR;
        for (int r = 0; r < s->nR; r++) {
            if (pulse == PULSE_SHORT && !(pc->range_m[r] < handover)) break;
            const int bin = tap_bin(p, pc->range_m[r]);
            if (bin < 0) continue;
            p->profile[bin] += yr[r] * yr[r] + yi[r] * yi[r];
            p->profileCnt[bin]++;
        }
    }
}

/* The cells of the selected beam of a sub-pulse's map into the joined range-Doppler map. */
static void rdmap_add(Processor *p, int pulse, const PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const PulseConfig *pc = &cfg->pulse[pulse];
    const double handover = cfg->handover_m;
    const int b = p->rdBeam, nfft = cfg->nfft;
    for (int r = 0; r < s->nR; r++) {
        if (pulse == PULSE_SHORT && !(pc->range_m[r] < handover)) break;
        const int bin = tap_bin(p, pc->range_m[r]);
        if (bin < 0) continue;
        for (int d = 0; d < nfft; d++) p->rd[(size_t)d * p->nBins + bin] = P_at(s, d, b, r);
    }
    for (int i = 0; i < s->nReports; i++) {                /* the reports of this beam */
        if (s->reports[i].beam != b) continue;
        if (p->nRdDets == p->capRdDets) {
            p->capRdDets = p->capRdDets ? 2 * p->capRdDets : 64;
            p->rdDets = xrealloc(p->rdDets, sizeof(Detection) * p->capRdDets);
        }
        p->rdDets[p->nRdDets++] = s->reports[i];
    }
    p->rdUpdated = 1;
}

static void rdmap_post(Processor *p, const PriData *pri, const ProcessorCallbacks *cb)
{
    const RadarConfig *cfg = p->cfg;
    const size_t n = (size_t)p->nBins * cfg->nfft;
    size_t m = 0;
    for (size_t i = 0; i < n; i++) if (p->rd[i] > 0) p->rdWork[m++] = p->rd[i];
    RdMap map;
    memset(&map, 0, sizeof map);
    map.beam = p->rdBeam;
    map.nBins = p->nBins;
    map.nfft = cfg->nfft;
    map.dR_m = p->dR;
    map.power = p->rd;
    map.vel_mps = p->rdVel;
    map.notchBins = cfg->notchBins;
    map.noise = m ? select_kth(p->rdWork, (int)m, (int)(m / 2)) : 0.0;
    map.time_s = pri->time_s;
    map.az_deg = pri->az_deg;
    map.dets = p->rdDets;
    map.nDets = p->nRdDets;
    cb->rdmap(cb->user, &map);
}

/* ------------------------------------------------------------------------ */
/* One PRI                                                                  */
/* ------------------------------------------------------------------------ */

static double elapsed_ms(double *t0)
{
    double t = now_seconds(), dt = (t - *t0) * 1e3;
    *t0 = t;
    return dt;
}

/* The extractor's plots are collected first, so that the extractor and the
 * tracker are timed separately (and the output callbacks are not timed). */
static void plot_callback(void *user, const Detection *plot)
{
    Processor *p = user;
    if (p->nNewPlots == p->capNewPlots) {
        p->capNewPlots = p->capNewPlots ? 2 * p->capNewPlots : 16;
        p->newPlots = xrealloc(p->newPlots, sizeof(Detection) * p->capNewPlots);
    }
    p->newPlots[p->nNewPlots++] = *plot;
}

/* Track updates and deletions are collected first and reported after the
 * tracker has run (so that the output callbacks are not timed). */
typedef struct {
    Track *list;
    int n, cap;
} TrackReports;

static void track_callback(void *user, const Track *track, int plotIndex)
{
    TrackReports *r = user;
    (void)plotIndex;
    if (r->n == r->cap) {
        r->cap = r->cap ? 2 * r->cap : 16;
        r->list = xrealloc(r->list, sizeof(Track) * r->cap);
    }
    r->list[r->n++] = *track;
}

/* Hands the collected plots to the output; the tracker gets them at the end of the PRI. */
static void plots_out(Processor *p, const ProcessorCallbacks *cb)
{
    for (int k = 0; k < p->nNewPlots; k++) {
        const Detection *plot = &p->newPlots[k];
        if (cb && cb->plot) cb->plot(cb->user, plot);
        if (!p->cfg->trackEnable) continue;
        if (p->nTrackPlots == p->capTrackPlots) {
            p->capTrackPlots = p->capTrackPlots ? 2 * p->capTrackPlots : 16;
            p->trackPlots = xrealloc(p->trackPlots, sizeof(Detection) * p->capTrackPlots);
        }
        p->trackPlots[p->nTrackPlots++] = *plot;
    }
    p->nNewPlots = 0;
}

/* The plots of this PRI to the tracker (GNN batches, deletions) at time t_s, as
 * track_update.m once per PRI; flush = 1 associates the open batch. */
static void run_tracker(Processor *p, double t_s, int flush, const ProcessorCallbacks *cb)
{
    if (!p->cfg->trackEnable) return;
    TrackReports r = { NULL, 0, 0 };
    double t0 = now_seconds();
    tracker_update(&p->tracker, p->trackPlots, p->nTrackPlots, t_s, flush, track_callback, &r);
    p->timing.track_ms += elapsed_ms(&t0);
    p->nTrackPlots = 0;
    for (int i = 0; i < r.n; i++)
        if (cb && cb->track) cb->track(cb->user, &r.list[i]);
    free(r.list);
}

/* The PRI has the layout of the configuration: channels, DDC samples, the
 * receive intervals of both sub-pulses (header, 1-based) and the DDC rate.
 * (None of them depends on the CFAR or on the pulse compression filters.) */
static int pri_matches_config(const RadarConfig *cfg, const PriData *pri)
{
    if (pri->nCh != cfg->nCh || pri->nSamples != cfg->nDdc || !pri->iq) return 0;
    const uint32_t *idx[N_PULSES] = { pri->shortIdx, pri->longIdx };
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        if (!pc->enable) continue;                          /* (its header interval is [0 0]) */
        if (idx[i][0] != (uint32_t)pc->rxFirst + 1 || idx[i][1] != (uint32_t)(pc->rxFirst + pc->rxCount))
            return 0;
    }
    return fabs(pri->ddcFs_Hz - cfg->fsOut_Hz) <= 1e-3 * cfg->fsOut_Hz;
}

void processor_flush(Processor *p, const ProcessorCallbacks *cb)
{
    for (int i = 0; i < N_PULSES; i++) {
        extractor_update(&p->ps[i].ext, NULL, 0, 1, plot_callback, p);
        plots_out(p, cb);
    }
    if (p->havePulse) run_tracker(p, p->lastTime, 1, cb);   /* at the newest PRI, as run_radar_sim.m */
}

void processor_process(Processor *p, const PriData *pri, const ProcessorCallbacks *cb)
{
    const RadarConfig *cfg = p->cfg;
    if (!pri_matches_config(cfg, pri)) {                    /* made for another configuration */
        p->timing.rejected++;
        return;
    }

    /* a gap in the pulse numbers (skipped sector): restart the CPIs, close the plots */
    if (p->havePulse && pri->pulseNumber != p->lastPulse + 1) {
        processor_flush(p, cb);
        for (int i = 0; i < N_PULSES; i++) reset_buffers(&p->ps[i]);
    }
    p->lastPulse = pri->pulseNumber;
    p->lastTime = pri->time_s;
    p->havePulse = 1;
    if (cb && cb->antenna) cb->antenna(cb->user, pri->az_deg, pri->time_s, pri->scanNumber);

    const double tStart = now_seconds();
    const int profileKind = cb && cb->profile ? p->profileKind : PROFILE_NONE;
    const int rdOn = cb && cb->rdmap && p->rdEnabled && pri->time_s - p->rdLast >= p->rdPeriod_s;
    if (profileKind) {
        memset(p->profile, 0, sizeof(float) * (size_t)p->nBins);
        memset(p->profileCnt, 0, sizeof(int) * (size_t)p->nBins);
    }
    if (p->rdEnabled && pri->time_s < p->rdLast) p->rdLast = -1e300;   /* a new run */
    if (rdOn) { p->nRdDets = 0; p->rdUpdated = 0; }
    const int firstOn = cfg->pulse[PULSE_SHORT].enable ? PULSE_SHORT : PULSE_LONG;   /* counts the updates */
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        PulseState *s = &p->ps[i];
        if (!pc->enable) continue;                          /* sub-pulse switched off */
        double t0 = now_seconds();

        /* fast time into the next CPI slot */
        load_segment(cfg, pc, pri, s);
        beamform(cfg, pc, s, p->nThreads);
        p->timing.dbf_ms += elapsed_ms(&t0);
        if (profileKind == PROFILE_DBF) profile_add(p, i, s);   /* (the display's time, not the block's) */
        const int slot = (s->last + 1) % s->nSlots;
        s->mfLast = (s->mfLast + 1) % cfg->nMtiTaps;
        matched_filter(cfg, pc, s, s->mfLast, p->nThreads);
        if (profileKind == PROFILE_MF) profile_add(p, i, s);
        s->last = slot;
        s->staggerIdx[slot] = pri->staggerIndex >= 1 && pri->staggerIndex <= cfg->nStagger ? pri->staggerIndex : 1;
        s->az[slot] = pri->az_deg;
        s->time[slot] = pri->time_s;
        s->newestPulse = pri->pulseNumber;
        if (s->count < s->nSlots) s->count++;
        s->nPushed = (s->nPushed + 1) % cfg->updateEvery;  /* pulses since the last update (mod N) */
        p->timing.mf_ms += elapsed_ms(&t0);

        /* slow time */
        if (s->count >= cfg->nMtiTaps) mti_newest(cfg, pc, s);
        if (s->count == s->nSlots) s->sinceUpdate++;
        p->timing.mti_ms += elapsed_ms(&t0);
        if (s->count < s->nSlots || s->nPushed != 0) continue;
        if (!p->slowTime) {                                 /* overload: this update is skipped */
            if (i == firstOn) p->timing.skipped++;
            continue;
        }
        /* a plot needs detections over the same azimuth extent also when
         * updates were skipped: an update counts for the pulses it covers */
        s->ext.weight = s->sinceUpdate / cfg->updateEvery > 1 ? s->sinceUpdate / cfg->updateEvery : 1;
        s->sinceUpdate = 0;
        doppler(p, pc, s);
        p->timing.doppler_ms += elapsed_ms(&t0);
        cfar(p, s, pc);
        prune_and_report(p, i, s);
        p->timing.cfar_ms += elapsed_ms(&t0);
        peak_interpolation(p, i, s);
        p->timing.interp_ms += elapsed_ms(&t0);
        if (rdOn) rdmap_add(p, i, s);
        if (cb && cb->detections) cb->detections(cb->user, i, s->reports, s->nReports);
        t0 = now_seconds();
        extractor_update(&s->ext, s->reports, s->nReports, 0, plot_callback, p);
        p->timing.extract_ms += elapsed_ms(&t0);
        plots_out(p, cb);
        if (i == firstOn) p->timing.updates++;
    }
    if (profileKind) {                                     /* mean of the samples in every bin */
        for (int b = 0; b < p->nBins; b++) p->profile[b] = p->profileCnt[b] ? p->profile[b] / p->profileCnt[b] : 0.0f;
        cb->profile(cb->user, profileKind, p->profile, p->nBins, p->dR, pri->az_deg, pri->time_s);
    }
    if (rdOn && p->rdUpdated) {
        rdmap_post(p, pri, cb);
        p->rdLast = pri->time_s;
    }
    run_tracker(p, pri->time_s, 0, cb);                    /* plots of this PRI, missed scans, deletion */
    p->timing.total_ms += (now_seconds() - tStart) * 1e3;
    p->timing.pris++;
}

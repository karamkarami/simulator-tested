/*
 * cfar.c - CFAR detectors (as cfar_detect.m, cfar_ca.m, cfar_goso.m and
 * cfar_os.m of the MATLAB simulator) and their false alarm probability
 * (cfar_pfa.m).
 *
 * Detector (cfg->cfarType): the statistic of the training cells of the
 * window of a cell that the threshold factor alpha multiplies:
 *   CA  their mean
 *   GO  the larger of the means of the two halves of the window
 *   SO  the smaller of them
 *       (halves: 2-D and 1-D range window: the cells at a smaller and at a
 *       larger range than the cell under test, the cells on its own range
 *       in neither; 1-D Doppler window: the bins below and above it; at a
 *       range edge an empty half is replaced by the other one)
 *   OS  the k-th smallest of them (k = cfg->cfarOsK for the full window,
 *       round(k * N / nTrain) for the N training cells at the range edges)
 * Window (cfg->cfarMethod): guardRange/trainRange and guardDoppler/
 * trainDoppler cells on each side (the 1-D windows have no window in the
 * other dimension). Doppler is circular. At the first and last range cells
 * the window is shifted inwards (its centre rc is the nearest cell whose
 * window lies inside the data) and the guard cells are cut at the edge.
 * Every detection must also be a local peak (the maximum of its 3 x 3 x 3
 * neighbourhood in beam x range x Doppler, the first step of the pruning,
 * prune_peaks.m): the OS statistic is computed only for those cells (the
 * result after the pruning is the same as with every cell tested). The
 * noise estimate of a detection (the SNR of the reports with
 * cfg.report.snrNoise = 'cfar') is the statistic / cfg->cfarNoiseScale (the
 * mean of the statistic for unit noise; 1 for CA).
 */
#include "cfar.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n ? n : 1, size);
    if (!p) abort();
    return p;
}

void cfar_scratch_init(CfarScratch *w, const RadarConfig *cfg, int maxNR)
{
    const int oR = cfg->guardRange + cfg->trainRange, oD = cfg->guardDoppler + cfg->trainDoppler;
    const size_t partialLen = (size_t)maxNR + 2 * oR + 32;
    memset(w, 0, sizeof *w);
    w->maxNR = maxNR;
    w->ringRows = 2 * oD + 2;
    w->ringO = xcalloc((size_t)w->ringRows * maxNR, sizeof(float));
    w->ringG = xcalloc((size_t)w->ringRows * maxNR, sizeof(float));
    w->partial = xcalloc(6 * partialLen, sizeof(float));
    w->accO = xcalloc(maxNR, sizeof(double));
    w->accG = xcalloc(maxNR, sizeof(double));
    w->colMax = xcalloc((size_t)maxNR + 2, sizeof(float));
    w->cO = xcalloc(maxNR, sizeof(double));
    w->cG = xcalloc(maxNR, sizeof(double));
    w->sO = xcalloc((size_t)maxNR + 1, sizeof(double));
    w->sG = xcalloc((size_t)maxNR + 1, sizeof(double));
    w->win = xcalloc((size_t)(2 * oR + 1) * (2 * oD + 1), sizeof(float));
}

void cfar_scratch_free(CfarScratch *w)
{
    free(w->ringO); free(w->ringG); free(w->partial); free(w->accO); free(w->accG); free(w->colMax);
    free(w->cO); free(w->cG); free(w->sO); free(w->sG); free(w->win);
    memset(w, 0, sizeof *w);
}

/* ------------------------------------------------------------------------ */
/* Local peaks                                                              */
/* ------------------------------------------------------------------------ */

/* The cell (b, r, d) is a local peak: no neighbour in beam (b-1 ... b+1),
 * range (r-1 ... r+1) and Doppler (d-1 ... d+1, circular) is stronger (as
 * prune_peaks.m). */
static int is_peak(const CfarMap *m, int b, int r, int d)
{
    const int nR = m->nR, nfft = m->nfft, nB = m->nBeams;
    const float v = m->P[(size_t)d * m->nCells + (size_t)b * nR + r];
    for (int dd = -1; dd <= 1; dd++) {
        const float *row = m->P + (size_t)((d + dd + nfft) % nfft) * m->nCells;
        for (int bb = b > 0 ? b - 1 : 0; bb <= b + 1 && bb < nB; bb++)
            for (int rr = r > 0 ? r - 1 : 0; rr <= r + 1 && rr < nR; rr++)
                if (row[(size_t)bb * nR + rr] > v) return 0;
    }
    return 1;
}

/* Keeps in det[f ... l] (one Doppler row d of beam b) only the local peaks.
 * Many detections in the row (a low threshold): one vectorised pass with the
 * maxima over the neighbouring beams and Doppler bins (cm, nR + 2 floats),
 * whose cost does not depend on the threshold; few: each one is checked. */
static void peaks_row(const CfarMap *m, int b, int d, int f, int l, float *cm, unsigned char *restrict det)
{
    const int nR = m->nR, nfft = m->nfft, nB = m->nBeams;
    int nDet = 0;
    for (int r = f; r <= l; r++) nDet += det[r];
    if (nDet == 0) return;
    if (8 * nDet < l - f + 1) {
        for (int r = f; r <= l; r++)
            if (det[r] && !is_peak(m, b, r, d)) det[r] = 0;
        return;
    }
    const int ra = f > 0 ? f - 1 : 0, rb = l < nR - 1 ? l + 1 : nR - 1;
    float *restrict c = cm + 1;                                  /* c[-1] ... c[nR] */
    c[-1] = c[nR] = -FLT_MAX;                                    /* no neighbour outside the range */
    for (int r = ra; r <= rb; r++) c[r] = -FLT_MAX;
    for (int dd = -1; dd <= 1; dd++) {
        const float *row = m->P + (size_t)((d + dd + nfft) % nfft) * m->nCells;
        for (int bb = b > 0 ? b - 1 : 0; bb <= b + 1 && bb < nB; bb++) {
            const float *restrict x = row + (size_t)bb * nR;
            for (int r = ra; r <= rb; r++) c[r] = c[r] > x[r] ? c[r] : x[r];
        }
    }
    const float *restrict x = m->P + (size_t)d * m->nCells + (size_t)b * nR;
    for (int r = f; r <= l; r++) {
        float mx = c[r - 1] > c[r] ? c[r - 1] : c[r];
        mx = mx > c[r + 1] ? mx : c[r + 1];
        det[r] = (unsigned char)(det[r] & (x[r] >= mx));
    }
}

/* Range cells the windows of the tested cells cf ... cl need: the outer
 * window of the centre rc of every tested cell (shifted inwards at the
 * edges) and the guard window of the cell use the cells e0 ... e1; the sums
 * are needed at a ... z. */
static void needed_span(int nR, int oR, int cf, int cl, int *a, int *z, int *e0, int *e1)
{
    const int rcA = cf < oR ? oR : cf > nR - 1 - oR ? nR - 1 - oR : cf;
    const int rcB = cl < oR ? oR : cl > nR - 1 - oR ? nR - 1 - oR : cl;
    *a = rcA < cf ? rcA : cf;
    *z = rcB > cl ? rcB : cl;
    *e0 = *a - oR > 0 ? *a - oR : 0;
    *e1 = *z + oR < nR - 1 ? *z + oR : nR - 1;
}

static int window_centre(int r, int oR, int nR)
{
    return r < oR ? oR : r > nR - 1 - oR ? nR - 1 - oR : r;
}

/* ------------------------------------------------------------------------ */
/* CA                                                                       */
/* ------------------------------------------------------------------------ */

/* Range box sums of the range segment x[0 ... nR-1] (zero outside it):
 * yO[r] = sum of x over r-hO ... r+hO, yG[r] over r-hG ... r+hG, for
 * r = o0 ... o0+no-1 (yO, yG are indexed like x). Partial sums of 1, 2, 4,
 * 8, ... neighbouring values are built with vectorised passes, and every
 * window is the sum of the partial sums given by the binary digits of its
 * length, added in one pass (no running sum, so no rounding drift; a sum
 * depends only on the values of its window, so any segment that holds the
 * window gives the same float result). buf: 6 x (nR + 2*hO + 32) floats. */
static void range_box_sums(const float *x, float *yO, float *yG, int nR, int o0, int no, int hO, int hG, float *buf)
{
    const int L = nR + 2 * hO, ld = L + 32;
    const int wO = 2 * hO + 1, wG = 2 * hG + 1;
    float *part[6];
    int top = 0;                                        /* part[k]: sums of 2^k values */
    while ((2 << top) <= wO && top < 5) top++;
    for (int k = 0; k <= top; k++) part[k] = buf + (size_t)k * ld;
    memset(part[0], 0, sizeof(float) * hO);
    memcpy(part[0] + hO, x, sizeof(float) * nR);
    memset(part[0] + hO + nR, 0, sizeof(float) * (hO + 32));
    for (int k = 1; k <= top; k++) {
        const int sh = 1 << (k - 1), n = L - (1 << k) + 1;
        const float *restrict a = part[k - 1];
        float *restrict b = part[k];
        for (int i = 0; i < n; i++) b[i] = a[i] + a[i + sh];
    }
    const int width[2] = { wO, wG }, start[2] = { 0, hO - hG };
    float *out[2] = { yO, yG };
    for (int w = 0; w < 2; w++) {
        const float *t[6];                              /* the terms, largest first */
        int nt = 0, off = start[w];
        for (int k = top; k >= 0; k--) {
            if (!(width[w] & (1 << k))) continue;
            t[nt++] = part[k] + off;
            off += 1 << k;
        }
        float *restrict y = out[w];
        const float *restrict t0 = t[0], *restrict t1 = t[nt > 1 ? 1 : 0], *restrict t2 = t[nt > 2 ? 2 : 0];
        const float *restrict t3 = t[nt > 3 ? 3 : 0], *restrict t4 = t[nt > 4 ? 4 : 0], *restrict t5 = t[nt > 5 ? 5 : 0];
        const int r1 = o0 + no;
        switch (nt) {                                   /* one pass, added in the same order */
        case 1:  for (int r = o0; r < r1; r++) y[r] = t0[r]; break;
        case 2:  for (int r = o0; r < r1; r++) y[r] = t0[r] + t1[r]; break;
        case 3:  for (int r = o0; r < r1; r++) y[r] = t0[r] + t1[r] + t2[r]; break;
        case 4:  for (int r = o0; r < r1; r++) y[r] = t0[r] + t1[r] + t2[r] + t3[r]; break;
        case 5:  for (int r = o0; r < r1; r++) y[r] = t0[r] + t1[r] + t2[r] + t3[r] + t4[r]; break;
        default: for (int r = o0; r < r1; r++) y[r] = t0[r] + t1[r] + t2[r] + t3[r] + t4[r] + t5[r]; break;
        }
    }
}

/* Threshold test of the cells f ... l of one Doppler row: det[r] = P[r] >
 * alpha/nTrain(r) * (outer sum - guard sum), with the outer sum of the centre
 * rc of the (shifted) range window. Three loops (left edge, interior, right
 * edge), so that the interior one vectorises. */
static void ca_test_row(const float *restrict P, const double *restrict aO, const double *restrict aG,
                        const double *restrict ts, unsigned char *restrict det, int f, int l, int oR, int nR)
{
    const int e0 = l < oR - 1 ? l : oR - 1, e1 = l < nR - 1 - oR ? l : nR - 1 - oR;
    int r = f;
    for (; r <= e0; r++) det[r] = P[r] > ts[r] * (aO[oR] - aG[r]);                   /* rc = oR       */
    for (; r <= e1; r++) det[r] = P[r] > ts[r] * (aO[r] - aG[r]);                    /* rc = r        */
    for (; r <= l; r++) det[r] = P[r] > ts[r] * (aO[nR - 1 - oR] - aG[r]);           /* rc = nR-1-oR  */
}

/* CA: one beam (or chunk) at a time, the range box sums of the Doppler rows
 * d-oD ... d+1+oD in a ring (a few rows, so they stay in the cache), the
 * running sums along Doppler in double (no drift); the noise is divided by
 * the actual number of training cells (m->cntInv, m->thrScale). */
static void ca_chunk(const RadarConfig *cfg, const CfarMap *m, int b, int cf, int cl, CfarScratch *w)
{
    const int nfft = cfg->nfft, nR = m->nR, R = w->ringRows;
    const int gR = cfg->guardRange, oR = gR + cfg->trainRange;
    const int gD = cfg->guardDoppler, oD = gD + cfg->trainDoppler;
    const int zeroBin = nfft / 2;
    const size_t nCells = (size_t)m->nCells;
    int a, z, e0, e1;
    needed_span(nR, oR, cf, cl, &a, &z, &e0, &e1);
    const float *Pb = m->P + (size_t)b * nR;
    /* range sums of the (unwrapped) Doppler row u in ring slot u mod R */
#define RING_ROW(u) do { const int u_ = (u), slot_ = ((u_ % R) + R) % R;                                     \
                         range_box_sums(Pb + (size_t)(((u_ % nfft) + nfft) % nfft) * nCells + e0,              \
                                        w->ringO + (size_t)slot_ * nR + e0, w->ringG + (size_t)slot_ * nR + e0, \
                                        e1 - e0 + 1, a - e0, z - a + 1, oR, gR, w->partial);                   \
                     } while (0)
#define RING(arr, u) (w->arr + (size_t)((((u) % R) + R) % R) * nR)
    for (int u = -oD; u <= oD; u++) RING_ROW(u);
    double *restrict aO = w->accO, *restrict aG = w->accG;
    for (int r = a; r <= z; r++) aO[r] = aG[r] = 0.0;
    for (int k = -oD; k <= oD; k++) {                          /* window of bin 0 */
        const float *restrict o = RING(ringO, k);
        for (int r = a; r <= z; r++) aO[r] += o[r];
    }
    for (int k = -gD; k <= gD; k++) {
        const float *restrict g = RING(ringG, k);
        for (int r = a; r <= z; r++) aG[r] += g[r];
    }
    for (int d = 0; d < nfft; d++) {
        {
            unsigned char *det = m->det + (size_t)d * nCells + (size_t)b * nR;
            ca_test_row(Pb + (size_t)d * nCells, aO, aG, m->thrScale, det, cf, cl, oR, nR);
            peaks_row(m, b, d, cf, cl, w->colMax, det);
            if (m->noise) {                                     /* SNR of the reports: CFAR noise */
                float *nz = m->noise + (size_t)d * nCells + (size_t)b * nR;
                const unsigned char *hit;
                for (int r = cf; r <= cl && (hit = memchr(det + r, 1, (size_t)(cl - r + 1))) != NULL; r++) {
                    r = (int)(hit - det);
                    nz[r] = (float)(m->cntInv[r] * (aO[window_centre(r, oR, nR)] - aG[r]));
                }
            }
        }
        if (d == nfft - 1) break;
        /* slide the Doppler windows to bin d + 1 */
        RING_ROW(d + 1 + oD);
        const float *restrict addO = RING(ringO, d + 1 + oD), *restrict subO = RING(ringO, d - oD);
        const float *restrict addG = RING(ringG, d + 1 + gD), *restrict subG = RING(ringG, d - gD);
        for (int r = a; r <= z; r++) {
            aO[r] += (double)addO[r] - subO[r];
            aG[r] += (double)addG[r] - subG[r];
        }
    }
#undef RING
#undef RING_ROW
}

/* ------------------------------------------------------------------------ */
/* GO, SO                                                                   */
/* ------------------------------------------------------------------------ */

/* GO / SO with the 2-D or the 1-D range window: Doppler sums of every range
 * cell over +-oD and +-gD (running, double), their prefix sums along range
 * (double), from which the two halves of every cell follow (as cfar_goso.m);
 * with the 1-D Doppler window: running sums of the bins below and above the
 * cell. det = P > alpha * (larger / smaller half mean), then the peak test. */
static void goso_chunk(const RadarConfig *cfg, const CfarMap *m, double alpha, int b, int cf, int cl, CfarScratch *w)
{
    const int nfft = cfg->nfft, nR = m->nR, go = cfg->cfarType == CFAR_GO;
    const int gR = cfg->guardRange, oR = gR + cfg->trainRange;
    const int gD = cfg->guardDoppler, oD = gD + cfg->trainDoppler;
    const int zeroBin = nfft / 2;
    const size_t nCells = (size_t)m->nCells;
    const double scale = cfg->cfarNoiseScale;
    const float *Pb = m->P + (size_t)b * nR;
#define ROW(u) (Pb + (size_t)((((u) % nfft) + nfft) % nfft) * nCells)
    double *cO = w->cO, *cG = w->cG, *sO = w->sO, *sG = w->sG;
    if (cfg->cfarMethod == CFAR_DOPPLER) {
        const double nh = oD - gD;
        for (int r = cf; r <= cl; r++) {                       /* halves of bin 0: cO below, cG above */
            double sa = 0.0, sb = 0.0;
            for (int k = -oD; k <= -gD - 1; k++) sa += ROW(k)[r];
            for (int k = gD + 1; k <= oD; k++) sb += ROW(k)[r];
            cO[r] = sa;  cG[r] = sb;
        }
        for (int d = 0; d < nfft; d++) {
            if (abs(d - zeroBin) > cfg->notchBins) {
                unsigned char *det = m->det + (size_t)d * nCells + (size_t)b * nR;
                const float *x = ROW(d);
                for (int r = cf; r <= cl; r++) {
                    const double ma = cO[r] / nh, mb = cG[r] / nh, s = go ? fmax(ma, mb) : fmin(ma, mb);
                    det[r] = (double)x[r] > alpha * s;
                    sO[r] = s;
                }
                peaks_row(m, b, d, cf, cl, w->colMax, det);
                if (m->noise) {
                    float *nz = m->noise + (size_t)d * nCells + (size_t)b * nR;
                    for (int r = cf; r <= cl; r++) if (det[r]) nz[r] = (float)(sO[r] / scale);
                }
            }
            const float *addA = ROW(d - gD), *subA = ROW(d - oD), *addB = ROW(d + 1 + oD), *subB = ROW(d + 1 + gD);
            for (int r = cf; r <= cl; r++) {
                cO[r] += (double)addA[r] - subA[r];
                cG[r] += (double)addB[r] - subB[r];
            }
        }
        return;
    }
    int a, z, e0, e1;
    needed_span(nR, oR, cf, cl, &a, &z, &e0, &e1);
    (void)a;  (void)z;
    for (int i = e0; i <= e1; i++) {                           /* Doppler sums of bin 0 */
        double so = 0.0, sg = 0.0;
        for (int k = -oD; k <= oD; k++) so += ROW(k)[i];
        for (int k = -gD; k <= gD; k++) sg += ROW(k)[i];
        cO[i] = so;  cG[i] = sg;
    }
    for (int d = 0; d < nfft; d++) {
        if (abs(d - zeroBin) > cfg->notchBins) {
            unsigned char *det = m->det + (size_t)d * nCells + (size_t)b * nR;
            const float *x = ROW(d);
            sO[0] = sG[0] = 0.0;                               /* sO[i - e0] = sum of the cells e0 ... i-1 */
            for (int i = e0; i <= e1; i++) {
                sO[i - e0 + 1] = sO[i - e0] + cO[i];
                sG[i - e0 + 1] = sG[i - e0] + cG[i];
            }
#define S(arr, i) arr[(i) - e0]
            /* interior cells (window not shifted, halves complete): one loop without branches */
            const int i0 = cf > oR ? cf : oR, i1 = cl < nR - 1 - oR ? cl : nR - 1 - oR;
            if (i0 <= i1) {
                const int nh = oR * (2 * oD + 1) - gR * (2 * gD + 1);
                const double *restrict pO = sO - e0, *restrict pG = sG - e0;
                double *restrict st = w->accO;
                unsigned char *restrict dt = det;
                for (int r = i0; r <= i1; r++) {
                    const double sa = (pO[r] - pO[r - oR]) - (pG[r] - pG[r - gR]);
                    const double sb = (pO[r + oR + 1] - pO[r + 1]) - (pG[r + gR + 1] - pG[r + 1]);
                    const double ma = sa / nh, mb = sb / nh;
                    const double s = go ? (ma > mb ? ma : mb) : (ma < mb ? ma : mb);
                    dt[r] = (double)x[r] > alpha * s;
                    st[r] = s;
                }
            }
            for (int r = cf; r <= cl; r++) {
                if (r == i0 && i0 <= i1) { r = i1; continue; }   /* done above */
                const int rc = window_centre(r, oR, nR), a1 = rc - oR, b2 = rc + oR;
                const int g1 = r - gR > 0 ? r - gR : 0, g2 = r + gR < nR - 1 ? r + gR : nR - 1;
                const double sa = (S(sO, r) - S(sO, a1)) - (S(sG, r) - S(sG, g1));
                const double sb = (S(sO, b2 + 1) - S(sO, r + 1)) - (S(sG, g2 + 1) - S(sG, r + 1));
                const int na = (r - a1) * (2 * oD + 1) - (r - g1) * (2 * gD + 1);
                const int nb = (b2 - r) * (2 * oD + 1) - (g2 - r) * (2 * gD + 1);
                double ma = sa / (na > 1 ? na : 1), mb = sb / (nb > 1 ? nb : 1);
                if (na == 0) ma = mb;
                if (nb == 0) mb = ma;
                const double s = go ? fmax(ma, mb) : fmin(ma, mb);
                det[r] = (double)x[r] > alpha * s;
                w->accO[r] = s;
            }
#undef S
            peaks_row(m, b, d, cf, cl, w->colMax, det);
            if (m->noise) {
                float *nz = m->noise + (size_t)d * nCells + (size_t)b * nR;
                for (int r = cf; r <= cl; r++) if (det[r]) nz[r] = (float)(w->accO[r] / scale);
            }
        }
        if (d == nfft - 1) break;
        const float *addO = ROW(d + 1 + oD), *subO = ROW(d - oD), *addG = ROW(d + 1 + gD), *subG = ROW(d - gD);
        for (int i = e0; i <= e1; i++) {                       /* slide to bin d + 1 */
            cO[i] += (double)addO[i] - subO[i];
            cG[i] += (double)addG[i] - subG[i];
        }
    }
#undef ROW
}

/* ------------------------------------------------------------------------ */
/* OS                                                                       */
/* ------------------------------------------------------------------------ */

/* k-th smallest (0-based) of v[0 ... n-1] (reorders v). */
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

/* OS: the local peaks of the tested cells are the candidates; a candidate
 * is detected when at least k of its N training cells t have alpha*t < P
 * (in double; i.e. P > alpha * X(k)). The count stops after the Doppler row
 * in which more than N - k cells are above (the cell cannot be detected). */
static void os_chunk(const RadarConfig *cfg, const CfarMap *m, double alpha, int b, int cf, int cl, CfarScratch *w)
{
    const int nfft = cfg->nfft, nR = m->nR;
    const int gR = cfg->guardRange, oR = gR + cfg->trainRange;
    const int gD = cfg->guardDoppler, oD = gD + cfg->trainDoppler;
    const int zeroBin = nfft / 2, nWin = (2 * oR + 1) * (2 * oD + 1);
    const size_t nCells = (size_t)m->nCells;
    const float *Pb = m->P + (size_t)b * nR;
    for (int d = 0; d < nfft; d++) {
        /* The range-limited notch is applied by cfar_chunk after detection. */
        unsigned char *det = m->det + (size_t)d * nCells + (size_t)b * nR;
        for (int r = cf; r <= cl; r++) det[r] = 1;
        peaks_row(m, b, d, cf, cl, w->colMax, det);              /* the candidates */
        const float *x = Pb + (size_t)d * nCells;
        for (int r = cf; r <= cl; r++) {
            if (!det[r]) continue;
            const int rc = window_centre(r, oR, nR);
            const double th = x[r];
            const int g1 = r - gR > 0 ? r - gR : 0, g2 = r + gR < nR - 1 ? r + gR : nR - 1;
            const int nr = nWin - (g2 - g1 + 1) * (2 * gD + 1);
            long kr = lround((double)cfg->cfarOsK * nr / cfg->nTrain);
            if (kr < 1) kr = 1;
            if (kr > nr) kr = nr;
            /* detected when at least kr training cells t have alpha*t < P, i.e. at most nr - kr
             * are not: stop as soon as more are not (most candidates, a fraction of the window) */
            const long maxAbove = nr - kr;
            long above = 0;
            for (int k = -oD; k <= oD && above <= maxAbove; k++) {
                const float *row = Pb + (size_t)(((d + k) % nfft + nfft) % nfft) * nCells;
                const int inG = abs(k) <= gD;
                int n = 0;
                if (inG) {                                   /* the row without its guard cells */
                    for (int rr = rc - oR; rr < g1; rr++) n += !(alpha * (double)row[rr] < th);
                    for (int rr = g2 + 1; rr <= rc + oR; rr++) n += !(alpha * (double)row[rr] < th);
                } else {
                    for (int rr = rc - oR; rr <= rc + oR; rr++) n += !(alpha * (double)row[rr] < th);
                }
                above += n;
            }
            det[r] = above <= maxAbove;
            if (det[r] && m->noise) {                            /* X(k): the noise estimate */
                int n = 0;
                for (int k = -oD; k <= oD; k++) {
                    const float *row = Pb + (size_t)(((d + k) % nfft + nfft) % nfft) * nCells;
                    const int inG = abs(k) <= gD;
                    for (int rr = rc - oR; rr <= rc + oR; rr++)
                        if (!(inG && abs(rr - r) <= gR)) w->win[n++] = row[rr];
                }
                m->noise[(size_t)d * nCells + (size_t)b * nR + r] =
                    (float)(select_kth(w->win, n, (int)kr - 1) / cfg->cfarNoiseScale);
            }
        }
    }
}

void cfar_chunk(const RadarConfig *cfg, const CfarMap *m, double alpha, int b, int cf, int cl,
                int notchLastCell, CfarScratch *w)
{
    switch (cfg->cfarType) {
    case CFAR_GO:
    case CFAR_SO: goso_chunk(cfg, m, alpha, b, cf, cl, w); break;
    case CFAR_OS: os_chunk(cfg, m, alpha, b, cf, cl, w); break;
    default:      ca_chunk(cfg, m, b, cf, cl, w); break;
    }
    /* Suppress only near-range cells in the MTI notch. Beyond notchLastCell
     * the CFAR result above is retained for every Doppler bin. */
    if (notchLastCell >= cf) {
        const int d0 = cfg->nfft / 2 - cfg->notchBins > 0 ? cfg->nfft / 2 - cfg->notchBins : 0;
        const int d1 = cfg->nfft / 2 + cfg->notchBins < cfg->nfft - 1 ? cfg->nfft / 2 + cfg->notchBins : cfg->nfft - 1;
        const int r1 = notchLastCell < cl ? notchLastCell : cl;
        const size_t nCells = (size_t)m->nCells;
        for (int d = d0; d <= d1; d++) {
            unsigned char *det = m->det + (size_t)d * nCells + (size_t)b * m->nR;
            memset(det + cf, 0, (size_t)(r1 - cf + 1));
            if (m->noise) {
                float *noise = m->noise + (size_t)d * nCells + (size_t)b * m->nR;
                memset(noise + cf, 0, sizeof(float) * (size_t)(r1 - cf + 1));
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/* False alarm probability                                                  */
/* ------------------------------------------------------------------------ */

/* Continued fraction of the incomplete beta function (modified Lentz). */
static double beta_cf(double a, double b, double x)
{
    const double tiny = 1e-300;
    double c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0);
    if (fabs(d) < tiny) d = tiny;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m <= 1000; m++) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
        d = 1.0 + aa * d;  if (fabs(d) < tiny) d = tiny;
        c = 1.0 + aa / c;  if (fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
        d = 1.0 + aa * d;  if (fabs(d) < tiny) d = tiny;
        c = 1.0 + aa / c;  if (fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < 1e-15) break;
    }
    return h;
}

/* Regularized incomplete beta function I_x(a, b) (MATLAB's betainc). */
static double beta_inc(double x, double a, double b)
{
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;
    const double lfront = lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log1p(-x);
    if (x < (a + 1.0) / (a + b + 2.0))
        return exp(lfront) * beta_cf(a, b, x) / a;
    return 1.0 - exp(lfront) * beta_cf(b, a, 1.0 - x) / b;
}

/* Regularized lower incomplete gamma function P(a, x) (MATLAB's gammainc(x, a)). */
static double gamma_p(double a, double x)
{
    if (!(x > 0.0)) return 0.0;
    const double lfront = -x + a * log(x) - lgamma(a);
    if (x < a + 1.0) {                                        /* series */
        double ap = a, del = 1.0 / a, sum = del;
        for (int n = 0; n < 100000; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * 1e-16) break;
        }
        return sum * exp(lfront);
    }
    const double tiny = 1e-300;                               /* continued fraction of Q(a, x) */
    double bb = x + 1.0 - a, c = 1.0 / tiny, d = 1.0 / bb, h = d;
    for (int i = 1; i < 100000; i++) {
        const double an = -i * (i - a);
        bb += 2.0;
        d = an * d + bb;  if (fabs(d) < tiny) d = tiny;
        c = bb + an / c;  if (fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < 1e-16) break;
    }
    return 1.0 - exp(lfront) * h;
}

/* Distribution function of the statistic for unit noise (cfar_statistic.m). */
static double statistic_cdf(const RadarConfig *cfg, double z, double L)
{
    if (z < 0.0) z = 0.0;
    switch (cfg->cfarType) {
    case CFAR_GO: { const double n = cfg->cfarNHalf * L, F = gamma_p(n, n * z); return F * F; }
    case CFAR_SO: { const double n = cfg->cfarNHalf * L, F = gamma_p(n, n * z); return F * (2.0 - F); }
    case CFAR_OS: {
        const int N = cfg->nTrain, k = cfg->cfarOsK;
        return beta_inc(gamma_p(L, L * z), k, N - k + 1);
    }
    default: { const double n = cfg->nTrain * L; return gamma_p(n, n * z); }
    }
}

double cfar_pfa(const RadarConfig *cfg, double alpha)
{
    const double L = cfg->nLooks;
    if (!(L > 0.0) || cfg->nTrain < 1 || !(alpha >= 0.0)) return NAN;
    if (cfg->cfarType == CFAR_CA) return beta_inc(1.0 / (1.0 + alpha / cfg->nTrain), cfg->nTrain * L, L);
    /* P(Z < S/alpha) over the density of S: 8-point Gauss-Legendre on 2000
     * panels of [0, sHi] (cfar_quadrature.m) */
    static const double x[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                  0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363 };
    static const double g[8] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    const int nP = 2000;
    const double sHi = (L + 10.0 * sqrt(L) + 75.0) / L, h = sHi / nP, lc = L * log(L) - lgamma(L);
    double sum = 0.0;
    for (int p = 0; p < nP; p++) {
        const double mid = (p + 0.5) * h;
        for (int i = 0; i < 8; i++) {
            const double s = mid + x[i] * h / 2.0;
            sum += g[i] * h / 2.0 * statistic_cdf(cfg, s / alpha, L) * exp(lc + (L - 1.0) * log(s) - L * s);
        }
    }
    return sum;
}

const char *cfar_type_name(int type)
{
    return type == CFAR_GO ? "GO" : type == CFAR_SO ? "SO" : type == CFAR_OS ? "OS" : "CA";
}

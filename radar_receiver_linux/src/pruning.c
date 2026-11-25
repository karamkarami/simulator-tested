/* Collect CFAR peaks and prune sidelobes before reporting. */
#include "detection_reports.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* MATLAB's find() order: Doppler, then range, then beam. */
static int cell_order(const void *a, const void *b)
{
    const DetectionCell *x = a, *y = b;
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
static void blank_dense(const Processor *p, PulseState *s, const DetectionCell *cells, int n, int pulseCells,
                        double ratio, unsigned char *keep)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams, nR = s->nR;
    const size_t nCells = (size_t)s->nCells, map = (size_t)nfft * nCells;
    const int rows = nfft > nB ? nfft : nB;
    if (!s->pkM) {                                  /* (room for the largest FFT of all PRFs) */
        const int nfftMax = cfg->maxNfft > nfft ? cfg->maxNfft : nfft, rowsMax = nfftMax > nB ? nfftMax : nB;
        s->pkM = proc_calloc((size_t)nfftMax * nCells, sizeof(float));
        s->pkA = proc_calloc((size_t)nfftMax * nCells, sizeof(float));
        s->pkE = proc_calloc((size_t)nfftMax * nR, sizeof(float));
        s->pkD = proc_calloc((size_t)nB * nR, sizeof(float));
        s->pkT = proc_calloc((size_t)rowsMax * nR, sizeof(float));
    }
    (void)map;  (void)rows;
    float *M = s->pkM, *A = s->pkA, *E = s->pkE, *D = s->pkD, *T = s->pkT;
    memset(M, 0, sizeof(float) * map);
    for (int i = 0; i < n; i++) M[(size_t)cells[i].d * nCells + (size_t)cells[i].b * nR + cells[i].r] = cells[i].pw;
    const size_t bufLen = 2 * ((size_t)nR + 4 * (size_t)pulseCells + 2) + 8;

    /* range sidelobes: beam +-1, Doppler +-1 (circular), range +-pulseCells */
    #pragma omp parallel num_threads(team_size(p->nThreads, nfft)) if(nfft > 1)
    {
        float *row = proc_calloc(nR, sizeof(float)), *buf = proc_calloc(bufLen, sizeof(float));
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
    float *buf = proc_calloc(bufLen, sizeof(float));
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
        const DetectionCell *c = &cells[i];
        const double x = ratio * c->pw;
        keep[i] = !((double)A[(size_t)c->d * nCells + (size_t)c->b * nR + c->r] > x ||
                    (double)E[(size_t)c->d * nR + c->r] > x || (double)D[(size_t)c->b * nR + c->r] > x);
    }
}

/* Keeps the CFAR detections that are not sidelobes of a much stronger peak
 * (processor_cfar() has kept only the local peaks: 3x3x3 beam x range x Doppler) and
 * converts them to reports. */
void prune_and_report(const Processor *p, int pulse, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams, nR = s->nR;
    const PulseConfig *pc = &cfg->pulse[pulse];
    DetectionCell *cells = NULL;
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
                    cells = proc_realloc(cells, sizeof(DetectionCell) * cap);
                }
                cells[n++] = (DetectionCell){ b, r, d, P_at(s, d, b, r) };
            }
        }
    }
    if (n > 1 && nB > 1) qsort(cells, n, sizeof(DetectionCell), cell_order);

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
    unsigned char *keep = proc_calloc(n > 0 ? n : 1, 1);
    Ranked *rank = proc_calloc(n > 0 ? n : 1, sizeof(Ranked));
    int *loK = proc_calloc(n > 0 ? n : 1, sizeof(int));
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

    build_detection_reports(p, pulse, s, cells, keep, n);
    free(keep);
    free(cells);
}

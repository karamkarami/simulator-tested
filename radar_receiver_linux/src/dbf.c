#include "processor_internal.h"

#include <string.h>

/* beam b = sum over channels of W[b][ch] * x[ch]  (or the single channel).
 * Three beams are updated per pass over a channel (fewer loads). */
void beamform(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int nThreads)
{
    const int n = pc->ddcCount, nB = cfg->nBeams, nCh = cfg->nChUsed > 0 ? cfg->nChUsed : cfg->nCh;
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

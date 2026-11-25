/* Moving target indication across successive pulses. */
#include "processor_internal.h"

#include <string.h>

/* MTI output of the newest L pulses: sum_l c[stagger of the oldest][l] * x[oldest + l]. */
void mti_newest(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s)
{
    const int L = cfg->nMtiTaps, N = s->nSlots;
    const int first = (s->last - (L - 1) + N) % N;               /* oldest pulse of the output */
    const double *c = cfg->mtiCoeffs + (size_t)(s->staggerIdx[first] - 1) * L;
    const int firstMf = (s->mfLast + 1) % L;                     /* its matched filter output  */
    const int nBlocks = (s->nCells + CELL_BLOCK - 1) / CELL_BLOCK;
    s->mtiLast = (s->mtiLast + 1) % s->nMtiOut;
    s->mtiStag[s->mtiLast] = s->staggerIdx[first];
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

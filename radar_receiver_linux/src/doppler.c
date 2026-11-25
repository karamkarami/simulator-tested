/* Slow-time Doppler processing for uniform and staggered PRFs. */
#include "processor_internal.h"

#include <string.h>

/* NUFFT (staggered PRF, cfg->dopplerMethod == DOPPLER_NUFFT): P[j][cell] = |sum_q A[p][j][q] x_q|^2
 * over the newest rowLen MTI outputs x_q (oldest first), with the matrix of the stagger position p
 * of the oldest one (its window, density compensation and noise scaling are in A); the cells
 * beyond cfg.mti.maxRange_m (no clutter filter) are divided by their gain. The cells of a block
 * are the fast index (32 lanes: the products vectorise). */
static void doppler_nufft(const Processor *p, const PulseConfig *pc, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const int K = cfg->nfft, N = cfg->rowLen;
    const int nBlocks = (s->nCells + CELL_BLOCK - 1) / CELL_BLOCK;
    const int oldest = (s->mtiLast - (N - 1) + 2 * s->nMtiOut) % s->nMtiOut;
    int ph = s->mtiStag[oldest] - 1;
    if (ph < 0 || ph >= cfg->nStagger) ph = 0;
    const float *Ar = cfg->nufftARe + (size_t)ph * K * N, *Ai = cfg->nufftAIm + (size_t)ph * K * N;
    const float *gB = cfg->nufftBypassGain + (size_t)ph * K;

    #pragma omp parallel num_threads(team_size(p->nThreads, nBlocks)) if(nBlocks > 1)
    {
        float *xr = p->scratch[thread_num()].xr;
        float *xi = p->scratch[thread_num()].xi;
        #pragma omp for schedule(static)
        for (int blk = 0; blk < nBlocks; blk++) {
            const int c0 = blk * CELL_BLOCK;
            const int nc = (s->nCells - c0 < CELL_BLOCK) ? s->nCells - c0 : CELL_BLOCK;
            for (int q = 0; q < N; q++) {                          /* the inputs, oldest first */
                const int slot = (oldest + q) % s->nMtiOut;
                memcpy(xr + (size_t)q * CELL_BLOCK, mti_row(s, s->mtiRe, blk, slot), sizeof(float) * nc);
                memcpy(xi + (size_t)q * CELL_BLOCK, mti_row(s, s->mtiIm, blk, slot), sizeof(float) * nc);
                if (nc < CELL_BLOCK) {
                    memset(xr + (size_t)q * CELL_BLOCK + nc, 0, sizeof(float) * (CELL_BLOCK - nc));
                    memset(xi + (size_t)q * CELL_BLOCK + nc, 0, sizeof(float) * (CELL_BLOCK - nc));
                }
            }
            for (int j = 0; j < K; j++) {
                float yr[CELL_BLOCK], yi[CELL_BLOCK];
                for (int c = 0; c < CELL_BLOCK; c++) yr[c] = yi[c] = 0.0f;
                const float *ar = Ar + (size_t)j * N, *ai = Ai + (size_t)j * N;
                for (int q = 0; q < N; q++) {
                    const float a = ar[q], b = ai[q];
                    const float *restrict ur = xr + (size_t)q * CELL_BLOCK, *restrict ui = xi + (size_t)q * CELL_BLOCK;
                    IVDEP
                    for (int c = 0; c < CELL_BLOCK; c++) {
                        yr[c] += a * ur[c] - b * ui[c];
                        yi[c] += a * ui[c] + b * ur[c];
                    }
                }
                float *restrict pd = s->P + (size_t)j * s->nCells + c0;
                for (int c = 0; c < nc; c++) pd[c] = yr[c] * yr[c] + yi[c] * yi[c];
                if (pc->mtiCells < s->nR) {                       /* unfiltered cells: their own gain */
                    const float g = 1.0f / gB[j];
                    for (int c = 0; c < nc; c++)
                        if ((c0 + c) % s->nR >= pc->mtiCells) pd[c] *= g;
                }
            }
        }
    }
}

/* Range-Doppler power map P[bin][cell] from the newest nRows*rowLen MTI outputs:
 * rows (fixed PRF: one row of consecutive outputs; staggered: de-interleaved by
 * stagger position), Hamming window, FFT, |X|^2 averaged over the rows, fftshift
 * and division by the noise gain of every bin. */
void doppler(const Processor *p, const PulseConfig *pc, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    if (cfg->dopplerMethod == DOPPLER_NUFFT) {
        doppler_nufft(p, pc, s);
        return;
    }
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

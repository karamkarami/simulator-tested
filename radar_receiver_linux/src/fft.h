/*
 * fft.h - radix-2 FFT of many signals at once (internal to the receiver).
 *
 * The signals are stored as rows: sample m of signal c is at re[m*stride + c].
 * Every butterfly works on whole rows, so the inner loop runs over the
 * signals and vectorises well (SSE/AVX).
 */
#ifndef FFT_H
#define FFT_H

#include <stddef.h>

typedef struct {
    int    n;          /* FFT length (power of 2)            */
    int   *bitrev;     /* bit-reversed index of every sample */
    float *twRe;       /* twiddles exp(-j*2*pi*k/n), k < n/2 */
    float *twIm;
} FftPlan;

int  fft_plan_init(FftPlan *plan, int n);   /* 0 on success */
void fft_plan_free(FftPlan *plan);

/* In-place forward FFT along the rows of re/im for signals 0 ... nSignals-1.
 * The input must already be in bit-reversed order: sample m of a signal is
 * stored in row plan->bitrev[m]. The output is in natural order. */
void fft_rows(const FftPlan *plan, float *re, float *im, size_t stride, int nSignals);

#endif

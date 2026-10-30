/*
 * fft.c - radix-2 / radix-2^2 decimation-in-time FFT of many signals at once.
 */
#include "fft.h"

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

int fft_plan_init(FftPlan *plan, int n)
{
    plan->n = n;
    plan->bitrev = malloc(sizeof(int) * n);
    plan->twRe = malloc(sizeof(float) * (n / 2 + 1));
    plan->twIm = malloc(sizeof(float) * (n / 2 + 1));
    if (!plan->bitrev || !plan->twRe || !plan->twIm || n < 1 || (n & (n - 1))) return -1;
    int bits = 0;
    while ((1 << bits) < n) bits++;
    for (int m = 0; m < n; m++) {
        int r = 0;
        for (int b = 0; b < bits; b++)
            if (m & (1 << b)) r |= 1 << (bits - 1 - b);
        plan->bitrev[m] = r;
    }
    for (int k = 0; k < n / 2; k++) {
        double a = -2.0 * PI * k / n;
        plan->twRe[k] = (float)cos(a);
        plan->twIm[k] = (float)sin(a);
    }
    return 0;
}

void fft_plan_free(FftPlan *plan)
{
    free(plan->bitrev);
    free(plan->twRe);
    free(plan->twIm);
    plan->bitrev = NULL;
    plan->twRe = plan->twIm = NULL;
}

/* Tell the compiler that the signal loops do not alias (vectorise them
 * without run-time overlap checks). */
#if defined(__clang__)
#  define IVDEP _Pragma("clang loop vectorize(assume_safety)")
#elif defined(__GNUC__)
#  define IVDEP _Pragma("GCC ivdep")
#elif defined(_MSC_VER)
#  define IVDEP __pragma(loop(ivdep))
#else
#  define IVDEP
#endif

/* Radix-2^2: the butterflies of two radix-2 stages are done in one pass over
 * the data (same operations and results as two radix-2 passes, half the
 * memory traffic). With a constant nS the compiler unrolls the signal loop. */
#if defined(__GNUC__)
__attribute__((always_inline))
#endif
static inline void fft_core(const FftPlan *plan, float *re, float *im, size_t stride, int nS)
{
    const int n = plan->n;
    int size = 2;
    if (n >= 4) {                    /* stages of size 2 and 4: twiddles 1 and -j, no multiplications */
        for (int s0 = 0; s0 < n; s0 += 4) {
            float *restrict r0 = re + (size_t)s0 * stride, *restrict r1 = r0 + stride;
            float *restrict r2 = r1 + stride, *restrict r3 = r2 + stride;
            float *restrict i0 = im + (size_t)s0 * stride, *restrict i1 = i0 + stride;
            float *restrict i2 = i1 + stride, *restrict i3 = i2 + stride;
            IVDEP
            for (int c = 0; c < nS; c++) {
                const float ar = r0[c] + r1[c], ai = i0[c] + i1[c], br = r0[c] - r1[c], bi = i0[c] - i1[c];
                const float cr = r2[c] + r3[c], ci = i2[c] + i3[c], dr = r2[c] - r3[c], di = i2[c] - i3[c];
                r0[c] = ar + cr;  i0[c] = ai + ci;  r2[c] = ar - cr;  i2[c] = ai - ci;
                r1[c] = br + di;  i1[c] = bi - dr;  r3[c] = br - di;  i3[c] = bi + dr;   /* b -+ j*d */
            }
        }
        size = 8;
    }
    for (; size <= n; size *= 4) {
        if (2 * size > n) {          /* one radix-2 stage left */
            const int half = size / 2, step = n / size;
            for (int s0 = 0; s0 < n; s0 += size)
                for (int k = 0; k < half; k++) {
                    const float wr = plan->twRe[k * step], wi = plan->twIm[k * step];
                    float *restrict ar = re + (size_t)(s0 + k) * stride, *restrict ai = im + (size_t)(s0 + k) * stride;
                    float *restrict br = re + (size_t)(s0 + k + half) * stride;
                    float *restrict bi = im + (size_t)(s0 + k + half) * stride;
                    IVDEP
                    for (int c = 0; c < nS; c++) {
                        const float tr = wr * br[c] - wi * bi[c], ti = wr * bi[c] + wi * br[c];
                        br[c] = ar[c] - tr;  bi[c] = ai[c] - ti;  ar[c] += tr;  ai[c] += ti;
                    }
                }
            break;
        }
        /* stages 'size' (half h) and '2*size' (half size) together */
        const int h = size / 2, M = 2 * size, st1 = n / size, st2 = n / M;
        for (int s0 = 0; s0 < n; s0 += M)
            for (int k = 0; k < h; k++) {
                const float w1r = plan->twRe[k * st1], w1i = plan->twIm[k * st1];
                const float w2r = plan->twRe[k * st2], w2i = plan->twIm[k * st2];
                const float w3r = plan->twRe[(k + h) * st2], w3i = plan->twIm[(k + h) * st2];
                float *restrict r0 = re + (size_t)(s0 + k) * stride, *restrict i0 = im + (size_t)(s0 + k) * stride;
                float *restrict r1 = re + (size_t)(s0 + k + h) * stride, *restrict i1 = im + (size_t)(s0 + k + h) * stride;
                float *restrict r2 = re + (size_t)(s0 + k + size) * stride;
                float *restrict i2 = im + (size_t)(s0 + k + size) * stride;
                float *restrict r3 = re + (size_t)(s0 + k + size + h) * stride;
                float *restrict i3 = im + (size_t)(s0 + k + size + h) * stride;
                IVDEP
                for (int c = 0; c < nS; c++) {
                    float tr = w1r * r1[c] - w1i * i1[c], ti = w1r * i1[c] + w1i * r1[c];   /* stage 1 */
                    const float a0r = r0[c] + tr, a0i = i0[c] + ti, a1r = r0[c] - tr, a1i = i0[c] - ti;
                    tr = w1r * r3[c] - w1i * i3[c];  ti = w1r * i3[c] + w1i * r3[c];
                    const float a2r = r2[c] + tr, a2i = i2[c] + ti, a3r = r2[c] - tr, a3i = i2[c] - ti;
                    tr = w2r * a2r - w2i * a2i;  ti = w2r * a2i + w2i * a2r;                   /* stage 2 */
                    r0[c] = a0r + tr;  i0[c] = a0i + ti;  r2[c] = a0r - tr;  i2[c] = a0i - ti;
                    tr = w3r * a3r - w3i * a3i;  ti = w3r * a3i + w3i * a3r;
                    r1[c] = a1r + tr;  i1[c] = a1i + ti;  r3[c] = a1r - tr;  i3[c] = a1i - ti;
                }
            }
    }
}

void fft_rows(const FftPlan *plan, float *re, float *im, size_t stride, int nSignals)
{
    switch (nSignals) {                                /* common widths: unrolled */
    case 1:  fft_core(plan, re, im, stride, 1);  break;
    case 24: fft_core(plan, re, im, stride, 24); break;
    case 32: fft_core(plan, re, im, stride, 32); break;
    case 64: fft_core(plan, re, im, stride, 64); break;
    default: fft_core(plan, re, im, stride, nSignals);
    }
}

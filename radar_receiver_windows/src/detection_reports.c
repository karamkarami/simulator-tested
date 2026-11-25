/* Detection reports and beam noise floor estimates. */
#include "detection_reports.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#if defined(__AVX2__)
#include <immintrin.h>
#endif

#define MEDIAN_SAMPLE 2048      /* noise floor: sample that brackets the median */

/* Noise floor median: lane numbers of the set bits of every 8-bit mask, in
 * order, and their count (left-packing of 8 values with AVX2). The same
 * values are written by every processor_create. */
static uint8_t g_packIdx[256][8], g_packCnt[256];

void init_pack_table(void)
{
    for (int m = 0; m < 256; m++) {
        int k = 0;
        for (int lane = 0; lane < 8; lane++)
            if (m >> lane & 1) g_packIdx[m][k++] = (uint8_t)lane;
        g_packCnt[m] = (uint8_t)k;
    }
}

/* k-th smallest value (0-based) of v[0..n-1]; v is reordered. */
float processor_select_kth(float *v, int n, int k)
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
double map_median(const float *base, int nRows, int rowLen, size_t stride, float *work, float *guess)
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
        const float lo = sLo <= 0 ? -FLT_MAX : processor_select_kth(work, S, sLo);
        const float hi = sHi >= S - 1 ? FLT_MAX : processor_select_kth(work, S, sHi);
        below = cnt = 0;
        for (int d = 0; d < nRows; d++) bracket_pass(base + (size_t)d * stride, rowLen, lo, hi, work, &cnt, &below);
        found = below <= kLo && kHi < below + cnt;
    }
    if (!found) {
        cnt = below = 0;
        for (int d = 0; d < nRows; d++, cnt += rowLen)
            memcpy(work + cnt, base + (size_t)d * stride, sizeof(float) * rowLen);
    }
    const float vHi = processor_select_kth(work, cnt, kHi - below);   /* work[0 ... kHi-below-1] <= vHi */
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
        s->reports = proc_realloc(s->reports, sizeof(Detection) * s->capReports);
    }
    s->reports[s->nReports++] = *d;
}

void build_detection_reports(const Processor *p, int pulse, PulseState *s,
                             const DetectionCell *cells, const unsigned char *keep, int n)
{
    const RadarConfig *cfg = p->cfg;
    const PulseConfig *pc = &cfg->pulse[pulse];
    const int nB = cfg->nBeams, nR = s->nR;

    s->nReports = 0;
    double sinSum = 0, cosSum = 0, tSum = 0;
    for (int k = 0; k < s->nSlots; k++) {
        sinSum += sin(s->az[k] * PROCESSOR_PI / 180);
        cosSum += cos(s->az[k] * PROCESSOR_PI / 180);
        tSum += s->time[k];
    }
    double azCenter = atan2(sinSum / s->nSlots, cosSum / s->nSlots) * 180 / PROCESSOR_PI;
    azCenter = fmod(azCenter + 360.0, 360.0);
    /* noise floors of the beams with detections, in parallel (not needed with
     * the CFAR noise of every cell) */
    double *floors = proc_calloc(nB, sizeof(double));
    int *beams = proc_calloc(nB, sizeof(int)), nNeed = 0;
    for (int i = 0; i < n && !s->noise; i++)
        if (keep[i] && floors[cells[i].b] == 0) { floors[cells[i].b] = -1; beams[nNeed++] = cells[i].b; }
    #pragma omp parallel for schedule(dynamic) num_threads(team_size(p->nThreads, nNeed)) if(nNeed > 1)
    for (int j = 0; j < nNeed; j++)
        floors[beams[j]] = beam_noise_floor(cfg, s, beams[j], &p->scratch[thread_num()]);
    free(beams);
    for (int i = 0; i < n; i++) {
        if (!keep[i]) continue;
        const DetectionCell *c = &cells[i];
        Detection d;
        memset(&d, 0, sizeof d);
        d.time_s = tSum / s->nSlots;
        d.az_deg = azCenter;
        d.range_m = pc->range_m[c->r];
        d.vr_mps = -cfg->dopFreq_Hz[c->d] * pc->lambda_m / 2;
        d.el_deg = cfg->dbfEnabled ? cfg->tilt_deg + asin(pc->beamU[c->b]) * 180 / PROCESSOR_PI : NAN;   /* beam centre (on its carrier) */
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
}

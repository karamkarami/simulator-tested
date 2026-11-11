/*
 * video.c - raw video of the PPI ("Video" display mode), as ppi_video_add.m.
 *
 * The DDC output of a PRI is averaged non-coherently over the receive
 * elements (mean of I^2 + Q^2, before the beamforming and the pulse
 * compression). The receive intervals of both sub-pulses are cut out and
 * joined: the short-pulse samples nearer than the first long-pulse range
 * cell (the handover), then all long-pulse samples. Receive sample k of a
 * sub-pulse is the range of an echo whose leading edge arrives at it,
 * range_m[0] + k * c / (2 fsOut). The result goes into range bins of that
 * size from 0 to the display range (bins without a sample stay 0).
 */
#include "radar_rx.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

int video_layout_init(VideoLayout *v, const RadarConfig *cfg, double maxRange_m)
{
    memset(v, 0, sizeof *v);
    if (!cfg || !(cfg->fsOut_Hz > 0) || !(maxRange_m > 0)) return -1;
    v->dR_m = 299792458.0 / (2.0 * cfg->fsOut_Hz);
    v->nBins = (int)ceil(maxRange_m / v->dR_m);
    if (v->nBins < 1) v->nBins = 1;
    const PulseConfig *ps = &cfg->pulse[PULSE_SHORT], *pl = &cfg->pulse[PULSE_LONG];
    const int cap = ps->rxCount + pl->rxCount;
    v->sample = malloc(sizeof(int) * (cap > 0 ? cap : 1));
    v->bin = malloc(sizeof(int) * (cap > 0 ? cap : 1));
    v->binInv = calloc((size_t)v->nBins, sizeof(float));
    v->acc = malloc(sizeof(float) * (cap > 0 ? cap : 1));
    if (!v->sample || !v->bin || !v->binInv || !v->acc) {
        video_layout_free(v);
        return -1;
    }
    const double handover_m = cfg->handover_m;       /* (Inf without the long pulse) */
    (void)pl;
    for (int p = 0; p < N_PULSES; p++) {
        const PulseConfig *pc = &cfg->pulse[p];
        if (!pc->enable) continue;
        for (int k = 0; k < pc->rxCount; k++) {
            const double r = pc->range_m[0] + k * v->dR_m;
            if (p == PULSE_SHORT && !(r < handover_m)) break;
            const int b = (int)floor(r / v->dR_m);
            const int n = pc->rxFirst + k;
            if (b < 0 || b >= v->nBins || n < 0 || n >= cfg->nDdc) continue;
            v->sample[v->nSamples] = n;
            v->bin[v->nSamples++] = b;
            v->binInv[b] += 1.0f;
        }
    }
    for (int b = 0; b < v->nBins; b++)                    /* 1 / (samples x channels) of every bin */
        v->binInv[b] = v->binInv[b] > 0 ? 1.0f / (v->binInv[b] * (float)cfg->nCh) : 0.0f;
    return 0;
}

void video_row(const VideoLayout *v, const PriData *pri, float *row)
{
    memset(row, 0, sizeof(float) * (size_t)v->nBins);
    if (!pri->iq || v->nSamples == 0) return;
    float *acc = v->acc;
    memset(acc, 0, sizeof(float) * (size_t)v->nSamples);
    /* the joined samples are runs of consecutive DDC samples (the receive
     * intervals): contiguous loops over every channel, which vectorise */
    for (int i0 = 0; i0 < v->nSamples;) {
        int n = 1;
        while (i0 + n < v->nSamples && v->sample[i0 + n] == v->sample[i0] + n) n++;
        const int first = v->sample[i0];
        if (first + n <= pri->nSamples) {
            for (int ch = 0; ch < pri->nCh; ch++) {
                const int16_t *restrict x = pri->iq + 2 * ((size_t)ch * (size_t)pri->nSamples + (size_t)first);
                float *restrict a = acc + i0;
                for (int k = 0; k < n; k++) {
                    const float re = x[2 * k], im = x[2 * k + 1];
                    a[k] += re * re + im * im;
                }
            }
        }
        i0 += n;
    }
    for (int i = 0; i < v->nSamples; i++) row[v->bin[i]] += acc[i];
    for (int b = 0; b < v->nBins; b++) row[b] *= v->binInv[b];
}

void video_layout_free(VideoLayout *v)
{
    free(v->sample);
    free(v->bin);
    free(v->binInv);
    free(v->acc);
    memset(v, 0, sizeof *v);
}

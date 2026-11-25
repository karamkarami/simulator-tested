/* Matched-filter history and pulse metadata for each coherent processing interval. */
#include "processor_internal.h"

void cpi_buffer_reset(PulseState *s)
{
    s->last = -1;
    s->mfLast = -1;
    s->sinceUpdate = 0;
    s->count = 0;
    s->nPushed = 0;
    s->mtiLast = -1;
    s->mtiCount = 0;
}

/* Reserve a matched-filter history row before its samples are written. */
int cpi_buffer_begin(PulseState *s, int nMtiTaps)
{
    const int slot = (s->last + 1) % s->nSlots;
    s->mfLast = (s->mfLast + 1) % nMtiTaps;
    return slot;
}

/* Publish the completed pulse and its metadata together. */
void cpi_buffer_push(PulseState *s, const RadarConfig *cfg, const PriData *pri, int slot)
{
    s->last = slot;
    s->staggerIdx[slot] = pri->staggerIndex >= 1 && pri->staggerIndex <= cfg->nStagger ? pri->staggerIndex : 1;
    s->az[slot] = pri->az_deg;
    s->time[slot] = pri->time_s;
    s->newestPulse = pri->pulseNumber;
    if (s->count < s->nSlots) s->count++;
    s->nPushed = (s->nPushed + 1) % cfg->updateEvery;  /* pulses since the last update (mod N) */
}

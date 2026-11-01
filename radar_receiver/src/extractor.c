/*
 * extractor.c - azimuth peak picking: the detections of the same target in
 * successive (sliding) CPIs are grouped (same beam +-1, range cell +-2,
 * Doppler bin +-1); when a group gets no detection for more than maxMiss CPIs
 * it is closed and its strongest detection becomes a plot (if it has at least
 * minHits detections).
 *
 * A detection joins the oldest open group within those distances, as in
 * plot_extractor.m. The candidates are found in a hash grid of the cells
 * (beam, range cell, Doppler bin / 8) instead of by a search through all open
 * groups, so the time per detection does not grow with the number of open
 * groups (a low CFAR threshold gives thousands of them).
 */
#include "extractor.h"

#include <stdlib.h>
#include <string.h>

#define DOP_SHIFT 3                 /* 8 Doppler bins per cell of the hash grid */

static void *grow(void *p, size_t size)
{
    p = realloc(p, size ? size : 1);
    if (!p) abort();                 /* out of memory (as the rest of the processing) */
    return p;
}

void extractor_init(Extractor *e, int nDop, int minHits, int maxMiss)
{
    memset(e, 0, sizeof(*e));
    e->nDop = nDop;
    e->minHits = minHits;
    e->maxMiss = maxMiss;
    e->weight = 1;
}

void extractor_free(Extractor *e)
{
    free(e->groups);
    free(e->freeSlots);
    free(e->open);
    free(e->bucket);
    memset(e, 0, sizeof(*e));
}

static int bucket_of(const Extractor *e, int beam, int rIdx, int dCell)
{
    unsigned h = (unsigned)beam * 0x9E3779B1u ^ (unsigned)rIdx * 0x85EBCA77u ^ (unsigned)dCell * 0xC2B2AE3Du;
    h ^= h >> 15;
    return (int)(h & (unsigned)(e->nBuckets - 1));
}

static int bucket_of_group(const Extractor *e, const PlotGroup *g)
{
    return bucket_of(e, g->beam, g->rIdx, g->dIdx >> DOP_SHIFT);
}

static void bucket_insert(Extractor *e, int slot)
{
    const int b = bucket_of_group(e, &e->groups[slot]);
    e->groups[slot].next = e->bucket[b];
    e->bucket[b] = slot;
}

static void bucket_remove(Extractor *e, int slot)
{
    int *link = &e->bucket[bucket_of_group(e, &e->groups[slot])];
    while (*link != slot) link = &e->groups[*link].next;
    *link = e->groups[slot].next;
}

/* Room for nGroups open groups: the pool, and at least 2 buckets per group. */
static void reserve(Extractor *e, int nGroups)
{
    if (nGroups > e->capacity) {
        int cap = e->capacity ? e->capacity : 16;
        while (cap < nGroups) cap *= 2;
        e->groups = grow(e->groups, sizeof(PlotGroup) * cap);
        e->freeSlots = grow(e->freeSlots, sizeof(int) * cap);
        e->open = grow(e->open, sizeof(int) * cap);
        for (int s = cap - 1; s >= e->capacity; s--) e->freeSlots[e->nFree++] = s;   /* lowest first */
        e->capacity = cap;
    }
    if (e->nBuckets < 2 * nGroups || e->nBuckets == 0) {
        int nb = e->nBuckets ? e->nBuckets : 256;
        while (nb < 2 * nGroups) nb *= 2;
        e->bucket = grow(e->bucket, sizeof(int) * nb);
        e->nBuckets = nb;
        for (int b = 0; b < nb; b++) e->bucket[b] = -1;
        for (int k = 0; k < e->n; k++) bucket_insert(e, e->open[k]);
    }
}

/* Slot of the oldest open group within beam +-1, range cell +-2 and Doppler
 * bin +-1 (circular) of the detection, -1 if none. */
static int find_group(const Extractor *e, const Detection *d)
{
    int dCell[3], nCell = 0, best = -1;
    long long bestId = 0;
    for (int dd = -1; dd <= 1; dd++) {                /* grid cells of the Doppler bins */
        const int c = (((d->dIdx + dd) % e->nDop + e->nDop) % e->nDop) >> DOP_SHIFT;
        int seen = 0;
        for (int i = 0; i < nCell; i++) seen |= dCell[i] == c;
        if (!seen) dCell[nCell++] = c;
    }
    for (int b = d->beam - 1; b <= d->beam + 1; b++) {
        if (b < 0 || (e->nBeams > 0 && b >= e->nBeams)) continue;
        for (int r = d->rIdx - 2; r <= d->rIdx + 2; r++) {
            if (r < 0) continue;
            for (int i = 0; i < nCell; i++)
                for (int s = e->bucket[bucket_of(e, b, r, dCell[i])]; s >= 0; s = e->groups[s].next) {
                    const PlotGroup *g = &e->groups[s];
                    if (best >= 0 && g->id >= bestId) continue;
                    int dDop = abs(g->dIdx - d->dIdx);
                    if (e->nDop - dDop < dDop) dDop = e->nDop - dDop;     /* circular Doppler distance */
                    if (abs(g->beam - d->beam) <= 1 && abs(g->rIdx - d->rIdx) <= 2 && dDop <= 1) {
                        best = s;
                        bestId = g->id;
                    }
                }
        }
    }
    return best;
}

void extractor_update(Extractor *e, const Detection *dets, int n, int flush,
                      void (*callback)(void *user, const Detection *plot), void *user)
{
    reserve(e, e->n + n);                     /* every detection may open a group */
    for (int i = 0; i < n; i++) {
        const Detection *d = &dets[i];
        int s = find_group(e, d);
        if (s < 0) {
            s = e->freeSlots[--e->nFree];
            PlotGroup *g = &e->groups[s];
            g->beam = d->beam;  g->rIdx = d->rIdx;  g->dIdx = d->dIdx;
            g->hits = e->weight;  g->miss = 0;  g->best = *d;
            g->id = e->nextId++;
            g->updated = 1;
            bucket_insert(e, s);
            e->open[e->n++] = s;
        } else {
            PlotGroup *g = &e->groups[s];
            if (bucket_of(e, d->beam, d->rIdx, d->dIdx >> DOP_SHIFT) != bucket_of_group(e, g)) {
                bucket_remove(e, s);
                g->beam = d->beam;  g->rIdx = d->rIdx;  g->dIdx = d->dIdx;
                bucket_insert(e, s);
            } else {
                g->beam = d->beam;  g->rIdx = d->rIdx;  g->dIdx = d->dIdx;
            }
            g->hits += e->weight;
            g->miss = 0;
            if (d->power > g->best.power) g->best = *d;
            g->updated = 1;
        }
    }

    /* groups without a detection get a miss; close those without detections for too long */
    int kept = 0;
    for (int k = 0; k < e->n; k++) {
        const int s = e->open[k];
        PlotGroup *g = &e->groups[s];
        if (!g->updated) g->miss++;
        g->updated = 0;
        if (g->miss > e->maxMiss || flush) {
            if (g->hits >= e->minHits && callback) {
                Detection plot = g->best;
                plot.hits = g->hits;
                callback(user, &plot);
            }
            bucket_remove(e, s);
            e->freeSlots[e->nFree++] = s;
        } else {
            e->open[kept++] = s;
        }
    }
    e->n = kept;
}

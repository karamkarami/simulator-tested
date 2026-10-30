/*
 * extractor.c - azimuth peak picking: the detections of the same target in
 * successive (sliding) CPIs are grouped (same beam +-1, range cell +-2,
 * Doppler bin +-1); when a group gets no detection for more than maxMiss CPIs
 * it is closed and its strongest detection becomes a plot (if it has at least
 * minHits detections).
 */
#include "extractor.h"

#include <stdlib.h>
#include <string.h>

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
    memset(e, 0, sizeof(*e));
}

static int find_group(const Extractor *e, const Detection *d)
{
    for (int k = 0; k < e->n; k++) {
        const PlotGroup *g = &e->groups[k];
        int dDop = abs(g->dIdx - d->dIdx);
        if (e->nDop - dDop < dDop) dDop = e->nDop - dDop;     /* circular Doppler distance */
        if (abs(g->beam - d->beam) <= 1 && abs(g->rIdx - d->rIdx) <= 2 && dDop <= 1) return k;
    }
    return -1;
}

void extractor_update(Extractor *e, const Detection *dets, int n, int flush,
                      void (*callback)(void *user, const Detection *plot), void *user)
{
    /* groups that exist before this update and were not updated get a miss */
    int nOld = e->n;
    char *updated = calloc((size_t)nOld + n + 1, 1);
    if (!updated) abort();                    /* out of memory (as the rest of the processing) */

    for (int i = 0; i < n; i++) {
        const Detection *d = &dets[i];
        int j = find_group(e, d);
        if (j < 0) {
            if (e->n == e->capacity) {
                e->capacity = e->capacity ? 2 * e->capacity : 16;
                e->groups = realloc(e->groups, sizeof(PlotGroup) * e->capacity);
                if (!e->groups) abort();
            }
            PlotGroup *g = &e->groups[e->n];
            g->beam = d->beam;  g->rIdx = d->rIdx;  g->dIdx = d->dIdx;
            g->hits = e->weight;  g->miss = 0;  g->best = *d;
            updated[e->n++] = 1;
        } else {
            PlotGroup *g = &e->groups[j];
            g->beam = d->beam;  g->rIdx = d->rIdx;  g->dIdx = d->dIdx;
            g->hits += e->weight;
            g->miss = 0;
            if (d->power > g->best.power) g->best = *d;
            updated[j] = 1;
        }
    }

    /* close the groups without detections for too long */
    int kept = 0;
    for (int k = 0; k < e->n; k++) {
        PlotGroup *g = &e->groups[k];
        if (!updated[k]) g->miss++;
        if (g->miss > e->maxMiss || flush) {
            if (g->hits >= e->minHits && callback) {
                Detection plot = g->best;
                plot.hits = g->hits;
                callback(user, &plot);
            }
        } else {
            e->groups[kept++] = *g;
        }
    }
    e->n = kept;
    free(updated);
}

/*
 * extractor.h - azimuth peak picking over successive CPIs (internal).
 * Same algorithm as src/plot_extractor.m of the MATLAB simulator.
 */
#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include "radar_rx.h"

typedef struct {
    int beam, rIdx, dIdx;   /* cell of the latest detection of the group */
    int hits, miss;
    Detection best;         /* strongest detection of the group           */
    long long id;           /* creation order (the first matching group is the oldest) */
    int next;               /* next group of the same hash bucket (-1: none) */
    int updated;            /* got a detection in the current update      */
} PlotGroup;

typedef struct {
    PlotGroup *groups;      /* pool of groups (slots)                     */
    int capacity;
    int *freeSlots, nFree;  /* unused slots                               */
    int *open, n;           /* slots of the open groups, oldest first     */
    int *bucket, nBuckets;  /* hash grid: first group of every bucket (-1: none) */
    long long nextId;
    int nDop;               /* Doppler bins (circular distance)           */
    int nBeams;             /* beams (0 = unknown): no neighbours beyond them are searched */
    int minHits, maxMiss;
    int weight;             /* CPI updates represented by the next update: 1, or
                               more after updates were skipped (overload)  */
} Extractor;

void extractor_init(Extractor *e, int nDop, int minHits, int maxMiss);
void extractor_free(Extractor *e);
/* Adds the detections of one CPI update. Closed groups with enough hits are
 * returned as plots through callback(user, plot). flush = 1 closes all groups. */
void extractor_update(Extractor *e, const Detection *dets, int n, int flush,
                      void (*callback)(void *user, const Detection *plot), void *user);

#endif

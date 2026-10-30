/*
 * tracker.h - EKF tracker of the plots (internal).
 * Same algorithm as src/track_update.m of the MATLAB simulator.
 */
#ifndef TRACKER_H
#define TRACKER_H

#include "radar_rx.h"

#define TRACK_MAX_STATE 9                /* 3-D, constant acceleration             */

typedef struct {
    uint32_t id;
    int      confirmed;
    double   x[TRACK_MAX_STATE];   /* [positions; velocities (; accelerations)], e.g. [x y vx vy] */
    double   P[TRACK_MAX_STATE * TRACK_MAX_STATE];   /* covariance, row-major n x n   */
    double   t;             /* time of the last update                          */
    int      hits, misses;
} TrackState;

typedef struct {            /* a plot waiting in the GNN batch                   */
    double   time_s, az_deg, range_m, el_deg;
    int      plotIndex;     /* counts the plots given to the tracker (1-based)   */
} PendingPlot;

typedef struct {
    int      dim, nd, n;                 /* measurement dimension, values per axis (2 CV, 3 CA), state size */
    double   R[3];                       /* measurement variances                   */
    double   q[3];                       /* white acceleration (CV) / jerk (CA) density per axis */
    double   sigmaV0, sigmaA0, gate, scanTime_s, height_m;
    double   batchGap_s, batchMax_s;     /* GNN batch                               */
    int      confirmHits, maxMissTentative, maxMiss;
    TrackState *tracks;
    int      nTracks, capacity;
    uint32_t nextId;
    PendingPlot *pending;                /* GNN batch                               */
    int      nPending, capPending, nPlots;
    double   pendingFirst, pendingLast;  /* arrival time of its first and newest plot */
} Tracker;

/* Track update (plotIndex = the plot, see PendingPlot) or deletion (deleted = 1, plotIndex 0). */
typedef void (*TrackerCallback)(void *user, const Track *track, int plotIndex);

void tracker_init(Tracker *t, const RadarConfig *cfg);
void tracker_free(Tracker *t);
/* The new plots of time t_s (call it every PRI, also with n = 0), as
 * track_update.m: the batch is associated (GNN) when the antenna has turned
 * batchGap_s without a new plot or batchMax_s after its first plot, the new
 * plots join the batch, flush = 1 associates it at once; then the tracks
 * with too many missed scans are deleted. Every update and deletion is
 * reported through callback (may be NULL). Plots without elevation are not
 * used in 3-D. */
void tracker_update(Tracker *t, const Detection *plots, int n, double t_s, int flush,
                    TrackerCallback callback, void *user);
void tracker_report(const Tracker *t, int index, Track *out);

#endif

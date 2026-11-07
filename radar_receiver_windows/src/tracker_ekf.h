/*
 * tracker_ekf.h - EKF tracker of the plots, GNN association (internal;
 * cfg.track.method = 'ekf'). Same algorithm as src/track_ekf.m of the MATLAB
 * simulator. Used through tracker.h.
 */
#ifndef TRACKER_EKF_H
#define TRACKER_EKF_H

#include "tracker_common.h"

#define TRACK_MAX_STATE 9                /* 3-D, constant acceleration             */

typedef struct {
    uint32_t id;
    int      confirmed;
    double   x[TRACK_MAX_STATE];   /* [positions; velocities (; accelerations)], e.g. [x y vx vy] */
    double   P[TRACK_MAX_STATE * TRACK_MAX_STATE];   /* covariance, row-major n x n   */
    double   t;             /* time of the last update                          */
    int      hits, misses;
} TrackState;

typedef struct {
    int      dim, nd, n;                 /* measurement dimension, values per axis (2 CV, 3 CA), state size */
    double   R[3];                       /* measurement variances                   */
    double   q[3];                       /* white acceleration (CV) / jerk (CA) density per axis */
    double   sigmaV0, sigmaA0, gate, scanTime_s, height_m;
    int      useVr, vrUnambiguous;       /* radial velocity measurement (as track_ekf.m) */
    double   Rvr, vrSpan, gateVr;
    double   batchGap_s, batchMax_s;     /* GNN batch                               */
    int      confirmHits, maxMissTentative, maxMiss;
    double   minSpeed;                   /* delete tracks slower than this [m/s], 0 = off */
    int      minSpeedHits;               /* updates before the speed is tested (2)  */
    TrackState *tracks;
    int      nTracks, capacity;
    uint32_t nextId;
    PendingPlot *pending;                /* GNN batch                               */
    int      nPending, capPending, nPlots;
    double   pendingFirst, pendingLast;  /* arrival time of its first and newest plot */
} EkfTracker;

void ekf_tracker_init(EkfTracker *t, const RadarConfig *cfg);
void ekf_tracker_free(EkfTracker *t);
/* The new plots of time t_s (call it every PRI, also with n = 0), as
 * track_ekf.m: the batch is associated (GNN) when the antenna has turned
 * batchGap_s without a new plot or batchMax_s after its first plot, the new
 * plots join the batch, flush = 1 associates it at once; then the tracks
 * with too many missed scans are deleted, and with minSpeed > 0 the tracks
 * with at least minSpeedHits updates whose speed |v| is below minSpeed. Every update and deletion is
 * reported through callback (may be NULL). Plots without elevation are not
 * used in 3-D. */
void ekf_tracker_update(EkfTracker *t, const Detection *plots, int n, double t_s, int flush,
                    TrackerCallback callback, void *user);
void ekf_tracker_report(const EkfTracker *t, int index, Track *out);

#endif

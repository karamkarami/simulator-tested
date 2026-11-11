/*
 * tracker_common.h - types and helpers shared by the trackers (internal):
 * tracker_ekf.c (EKF, GNN) and tracker_pmbm.c (PMBM, IMM-UKF), selected by
 * tracker.c (cfg.track.method).
 */
#ifndef TRACKER_COMMON_H
#define TRACKER_COMMON_H

#include "radar_rx.h"

typedef struct {            /* a plot waiting in the batch                       */
    double   time_s, az_deg, range_m, el_deg, vr_mps;
    double   vrSpan;        /* its radial velocity is known modulo this (the Doppler span on
                               the carrier of its sub-pulse) [m/s]                */
    int      plotIndex;     /* counts the plots given to the tracker (1-based)   */
} PendingPlot;

/* Track update (plotIndex = the plot, see PendingPlot) or deletion (deleted = 1, plotIndex 0). */
typedef void (*TrackerCallback)(void *user, const Track *track, int plotIndex);

/* Minimum-cost assignment of every row of A (n x m, n <= m, row-major) to its
 * own column (shortest augmenting paths with potentials, as hungarian() in
 * track_ekf.m and track_pmbm.m). rowOf[j] = row (1-based) of column j, 0 = none. */
void tracker_hungarian(const double *A, int n, int m, int *rowOf);

/* The Doppler span of a plot: on the carrier of its sub-pulse (cfg->pulse[].vrSpan_mps), else
 * cfg->trackVrSpan. */
static inline double tracker_plot_vr_span(const RadarConfig *cfg, const Detection *d)
{
    const double s = d && d->pulse >= 0 && d->pulse < N_PULSES ? cfg->pulse[d->pulse].vrSpan_mps : 0.0;
    return s > 0.0 ? s : cfg->trackVrSpan;
}

#endif

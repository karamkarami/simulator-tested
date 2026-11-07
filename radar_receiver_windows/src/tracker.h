/*
 * tracker.h - the tracking block of the receiver (internal): the tracker of
 * cfg.track.method, as track_update.m of the MATLAB simulator:
 *   TRACK_METHOD_EKF   tracker_ekf.c   extended Kalman filter, GNN association
 *   TRACK_METHOD_PMBM  tracker_pmbm.c  PMBM multi-target filter, IMM-UKF
 */
#ifndef TRACKER_H
#define TRACKER_H

#include "tracker_common.h"
#include "tracker_ekf.h"
#include "tracker_pmbm.h"

typedef struct {
    int         method;     /* TRACK_METHOD_EKF or TRACK_METHOD_PMBM            */
    EkfTracker  ekf;
    PmbmTracker pmbm;
} Tracker;

void tracker_init(Tracker *t, const RadarConfig *cfg);
void tracker_free(Tracker *t);
/* The new plots of time t_s (call it every PRI, also with n = 0); flush = 1
 * also associates the plots waiting in the batch. Every track update (with
 * the plot) and deletion is reported through callback (may be NULL). */
void tracker_update(Tracker *t, const Detection *plots, int n, double t_s, int flush,
                    TrackerCallback callback, void *user);
int  tracker_count(const Tracker *t);                    /* current tracks           */
void tracker_report(const Tracker *t, int index, Track *out);

#endif

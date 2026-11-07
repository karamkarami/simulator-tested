/*
 * tracker_pmbm.h - multi-target Bayesian tracker (internal;
 * cfg.track.method = 'pmbm'): Poisson multi-Bernoulli mixture (PMBM) filter
 * with IMM-UKF state estimation. Same algorithm as src/track_pmbm.m of the
 * MATLAB simulator. Used through tracker.h.
 */
#ifndef TRACKER_PMBM_H
#define TRACKER_PMBM_H

#include "tracker_common.h"

#define PMBM_N    10        /* largest state: 3-D [position; velocity; acceleration; turn rate] */
#define PMBM_M    3         /* largest number of IMM models                               */
#define PMBM_HIST 8         /* newest plots kept per single-target hypothesis             */

typedef struct {            /* single-target hypothesis (STH)                             */
    double r;               /* existence probability (0: no target)                       */
    double x[PMBM_M][PMBM_N];               /* state of every model                       */
    double P[PMBM_M][PMBM_N * PMBM_N];      /* covariance of every model (row-major)      */
    double mu[PMBM_M];      /* model probabilities                                        */
    double t, tOpp, tNext;  /* time of the state, of the last look, of the next look      */
    int    hits, misses, plot, batch;       /* plots, missed looks, newest plot (index), its batch */
    int    hist[PMBM_HIST], nHist;          /* newest plots                               */
} PmbmSth;

typedef struct {            /* Bernoulli component = track                                */
    uint32_t id;            /* 0 until it survives its first batch                        */
    int      confirmed;
    double   rMarg;         /* existence over the global hypotheses                       */
    double   rTarget;       /* probability that a target made its newest plot             */
    int      plot;          /* its first plot                                             */
    PmbmSth *h;             /* its single-target hypotheses                               */
    int      nh, cap;
} PmbmTrack;

typedef struct {
    /* parameters (track_setup.m, export_receiver_config.m) */
    int    dim, n, M, models[PMBM_M];
    double R[3], Rvr, vrSpan, sigmaV0, sigmaA0, gate, gateVr;
    int    useVr, vrUnambiguous;
    double scanTime_s, height_m, az0_deg, rate_dps, batchGap_s, batchMax_s;
    int    confirmHits, maxMissTentative, maxMiss, minSpeedHits;
    double minSpeed;
    double pd, ps, logLambdaC, logLambdaU, logVrDensity, logHypPrune, rPrune, rConfirm, rDelete;
    int    maxHyp;
    double mu0[PMBM_M], switchLambda, qCv[3], qCa[3], qCt[3], qTurn, sigmaTurn0, sigmaClimb0;
    double varUnusedA, varUnusedW, ukfScale, Wm0, Wc0, Wi, latency_s, maintStep_s;
    /* multi-Bernoulli mixture */
    PmbmTrack *T;           /* tracks [nT]                                                */
    int    nT, capT;
    int   *G;               /* global hypotheses: STH of every track [nG x nT], row-major */
    double *logw;           /* their log weights [nG] (normalised)                        */
    int    nG;
    uint32_t nextId;
    int    batch, logged, nPlots;
    double tMaint;
    PendingPlot *pending;   /* the batch                                                  */
    int    nPending, capPending;
    double pendingFirst, pendingLast;
    uint32_t *confId;       /* confirmed tracks shown: number and newest plot             */
    int   *confPlot, nConf, capConf;
    Track *shown;           /* the tracks of the most likely global hypothesis [nShown]   */
    int    nShown, capShown;
} PmbmTracker;

void pmbm_tracker_init(PmbmTracker *t, const RadarConfig *cfg);
void pmbm_tracker_free(PmbmTracker *t);
/* As track_pmbm.m: the new plots of time t_s (every PRI, also n = 0); flush
 * associates the waiting batch. callback gets every update of a track of the
 * most likely global hypothesis by a plot (as the MATLAB track log) and every
 * track that leaves it (deleted = 1). */
void pmbm_tracker_update(PmbmTracker *t, const Detection *plots, int n, double t_s, int flush,
                         TrackerCallback callback, void *user);
void pmbm_tracker_report(const PmbmTracker *t, int index, Track *out);

#endif

/*
 * tracker_tomht.h - track-oriented multiple hypothesis tracker (internal;
 * cfg.track.method = 'tomht'): track trees scored by their log-likelihood
 * ratio, k-best global hypotheses (Murty) and N-scan pruning, IMM-UKF state
 * estimation of every track hypothesis. Same algorithm as src/track_tomht.m of
 * the MATLAB simulator. Used through tracker.h.
 */
#ifndef TRACKER_TOMHT_H
#define TRACKER_TOMHT_H

#include "tracker_common.h"
#include "tracker_pmbm.h"            /* PMBM_N, PMBM_M, PMBM_HIST: the IMM-UKF sizes */

#define TOMHT_DEC      24            /* newest decisions (looks) kept per track hypothesis */
#define TOMHT_MAX_NSCAN 16           /* largest N of the N-scan pruning (< TOMHT_DEC)   */

typedef struct {            /* track hypothesis (a leaf of a track tree)                   */
    double r;               /* 1: a target, 0: "no target" (the plots are false)           */
    double x[PMBM_M][PMBM_N];               /* state of every model                       */
    double P[PMBM_M][PMBM_N * PMBM_N];      /* covariance of every model (row-major)      */
    double mu[PMBM_M];      /* model probabilities                                        */
    double t, tOpp, tNext;  /* time of the state, of the last look, of the next look      */
    int    hits, misses, plot, batch;       /* plots, missed looks, newest plot, its batch */
    int    hist[PMBM_HIST], nHist;          /* newest plots                               */
    double score, maxScore; /* log-likelihood ratio, its largest value                     */
    double decT[TOMHT_DEC]; /* decisions of the newest looks: time and plot (0 = missed)  */
    int    decP[TOMHT_DEC], nDec;
    uint64_t decHash;       /* the older decisions (FNV-1a in their order)                */
    int    nDecOld;
} TomhtSth;

typedef struct {            /* track tree                                                 */
    uint32_t id;
    int      confirmed;
    double   rMarg;         /* probability over the global hypotheses                     */
    double   rTarget;       /* probability that a target made its newest plot             */
    int      plot;          /* its first plot                                             */
    double   t0;            /* time of its first plot                                     */
    TomhtSth *h;            /* its leaves                                                 */
    int      nh, cap;
} TomhtTrack;

typedef struct {
    /* parameters (track_setup.m, export_receiver_config.m) */
    int    dim, n, M, models[PMBM_M];
    double R[3], Rvr, vrSpan, sigmaV0, sigmaA0, gate, gateVr;
    double vrSpanPulse[N_PULSES];
    int    useVr, vrUnambiguous;
    double scanTime_s, height_m, az0_deg, rate_dps, batchGap_s, batchMax_s;
    int    confirmHits, maxMissTentative, maxMiss, minSpeedHits;
    double minSpeed;
    double pd, logLambdaC, logLambdaU, logVrDensity, logHypPrune;
    int    maxHyp, kBest, nScan;
    double confirmScore, deleteScore, pDelete;
    double mu0[PMBM_M], switchLambda, qCv[3], qCa[3], qCt[3], qTurn, sigmaTurn0, sigmaClimb0;
    double varUnusedA, varUnusedW, ukfScale, Wm0, Wc0, Wi, latency_s, maintStep_s;
    /* track trees and global hypotheses */
    TomhtTrack *T;
    int    nT, capT;
    int   *G;               /* global hypotheses: leaf of every tree [nG x nT], row-major */
    double *logw;           /* their scores [nG]                                          */
    int    nG;
    uint32_t nextId;
    int    batch, logged, nPlots;
    double tMaint, tNow;
    PendingPlot *pending;
    int    nPending, capPending;
    double pendingFirst, pendingLast;
    uint32_t *confId;
    int   *confPlot, nConf, capConf;
    Track *shown;           /* the tracks of the best global hypothesis [nShown]          */
    int    nShown, capShown;
} TomhtTracker;

void tomht_tracker_init(TomhtTracker *t, const RadarConfig *cfg);
void tomht_tracker_free(TomhtTracker *t);
void tomht_tracker_update(TomhtTracker *t, const Detection *plots, int n, double t_s, int flush,
                          TrackerCallback callback, void *user);
void tomht_tracker_report(const TomhtTracker *t, int index, Track *out);

#endif

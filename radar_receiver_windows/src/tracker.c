/*
 * tracker.c - the tracking block: the tracker of cfg.track.method (as
 * track_update.m), and the assignment shared by both trackers.
 */
#include "tracker.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void tracker_init(Tracker *t, const RadarConfig *cfg)
{
    memset(t, 0, sizeof *t);
    t->method = cfg->trackMethod;
    if (t->method == TRACK_METHOD_PMBM) pmbm_tracker_init(&t->pmbm, cfg);
    else if (t->method == TRACK_METHOD_TOMHT) tomht_tracker_init(&t->tomht, cfg);
    else ekf_tracker_init(&t->ekf, cfg);
}

void tracker_free(Tracker *t)
{
    if (t->method == TRACK_METHOD_PMBM) pmbm_tracker_free(&t->pmbm);
    else if (t->method == TRACK_METHOD_TOMHT) tomht_tracker_free(&t->tomht);
    else ekf_tracker_free(&t->ekf);
}

void tracker_update(Tracker *t, const Detection *plots, int n, double t_s, int flush,
                    TrackerCallback callback, void *user)
{
    if (t->method == TRACK_METHOD_PMBM) pmbm_tracker_update(&t->pmbm, plots, n, t_s, flush, callback, user);
    else if (t->method == TRACK_METHOD_TOMHT) tomht_tracker_update(&t->tomht, plots, n, t_s, flush, callback, user);
    else ekf_tracker_update(&t->ekf, plots, n, t_s, flush, callback, user);
}

int tracker_count(const Tracker *t)
{
    return t->method == TRACK_METHOD_PMBM ? t->pmbm.nShown : t->method == TRACK_METHOD_TOMHT ? t->tomht.nShown
                                                                                       : t->ekf.nTracks;
}

void tracker_report(const Tracker *t, int index, Track *out)
{
    if (t->method == TRACK_METHOD_PMBM) pmbm_tracker_report(&t->pmbm, index, out);
    else if (t->method == TRACK_METHOD_TOMHT) tomht_tracker_report(&t->tomht, index, out);
    else ekf_tracker_report(&t->ekf, index, out);
}

void tracker_name(const RadarConfig *cfg, char *buf, size_t len)
{
    if (cfg->trackMethod == TRACK_METHOD_PMBM || cfg->trackMethod == TRACK_METHOD_TOMHT) {
        static const char *names[3] = { "CV", "CA", "CT" };
        char models[16] = "";
        for (int j = 0; j < cfg->pmbmNModels; j++) {
            if (j) strcat(models, "/");
            strcat(models, names[cfg->pmbmModels[j]]);
        }
        snprintf(buf, len, "%s %d-D IMM %s", cfg->trackMethod == TRACK_METHOD_TOMHT ? "TOMHT" : "PMBM", cfg->trackDim,
                 models);
    } else {
        snprintf(buf, len, "EKF %d-D %s", cfg->trackDim, cfg->trackModel == TRACK_CA ? "CA" : "CV");
    }
}

void tracker_hungarian(const double *A, int n, int m, int *rowOf)
{
    double *u = calloc(n + 1, sizeof(double)), *v = calloc(m + 1, sizeof(double));
    double *minv = malloc(sizeof(double) * (m + 1));
    int *p = calloc(m + 1, sizeof(int)), *way = calloc(m + 1, sizeof(int));
    char *used = malloc(m + 1);
    if (!u || !v || !minv || !p || !way || !used) abort();
    for (int i = 1; i <= n; i++) {
        p[0] = i;
        int j0 = 0;
        for (int j = 0; j <= m; j++) { minv[j] = INFINITY; used[j] = 0; }
        do {
            used[j0] = 1;
            const int i0 = p[j0];
            double delta = INFINITY;
            int j1 = 0;
            for (int j = 1; j <= m; j++)
                if (!used[j]) {
                    const double cur = A[(size_t)(i0 - 1) * m + (j - 1)] - u[i0] - v[j];
                    if (cur < minv[j]) { minv[j] = cur; way[j] = j0; }
                    if (minv[j] < delta) { delta = minv[j]; j1 = j; }
                }
            for (int j = 0; j <= m; j++)
                if (used[j]) { u[p[j]] += delta; v[j] -= delta; }
                else minv[j] -= delta;
            j0 = j1;
        } while (p[j0] != 0);
        do {                                               /* augment along the path */
            const int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    for (int j = 1; j <= m; j++) rowOf[j - 1] = p[j];
    free(u); free(v); free(minv); free(p); free(way); free(used);
}

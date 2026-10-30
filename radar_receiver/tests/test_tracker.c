/*
 * test_tracker.c - the C EKF tracker (tracker.c) against the MATLAB tracker.
 *
 *   test_tracker CONFIG                      self tests: two crossing targets, 3-D and 2-D;
 *                                            GNN with two targets 200 m apart
 *   test_tracker CONFIG PLOTS.csv LOG.csv TEND
 *       CONFIG:    exported for the case (export_receiver_config.m): dimension,
 *                  motion model, gate and GNN batch as in tests/test_tracker.m
 *       PLOTS.csv: time_s,az_deg,range_m,el_deg of every plot in processing order
 *       LOG.csv:   the MATLAB track log (track_update.m): time,id,confirmed,x,y,z,vx,vy,vz,hits,plotIndex
 *       TEND:      time of the last call (the open batch is associated)
 *       Every plot is given to the tracker in a call of its own at the plot
 *       time (as tests/test_tracker.m does), then a call at TEND with flush;
 *       every track update must match the log in the same order: same plot,
 *       track number, state and number of updates, positions and velocities
 *       within 1e-6 relative.
 *
 * Prints PASS or FAIL and returns 0 on success.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "radar_rx.h"
#include "tracker.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
static int failures = 0;

#define CHECK(cond, ...)                                                     \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL: " __VA_ARGS__);                                    \
            printf("\n");                                                    \
            failures++;                                                      \
        }                                                                    \
    } while (0)

static double gauss(void)
{
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0), u2 = rand() / (RAND_MAX + 1.0);
    return sqrt(-2 * log(u1)) * cos(2 * PI * u2);
}

/* Two targets crossing each other, 6 scans, plot errors as configured. */
static void self_test(RadarConfig *cfg, int dim)
{
    cfg->trackDim = dim;
    Tracker t;
    tracker_init(&t, cfg);
    const double p0[2][3] = { { 20e3, 5e3, 2000 }, { 20e3, -5e3, 3000 } };
    const double v[2][3] = { { -50, -150, 0 }, { -50, 150, 0 } };
    const double T = cfg->scanTime_s, h = cfg->height_m;
    srand(3);
    for (int scan = 0; scan < 6; scan++) {
        for (int i = 0; i < 2; i++) {
            const double tt = scan * T + i * 0.1;
            const double x = p0[i][0] + v[i][0] * tt, y = p0[i][1] + v[i][1] * tt, z = p0[i][2] - h;
            const double R = sqrt(x * x + y * y + z * z);
            Detection d;
            memset(&d, 0, sizeof d);
            d.time_s = tt;
            d.range_m = R + cfg->trackSigma[0] * gauss();
            d.az_deg = (atan2(x, y) + cfg->trackSigma[1] * gauss()) * 180.0 / PI;
            d.el_deg = dim == 3 ? (asin(z / R) + cfg->trackSigma[2] * gauss()) * 180.0 / PI : NAN;
            tracker_update(&t, &d, 1, tt, 0, NULL, NULL);
        }
    }
    tracker_update(&t, NULL, 0, 6 * T, 1, NULL, NULL);
    int confirmed = 0;
    for (int j = 0; j < t.nTracks; j++) {
        Track r;
        tracker_report(&t, j, &r);
        if (!r.confirmed) continue;
        confirmed++;
        const double tt = r.time_s;
        double dist[2];
        for (int k = 0; k < 2; k++)
            dist[k] = hypot(r.pos_m[0] - (p0[k][0] + v[k][0] * tt), r.pos_m[1] - (p0[k][1] + v[k][1] * tt));
        const int i = dist[0] <= dist[1] ? 0 : 1;          /* nearest target (they cross at t = 33 s) */
        const double ex = r.pos_m[0] - (p0[i][0] + v[i][0] * tt), ey = r.pos_m[1] - (p0[i][1] + v[i][1] * tt);
        const double ev = hypot(r.vel_mps[0] - v[i][0], r.vel_mps[1] - v[i][1]);
        printf("  %d-D: track %u, %d updates, position error %.0f m, velocity error %.1f m/s", dim, r.id, r.hits,
               hypot(ex, ey), ev);
        if (dim == 3) printf(", height %.0f m (true %.0f)", r.alt_m, p0[i][2]);
        printf("\n");
        CHECK(r.hits == 6, "%d-D: track %u has %d updates (expected 6)", dim, r.id, r.hits);
        if (dim == 3) CHECK(hypot(ex, ey) < 300 && ev < 30, "%d-D: track %u too far from its target", dim, r.id);
        if (dim == 3) CHECK(fabs(r.alt_m - p0[i][2]) < 400, "%d-D: height of track %u", dim, r.id);
    }
    CHECK(confirmed == 2, "%d-D: %d confirmed tracks (expected 2)", dim, confirmed);
    tracker_free(&t);
}

/* GNN: two targets 200 m apart in range at azimuth 30 deg, receding at 100 m/s,
 * plots without errors (as gnn_case in tests/test_tracker.m). At the second
 * scan plot 1 is closer to track 2: only the global assignment keeps the tracks. */
static void gnn_test(RadarConfig *cfg)
{
    cfg->trackDim = 2;
    Tracker t;
    tracker_init(&t, cfg);
    const double R0[2] = { 20e3, 20.2e3 }, v = 100.0, az = 30.0, T = cfg->scanTime_s;
    for (int scan = 0; scan < 6; scan++) {
        Detection d[2];
        memset(d, 0, sizeof d);
        for (int i = 0; i < 2; i++) {
            d[i].time_s = scan * T;
            d[i].range_m = R0[i] + v * scan * T;
            d[i].az_deg = az;
            d[i].el_deg = NAN;
        }
        tracker_update(&t, d, 2, scan * T, 1, NULL, NULL);
    }
    CHECK(t.nTracks == 2, "GNN: %d tracks (expected 2)", t.nTracks);
    for (int j = 0; j < t.nTracks; j++) {
        Track r;
        tracker_report(&t, j, &r);
        const double R = hypot(r.pos_m[0], r.pos_m[1]);
        const int i = fabs(R - R0[1] - v * r.time_s) < fabs(R - R0[0] - v * r.time_s);
        const double dv = hypot(r.vel_mps[0] - v * sin(az * PI / 180), r.vel_mps[1] - v * cos(az * PI / 180));
        printf("  GNN: track %u on target %d, %d updates, range error %.2f m, velocity error %.3f m/s\n", r.id, i + 1,
               r.hits, R - R0[i] - v * r.time_s, dv);
        CHECK(r.confirmed && r.hits == 6 && (int)r.id == i + 1 && fabs(R - R0[i] - v * r.time_s) < 5 && dv < 1,
              "GNN: track %u does not follow target %d", r.id, i + 1);
    }
    tracker_free(&t);
}

typedef struct {
    FILE *log;
    double ref[11];
    int haveRef, updates, bad;
    double maxRel;
} Compare;

static int read_ref(Compare *c)
{
    double *r = c->ref;
    return fscanf(c->log, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &r[0], &r[1], &r[2], &r[3], &r[4], &r[5],
                  &r[6], &r[7], &r[8], &r[9], &r[10]) == 11;
}

static void on_update(void *user, const Track *r, int plotIndex)
{
    Compare *c = user;
    if (r->deleted) return;
    c->updates++;
    if (!c->haveRef || (int)c->ref[10] != plotIndex) {
        c->bad++;
        printf("FAIL: plot %d updated track %u, MATLAB updated with plot %.0f\n", plotIndex, r->id,
               c->haveRef ? c->ref[10] : -1.0);
        return;
    }
    const double v[9] = { r->time_s, (double)r->id, (double)r->confirmed, r->pos_m[0], r->pos_m[1], r->pos_m[2],
                          r->vel_mps[0], r->vel_mps[1], r->vel_mps[2] };
    int same = v[1] == c->ref[1] && v[2] == c->ref[2] && r->hits == (int)c->ref[9];
    for (int k = 3; k < 9; k++) {
        const double rel = fabs(v[k] - c->ref[k]) / (fabs(c->ref[k]) + 1.0);
        if (rel > c->maxRel) c->maxRel = rel;
        same &= rel < 1e-6;
    }
    if (!same) {
        c->bad++;
        if (c->bad < 5) printf("FAIL: plot %d: C track %u (conf %d, hits %d), MATLAB track %.0f (conf %.0f, hits %.0f)\n",
                               plotIndex, r->id, r->confirmed, r->hits, c->ref[1], c->ref[2], c->ref[9]);
    }
    c->haveRef = read_ref(c);
}

static int compare_with_matlab(RadarConfig *cfg, const char *plotsPath, const char *logPath, double tEnd)
{
    FILE *fp = fopen(plotsPath, "r");
    Compare c;
    memset(&c, 0, sizeof c);
    c.log = fopen(logPath, "r");
    if (!fp || !c.log) {
        printf("FAIL: cannot open %s or %s\n", plotsPath, logPath);
        if (fp) fclose(fp);
        if (c.log) fclose(c.log);
        return 0;
    }
    Tracker t;
    tracker_init(&t, cfg);
    c.haveRef = read_ref(&c);
    int nPlots = 0;
    double tp, az, R, el;
    while (fscanf(fp, "%lf,%lf,%lf,%lf", &tp, &az, &R, &el) == 4) {
        nPlots++;
        Detection d;
        memset(&d, 0, sizeof d);
        d.time_s = tp;  d.az_deg = az;  d.range_m = R;  d.el_deg = el;
        tracker_update(&t, &d, 1, tp, 0, on_update, &c);
    }
    tracker_update(&t, NULL, 0, tEnd, 1, on_update, &c);
    if (c.haveRef) { c.bad++; printf("FAIL: MATLAB has more track updates than C\n"); }
    printf("  %d-D %s: %d plots, %d track updates, %d different from MATLAB, max relative state difference %.1e\n",
           t.dim, t.nd == 3 ? "CA" : "CV", nPlots, c.updates, c.bad, c.maxRel);
    fclose(fp);
    fclose(c.log);
    tracker_free(&t);
    failures += c.bad > 0;
    return c.bad == 0;
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("usage: test_tracker CONFIG [PLOTS.csv LOG.csv TEND]\n");
        return 2;
    }
    RadarConfig cfg;
    char err[256];
    if (radar_config_load(argv[1], &cfg, err, sizeof err) != 0) {
        printf("FAIL: configuration %s: %s\n", argv[1], err);
        return 1;
    }
    if (!cfg.trackEnable) printf("(tracking is off in the configuration; its parameters are used anyway)\n");
    if (argc >= 5) {
        compare_with_matlab(&cfg, argv[2], argv[3], atof(argv[4]));
    } else {
        self_test(&cfg, 3);
        self_test(&cfg, 2);
        gnn_test(&cfg);
    }
    radar_config_free(&cfg);
    printf(failures ? "FAIL\n" : "PASS\n");
    return failures ? 1 : 0;
}

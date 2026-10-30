/*
 * tracker.c - extended Kalman filter tracking of the plots (as track_update.m).
 *
 * State relative to the antenna (x East, y North, z up), per axis the
 * position and its derivatives, ordered [positions; velocities (;
 * accelerations)]: 2-D [x y vx vy], 3-D [x y z vx vy vz], constant
 * acceleration [... ax ay (az)]. Prediction per axis: constant velocity
 * F = [1 dt; 0 1], Q = q*[dt^3/3 dt^2/2; dt^2/2 dt] (white acceleration), or
 * constant acceleration F = [1 dt dt^2/2; 0 1 dt; 0 0 1],
 * Q = q*[dt^5/20 dt^4/8 dt^3/6; dt^4/8 dt^3/3 dt^2/2; dt^3/6 dt^2/2 dt]
 * (white jerk). Measurement of a plot: 2-D [ground range, azimuth] (slant
 * range when the elevation is unknown), 3-D [range, azimuth, elevation];
 * h(x) is linearised at the predicted state.
 * Association: global nearest neighbour on batches of plots. A batch is
 * complete when the antenna has turned batchGap_s without a new plot, or
 * batchMax_s after its first plot. Cost of a plot-track pair d^2 + ln|S| when
 * the plot lies in the chi-square gate of the track; the assignment (every
 * track at most one plot) has the most pairs and, among those, the smallest
 * total cost (Hungarian algorithm). The plots are then taken in their order:
 * an assigned plot updates its track, the others start tentative tracks
 * (zero velocity and acceleration, std sigmaV0, sigmaA0). Confirmation after
 * confirmHits updates; deletion after more than maxMiss (tentative:
 * maxMissTentative) missed scans, a scan being missed when the antenna passed
 * the track half a scan ago without an update.
 */
#include "tracker.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NMAX TRACK_MAX_STATE

void tracker_init(Tracker *t, const RadarConfig *cfg)
{
    memset(t, 0, sizeof *t);
    t->dim = cfg->trackDim == 3 ? 3 : 2;
    t->nd = cfg->trackModel == TRACK_CA ? 3 : 2;
    t->n = t->nd * t->dim;
    for (int i = 0; i < 3; i++) {
        t->R[i] = cfg->trackSigma[i] * cfg->trackSigma[i];
        t->q[i] = cfg->trackQ[i];
    }
    t->sigmaV0 = cfg->trackSigmaV0;
    t->sigmaA0 = cfg->trackSigmaA0;
    t->gate = cfg->trackGate;
    t->scanTime_s = cfg->scanTime_s;
    t->batchGap_s = cfg->trackBatchGap_s;
    t->batchMax_s = cfg->trackBatchMax_s;
    t->height_m = cfg->height_m;
    t->confirmHits = cfg->trackConfirmHits;
    t->maxMissTentative = cfg->trackMaxMissTentative;
    t->maxMiss = cfg->trackMaxMiss;
    t->nextId = 1;
    t->pendingFirst = t->pendingLast = -INFINITY;
}

void tracker_free(Tracker *t)
{
    free(t->tracks);
    free(t->pending);
    memset(t, 0, sizeof *t);
}

/* ---- small dense matrices (row-major) ----------------------------------- */

/* C = A (r x k) * B (k x c) */
static void mul(const double *A, const double *B, double *C, int r, int k, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            double s = 0.0;
            for (int m = 0; m < k; m++) s += A[i * k + m] * B[m * c + j];
            C[i * c + j] = s;
        }
}

/* C = A (r x k) * B' (B: c x k) */
static void mul_bt(const double *A, const double *B, double *C, int r, int k, int c)
{
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++) {
            double s = 0.0;
            for (int m = 0; m < k; m++) s += A[i * k + m] * B[j * k + m];
            C[i * c + j] = s;
        }
}

/* Inverse and determinant of a symmetric 2x2 or 3x3 matrix. Returns 0 if singular. */
static int inv_small(const double *S, int m, double *Si, double *det)
{
    if (m == 2) {
        *det = S[0] * S[3] - S[1] * S[2];
        if (!(*det > 0.0)) return 0;
        Si[0] = S[3] / *det;  Si[1] = -S[1] / *det;
        Si[2] = -S[2] / *det; Si[3] = S[0] / *det;
        return 1;
    }
    const double a = S[0], b = S[1], c = S[2], d = S[3], e = S[4], f = S[5], g = S[6], h = S[7], k = S[8];
    const double A = e * k - f * h, B = -(d * k - f * g), C = d * h - e * g;
    *det = a * A + b * B + c * C;
    if (!(*det > 0.0)) return 0;
    Si[0] = A / *det;  Si[1] = -(b * k - c * h) / *det; Si[2] = (b * f - c * e) / *det;
    Si[3] = B / *det;  Si[4] = (a * k - c * g) / *det;  Si[5] = -(a * f - c * d) / *det;
    Si[6] = C / *det;  Si[7] = -(a * h - b * g) / *det; Si[8] = (a * e - b * d) / *det;
    return 1;
}

/* ---- EKF ------------------------------------------------------------------ */

/* F = kron(F1, I_dim), Q = kron(Q1, diag(q)) over dt */
static void predict(const Tracker *t, const TrackState *s, double dt, double *x, double *P)
{
    const int d = t->dim, n = t->n, nd = t->nd;
    const double a = fabs(dt);
    double F1[9], Q1[9];
    if (nd == 2) {
        const double f[4] = { 1.0, dt, 0.0, 1.0 };
        const double q[4] = { a * a * a / 3.0, a * a / 2.0, a * a / 2.0, a };
        memcpy(F1, f, sizeof f);
        memcpy(Q1, q, sizeof q);
    } else {
        const double f[9] = { 1.0, dt, dt * dt / 2.0, 0.0, 1.0, dt, 0.0, 0.0, 1.0 };
        const double a2 = a * a, a3 = a2 * a, a4 = a3 * a, a5 = a4 * a;
        const double q[9] = { a5 / 20.0, a4 / 8.0, a3 / 6.0, a4 / 8.0, a3 / 3.0, a2 / 2.0, a3 / 6.0, a2 / 2.0, a };
        memcpy(F1, f, sizeof f);
        memcpy(Q1, q, sizeof q);
    }
    double F[NMAX * NMAX], FP[NMAX * NMAX];
    memset(F, 0, sizeof F);
    for (int bi = 0; bi < nd; bi++)
        for (int bj = 0; bj < nd; bj++)
            for (int i = 0; i < d; i++) F[(bi * d + i) * n + bj * d + i] = F1[bi * nd + bj];
    mul(F, s->x, x, n, n, 1);
    mul(F, s->P, FP, n, n, n);
    mul_bt(FP, F, P, n, n, n);
    for (int bi = 0; bi < nd; bi++)
        for (int bj = 0; bj < nd; bj++)
            for (int i = 0; i < d; i++) P[(bi * d + i) * n + bj * d + i] += Q1[bi * nd + bj] * t->q[i];
}

/* Innovation z - h(x) (azimuth wrapped) and Jacobian H (dim x n). */
static void innovation(const Tracker *t, const double *z, const double *x, double *nu, double *H)
{
    const int n = t->n;
    const double px = x[0], py = x[1], rh2 = px * px + py * py, rh = sqrt(rh2);
    memset(H, 0, sizeof(double) * t->dim * n);
    if (t->dim == 3) {
        const double pz = x[2], R2 = rh2 + pz * pz, R = sqrt(R2);
        nu[0] = z[0] - R;
        nu[1] = z[1] - atan2(px, py);
        nu[2] = z[2] - atan2(pz, rh);
        H[0] = px / R;                    H[1] = py / R;                    H[2] = pz / R;
        H[n + 0] = py / rh2;              H[n + 1] = -px / rh2;
        H[2 * n + 0] = -px * pz / (R2 * rh); H[2 * n + 1] = -py * pz / (R2 * rh); H[2 * n + 2] = rh / R2;
    } else {
        nu[0] = z[0] - rh;
        nu[1] = z[1] - atan2(px, py);
        H[0] = px / rh;      H[1] = py / rh;
        H[n + 0] = py / rh2; H[n + 1] = -px / rh2;
    }
    nu[1] = fmod(nu[1] + PI, 2.0 * PI);                   /* MATLAB mod: result in [0, 2*pi) */
    if (nu[1] < 0.0) nu[1] += 2.0 * PI;
    nu[1] -= PI;
}

/* Track predicted to time tp, innovation, Jacobian and S^-1. Returns 0 if S is singular. */
static int predict_measure(const Tracker *t, const TrackState *s, const double *z, double tp, double *x, double *P,
                           double *nu, double *H, double *Si, double *det)
{
    const int d = t->dim, n = t->n;
    double HP[3 * NMAX], S[9];
    memset(x, 0, sizeof(double) * NMAX);
    predict(t, s, tp - s->t, x, P);
    innovation(t, z, x, nu, H);
    mul(H, P, HP, d, n, n);
    mul_bt(HP, H, S, d, n, d);
    for (int i = 0; i < d; i++) S[i * d + i] += t->R[i];
    return inv_small(S, d, Si, det);
}

/* Measurement vector of a plot; 0 if it cannot be used (3-D without elevation). */
static int measurement(const Tracker *t, const PendingPlot *b, double *z)
{
    const double az = b->az_deg * PI / 180.0, el = b->el_deg * PI / 180.0;
    z[0] = z[1] = z[2] = 0.0;
    if (t->dim == 3) {
        if (!isfinite(el)) return 0;
        z[0] = b->range_m;  z[1] = az;  z[2] = el;
    } else {
        z[0] = isfinite(el) ? b->range_m * cos(el) : b->range_m;
        z[1] = az;
    }
    return 1;
}

/* ---- GNN assignment --------------------------------------------------------- */

/* Minimum-cost assignment of every row of A (n x m, n <= m, row-major) to its
 * own column (shortest augmenting paths with potentials, as hungarian() in
 * track_update.m). rowOf[j] = row (1-based) of column j, 0 = none. */
static void hungarian(const double *A, int n, int m, int *rowOf)
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

/* cost [nT x m] (INFINITY outside the gate) -> trackOf[i] = track of plot i, -1 = none.
 * Most finite pairs first, then the smallest total cost (as gnn_assign() in track_update.m). */
static void gnn_assign(const double *cost, int nT, int m, int *trackOf)
{
    for (int i = 0; i < m; i++) trackOf[i] = -1;
    int *rows = malloc(sizeof(int) * (nT + 1)), *cols = malloc(sizeof(int) * (m + 1)), nr = 0, nc = 0;
    if (!rows || !cols) abort();
    for (int j = 0; j < nT; j++)
        for (int i = 0; i < m; i++)
            if (isfinite(cost[(size_t)j * m + i])) { rows[nr++] = j; break; }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < nT; j++)
            if (isfinite(cost[(size_t)j * m + i])) { cols[nc++] = i; break; }
    if (nr > 0) {
        double lo = INFINITY, hi = -INFINITY;
        for (int r = 0; r < nr; r++)
            for (int c = 0; c < nc; c++) {
                const double v = cost[(size_t)rows[r] * m + cols[c]];
                if (isfinite(v)) { if (v < lo) lo = v; if (v > hi) hi = v; }
            }
        /* a forbidden pair costs more than any assignment with one more allowed pair */
        const double big = (hi - lo + 1.0) * ((nr < nc ? nr : nc) + 1);
        const int tr = nr > nc;                            /* rows must not outnumber the columns */
        const int n1 = tr ? nc : nr, m1 = tr ? nr : nc;
        double *C = malloc(sizeof(double) * (size_t)n1 * m1);
        int *rowOf = malloc(sizeof(int) * m1);
        if (!C || !rowOf) abort();
        for (int r = 0; r < nr; r++)
            for (int c = 0; c < nc; c++) {
                const double v = cost[(size_t)rows[r] * m + cols[c]];
                C[tr ? (size_t)c * m1 + r : (size_t)r * m1 + c] = isfinite(v) ? v - lo : big;
            }
        hungarian(C, n1, m1, rowOf);
        for (int k = 0; k < m1; k++) {
            if (rowOf[k] == 0) continue;
            const int r = tr ? k : rowOf[k] - 1, c = tr ? rowOf[k] - 1 : k;
            if (isfinite(cost[(size_t)rows[r] * m + cols[c]])) trackOf[cols[c]] = rows[r];
        }
        free(C);
        free(rowOf);
    }
    free(rows);
    free(cols);
}

/* ---- tracks ----------------------------------------------------------------- */

static void ekf_update(const Tracker *t, TrackState *s, const double *xb, const double *Pb, const double *nub,
                       const double *Hb, const double *Sib, double tp)
{
    const int d = t->dim, n = t->n;
    double PHt[NMAX * 3], K[NMAX * 3], IKH[NMAX * NMAX], T1[NMAX * NMAX], P1[NMAX * NMAX], KR[NMAX * 3], P2[NMAX * NMAX];
    mul_bt(Pb, Hb, PHt, n, n, d);                          /* P H' (n x d) */
    mul(PHt, Sib, K, n, d, d);                             /* K = P H' S^-1 */
    for (int i = 0; i < n; i++) {
        double dx = 0.0;
        for (int m = 0; m < d; m++) dx += K[i * d + m] * nub[m];
        s->x[i] = xb[i] + dx;
    }
    mul(K, Hb, IKH, n, d, n);                              /* Joseph form */
    for (int i = 0; i < n * n; i++) IKH[i] = -IKH[i];
    for (int i = 0; i < n; i++) IKH[i * n + i] += 1.0;
    mul(IKH, Pb, T1, n, n, n);
    mul_bt(T1, IKH, P1, n, n, n);
    for (int i = 0; i < n; i++)
        for (int m = 0; m < d; m++) KR[i * d + m] = K[i * d + m] * t->R[m];
    mul_bt(KR, K, P2, n, d, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) s->P[i * n + j] = P1[i * n + j] + P2[i * n + j];
    for (int i = 0; i < n; i++)                            /* symmetrise */
        for (int j = i + 1; j < n; j++) {
            const double v = 0.5 * (s->P[i * n + j] + s->P[j * n + i]);
            s->P[i * n + j] = s->P[j * n + i] = v;
        }
    s->t = tp;
    s->hits++;
    s->misses = 0;
    s->confirmed = s->confirmed || s->hits >= t->confirmHits;
}

/* New tentative track at the measured position; returns its index or -1. */
static int new_track(Tracker *t, const double *z, double tp)
{
    const int d = t->dim, n = t->n;
    if (t->nTracks == t->capacity) {
        const int cap = t->capacity ? 2 * t->capacity : 32;
        TrackState *p = realloc(t->tracks, sizeof(TrackState) * cap);
        if (!p) return -1;
        t->tracks = p;
        t->capacity = cap;
    }
    TrackState *s = &t->tracks[t->nTracks];
    memset(s, 0, sizeof *s);
    const double az = z[1];
    double J[9], JR[9], Ppos[9];
    if (d == 3) {
        const double R = z[0], el = z[2], ce = cos(el), se = sin(el), sa = sin(az), ca = cos(az);
        s->x[0] = R * ce * sa;  s->x[1] = R * ce * ca;  s->x[2] = R * se;
        const double Jm[9] = { ce * sa, R * ce * ca, -R * se * sa,
                               ce * ca, -R * ce * sa, -R * se * ca,
                               se, 0.0, R * ce };
        memcpy(J, Jm, sizeof Jm);
    } else {
        const double rho = z[0], sa = sin(az), ca = cos(az);
        s->x[0] = rho * sa;  s->x[1] = rho * ca;
        const double Jm[4] = { sa, rho * ca, ca, -rho * sa };
        memcpy(J, Jm, sizeof Jm);
    }
    for (int i = 0; i < d; i++)
        for (int m = 0; m < d; m++) JR[i * d + m] = J[i * d + m] * t->R[m];
    mul_bt(JR, J, Ppos, d, d, d);                          /* J R J' */
    for (int i = 0; i < d; i++)
        for (int m = 0; m < d; m++) s->P[i * n + m] = Ppos[i * d + m];
    for (int i = d; i < 2 * d; i++) s->P[i * n + i] = t->sigmaV0 * t->sigmaV0;
    for (int i = 2 * d; i < n; i++) s->P[i * n + i] = t->sigmaA0 * t->sigmaA0;
    s->id = t->nextId++;
    s->t = tp;
    s->hits = 1;
    s->confirmed = t->confirmHits <= 1;
    return t->nTracks++;
}

static void report_update(const Tracker *t, int j, int plotIndex, TrackerCallback callback, void *user)
{
    if (!callback) return;
    Track r;
    tracker_report(t, j, &r);
    callback(user, &r, plotIndex);
}

/* GNN association of the batch, then the updates and the new tracks. */
static void associate(Tracker *t, TrackerCallback callback, void *user)
{
    const int m = t->nPending, nT = t->nTracks, d = t->dim;
    double *Z = malloc(sizeof(double) * 3 * (m + 1)), *cost = malloc(sizeof(double) * ((size_t)nT * m + 1));
    int *usable = malloc(sizeof(int) * (m + 1)), *trackOf = malloc(sizeof(int) * (m + 1));
    if (!Z || !cost || !usable || !trackOf) abort();
    for (int i = 0; i < m; i++) usable[i] = measurement(t, &t->pending[i], Z + 3 * i);
    for (int j = 0; j < nT; j++)
        for (int i = 0; i < m; i++) {
            double x[NMAX], P[NMAX * NMAX], nu[3], H[3 * NMAX], Si[9], det, Sinu[3], d2 = 0.0;
            cost[(size_t)j * m + i] = INFINITY;
            if (!usable[i] || !predict_measure(t, &t->tracks[j], Z + 3 * i, t->pending[i].time_s, x, P, nu, H, Si, &det))
                continue;
            mul(Si, nu, Sinu, d, d, 1);
            for (int k = 0; k < d; k++) d2 += nu[k] * Sinu[k];
            if (d2 <= t->gate) cost[(size_t)j * m + i] = d2 + log(det);
        }
    gnn_assign(cost, nT, m, trackOf);
    for (int i = 0; i < m; i++) {
        if (!usable[i]) continue;
        const double tp = t->pending[i].time_s;
        int j = trackOf[i];
        if (j >= 0) {                                      /* EKF update */
            double x[NMAX], P[NMAX * NMAX], nu[3], H[3 * NMAX], Si[9], det;
            predict_measure(t, &t->tracks[j], Z + 3 * i, tp, x, P, nu, H, Si, &det);
            ekf_update(t, &t->tracks[j], x, P, nu, H, Si, tp);
        } else {
            j = new_track(t, Z + 3 * i, tp);
            if (j < 0) continue;
        }
        report_update(t, j, t->pending[i].plotIndex, callback, user);
    }
    t->nPending = 0;
    free(Z); free(cost); free(usable); free(trackOf);
}

void tracker_update(Tracker *t, const Detection *plots, int n, double t_s, int flush,
                    TrackerCallback callback, void *user)
{
    /* 1) the batch is complete */
    if (t->nPending > 0 && (t_s - t->pendingLast >= t->batchGap_s || t_s - t->pendingFirst >= t->batchMax_s))
        associate(t, callback, user);
    /* 2) the new plots join the batch */
    for (int k = 0; k < n; k++) {
        if (t->nPending == t->capPending) {
            const int cap = t->capPending ? 2 * t->capPending : 16;
            PendingPlot *p = realloc(t->pending, sizeof(PendingPlot) * cap);
            if (!p) return;
            t->pending = p;
            t->capPending = cap;
        }
        if (t->nPending == 0) t->pendingFirst = t_s;
        PendingPlot *b = &t->pending[t->nPending++];
        b->time_s = plots[k].time_s;
        b->az_deg = plots[k].az_deg;
        b->range_m = plots[k].range_m;
        b->el_deg = plots[k].el_deg;
        b->plotIndex = ++t->nPlots;
        t->pendingLast = t_s;
    }
    if (flush && t->nPending > 0) associate(t, callback, user);

    /* 3) missed scans and deletion */
    int keep = 0;
    for (int j = 0; j < t->nTracks; j++) {
        TrackState *s = &t->tracks[j];
        const double m = floor((t_s - s->t) / t->scanTime_s - 0.5);
        s->misses = m > 0.0 ? (int)m : 0;
        const int limit = s->confirmed ? t->maxMiss : t->maxMissTentative;
        if (s->misses > limit) {
            if (callback) {
                Track r;
                tracker_report(t, j, &r);
                r.deleted = 1;
                callback(user, &r, 0);
            }
            continue;
        }
        if (keep != j) t->tracks[keep] = *s;
        keep++;
    }
    t->nTracks = keep;
}

void tracker_report(const Tracker *t, int index, Track *out)
{
    const TrackState *s = &t->tracks[index];
    const int d = t->dim;
    memset(out, 0, sizeof *out);
    out->id = s->id;
    out->confirmed = s->confirmed;
    out->dim = d;
    out->time_s = s->t;
    for (int i = 0; i < d; i++) {
        out->pos_m[i] = s->x[i];
        out->vel_mps[i] = s->x[d + i];
        out->acc_mps2[i] = t->nd == 3 ? s->x[2 * d + i] : 0.0;
    }
    out->alt_m = d == 3 ? s->x[2] + t->height_m : NAN;
    out->hits = s->hits;
    out->misses = s->misses;
}

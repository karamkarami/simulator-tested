/*
 * tracker_tomht.c - track-oriented multiple hypothesis tracker (TOMHT) with
 * IMM-UKF state estimation (cfg.track.method = 'tomht'; as track_tomht.m,
 * which describes the algorithm; every step here follows it in the same order).
 *
 * Every track is a tree of track hypotheses (leaves: a history of plots or
 * missed looks, an IMM-UKF density and its score, the log-likelihood ratio
 * target / false plots), with a "no target" leaf. The global hypotheses G pick
 * one leaf per tree (no plot in two tracks); logw is the sum of their scores.
 * A batch of plots updates every global hypothesis with its k best assignments
 * (Murty's algorithm: plot to track, or a new tentative track); hypothesis
 * management: k-best (maxHyp global hypotheses, logHypPrune), N-scan pruning
 * (the decisions older than nScan scans are those of the best global
 * hypothesis), deletion by the score (deleteScore below its maximum), by the
 * missed scans and by the probability of the tree (pDelete). The state
 * estimation (IMM of CV, CA, CT, unscented transform), the assignment (Murty,
 * Hungarian) and the antenna looks are those of the PMBM (tracker_pmbm.c).
 *
 * The decision history of a leaf (time and plot, 0 = missed look, of every
 * look) is kept in a ring of the newest TOMHT_DEC decisions; older ones are
 * folded into a 64-bit hash (FNV-1a in order), so two histories before the
 * N-scan time compare by (hash, count) of that prefix, as the exact comparison
 * of track_tomht.m (the ring holds more than the nScan newest looks).
 */
#include "tracker_tomht.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define NX PMBM_N
#define NZ 4

static void *xmalloc(size_t size)
{
    void *p = malloc(size ? size : 1);
    if (!p) abort();
    return p;
}

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n ? n : 1, size);
    if (!p) abort();
    return p;
}

static void *xrealloc(void *p, size_t size)
{
    p = realloc(p, size ? size : 1);
    if (!p) abort();
    return p;
}

/* ---- small dense linear algebra (row-major) ----------------------------- */

/* Cholesky factor L (lower) of the n x n matrix A; 0 if A is not positive definite. */
static int chol_lower(const double *A, int n, double *L)
{
    memset(L, 0, sizeof(double) * n * n);
    for (int j = 0; j < n; j++) {
        double s = A[j * n + j];
        for (int k = 0; k < j; k++) s -= L[j * n + k] * L[j * n + k];
        if (!(s > 0.0)) return 0;
        const double d = sqrt(s);
        L[j * n + j] = d;
        for (int i = j + 1; i < n; i++) {
            double v = A[i * n + j];
            for (int k = 0; k < j; k++) v -= L[i * n + k] * L[j * n + k];
            L[i * n + j] = v / d;
        }
    }
    return 1;
}

/* x = A^-1 b with the Cholesky factor L of A. */
static void chol_solve(const double *L, int n, const double *b, double *x)
{
    double y[NX];
    for (int i = 0; i < n; i++) {
        double v = b[i];
        for (int k = 0; k < i; k++) v -= L[i * n + k] * y[k];
        y[i] = v / L[i * n + i];
    }
    for (int i = n - 1; i >= 0; i--) {
        double v = y[i];
        for (int k = i + 1; k < n; k++) v -= L[k * n + i] * x[k];
        x[i] = v / L[i * n + i];
    }
}

static void symmetrise(double *P, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++) {
            const double v = (P[i * n + j] + P[j * n + i]) / 2.0;
            P[i * n + j] = P[j * n + i] = v;
        }
}

static double wrap_pi(double a)
{
    double r = fmod(a + PI, 2.0 * PI);                 /* MATLAB mod: result in [0, 2*pi) */
    if (r < 0.0) r += 2.0 * PI;
    return r - PI;
}

/* 2n+1 sigma points of N(x, P), X[s * n + i] (Cholesky factor; a small
 * diagonal loading when P is only semi-definite, as track_tomht.m / track_pmbm.m). */
static void sigma_points(const TomhtTracker *t, const double *x, const double *P, int n, double *X)
{
    double L[NX * NX], A[NX * NX];
    if (!chol_lower(P, n, L)) {
        double mx = 0.0;
        for (int i = 0; i < n; i++) if (fabs(P[i * n + i]) > mx) mx = fabs(P[i * n + i]);
        double load = 1e-12 * mx;
        for (;;) {
            memcpy(A, P, sizeof(double) * n * n);
            for (int i = 0; i < n; i++) A[i * n + i] += load;
            if (chol_lower(A, n, L) || !(load > 0.0) || load > 1e300) break;
            load *= 10.0;
        }
    }
    for (int i = 0; i < n; i++) X[i] = x[i];
    for (int c = 0; c < n; c++)
        for (int i = 0; i < n; i++) {
            const double v = t->ukfScale * L[i * n + c];
            X[(1 + c) * n + i] = x[i] + v;
            X[(1 + n + c) * n + i] = x[i] - v;
        }
}

static double weight_m(const TomhtTracker *t, int s) { return s == 0 ? t->Wm0 : t->Wi; }
static double weight_c(const TomhtTracker *t, int s) { return s == 0 ? t->Wc0 : t->Wi; }

/* ---- IMM-UKF ------------------------------------------------------------- */

/* Combined estimate of the IMM (moment matching; P may be NULL). */
static void imm_estimate(const TomhtTracker *t, const TomhtSth *s, double *x, double *P)
{
    const int n = t->n;
    for (int i = 0; i < n; i++) {
        double v = 0.0;
        for (int j = 0; j < t->M; j++) v += s->x[j][i] * s->mu[j];
        x[i] = v;
    }
    if (!P) return;
    memset(P, 0, sizeof(double) * n * n);
    for (int j = 0; j < t->M; j++) {
        double dx[NX];
        for (int i = 0; i < n; i++) dx[i] = s->x[j][i] - x[i];
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++) P[a * n + b] += s->mu[j] * (s->P[j][a * n + b] + dx[a] * dx[b]);
    }
}

/* Coordinated turn in the horizontal plane (turn rate x[n-1]), constant vertical velocity, acceleration 0. */
static void f_ct(double *x, double dt, int d, int n)
{
    const double w = x[n - 1], vx = x[d], vy = x[d + 1], wt = w * dt;
    double sw, cw;
    if (fabs(wt) < 1e-6) { sw = dt - w * w * dt * dt * dt / 6.0;  cw = w * dt * dt / 2.0; }
    else { sw = sin(wt) / w;  cw = (1.0 - cos(wt)) / w; }
    x[0] += sw * vx - cw * vy;
    x[1] += cw * vx + sw * vy;
    x[d] = cos(wt) * vx - sin(wt) * vy;
    x[d + 1] = sin(wt) * vx + cos(wt) * vy;
    if (d == 3) x[2] += dt * x[5];
    for (int i = 2 * d; i < 3 * d; i++) x[i] = 0.0;
}

/* Q of the white-acceleration (nd = 2) or white-jerk (nd = 3) model per axis, added to P. */
static void add_q(double *P, int n, int d, int nd, double a, const double *q)
{
    double Q1[9];
    if (nd == 2) {
        const double v[4] = { a * a * a / 3.0, a * a / 2.0, a * a / 2.0, a };
        memcpy(Q1, v, sizeof v);
    } else {
        const double a2 = a * a, a3 = a2 * a, a4 = a3 * a, a5 = a4 * a;
        const double v[9] = { a5 / 20.0, a4 / 8.0, a3 / 6.0, a4 / 8.0, a3 / 3.0, a2 / 2.0, a3 / 6.0, a2 / 2.0, a };
        memcpy(Q1, v, sizeof v);
    }
    for (int bi = 0; bi < nd; bi++)
        for (int bj = 0; bj < nd; bj++)
            for (int i = 0; i < d; i++) P[(bi * d + i) * n + bj * d + i] += Q1[bi * nd + bj] * q[i];
}

/* Prediction of one model over dt (track_tomht.m / track_pmbm.m predict_model). */
static void predict_model(const TomhtTracker *t, const double *x0, const double *P0, double dt, int model,
                          double *x, double *P)
{
    const int d = t->dim, n = t->n, iw = n - 1;
    const double a = fabs(dt);
    if (model == PMBM_CT) {
        double X[(2 * NX + 1) * NX];
        sigma_points(t, x0, P0, n, X);
        for (int s = 0; s < 2 * n + 1; s++) f_ct(X + s * n, dt, d, n);
        for (int i = 0; i < n; i++) {
            double v = 0.0;
            for (int s = 0; s < 2 * n + 1; s++) v += weight_m(t, s) * X[s * n + i];
            x[i] = v;
        }
        memset(P, 0, sizeof(double) * n * n);
        for (int s = 0; s < 2 * n + 1; s++) {
            const double w = weight_c(t, s);
            for (int i = 0; i < n; i++)
                for (int k = 0; k < n; k++) P[i * n + k] += w * (X[s * n + i] - x[i]) * (X[s * n + k] - x[k]);
        }
        add_q(P, n, d, 2, a, t->qCt);
        for (int i = 2 * d; i < 3 * d; i++) P[i * n + i] += t->varUnusedA;
        P[iw * n + iw] += t->qTurn * a;
    } else {
        double F[NX * NX], FP[NX * NX];
        memset(F, 0, sizeof F);
        for (int i = 0; i < d; i++) {
            F[i * n + i] = 1.0;
            F[i * n + d + i] = dt;
            F[(d + i) * n + d + i] = 1.0;
            if (model == PMBM_CA) {
                F[i * n + 2 * d + i] = dt * dt / 2.0;
                F[(d + i) * n + 2 * d + i] = dt;
                F[(2 * d + i) * n + 2 * d + i] = 1.0;
            }
        }
        for (int i = 0; i < n; i++) {
            double v = 0.0;
            for (int k = 0; k < n; k++) v += F[i * n + k] * x0[k];
            x[i] = v;
        }
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++) {
                double v = 0.0;
                for (int m = 0; m < n; m++) v += F[i * n + m] * P0[m * n + k];
                FP[i * n + k] = v;
            }
        for (int i = 0; i < n; i++)
            for (int k = 0; k < n; k++) {
                double v = 0.0;
                for (int m = 0; m < n; m++) v += FP[i * n + m] * F[k * n + m];
                P[i * n + k] = v;
            }
        if (model == PMBM_CA) {
            add_q(P, n, d, 3, a, t->qCa);
        } else {
            add_q(P, n, d, 2, a, t->qCv);
            for (int i = 2 * d; i < 3 * d; i++) P[i * n + i] += t->varUnusedA;
        }
        P[iw * n + iw] += t->varUnusedW;
    }
    symmetrise(P, n);
}

/* Measurement of a state: 3-D [range; azimuth; elevation], 2-D [ground range;
 * azimuth], and the range rate (2-D: of the ground range) when withVr. */
static void h_meas(const double *x, int d, int withVr, double *z)
{
    const double px = x[0], py = x[1], rh = sqrt(px * px + py * py);
    double R;
    if (d == 3) {
        const double pz = x[2];
        R = sqrt(rh * rh + pz * pz);
        z[0] = R;  z[1] = atan2(px, py);  z[2] = atan2(pz, rh);
    } else {
        R = rh;
        z[0] = rh;  z[1] = atan2(px, py);
    }
    if (withVr) {
        double pv = 0.0;
        for (int i = 0; i < d; i++) pv += x[i] * x[d + i];
        z[d] = pv / R;
    }
}

/* Unscented transform of the measurement: mean zh, covariance S (+ Rm on the
 * diagonal), cross covariance C [n x nz]. */
static void ukf_measure(const TomhtTracker *t, const double *x, const double *P, const double *Rm, int withVr,
                        double *zh, double *S, double *C)
{
    const int n = t->n, d = t->dim, nz = d + withVr, ns = 2 * n + 1;
    double X[(2 * NX + 1) * NX], Zs[(2 * NX + 1) * NZ];
    sigma_points(t, x, P, n, X);
    for (int s = 0; s < ns; s++) h_meas(X + s * n, d, withVr, Zs + s * nz);
    for (int s = 0; s < ns; s++) Zs[s * nz + 1] = Zs[1] + wrap_pi(Zs[s * nz + 1] - Zs[1]);   /* around the centre */
    for (int k = 0; k < nz; k++) {
        double v = 0.0;
        for (int s = 0; s < ns; s++) v += weight_m(t, s) * Zs[s * nz + k];
        zh[k] = v;
    }
    memset(S, 0, sizeof(double) * nz * nz);
    memset(C, 0, sizeof(double) * n * nz);
    for (int s = 0; s < ns; s++) {
        const double w = weight_c(t, s);
        double dz[NZ];
        for (int k = 0; k < nz; k++) dz[k] = Zs[s * nz + k] - zh[k];
        for (int a = 0; a < nz; a++)
            for (int b = 0; b < nz; b++) S[a * nz + b] += w * dz[a] * dz[b];
        for (int i = 0; i < n; i++) {
            const double dx = X[s * n + i] - x[i];
            for (int k = 0; k < nz; k++) C[i * nz + k] += w * dx * dz[k];
        }
    }
    for (int k = 0; k < nz; k++) S[k * nz + k] += Rm[k];
    symmetrise(S, nz);
}

/* IMM-UKF update of s with the plot z (radial velocity zvr) of time tp (track_tomht.m / track_pmbm.m
 * imm_ukf_update): returns 1 inside the gate of a model, *logl the log likelihood,
 * *out the updated STH (r = 1). */
static int imm_ukf_update(const TomhtTracker *t, const TomhtSth *s, const double *z, double zvr, double span, double tp, int hasVr,
                          double *logl, TomhtSth *out)
{
    const int M = t->M, d = t->dim, n = t->n;
    const double dt = tp - s->t;
    double Pi[PMBM_M][PMBM_M], c[PMBM_M], W[PMBM_M][PMBM_M];
    const double lk = pow(t->switchLambda, fabs(dt) / t->scanTime_s);
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) Pi[i][j] = (1.0 - lk) / M + (i == j ? lk : 0.0);
    for (int j = 0; j < M; j++) {
        double v = 0.0;
        for (int i = 0; i < M; i++) v += s->mu[i] * Pi[i][j];
        c[j] = v;
    }
    for (int i = 0; i < M; i++)
        for (int j = 0; j < M; j++) W[i][j] = Pi[i][j] * s->mu[i] / (c[j] > DBL_MIN ? c[j] : DBL_MIN);
    *out = *s;
    double xp[PMBM_M][NX], Pp[PMBM_M][NX * NX];
    for (int j = 0; j < M; j++) {
        double x0[NX], P0[NX * NX];
        for (int a = 0; a < n; a++) {
            double v = 0.0;
            for (int i = 0; i < M; i++) v += s->x[i][a] * W[i][j];
            x0[a] = v;
        }
        memset(P0, 0, sizeof P0);
        for (int i = 0; i < M; i++) {
            double dx[NX];
            for (int a = 0; a < n; a++) dx[a] = s->x[i][a] - x0[a];
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++) P0[a * n + b] += W[i][j] * (s->P[i][a * n + b] + dx[a] * dx[b]);
        }
        symmetrise(P0, n);
        predict_model(t, x0, P0, dt, t->models[j], xp[j], Pp[j]);
    }
    /* measurement prediction of every model (with the radial velocity when it is in the plot) */
    const int nzFull = d + hasVr;
    double Zh[PMBM_M][NZ], S[PMBM_M][NZ * NZ], Cx[PMBM_M][NX * NZ], Rm[NZ];
    for (int k = 0; k < d; k++) Rm[k] = t->R[k];
    if (hasVr) Rm[d] = t->Rvr;
    for (int j = 0; j < M; j++) ukf_measure(t, xp[j], Pp[j], Rm, hasVr, Zh[j], S[j], Cx[j]);
    int useVr = hasVr, nz = nzFull;
    double zz[NZ];
    for (int k = 0; k < d; k++) zz[k] = z[k];
    if (hasVr) {
        double k = 0.0;
        if (!t->vrUnambiguous) {
            double sdMax = 0.0, hv = 0.0;
            for (int j = 0; j < M; j++) {
                const double v = S[j][nzFull * nzFull - 1] - t->Rvr, sd = sqrt(v > 0.0 ? v : 0.0);
                if (sd > sdMax) sdMax = sd;
                hv += Zh[j][d] * c[j];
            }
            useVr = sdMax < span / 6.0;
            k = round((hv - zvr) / span);
        }
        if (useVr) {
            zz[d] = zvr + k * span;
        } else {                                       /* drop the radial velocity row and column */
            nz = d;
            for (int j = 0; j < M; j++) {
                double S2[NZ * NZ], C2[NX * NZ];
                for (int a = 0; a < d; a++)
                    for (int b = 0; b < d; b++) S2[a * d + b] = S[j][a * nzFull + b];
                for (int a = 0; a < n; a++)
                    for (int b = 0; b < d; b++) C2[a * d + b] = Cx[j][a * nzFull + b];
                memcpy(S[j], S2, sizeof S2);
                memcpy(Cx[j], C2, sizeof C2);
            }
        }
    }
    const double gate = useVr ? t->gateVr : t->gate;
    double logLam[PMBM_M], nu[PMBM_M][NZ], Lc[PMBM_M][NZ * NZ], d2min = INFINITY;
    for (int j = 0; j < M; j++) {
        for (int k = 0; k < nz; k++) nu[j][k] = zz[k] - Zh[j][k];
        nu[j][1] = wrap_pi(nu[j][1]);
        double sv[NZ], logdet = 0.0, d2 = 0.0;
        if (!chol_lower(S[j], nz, Lc[j])) return 0;
        chol_solve(Lc[j], nz, nu[j], sv);
        for (int k = 0; k < nz; k++) { d2 += nu[j][k] * sv[k]; logdet += 2.0 * log(Lc[j][k * nz + k]); }
        logLam[j] = -0.5 * (d2 + nz * log(2.0 * PI) + logdet);
        if (d2 < d2min) d2min = d2;
    }
    double lw[PMBM_M], mx = -INFINITY;
    for (int j = 0; j < M; j++) {
        lw[j] = log(c[j] > DBL_MIN ? c[j] : DBL_MIN) + logLam[j];
        if (lw[j] > mx) mx = lw[j];
    }
    double sum = 0.0;
    for (int j = 0; j < M; j++) sum += exp(lw[j] - mx);
    *logl = isfinite(mx) ? mx + log(sum) : mx;
    if (hasVr && !useVr) *logl += t->logVrDensity;      /* the plot velocity is not explained: uniform */
    if (!(d2min <= gate)) return 0;
    for (int j = 0; j < M; j++) {                       /* UKF update of every model */
        double K[NX * NZ], KS[NX * NZ];
        for (int a = 0; a < n; a++) {                   /* K = C S^-1 (S symmetric): rows solve S k = c */
            double row[NZ], kr[NZ];
            for (int b = 0; b < nz; b++) row[b] = Cx[j][a * nz + b];
            chol_solve(Lc[j], nz, row, kr);
            for (int b = 0; b < nz; b++) K[a * nz + b] = kr[b];
        }
        for (int a = 0; a < n; a++) {
            double v = 0.0;
            for (int b = 0; b < nz; b++) v += K[a * nz + b] * nu[j][b];
            out->x[j][a] = xp[j][a] + v;
        }
        for (int a = 0; a < n; a++)
            for (int b = 0; b < nz; b++) {
                double v = 0.0;
                for (int m = 0; m < nz; m++) v += K[a * nz + m] * S[j][m * nz + b];
                KS[a * nz + b] = v;
            }
        for (int a = 0; a < n; a++)
            for (int b = 0; b < n; b++) {
                double v = 0.0;
                for (int m = 0; m < nz; m++) v += KS[a * nz + m] * K[b * nz + m];
                out->P[j][a * n + b] = Pp[j][a * n + b] - v;
            }
        symmetrise(out->P[j], n);
    }
    double msum = 0.0;
    for (int j = 0; j < M; j++) { out->mu[j] = exp(lw[j] - mx); msum += out->mu[j]; }
    for (int j = 0; j < M; j++) out->mu[j] /= msum;
    out->t = tp;
    out->tOpp = tp;
    out->r = 1.0;
    out->hits = s->hits + 1;
    out->misses = 0;
    return 1;
}

/* ---- births, looks, gates ------------------------------------------------- */

static int measurement(const TomhtTracker *t, const PendingPlot *b, double *z, double *zvr)
{
    const double az = b->az_deg * PI / 180.0, el = b->el_deg * PI / 180.0;
    z[0] = z[1] = z[2] = 0.0;
    *zvr = b->vr_mps;
    if (t->dim == 3) {
        z[0] = b->range_m;  z[1] = az;  z[2] = el;
        return isfinite(el);
    }
    z[0] = isfinite(el) ? b->range_m * cos(el) : b->range_m;
    z[1] = az;
    if (isfinite(el)) *zvr /= cos(el);
    return 1;
}

static double next_look(const TomhtTracker *t, const TomhtSth *s)
{
    const double T = t->scanTime_s;
    double xm[NX];
    imm_estimate(t, s, xm, NULL);
    const int d = t->dim;
    const double t0 = s->tOpp + T / 2.0, antAz = t->az0_deg + t->rate_dps * t0;
    double tn = t0 + T / 2.0;
    for (int it = 0; it < 4; it++) {
        const double qx = xm[0] + xm[d] * (tn - s->t), qy = xm[1] + xm[d + 1] * (tn - s->t);
        const double az = atan2(qx, qy) * 180.0 / PI;
        double dAz = fmod(az - antAz, 360.0);
        if (dAz < 0.0) dAz += 360.0;
        tn = t0 + dAz / t->rate_dps;
    }
    return tn;
}

static int coarse_gate(const TomhtTracker *t, const double *xm, const double *Pm, double ts, const double *z, double tp)
{
    const int d = t->dim, n = t->n;
    const double dt = fabs(tp - ts);
    const double qx = xm[0] + xm[d] * (tp - ts), qy = xm[1] + xm[d + 1] * (tp - ts);
    const double rh = d == 3 ? z[0] * cos(z[2]) : z[0];
    const double zx = rh * sin(z[1]), zy = rh * cos(z[1]);
    const double s2 = Pm[0] + Pm[n + 1] + dt * dt * (Pm[d * n + d] + Pm[(d + 1) * n + d + 1]) +
                      dt * dt * dt * dt / 4.0 * (Pm[2 * d * n + 2 * d] + Pm[(2 * d + 1) * n + 2 * d + 1]);
    const double sm = sqrt(t->R[0]) + rh * sqrt(t->R[1]);
    return hypot(zx - qx, zy - qy) <= 6.0 * sqrt(s2) + 25.0 * dt * dt + 6.0 * sm + 1000.0;
}

static TomhtSth empty_sth(void)
{
    TomhtSth s;
    memset(&s, 0, sizeof s);
    s.t = s.tOpp = NAN;
    s.tNext = INFINITY;
    return s;
}

static TomhtSth new_sth(const TomhtTracker *t, const double *z, double zvr, const PendingPlot *b, double r)
{
    const int d = t->dim, n = t->n;
    double x[NX] = { 0 }, P[NX * NX] = { 0 }, J[9], pos[3];
    const double az = z[1];
    if (d == 3) {
        const double R = z[0], el = z[2], ce = cos(el), se = sin(el), sa = sin(az), ca = cos(az);
        pos[0] = R * ce * sa;  pos[1] = R * ce * ca;  pos[2] = R * se;
        const double Jm[9] = { ce * sa, R * ce * ca, -R * se * sa, ce * ca, -R * ce * sa, -R * se * ca, se, 0.0, R * ce };
        memcpy(J, Jm, sizeof Jm);
    } else {
        const double rho = z[0], sa = sin(az), ca = cos(az);
        pos[0] = rho * sa;  pos[1] = rho * ca;
        const double Jm[4] = { sa, rho * ca, ca, -rho * sa };
        memcpy(J, Jm, sizeof Jm);
    }
    for (int i = 0; i < d; i++) x[i] = pos[i];
    for (int i = 0; i < d; i++)                         /* J R J' */
        for (int k = 0; k < d; k++) {
            double v = 0.0;
            for (int m = 0; m < d; m++) v += J[i * d + m] * t->R[m] * J[k * d + m];
            P[i * n + k] = v;
        }
    double D[9] = { 0 };
    for (int i = 0; i < d; i++) D[i * d + i] = t->sigmaV0 * t->sigmaV0;
    if (d == 3) D[8] = t->sigmaClimb0 * t->sigmaClimb0;
    for (int i = 0; i < d; i++)
        for (int k = 0; k < d; k++) P[(d + i) * n + d + k] = D[i * d + k];
    for (int i = 0; i < d; i++) P[(2 * d + i) * n + 2 * d + i] = t->sigmaA0 * t->sigmaA0;
    P[(n - 1) * n + n - 1] = t->sigmaTurn0 * t->sigmaTurn0;
    if (t->useVr && t->vrUnambiguous && isfinite(zvr)) {   /* the prior conditioned on the radial velocity */
        double r2 = 0.0, u[3], Du[3], s2 = t->Rvr;
        for (int i = 0; i < d; i++) r2 += pos[i] * pos[i];
        for (int i = 0; i < d; i++) u[i] = pos[i] / sqrt(r2);
        for (int i = 0; i < d; i++) {
            double v = 0.0;
            for (int k = 0; k < d; k++) v += D[i * d + k] * u[k];
            Du[i] = v;
        }
        double uDu = 0.0;
        for (int i = 0; i < d; i++) uDu += u[i] * Du[i];
        s2 = uDu + t->Rvr;
        for (int i = 0; i < d; i++) x[d + i] = Du[i] * zvr / s2;
        double Pv[9];
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++) Pv[i * d + k] = D[i * d + k] - Du[i] * Du[k] / s2;
        for (int i = 0; i < d; i++)
            for (int k = 0; k < d; k++) P[(d + i) * n + d + k] = (Pv[i * d + k] + Pv[k * d + i]) / 2.0;
    }
    TomhtSth s = empty_sth();
    s.r = r;
    for (int j = 0; j < t->M; j++) {
        memcpy(s.x[j], x, sizeof(double) * n);
        memcpy(s.P[j], P, sizeof(double) * n * n);
        s.mu[j] = t->mu0[j];
    }
    s.t = s.tOpp = b->time_s;
    s.hits = 1;
    s.plot = b->plotIndex;
    s.hist[0] = b->plotIndex;
    s.nHist = 1;
    return s;
}

/* ---- assignment: Murty's k best ----------------------------------------------- */

/* Minimum-cost assignment of every row of C (n x m, n <= m, INFINITY =
 * forbidden) to its own column: a[i] = column of row i; returns 0 when no
 * assignment avoids the forbidden pairs (track_tomht.m / track_pmbm.m lap). */
static int lap(const double *C, int n, int m, int *a, double *cost)
{
    if (n == 0) { *cost = 0.0; return 1; }
    double lo = INFINITY, hi = -INFINITY;
    for (int i = 0; i < n; i++) {
        int any = 0;
        for (int j = 0; j < m; j++) {
            const double v = C[i * m + j];
            if (isfinite(v)) { any = 1; if (v < lo) lo = v; if (v > hi) hi = v; }
        }
        if (!any) return 0;
    }
    double *A = xmalloc(sizeof(double) * n * m);
    int *rowOf = xmalloc(sizeof(int) * m);
    const double big = (hi - lo + 1.0) * (n + 1);
    for (int k = 0; k < n * m; k++) A[k] = isfinite(C[k]) ? C[k] - lo : big;
    tracker_hungarian(A, n, m, rowOf);
    for (int i = 0; i < n; i++) a[i] = -1;
    for (int j = 0; j < m; j++) if (rowOf[j] > 0) a[rowOf[j] - 1] = j;
    free(A);
    free(rowOf);
    double c = 0.0;
    for (int i = 0; i < n; i++) {
        if (a[i] < 0 || !isfinite(C[i * m + a[i]])) return 0;
        c += C[i * m + a[i]];
    }
    *cost = c;
    return 1;
}

typedef struct { double *C; int *a; double c; } MurtyNode;

/* The k best assignments (track_tomht.m / track_pmbm.m murty): A[r * n + i] = column of row i
 * in the r-th best, cost[r]; returns their number (A, cost allocated). */
static int murty(const double *C, int n, int m, int k, int **Aout, double **costOut)
{
    int nOut = 0, capOut = 8;
    int *A = xmalloc(sizeof(int) * capOut * (n ? n : 1));
    double *cost = xmalloc(sizeof(double) * capOut);
    MurtyNode *Q = NULL;
    int nQ = 0, capQ = 0;
    {
        MurtyNode first = { xmalloc(sizeof(double) * n * m), xmalloc(sizeof(int) * (n ? n : 1)), 0.0 };
        memcpy(first.C, C, sizeof(double) * n * m);
        if (lap(first.C, n, m, first.a, &first.c)) {
            Q = xrealloc(Q, sizeof(MurtyNode) * (capQ = 16));
            Q[nQ++] = first;
        } else {
            free(first.C);
            free(first.a);
        }
    }
    double *Cw = xmalloc(sizeof(double) * n * m);
    while (nQ > 0 && nOut < k) {
        int b = 0;
        for (int q = 1; q < nQ; q++) if (Q[q].c < Q[b].c) b = q;
        MurtyNode node = Q[b];
        memmove(Q + b, Q + b + 1, sizeof(MurtyNode) * (nQ - b - 1));
        nQ--;
        if (nOut == capOut) {
            capOut *= 2;
            A = xrealloc(A, sizeof(int) * capOut * (n ? n : 1));
            cost = xrealloc(cost, sizeof(double) * capOut);
        }
        memcpy(A + nOut * n, node.a, sizeof(int) * n);
        cost[nOut++] = node.c;
        memcpy(Cw, node.C, sizeof(double) * n * m);
        for (int r = 0; r < n; r++) {
            MurtyNode sub = { xmalloc(sizeof(double) * n * m), xmalloc(sizeof(int) * n), 0.0 };
            memcpy(sub.C, Cw, sizeof(double) * n * m);
            sub.C[r * m + node.a[r]] = INFINITY;        /* without the pair (r, a(r)) */
            if (lap(sub.C, n, m, sub.a, &sub.c)) {
                if (nQ == capQ) Q = xrealloc(Q, sizeof(MurtyNode) * (capQ = capQ ? 2 * capQ : 16));
                Q[nQ++] = sub;
            } else {
                free(sub.C);
                free(sub.a);
            }
            const double keep = Cw[r * m + node.a[r]];   /* following subproblems keep (r, a(r)) */
            for (int j = 0; j < m; j++) Cw[r * m + j] = INFINITY;
            for (int i = 0; i < n; i++) Cw[i * m + node.a[r]] = INFINITY;
            Cw[r * m + node.a[r]] = keep;
        }
        free(node.C);
        free(node.a);
    }
    for (int q = 0; q < nQ; q++) { free(Q[q].C); free(Q[q].a); }
    free(Q);
    free(Cw);
    *Aout = A;
    *costOut = cost;
    return nOut;
}


/* ---- init ---------------------------------------------------------------- */

void tomht_tracker_init(TomhtTracker *t, const RadarConfig *cfg)
{
    memset(t, 0, sizeof *t);
    t->dim = cfg->trackDim == 3 ? 3 : 2;
    t->n = 3 * t->dim + 1;
    t->M = cfg->pmbmNModels;
    for (int j = 0; j < t->M; j++) {
        t->models[j] = cfg->pmbmModels[j];
        t->mu0[j] = cfg->pmbmMu0[j];
    }
    for (int i = 0; i < 3; i++) {
        t->R[i] = cfg->trackSigma[i] * cfg->trackSigma[i];
        t->qCv[i] = cfg->pmbmQCv[i];
        t->qCa[i] = cfg->pmbmQCa[i];
        t->qCt[i] = cfg->pmbmQCt[i];
    }
    t->Rvr = cfg->trackSigmaVr * cfg->trackSigmaVr;
    t->vrSpan = cfg->trackVrSpan;
    for (int p = 0; p < N_PULSES; p++) {
        Detection d;
        memset(&d, 0, sizeof d);
        d.pulse = p;
        t->vrSpanPulse[p] = tracker_plot_vr_span(cfg, &d);
    }
    t->sigmaV0 = cfg->trackSigmaV0;
    t->sigmaA0 = cfg->trackSigmaA0;
    t->gate = cfg->trackGate;
    t->gateVr = cfg->trackGateVr;
    t->useVr = cfg->trackUseVr;
    t->vrUnambiguous = cfg->trackVrUnambiguous;
    t->scanTime_s = cfg->scanTime_s;
    t->height_m = cfg->height_m;
    t->az0_deg = cfg->trackAz0_deg;
    t->rate_dps = cfg->rate_dps;
    t->batchGap_s = cfg->trackBatchGap_s;
    t->batchMax_s = cfg->trackBatchMax_s;
    t->confirmHits = cfg->trackConfirmHits;
    t->maxMissTentative = cfg->trackMaxMissTentative;
    t->maxMiss = cfg->trackMaxMiss;
    t->minSpeed = cfg->trackMinSpeed;
    t->minSpeedHits = 2;
    t->pd = cfg->tomhtPd;
    t->logLambdaC = cfg->pmbmLogLambdaC;
    t->logLambdaU = cfg->pmbmLogLambdaU;
    t->logVrDensity = cfg->pmbmLogVrDensity;
    t->logHypPrune = cfg->tomhtLogHypPrune;
    t->maxHyp = cfg->tomhtMaxHyp;
    t->kBest = cfg->tomhtKBest;
    t->nScan = cfg->tomhtNScan;
    t->confirmScore = cfg->tomhtConfirmScore;
    t->deleteScore = cfg->tomhtDeleteScore;
    t->pDelete = cfg->tomhtPDelete;
    t->switchLambda = cfg->pmbmSwitchLambda;
    t->qTurn = cfg->pmbmQTurn;
    t->sigmaTurn0 = cfg->pmbmSigmaTurn0;
    t->sigmaClimb0 = cfg->pmbmSigmaClimb0;
    t->varUnusedA = cfg->pmbmVarUnusedA;
    t->varUnusedW = cfg->pmbmVarUnusedW;
    t->ukfScale = cfg->pmbmUkfScale;
    t->Wm0 = cfg->pmbmWm0;
    t->Wc0 = cfg->pmbmWc0;
    t->Wi = cfg->pmbmWi;
    t->latency_s = cfg->pmbmLatency_s;
    t->maintStep_s = cfg->pmbmMaintStep_s;
    t->nextId = 1;
    t->tMaint = -INFINITY;
    t->tNow = -INFINITY;
    t->pendingFirst = t->pendingLast = -INFINITY;
    t->nG = 1;                                        /* one global hypothesis: no track */
    t->G = xcalloc(1, sizeof(int));
    t->logw = xcalloc(1, sizeof(double));
}

void tomht_tracker_free(TomhtTracker *t)
{
    for (int i = 0; i < t->nT; i++) free(t->T[i].h);
    free(t->T); free(t->G); free(t->logw); free(t->pending);
    free(t->confId); free(t->confPlot); free(t->shown);
    memset(t, 0, sizeof *t);
}


/* ---- decision history (N-scan) -------------------------------------------- */

static uint64_t fnv_mix(uint64_t h, const void *data, size_t n)
{
    const unsigned char *p = data;
    for (size_t i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211ULL; }
    return h;
}

static uint64_t fnv_row(uint64_t h, double time, double plot)
{
    h = fnv_mix(h, &time, sizeof time);
    return fnv_mix(h, &plot, sizeof plot);
}

#define FNV_START 14695981039346656037ULL

static void dec_push(TomhtSth *s, double time, int plot)
{
    if (s->nDec == TOMHT_DEC) {                          /* the oldest decision into the hash */
        s->decHash = fnv_row(s->decHash, s->decT[0], (double)s->decP[0]);
        s->nDecOld++;
        memmove(s->decT, s->decT + 1, sizeof(double) * (TOMHT_DEC - 1));
        memmove(s->decP, s->decP + 1, sizeof(int) * (TOMHT_DEC - 1));
        s->nDec--;
    }
    s->decT[s->nDec] = time;
    s->decP[s->nDec++] = plot;
}

typedef struct { uint64_t hash; int count; } Prefix;

/* The decisions before tCut (track_tomht.m prefix): the "no target" leaf of a
 * tree born before tCut is the row [-1 -1], a deleted track ends with [-2 -2]. */
static Prefix prefix(const TomhtSth *s, double t0, double tCut)
{
    Prefix p = { FNV_START, 0 };
    const int total = s->nDecOld + s->nDec;
    if (!(s->r > 0.0) && total == 0) {
        if (t0 < tCut) { p.hash = fnv_row(p.hash, -1.0, -1.0); p.count = 1; }
        return p;
    }
    p.hash = s->nDecOld ? s->decHash : FNV_START;
    p.count = s->nDecOld;
    for (int k = 0; k < s->nDec; k++)
        if (s->decT[k] < tCut) { p.hash = fnv_row(p.hash, s->decT[k], (double)s->decP[k]); p.count++; }
    if (!(s->r > 0.0)) { p.hash = fnv_row(p.hash, -2.0, -2.0); p.count++; }
    return p;
}

/* No plot at the look of the track hypothesis: the score falls by log(1 - pd);
 * returns the change of the hypothesis weight (a deleted track: no target, 0). */
static double miss_update(const TomhtTracker *t, TomhtSth *s, int confirmed)
{
    double dl = log(1.0 - t->pd);
    s->score += dl;
    dec_push(s, s->tNext, 0);
    s->tOpp = s->tNext;
    s->misses++;
    const int limit = confirmed ? t->maxMiss : t->maxMissTentative;
    if (s->misses > limit || s->score < s->maxScore - t->deleteScore) {
        s->r = 0.0;
        dl -= s->score;
    }
    s->tNext = s->r > 0.0 ? next_look(t, s) : INFINITY;
    return dl;
}

/* ---- pruning --------------------------------------------------------------- */

static int G_at(const TomhtTracker *t, int g, int i) { return t->G[(size_t)g * t->nT + i]; }

static void remove_tracks(TomhtTracker *t, const int *del)
{
    int nNew = 0;
    for (int i = 0; i < t->nT; i++) if (!del[i]) nNew++;
    int *G = xmalloc(sizeof(int) * (size_t)t->nG * (nNew ? nNew : 1));
    for (int g = 0; g < t->nG; g++) {
        int c = 0;
        for (int i = 0; i < t->nT; i++) if (!del[i]) G[(size_t)g * nNew + c++] = G_at(t, g, i);
    }
    int c = 0;
    for (int i = 0; i < t->nT; i++) {
        if (del[i]) { free(t->T[i].h); continue; }
        t->T[c++] = t->T[i];
    }
    free(t->G);
    t->G = G;
    t->nT = nNew;
}

static const TomhtTracker *g_sortTracker;
static int cmp_rows(const void *pa, const void *pb)
{
    const int a = *(const int *)pa, b = *(const int *)pb;
    const TomhtTracker *t = g_sortTracker;
    for (int i = 0; i < t->nT; i++) {
        const int x = G_at(t, a, i), y = G_at(t, b, i);
        if (x != y) return x < y ? -1 : 1;
    }
    return a < b ? -1 : a > b;
}

/* Equal global hypotheses are merged, the best score kept (MATLAB unique(G, 'rows')). */
static void merge_rows(TomhtTracker *t)
{
    if (t->nT == 0) {
        t->nG = 1;
        t->logw[0] = 0.0;
        return;
    }
    int *ord = xmalloc(sizeof(int) * t->nG);
    for (int g = 0; g < t->nG; g++) ord[g] = g;
    g_sortTracker = t;
    qsort(ord, t->nG, sizeof(int), cmp_rows);
    int *grp = xmalloc(sizeof(int) * t->nG), nU = 0;
    for (int k = 0; k < t->nG; k++) {
        if (k > 0 && memcmp(t->G + (size_t)ord[k] * t->nT, t->G + (size_t)ord[k - 1] * t->nT, sizeof(int) * t->nT) == 0)
            grp[ord[k]] = nU - 1;
        else
            grp[ord[k]] = nU++;
    }
    if (nU == t->nG) { free(ord); free(grp); return; }
    int *G = xmalloc(sizeof(int) * (size_t)nU * t->nT);
    double *lw = xmalloc(sizeof(double) * nU);
    for (int u = 0; u < nU; u++) lw[u] = -INFINITY;
    for (int g = 0; g < t->nG; g++) {
        memcpy(G + (size_t)grp[g] * t->nT, t->G + (size_t)g * t->nT, sizeof(int) * t->nT);
        if (t->logw[g] > lw[grp[g]]) lw[grp[g]] = t->logw[g];
    }
    free(t->G); free(t->logw);
    t->G = G;  t->logw = lw;  t->nG = nU;
    free(ord); free(grp);
}

static const double *g_sortW;
static int cmp_desc(const void *pa, const void *pb)
{
    const int a = *(const int *)pa, b = *(const int *)pb;
    if (g_sortW[a] != g_sortW[b]) return g_sortW[a] > g_sortW[b] ? -1 : 1;
    return a < b ? -1 : a > b;
}

static int cmp_plot(const void *pa, const void *pb)
{
    const int a = *(const int *)pa, b = *(const int *)pb;
    const int x = g_sortTracker->T[a].plot, y = g_sortTracker->T[b].plot;
    if (x != y) return x < y ? -1 : 1;
    return a < b ? -1 : a > b;
}

static int best_hyp(const TomhtTracker *t)
{
    int best = 0;
    for (int g = 1; g < t->nG; g++) if (t->logw[g] > t->logw[best]) best = g;
    return best;
}

static double log_sum_w(const TomhtTracker *t)
{
    double mx = -INFINITY, s = 0.0;
    for (int g = 0; g < t->nG; g++) if (t->logw[g] > mx) mx = t->logw[g];
    if (!isfinite(mx)) return mx;
    for (int g = 0; g < t->nG; g++) s += exp(t->logw[g] - mx);
    return mx + log(s);
}

/* N-scan pruning (track_tomht.m nscan_prune). */
static void nscan_prune(TomhtTracker *t)
{
    if (t->nT == 0 || t->nG == 0) return;
    const int best = best_hyp(t);
    const double tCut = t->tNow - (t->nScan - 0.5) * t->scanTime_s;
    char *keep = xmalloc(t->nG);
    memset(keep, 1, t->nG);
    for (int i = 0; i < t->nT; i++) {
        const TomhtTrack *tk = &t->T[i];
        const int hb = G_at(t, best, i);
        const Prefix ref = prefix(&tk->h[hb], tk->t0, tCut);
        char *done = xcalloc(tk->nh, 1);
        for (int g = 0; g < t->nG; g++) {
            const int h = G_at(t, g, i);
            if (h == hb || done[h]) continue;
            done[h] = 1;
            const Prefix q = prefix(&tk->h[h], tk->t0, tCut);
            if (q.hash != ref.hash || q.count != ref.count)
                for (int g2 = 0; g2 < t->nG; g2++) if (G_at(t, g2, i) == h) keep[g2] = 0;
        }
        free(done);
    }
    int nK = 0;
    for (int g = 0; g < t->nG; g++) {
        if (!keep[g]) continue;
        if (nK != g) {
            memmove(t->G + (size_t)nK * t->nT, t->G + (size_t)g * t->nT, sizeof(int) * t->nT);
            t->logw[nK] = t->logw[g];
        }
        nK++;
    }
    t->nG = nK;
    free(keep);
}

/* Pruning (k-best, N-scan, deletion), track numbers, confirmation (track_tomht.m cleanup). */
static void cleanup(TomhtTracker *t)
{
    const int d = t->dim;
    for (int i = 0; i < t->nT; i++) {
        TomhtTrack *tk = &t->T[i];
        int zero0 = -1;
        for (int h = 0; h < tk->nh; h++) {
            TomhtSth *s = &tk->h[h];
            if (s->r > 0.0 && t->minSpeed > 0.0 && s->hits >= t->minSpeedHits) {
                double xm[NX], v2 = 0.0;
                imm_estimate(t, s, xm, NULL);
                for (int k = d; k < 2 * d; k++) v2 += xm[k] * xm[k];
                if (sqrt(v2) < t->minSpeed) {
                    for (int g = 0; g < t->nG; g++) if (G_at(t, g, i) == h) t->logw[g] -= s->score;
                    s->r = 0.0;
                }
            }
        }
        int nZero = 0;
        for (int h = 0; h < tk->nh; h++) if (tk->h[h].r <= 0.0) { if (zero0 < 0) zero0 = h; nZero++; }
        if (nZero > 1)
            for (int g = 0; g < t->nG; g++) {
                int *v = &t->G[(size_t)g * t->nT + i];
                if (tk->h[*v].r <= 0.0) *v = zero0;
            }
    }
    for (int pass = 1; pass <= 2; pass++) {
        merge_rows(t);
        /* global hypotheses: weight and number (k-best) */
        double mx = -INFINITY;
        for (int g = 0; g < t->nG; g++) if (t->logw[g] > mx) mx = t->logw[g];
        int *ord = xmalloc(sizeof(int) * t->nG), nK = 0;
        for (int g = 0; g < t->nG; g++) if (t->logw[g] >= mx + t->logHypPrune) ord[nK++] = g;
        g_sortW = t->logw;
        qsort(ord, nK, sizeof(int), cmp_desc);
        if (nK > t->maxHyp) nK = t->maxHyp;
        int *G = xmalloc(sizeof(int) * (size_t)nK * (t->nT ? t->nT : 1));
        double *lw = xmalloc(sizeof(double) * (nK ? nK : 1));
        for (int k = 0; k < nK; k++) {
            memcpy(G + (size_t)k * t->nT, t->G + (size_t)ord[k] * t->nT, sizeof(int) * t->nT);
            lw[k] = t->logw[ord[k]];
        }
        free(ord); free(t->G); free(t->logw);
        t->G = G;  t->logw = lw;  t->nG = nK;
        nscan_prune(t);
        /* trees: probability over the global hypotheses */
        const double lse = log_sum_w(t);
        int *del = xcalloc(t->nT ? t->nT : 1, sizeof(int)), any = 0;
        for (int i = 0; i < t->nT; i++) {
            double rm = 0.0;
            for (int g = 0; g < t->nG; g++) rm += exp(t->logw[g] - lse) * (t->T[i].h[G_at(t, g, i)].r > 0.0);
            t->T[i].rMarg = rm;
            del[i] = rm < t->pDelete;
            any |= del[i];
        }
        if (!any && pass > 1) { free(del); break; }
        for (int i = 0; i < t->nT; i++) {               /* the rows holding a deleted track: its plots are false */
            if (!del[i]) continue;
            for (int g = 0; g < t->nG; g++) {
                const TomhtSth *s = &t->T[i].h[G_at(t, g, i)];
                if (s->r > 0.0) t->logw[g] -= s->score;
            }
        }
        remove_tracks(t, del);
        free(del);
    }
    /* unused leaves are removed */
    for (int i = 0; i < t->nT; i++) {
        TomhtTrack *tk = &t->T[i];
        int *map = xmalloc(sizeof(int) * tk->nh), nUsed = 0;
        for (int h = 0; h < tk->nh; h++) map[h] = -1;
        for (int g = 0; g < t->nG; g++) map[G_at(t, g, i)] = 1;
        for (int h = 0; h < tk->nh; h++) if (map[h] > 0) { map[h] = nUsed; tk->h[nUsed++] = tk->h[h]; }
        tk->nh = nUsed;
        for (int g = 0; g < t->nG; g++) t->G[(size_t)g * t->nT + i] = map[G_at(t, g, i)];
        free(map);
    }
    /* probability that a target made a plot (over the global hypotheses) */
    const double lse = log_sum_w(t);
    int nPl = 0, capPl = 16;
    int *pl = xmalloc(sizeof(int) * capPl);
    double *pw = xmalloc(sizeof(double) * capPl);
    for (int i = 0; i < t->nT; i++)
        for (int g = 0; g < t->nG; g++) {
            const TomhtSth *s = &t->T[i].h[G_at(t, g, i)];
            if (!(s->r > 0.0)) continue;
            int k = 0;
            while (k < nPl && pl[k] != s->plot) k++;
            if (k == nPl) {
                if (nPl == capPl) {
                    capPl *= 2;
                    pl = xrealloc(pl, sizeof(int) * capPl);
                    pw = xrealloc(pw, sizeof(double) * capPl);
                }
                pl[nPl] = s->plot;
                pw[nPl++] = 0.0;
            }
            pw[k] += exp(t->logw[g] - lse);
        }
    /* track numbers (in the order of their first plot) and confirmation */
    const int best = best_hyp(t);
    int *ord = xmalloc(sizeof(int) * (t->nT ? t->nT : 1));
    for (int i = 0; i < t->nT; i++) ord[i] = i;
    g_sortTracker = t;
    qsort(ord, t->nT, sizeof(int), cmp_plot);
    for (int k = 0; k < t->nT; k++) {
        TomhtTrack *tk = &t->T[ord[k]];
        if (tk->id == 0) tk->id = t->nextId++;
        tk->rTarget = 0.0;
        if (t->nG < 1) continue;
        const TomhtSth *s = &tk->h[G_at(t, best, ord[k])];
        if (!(s->r > 0.0)) continue;
        for (int q = 0; q < nPl; q++) if (pl[q] == s->plot) tk->rTarget = pw[q] < 1.0 ? pw[q] : 1.0;
        if (s->score >= t->confirmScore && s->hits >= t->confirmHits) tk->confirmed = 1;
    }
    free(ord); free(pl); free(pw);
}

/* ---- association ------------------------------------------------------------- */

static int cmp_time(const void *pa, const void *pb)
{
    const PendingPlot *a = pa, *b = pb;
    if (a->time_s != b->time_s) return a->time_s < b->time_s ? -1 : 1;
    return a->plotIndex < b->plotIndex ? -1 : a->plotIndex > b->plotIndex;
}

static void push_sth(TomhtTrack *tk, const TomhtSth *s)
{
    if (tk->nh == tk->cap) {
        tk->cap = tk->cap ? 2 * tk->cap : 4;
        tk->h = xrealloc(tk->h, sizeof(TomhtSth) * tk->cap);
    }
    tk->h[tk->nh++] = *s;
}

/* Update of the trees and the global hypotheses with the plots of the batch (track_tomht.m associate). */
static void associate(TomhtTracker *t)
{
    const int nT = t->nT;
    const double T = t->scanTime_s;
    PendingPlot *B = xmalloc(sizeof(PendingPlot) * (t->nPending ? t->nPending : 1));
    double *Z = xmalloc(sizeof(double) * 3 * (t->nPending + 1)), *Zvr = xmalloc(sizeof(double) * (t->nPending + 1));
    int nPl = 0;
    t->batch++;
    for (int k = 0; k < t->nPending; k++) {
        double z[3], zvr;
        if (!measurement(t, &t->pending[k], z, &zvr)) continue;
        B[nPl++] = t->pending[k];
    }
    t->nPending = 0;
    if (nPl == 0) { free(B); free(Z); free(Zvr); return; }
    qsort(B, nPl, sizeof(PendingPlot), cmp_time);
    double tMax = -INFINITY;
    for (int j = 0; j < nPl; j++) {
        measurement(t, &B[j], Z + 3 * j, &Zvr[j]);
        if (B[j].time_s > tMax) tMax = B[j].time_s;
    }
    const double tMissLimit = tMax + t->batchGap_s;
    int *hasVr = xmalloc(sizeof(int) * nPl);
    double *logC = xmalloc(sizeof(double) * nPl), *logNew = xmalloc(sizeof(double) * nPl);
    for (int j = 0; j < nPl; j++) {
        hasVr[j] = t->useVr && isfinite(Zvr[j]);
        logC[j] = t->logLambdaC + hasVr[j] * t->logVrDensity;
        logNew[j] = log(t->pd) + t->logLambdaU - t->logLambdaC;   /* (the vr density cancels) */
    }
    /* track hypotheses of the trees for the plots: score increments and updated densities */
    double **L = xcalloc(nT ? nT : 1, sizeof(double *)), **L0 = xcalloc(nT ? nT : 1, sizeof(double *));
    int **upd = xcalloc(nT ? nT : 1, sizeof(int *)), **miss = xcalloc(nT ? nT : 1, sizeof(int *));
    TomhtSth *pool = NULL;
    int nPool = 0, capPool = 0;
    for (int i = 0; i < nT; i++) {
        const TomhtTrack *tk = &t->T[i];
        const int nh = tk->nh;
        L[i] = xmalloc(sizeof(double) * nh * nPl);
        L0[i] = xcalloc(nh, sizeof(double));
        upd[i] = xmalloc(sizeof(int) * nh * nPl);
        miss[i] = xmalloc(sizeof(int) * nh);
        for (int k = 0; k < nh * nPl; k++) { L[i][k] = -INFINITY; upd[i][k] = -1; }
        char *used = xcalloc(nh, 1);
        for (int g = 0; g < t->nG; g++) used[G_at(t, g, i)] = 1;
        for (int h = 0; h < nh; h++) {
            miss[i][h] = -1;
            const TomhtSth *s = &tk->h[h];
            if (!used[h] || !(s->r > 0.0)) continue;
            double xm[NX], Pm[NX * NX];
            imm_estimate(t, s, xm, Pm);
            int anyFinite = 0;
            for (int j = 0; j < nPl; j++) {
                const double tp = B[j].time_s;
                if (tp < s->tOpp + T / 2.0) continue;
                if (!coarse_gate(t, xm, Pm, s->t, Z + 3 * j, tp)) continue;
                TomhtSth hn;
                double logl;
                if (!imm_ukf_update(t, s, Z + 3 * j, Zvr[j], B[j].vrSpan, tp, hasVr[j], &logl, &hn)) continue;
                const double dl = log(t->pd) + logl - logC[j];
                L[i][h * nPl + j] = dl;
                hn.score = s->score + dl;
                hn.maxScore = s->maxScore > hn.score ? s->maxScore : hn.score;
                hn.plot = B[j].plotIndex;
                hn.batch = t->batch;
                if (hn.nHist >= PMBM_HIST) {
                    memmove(hn.hist, hn.hist + hn.nHist - (PMBM_HIST - 1), sizeof(int) * (PMBM_HIST - 1));
                    hn.nHist = PMBM_HIST - 1;
                }
                hn.hist[hn.nHist++] = B[j].plotIndex;
                dec_push(&hn, tp, B[j].plotIndex);
                hn.tNext = next_look(t, &hn);
                if (nPool == capPool) pool = xrealloc(pool, sizeof(TomhtSth) * (capPool = capPool ? 2 * capPool : 16));
                pool[nPool] = hn;
                upd[i][h * nPl + j] = nPool++;
                anyFinite = 1;
            }
            if (anyFinite && s->tNext <= tMissLimit) {
                TomhtSth ms = *s;
                L0[i][h] = miss_update(t, &ms, tk->confirmed);
                if (nPool == capPool) pool = xrealloc(pool, sizeof(TomhtSth) * (capPool = capPool ? 2 * capPool : 16));
                pool[nPool] = ms;
                miss[i][h] = nPool++;
            }
        }
        free(used);
    }
    /* global hypotheses: k best assignments of each */
    const int nC = nT + nPl;
    double mx = -INFINITY, wsum = 0.0;
    for (int g = 0; g < t->nG; g++) if (t->logw[g] > mx) mx = t->logw[g];
    double *w = xmalloc(sizeof(double) * t->nG);
    for (int g = 0; g < t->nG; g++) { w[g] = exp(t->logw[g] - mx); wsum += w[g]; }
    for (int g = 0; g < t->nG; g++) w[g] /= wsum;
    int nRows = 0, capRows = 16;
    int *rows = xmalloc(sizeof(int) * capRows * nC);
    double *rowW = xmalloc(sizeof(double) * capRows);
    int **keys = xcalloc(nT ? nT : 1, sizeof(int *)), *nKeys = xcalloc(nT ? nT : 1, sizeof(int));
    int *capKeys = xcalloc(nT ? nT : 1, sizeof(int));
    int *part = xmalloc(sizeof(int) * (nT ? nT : 1)), *plotOf = xmalloc(sizeof(int) * (nT ? nT : 1));
    for (int g = 0; g < t->nG; g++) {
        int nP = 0;
        for (int i = 0; i < nT; i++) {
            const int h = G_at(t, g, i);
            for (int j = 0; j < nPl; j++) if (isfinite(L[i][h * nPl + j])) { part[nP++] = i; break; }
        }
        const int m = nP + nPl;
        double *C = xmalloc(sizeof(double) * nPl * m), base = 0.0;
        for (int k = 0; k < nPl * m; k++) C[k] = INFINITY;
        for (int c = 0; c < nP; c++) {
            const int i = part[c], h = G_at(t, g, i);
            for (int j = 0; j < nPl; j++) C[j * m + c] = -(L[i][h * nPl + j] - L0[i][h]);
            base += L0[i][h];
        }
        for (int j = 0; j < nPl; j++) C[j * m + nP + j] = -logNew[j];
        int k = (int)ceil(t->maxHyp * w[g]);
        if (k > t->kBest) k = t->kBest;
        if (k < 1) k = 1;
        int *A;
        double *cost;
        const int nA = murty(C, nPl, m, k, &A, &cost);
        for (int a = 0; a < nA; a++) {
            if (nRows == capRows) {
                capRows *= 2;
                rows = xrealloc(rows, sizeof(int) * capRows * nC);
                rowW = xrealloc(rowW, sizeof(double) * capRows);
            }
            int *row = rows + (size_t)nRows * nC;
            for (int i = 0; i < nT; i++) plotOf[i] = 0;
            for (int j = 0; j < nPl; j++) if (A[a * nPl + j] < nP) plotOf[part[A[a * nPl + j]]] = j + 1;
            for (int i = 0; i < nT; i++) {
                const int kh = G_at(t, g, i), kj = plotOf[i];
                int idx = -1;
                for (int q = 0; q < nKeys[i]; q++) if (keys[i][2 * q] == kh && keys[i][2 * q + 1] == kj) { idx = q; break; }
                if (idx < 0) {
                    if (nKeys[i] == capKeys[i]) keys[i] = xrealloc(keys[i], sizeof(int) * 2 * (capKeys[i] = capKeys[i] ? 2 * capKeys[i] : 4));
                    keys[i][2 * nKeys[i]] = kh;
                    keys[i][2 * nKeys[i] + 1] = kj;
                    idx = nKeys[i]++;
                }
                row[i] = idx;
            }
            for (int j = 0; j < nPl; j++) row[nT + j] = A[a * nPl + j] == nP + j ? 1 : 0;   /* 1 new track, 0 in an older one */
            rowW[nRows++] = t->logw[g] + base - cost[a];
        }
        free(A); free(cost); free(C);
    }
    /* the children become the leaves of the trees */
    for (int i = 0; i < nT; i++) {
        TomhtTrack *tk = &t->T[i];
        TomhtSth *Hn = xmalloc(sizeof(TomhtSth) * (nKeys[i] ? nKeys[i] : 1));
        for (int q = 0; q < nKeys[i]; q++) {
            const int h = keys[i][2 * q], j = keys[i][2 * q + 1];
            if (j > 0) Hn[q] = pool[upd[i][h * nPl + j - 1]];
            else if (miss[i][h] >= 0) Hn[q] = pool[miss[i][h]];
            else Hn[q] = tk->h[h];
        }
        free(tk->h);
        tk->h = Hn;
        tk->nh = tk->cap = nKeys[i];
    }
    if (t->nT + nPl > t->capT) {
        t->capT = 2 * (t->nT + nPl);
        t->T = xrealloc(t->T, sizeof(TomhtTrack) * t->capT);
    }
    for (int j = 0; j < nPl; j++) {                      /* new trees, one per plot */
        TomhtTrack *tk = &t->T[t->nT++];
        memset(tk, 0, sizeof *tk);
        TomhtSth e = empty_sth(), born = new_sth(t, Z + 3 * j, Zvr[j], &B[j], 1.0);
        born.batch = t->batch;
        born.score = born.maxScore = logNew[j];
        dec_push(&born, B[j].time_s, B[j].plotIndex);
        born.tNext = next_look(t, &born);
        push_sth(tk, &e);
        push_sth(tk, &born);
        tk->rMarg = 1.0;
        tk->plot = B[j].plotIndex;
        tk->t0 = B[j].time_s;
    }
    free(t->G); free(t->logw);
    t->G = rows;  t->logw = rowW;  t->nG = nRows;
    if (tMax > t->tNow) t->tNow = tMax;
    for (int i = 0; i < nT; i++) { free(L[i]); free(L0[i]); free(upd[i]); free(miss[i]); free(keys[i]); }
    free(L); free(L0); free(upd); free(miss); free(keys); free(nKeys); free(capKeys); free(part); free(plotOf);
    free(pool); free(w); free(B); free(Z); free(Zvr); free(hasVr); free(logC); free(logNew);
    cleanup(t);
}

/* Missed looks before tSafe (track_tomht.m maintain). Returns 1 when something changed. */
static int maintain(TomhtTracker *t, double tSafe)
{
    int changed = 0;
    for (int i = 0; i < t->nT; i++)
        for (int h = 0; h < t->T[i].nh; h++) {
            TomhtSth *s = &t->T[i].h[h];
            while (s->r > 0.0 && s->tNext <= tSafe) {
                const double l0 = miss_update(t, s, t->T[i].confirmed);
                for (int g = 0; g < t->nG; g++) if (G_at(t, g, i) == h) t->logw[g] += l0;
                changed = 1;
            }
        }
    if (changed) {
        if (tSafe > t->tNow) t->tNow = tSafe;
        cleanup(t);
    }
    return changed;
}

/* ---- output ------------------------------------------------------------------- */

static void make_track(const TomhtTracker *t, const TomhtTrack *tk, const TomhtSth *s, Track *r)
{
    const int d = t->dim;
    double x[NX];
    imm_estimate(t, s, x, NULL);
    memset(r, 0, sizeof *r);
    r->id = tk->id;
    r->confirmed = tk->confirmed;
    r->dim = d;
    r->time_s = s->t;
    for (int i = 0; i < d; i++) {
        r->pos_m[i] = x[i];
        r->vel_mps[i] = x[d + i];
        r->acc_mps2[i] = x[2 * d + i];
    }
    r->alt_m = d == 3 ? x[2] + t->height_m : NAN;
    r->hits = s->hits;
    r->misses = s->misses;
    r->method = TRACK_METHOD_TOMHT;
    r->existence = tk->rTarget;
    r->nModes = t->M;
    for (int j = 0; j < t->M; j++) { r->modes[j] = t->models[j]; r->modeProb[j] = s->mu[j]; }
}

/* A confirmed track shown before that is not in the most likely hypothesis
 * any more hands its number to the track of this hypothesis that holds its
 * newest plot (track_tomht.m keep_numbers). */
static void keep_numbers(TomhtTracker *t, int best)
{
    if (t->nConf == 0) return;
    int *exists = xmalloc(sizeof(int) * (t->nT ? t->nT : 1));
    for (int i = 0; i < t->nT; i++) exists[i] = t->T[i].h[G_at(t, best, i)].r > 0.0;
    for (int s = 0; s < t->nConf; s++) {
        int iP = -1;
        for (int i = 0; i < t->nT; i++) if (t->T[i].id == t->confId[s]) { iP = i; break; }
        if (iP >= 0 && exists[iP]) continue;             /* still shown */
        for (int iQ = 0; iQ < t->nT; iQ++) {
            if (!exists[iQ]) continue;
            int shownConf = 0;
            for (int q = 0; q < t->nConf; q++) shownConf |= t->confId[q] == t->T[iQ].id;
            if (shownConf) continue;
            const TomhtSth *hq = &t->T[iQ].h[G_at(t, best, iQ)];
            int holds = 0;
            for (int k = 0; k < hq->nHist; k++) holds |= hq->hist[k] == t->confPlot[s];
            if (!holds) continue;
            if (iP >= 0) {
                t->T[iP].id = t->T[iQ].id;
                t->T[iP].confirmed = t->T[iQ].confirmed;
            }
            t->T[iQ].id = t->confId[s];
            t->T[iQ].confirmed = 1;
            break;
        }
    }
    free(exists);
}

static int cmp_track_id(const void *a, const void *b)
{
    const Track *x = a, *y = b;
    return x->id < y->id ? -1 : x->id > y->id;
}

typedef struct { int plot, idx; } LogRow;
static int cmp_logrow(const void *a, const void *b)
{
    const LogRow *x = a, *y = b;
    return x->plot < y->plot ? -1 : x->plot > y->plot;
}

/* The tracks of the most likely global hypothesis, the updates by the plots of
 * the newest batch (callback, as the MATLAB log) and the deleted tracks. */
static void publish(TomhtTracker *t, TrackerCallback callback, void *user)
{
    Track *prev = t->shown;
    const int nPrev = t->nShown;
    t->shown = xmalloc(sizeof(Track) * (t->nT ? t->nT : 1));
    t->nShown = 0;
    t->capShown = t->nT;
    LogRow *rows = xmalloc(sizeof(LogRow) * (t->nT ? t->nT : 1));
    int nRows = 0, best = 0;
    if (t->nG > 0 && t->nT > 0) {
        best = best_hyp(t);
        keep_numbers(t, best);
        for (int i = 0; i < t->nT; i++) {
            const TomhtSth *s = &t->T[i].h[G_at(t, best, i)];
            if (!(s->r > 0.0)) continue;
            make_track(t, &t->T[i], s, &t->shown[t->nShown]);
            if (s->batch == t->batch && s->plot > t->logged) {   /* updated by a plot of this batch */
                rows[nRows].plot = s->plot;
                rows[nRows++].idx = t->nShown;
            }
            t->nShown++;
        }
    }
    if (nRows > 0) {
        qsort(rows, nRows, sizeof(LogRow), cmp_logrow);
        for (int k = 0; k < nRows; k++)
            if (callback) callback(user, &t->shown[rows[k].idx], rows[k].plot);
        t->logged = rows[nRows - 1].plot;
    }
    free(rows);
    qsort(t->shown, t->nShown, sizeof(Track), cmp_track_id);
    for (int k = 0; k < nPrev; k++) {                    /* tracks that left the most likely hypothesis */
        int found = 0;
        for (int q = 0; q < t->nShown && !found; q++) found = t->shown[q].id == prev[k].id;
        if (!found && callback) {
            Track r = prev[k];
            r.deleted = 1;
            callback(user, &r, 0);
        }
    }
    free(prev);
    /* the confirmed tracks shown, with their newest plot (keep_numbers) */
    t->nConf = 0;
    if (t->nG > 0 && t->nT > 0) {
        if (t->capConf < t->nT) {
            t->capConf = t->nT;
            t->confId = xrealloc(t->confId, sizeof(uint32_t) * t->capConf);
            t->confPlot = xrealloc(t->confPlot, sizeof(int) * t->capConf);
        }
        for (int i = 0; i < t->nT; i++) {
            const TomhtSth *s = &t->T[i].h[G_at(t, best, i)];
            if (s->r > 0.0 && t->T[i].confirmed) {
                t->confId[t->nConf] = t->T[i].id;
                t->confPlot[t->nConf++] = s->plot;
            }
        }
    }
}

void tomht_tracker_update(TomhtTracker *t, const Detection *plots, int n, double t_s, int flush,
                         TrackerCallback callback, void *user)
{
    int changed = 0;
    /* 1) the batch is complete */
    if (t->nPending > 0 && (t_s - t->pendingLast >= t->batchGap_s || t_s - t->pendingFirst >= t->batchMax_s)) {
        associate(t);
        maintain(t, t_s - t->latency_s);
        t->tMaint = t_s;
        changed = 1;
    }
    /* 2) the new plots join the batch */
    for (int k = 0; k < n; k++) {
        if (t->nPending == t->capPending) {
            t->capPending = t->capPending ? 2 * t->capPending : 16;
            t->pending = xrealloc(t->pending, sizeof(PendingPlot) * t->capPending);
        }
        if (t->nPending == 0) t->pendingFirst = t_s;
        PendingPlot *b = &t->pending[t->nPending++];
        b->time_s = plots[k].time_s;
        b->az_deg = plots[k].az_deg;
        b->range_m = plots[k].range_m;
        b->el_deg = plots[k].el_deg;
        b->vr_mps = plots[k].vr_mps;
        b->vrSpan = plots[k].pulse >= 0 && plots[k].pulse < N_PULSES ? t->vrSpanPulse[plots[k].pulse] : t->vrSpan;
        b->plotIndex = ++t->nPlots;
        t->pendingLast = t_s;
    }
    if (flush && t->nPending > 0) {
        associate(t);
        changed = 1;
    }
    /* 3) missed looks (only when no plot is waiting) */
    if (t->nPending == 0 && (flush || t_s - t->tMaint >= t->maintStep_s)) {
        changed |= maintain(t, t_s - t->latency_s);
        t->tMaint = t_s;
    }
    if (changed) publish(t, callback, user);
}

void tomht_tracker_report(const TomhtTracker *t, int index, Track *out)
{
    *out = t->shown[index];
}

/*
 * export_config.c - radar_export_config: writes the receiver configuration
 * (config/radar_config.bin) from the MATLAB configuration file radar_config.m,
 * without MATLAB. The same file as export_receiver_config.m writes:
 *
 *   radar_export_config [--config-m FILE] [--out FILE] [--set STATEMENT]... [--quiet]
 *     --config-m FILE   the configuration function (default ../radar_config.m of this project)
 *     --out FILE        the receiver configuration (default config/radar_config.bin of this project)
 *     --set STATEMENT   a MATLAB assignment run after the file, e.g. --set "cfg.prf.mode = 'staggered';"
 *                       (repeatable; like changing cfg in MATLAB before export_receiver_config(cfg))
 *     --quiet           print only errors
 *
 * radar_config.m is evaluated by a small interpreter of MATLAB assignments
 * (mfile.c); the input files it names (antenna pattern, DDC filters, DBF
 * matrix, stagger PRIs, phase codes, FIR clutter filter) are read with
 * mat_file.c. Then everything export_receiver_config.m writes is derived as
 * in setup_radar.m: 3 dB azimuth beamwidth, DDC rate and record timing,
 * receive intervals, pulse compression filters (Taylor-windowed LFM, phase
 * codes with the matched or the mismatch filter and its peak shift), range
 * axes and the handover, CFAR-tested cells, CPI, MTI weights (canceller or
 * FIR), Doppler window, noise gain and looks of every bin, median/mean of the
 * noise, the CFAR detector (CA, GO, SO, OS) and its factor from the Pfa,
 * extraction settings, EKF or PMBM tracker settings, peak interpolation. The results agree with MATLAB to rounding (the CFAR factor
 * and the looks come from root finding: ~1e-12).
 *
 * Unlike MATLAB, missing input files or default files made for other
 * configuration values (their 'madeFor') are not rebuilt: the program says
 * which file to rebuild (make_input_files.m) or replace.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "mat_file.h"
#include "mfile.h"

#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef RADAR_SOURCE_DIR
#  define RADAR_SOURCE_DIR "."
#endif

#define PI 3.14159265358979323846

static jmp_buf g_jb;
static char g_err[1024];
static int g_quiet = 0;

static void fail(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(g_err, sizeof g_err, fmt, ap);
    va_end(ap);
    longjmp(g_jb, 1);
}

static void warn(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "warning: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
}

static void *xcalloc(size_t n, size_t size)
{
    void *p = calloc(n ? n : 1, size);
    if (!p) fail("out of memory");
    return p;
}

/* ------------------------------------------------------------------------ */
/* Reading cfg                                                              */
/* ------------------------------------------------------------------------ */

static const MValue *g_cfg;

static const MValue *field(const char *path, int required)
{
    const MValue *v = mvalue_get(g_cfg, path);
    if (!v) {
        if (required) fail("cfg.%s is missing in the configuration", path);
        return NULL;
    }
    if (v->type == MV_UNEVAL) fail("cfg.%s: %s", path, v->str);
    return v;
}

static int has(const char *path) { return mvalue_get(g_cfg, path) != NULL; }

static int empty(const char *path)
{
    const MValue *v = field(path, 1);
    return v->type == MV_NUM && v->rows * v->cols == 0;
}

static double num(const char *path)
{
    const MValue *v = field(path, 1);
    if (v->type != MV_NUM || v->rows * v->cols != 1) fail("cfg.%s must be a number", path);
    return v->re[0];
}

static double num_or(const char *path, double def)
{
    return has(path) && !empty(path) ? num(path) : def;
}

static const char *str(const char *path)
{
    const MValue *v = field(path, 1);
    if (v->type != MV_STR) fail("cfg.%s must be a string", path);
    return v->str;
}

static int same_text(const char *a, const char *b)            /* strcmpi */
{
    for (; *a && *b; a++, b++)
        if ((*a | 32) != (*b | 32)) return 0;
    return *a == *b;
}

/* ------------------------------------------------------------------------ */
/* Records of the configuration file                                        */
/* ------------------------------------------------------------------------ */

typedef struct {
    char    name[32];
    int     cplx, rows, cols;
    double *re, *im;
} Rec;

static Rec *g_rec;
static int g_nRec, g_capRec;

static Rec *add_rec(const char *name, int rows, int cols, int cplx)
{
    if (strlen(name) >= 32) fail("record name %s too long", name);
    if (g_nRec == g_capRec) {
        g_capRec = g_capRec ? 2 * g_capRec : 128;
        g_rec = realloc(g_rec, sizeof(Rec) * (size_t)g_capRec);
        if (!g_rec) fail("out of memory");
    }
    Rec *r = &g_rec[g_nRec++];
    memset(r, 0, sizeof *r);
    strcpy(r->name, name);
    r->rows = rows;
    r->cols = cols;
    r->cplx = cplx;
    r->re = xcalloc((size_t)(rows * cols), sizeof(double));
    if (cplx) r->im = xcalloc((size_t)(rows * cols), sizeof(double));
    return r;
}

static void rec_scalar(const char *name, double x) { add_rec(name, 1, 1, 0)->re[0] = x; }

static void rec_vec(const char *name, const double *x, int n, int column)
{
    Rec *r = add_rec(name, column ? n : 1, column ? 1 : n, 0);
    memcpy(r->re, x, sizeof(double) * (size_t)n);
}

static void rec_cvec(const char *name, const double *re, const double *im, int n)
{
    Rec *r = add_rec(name, 1, n, 1);
    memcpy(r->re, re, sizeof(double) * (size_t)n);
    memcpy(r->im, im, sizeof(double) * (size_t)n);
}

static void pfx(char *buf, const char *a, const char *b) { snprintf(buf, 32, "%s%s", a, b); }

static int write_file(const char *path)
{
    FILE *f = fopen(path, "wb");
    if (!f) return -1;
    uint32_t h[2] = { 2, (uint32_t)g_nRec };
    fwrite("RADRXCFG", 1, 8, f);
    fwrite(h, 4, 2, f);                                  /* little-endian hosts (x86, ARM) */
    for (int i = 0; i < g_nRec; i++) {
        const Rec *r = &g_rec[i];
        char name[32];
        memset(name, 0, sizeof name);
        memcpy(name, r->name, strlen(r->name));
        fwrite(name, 1, 32, f);
        uint32_t t[3] = { (uint32_t)r->cplx, (uint32_t)r->rows, (uint32_t)r->cols };
        fwrite(t, 4, 3, f);
        for (int k = 0; k < r->rows * r->cols; k++) {
            fwrite(&r->re[k], 8, 1, f);
            if (r->cplx) fwrite(&r->im[k], 8, 1, f);
        }
    }
    return fclose(f) == 0 ? 0 : -1;
}

/* ------------------------------------------------------------------------ */
/* Numerics: incomplete gamma and beta functions, root finding, eigenvalues */
/* ------------------------------------------------------------------------ */

/* Regularised lower incomplete gamma P(a, x) (MATLAB gammainc(x, a)). */
static double gammainc_p(double x, double a)
{
    if (x <= 0) return 0.0;
    const double lg = a * log(x) - x - lgamma(a);
    if (x < a + 1.0) {                                    /* series */
        double ap = a, sum = 1.0 / a, del = sum;
        for (int n = 0; n < 10000; n++) {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * 1e-17) break;
        }
        return sum * exp(lg);
    }
    double b = x + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;   /* continued fraction (Lentz) */
    for (int i = 1; i < 10000; i++) {
        const double an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < 1e-300) d = 1e-300;
        c = b + an / c;
        if (fabs(c) < 1e-300) c = 1e-300;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < 1e-17) break;
    }
    return 1.0 - exp(lg) * h;
}

static double beta_cf(double a, double b, double x)
{
    const double tiny = 1e-300;
    double c = 1.0, d = 1.0 - (a + b) * x / (a + 1.0);
    if (fabs(d) < tiny) d = tiny;
    d = 1.0 / d;
    double h = d;
    for (int m = 1; m < 100000; m++) {
        const int m2 = 2 * m;
        double aa = m * (b - m) * x / ((a + m2 - 1.0) * (a + m2));
        d = 1.0 + aa * d; if (fabs(d) < tiny) d = tiny;
        c = 1.0 + aa / c; if (fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0));
        d = 1.0 + aa * d; if (fabs(d) < tiny) d = tiny;
        c = 1.0 + aa / c; if (fabs(c) < tiny) c = tiny;
        d = 1.0 / d;
        const double del = d * c;
        h *= del;
        if (fabs(del - 1.0) < 1e-16) break;
    }
    return h;
}

/* Regularised incomplete beta I_x(a, b) (MATLAB betainc(x, a, b)). */
static double betainc_i(double x, double a, double b)
{
    if (x <= 0) return 0.0;
    if (x >= 1) return 1.0;
    const double lbt = lgamma(a + b) - lgamma(a) - lgamma(b) + a * log(x) + b * log1p(-x);
    if (x < (a + 1.0) / (a + b + 2.0)) return exp(lbt) * beta_cf(a, b, x) / a;
    return 1.0 - exp(lbt) * beta_cf(b, a, 1.0 - x) / b;
}

/* Root of f in [a, b] (f(a), f(b) of opposite signs): Brent's method, as fzero, to machine precision. */
static double fzero_brent(double (*f)(double, const void *), const void *ctx, double a, double b)
{
    double fa = f(a, ctx), fb = f(b, ctx);
    if (fa == 0) return a;
    if (fb == 0) return b;
    if ((fa > 0) == (fb > 0)) fail("internal: no sign change for the root finder");
    double c = a, fc = fa, d = b - a, e = d;
    for (int it = 0; it < 500; it++) {
        if ((fb > 0) == (fc > 0)) { c = a; fc = fa; d = e = b - a; }
        if (fabs(fc) < fabs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
        const double tol = 2.0 * 2.220446049250313e-16 * fabs(b), m = 0.5 * (c - b);
        if (fabs(m) <= tol || fb == 0) break;
        if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
            d = e = m;
        } else {
            double s = fb / fa, p, q;
            if (a == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                const double qq = fa / fc, r = fb / fc;
                p = s * (2.0 * m * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0) q = -q; else p = -p;
            if (2.0 * p < 3.0 * m * q - fabs(tol * q) && p < fabs(0.5 * e * q)) { e = d; d = p / q; }
            else d = e = m;
        }
        a = b;
        fa = fb;
        b += fabs(d) > tol ? d : (m > 0 ? tol : -tol);
        fb = f(b, ctx);
    }
    return b;
}

typedef struct { double a, p; } GammaInv;
static double gammainv_f(double x, const void *ctx)
{
    const GammaInv *g = ctx;
    return gammainc_p(x, g->a) - g->p;
}

/* MATLAB gammaincinv(p, a): x with P(a, x) = p. */
static double gammaincinv(double p, double a)
{
    GammaInv g = { a, p };
    double hi = a + 10.0;
    while (gammainv_f(hi, &g) < 0) hi *= 2.0;
    return fzero_brent(gammainv_f, &g, 0.0, hi);
}

typedef struct { double L; } MedianCtx;
static double median_f(double x, const void *ctx)
{
    const MedianCtx *m = ctx;
    return gammainc_p(m->L * x, m->L) - 0.5;
}

/* The CFAR detector (sys.cfar): type 0 CA, 1 GO, 2 SO, 3 OS. */
typedef struct { int type, nTrain, nHalf, osK; double L; } Detector;

/* cfar_statistic.m: P(Z <= z) of the noise statistic for Gamma(L) noise of mean 1. */
static double cfar_statistic(double z, const Detector *c)
{
    const double L = c->L;
    if (z < 0) z = 0;
    switch (c->type) {
    case 0: return gammainc_p(c->nTrain * L * z, c->nTrain * L);
    case 1: { const double F = gammainc_p(c->nHalf * L * z, c->nHalf * L); return F * F; }
    case 2: { const double F = gammainc_p(c->nHalf * L * z, c->nHalf * L); return F * (2 - F); }
    default: return betainc_i(gammainc_p(L * z, L), c->osK, c->nTrain - c->osK + 1);
    }
}

/* cfar_quadrature.m: 8-point Gauss-Legendre on 2000 panels of [0, sHi]. */
static int g_nq;
static double *g_qs, *g_qw;
static void cfar_quadrature(double L)
{
    static const double x[8] = { -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                                  0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363 };
    static const double g[8] = { 0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                                 0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763 };
    const double sHi = (L + 10 * sqrt(L) + 75) / L;
    const int nP = 2000;
    const double h = sHi / nP;
    free(g_qs);
    free(g_qw);
    g_nq = 8 * nP;
    g_qs = xcalloc((size_t)g_nq, sizeof(double));
    g_qw = xcalloc((size_t)g_nq, sizeof(double));
    for (int p = 0; p < nP; p++) {
        const double mid = (p + 0.5) * h;
        for (int i = 0; i < 8; i++) {
            g_qs[p * 8 + i] = mid + x[i] * h / 2;
            g_qw[p * 8 + i] = g[i] * h / 2;
        }
    }
}

/* cfar_pfa.m */
static double cfar_pfa(double alpha, const Detector *c)
{
    const double L = c->L;
    if (c->type == 0) {
        const double N = c->nTrain;
        return betainc_i(1.0 / (1.0 + alpha / N), N * L, L);
    }
    double pfa = 0;
    for (int i = 0; i < g_nq; i++) {
        const double s = g_qs[i];
        const double lf = L * log(L) + (L - 1) * log(s) - L * s - lgamma(L);
        pfa += g_qw[i] * cfar_statistic(s / alpha, c) * exp(lf);
    }
    return pfa;
}

typedef struct { const Detector *c; double logPfa; } AlphaCtx;
static double alpha_f(double a, const void *ctx)
{
    const AlphaCtx *x = ctx;
    return log(cfar_pfa(a, x->c)) - x->logPfa;
}

/* cfar_alpha.m */
static double cfar_alpha(double pfa, const Detector *c)
{
    AlphaCtx x = { c, log(pfa) };
    const double N = c->nTrain, alpha0 = N * (pow(pfa, -1.0 / N) - 1.0);
    double lo = 1e-3 * alpha0, hi = 10 * alpha0;
    while (alpha_f(hi, &x) > 0) hi *= 4;                 /* SO, OS with a high rank: larger factors */
    while (alpha_f(lo, &x) < 0) lo /= 4;
    return fzero_brent(alpha_f, &x, lo, hi);
}

/* Largest eigenvalue of an n x n real symmetric matrix A (row-major, destroyed): cyclic Jacobi. */
static double max_eig_sym(double *A, int n)
{
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
        if (off < 1e-30) break;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
                const double apq = A[p * n + q];
                if (fabs(apq) < 1e-300) continue;
                const double theta = (A[q * n + q] - A[p * n + p]) / (2.0 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
                const double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; k++) {                   /* columns p, q */
                    const double akp = A[k * n + p], akq = A[k * n + q];
                    A[k * n + p] = c * akp - s * akq;
                    A[k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) {                   /* rows p, q */
                    const double apk = A[p * n + k], aqk = A[q * n + k];
                    A[p * n + k] = c * apk - s * aqk;
                    A[q * n + k] = s * apk + c * aqk;
                }
            }
    }
    double m = A[0];
    for (int i = 1; i < n; i++) m = A[i * n + i] > m ? A[i * n + i] : m;
    return m;
}

static int nextpow2(double x)                             /* 2^p >= |x| */
{
    int p = 0;
    while (ldexp(1.0, p) < fabs(x)) p++;
    return p;
}

static double sind(double x)                              /* as Octave: exact zeros at multiples of 180 */
{
    const double i = x / 180.0;
    return i == floor(i) && isfinite(i) ? 0.0 : sin(i * PI);
}

/* ------------------------------------------------------------------------ */
/* Input files                                                              */
/* ------------------------------------------------------------------------ */

typedef struct { const char *name; double v[4]; int n; } MadeFor;

/* Loads cfg.files.<key>. A 'madeFor' in the file (a default file of
 * make_input_files.m) must match the values the configuration needs: if it
 * does not and canBuild is set, *stale is set (the caller then computes the
 * default data as make_input_files.m; MATLAB would rewrite the file),
 * otherwise this is an error. */
static MatFile load_input_ex(const char *key, const MadeFor *mf, int nMf, int canBuild, int *stale)
{
    if (stale) *stale = 0;
    char path[64];
    snprintf(path, sizeof path, "files.%s", key);
    const char *file = str(path);
    MatFile m;
    char err[256];
    if (mat_load(file, &m, err, sizeof err) != 0)
        fail("cfg.files.%s: %s (create the default input files in MATLAB with make_input_files(cfg), "
             "or give your own file)", key, err);
    const MatVar *made = mat_get(&m, "madeFor");
    if (made && nMf > 0) {
        int ok = made->mxClass == MAT_CLASS_STRUCT && made->nFields == nMf;
        const char *bad = ok ? NULL : "its fields";
        for (int i = 0; i < nMf && ok; i++) {
            const MatVar *f = mat_field(made, mf[i].name);
            ok = f && (int)f->numel == mf[i].n;
            for (int k = 0; ok && k < mf[i].n; k++) ok = f->re && f->re[k] == mf[i].v[k];
            if (!ok) bad = mf[i].name;
        }
        if (!ok && canBuild) {
            if (!g_quiet)
                printf("note: %s was made for other configuration values (madeFor.%s): its default data are computed "
                       "as make_input_files.m does (MATLAB would also rewrite the file)\n", file, bad);
            *stale = 1;
        } else if (!ok)
            fail("%s is a default input file made for other configuration values (madeFor.%s). MATLAB "
                 "rebuilds it: run make_input_files(cfg, '%s') in MATLAB (or export_receiver_config once), "
                 "or replace it with your own file without 'madeFor'", file, bad, key);
    }
    return m;
}

static MatFile load_input(const char *key, const MadeFor *mf, int nMf)
{
    return load_input_ex(key, mf, nMf, 0, NULL);
}

static const MatVar *mat_var(const MatFile *m, const char *name, const char *key)
{
    const MatVar *v = mat_get(m, name);
    if (!v || !v->re) fail("cfg.files.%s: the variable '%s' is missing", key, name);
    return v;
}

/* ------------------------------------------------------------------------ */
/* Waveforms and pulse compression (setup_radar.m: compression_filter)       */
/* ------------------------------------------------------------------------ */

typedef struct {
    const char *name;
    int    enable;                /* cfg.wave.<pulse>.enable                   */
    double fc_Hz, lambda_m;       /* its carrier (cfg.wave.<pulse>.fc_Hz)      */
    double *beamU;                /* beam directions on its carrier            */
    int    phaseCode, mismatch;
    double width_s, bandwidth_Hz, chip_s, tTx_s;
    int    nCode;
    double *codeRe, *codeIm;      /* phase code chips                          */
    int    nMis;
    double *misRe, *misIm;        /* mismatch filter, one coefficient per period */
    double misRate_Hz;
    /* results */
    int    nTaps, mfShift;
    double *tapsRe, *tapsIm;
    int    rxFirst, rxCount;      /* 1-based first receive sample (as MATLAB) */
    int    nRange, ddcFirst1, ddcCount;
    double *delay_s, *range_m;
    int    testFirst, testLast;   /* 1-based */
} Wave;

static void envelope(const Wave *w, double t, double *re, double *im)
{
    *re = *im = 0.0;
    if (t < 0 || t >= w->width_s) return;
    if (!w->phaseCode) {
        const double k = w->bandwidth_Hz / w->width_s, d = t - w->width_s / 2;
        const double ph = PI * k * (d * d);
        *re = cos(ph);
        *im = sin(ph);
    } else {
        double m = floor(t / w->chip_s);
        if (m < 0) m = 0;
        if (m > w->nCode - 1) m = w->nCode - 1;
        *re = w->codeRe[(int)m];
        *im = w->codeIm[(int)m];
    }
}

/* taylor_window.m */
static void taylor_window(int N, int nbar, double sll_dB, double *w)
{
    const double A = acosh(pow(10.0, -sll_dB / 20.0)) / PI;
    const double sp2 = nbar * nbar / (A * A + (nbar - 0.5) * (nbar - 0.5));
    for (int k = 0; k < N; k++) w[k] = 1.0;
    for (int m = 1; m <= nbar - 1; m++) {
        double num = 1.0, den = 1.0;
        for (int n = 1; n <= nbar - 1; n++) {
            num *= 1.0 - (m * (double)m / sp2) / (A * A + (n - 0.5) * (n - 0.5));
            if (n != m) den *= 1.0 - (double)m * m / ((double)n * n);
        }
        const double Fm = ((m + 1) % 2 ? -1.0 : 1.0) * num / (2.0 * den);
        for (int k = 0; k < N; k++) {
            const double xi = (k - 0.5 * N + 0.5) / N;
            w[k] += 2.0 * Fm * cos(2.0 * PI * m * xi);
        }
    }
}

typedef struct {
    double fs, fsOut;
    int    D1, D2, L1, L2;
    double *h1, *h2;              /* real filters (unit DC gain shape)          */
} Ddc;

/* fir_decimate_valid.m (complex input, real filter). */
static int fir_decimate(const double *xr, const double *xi, int N, const double *h, int L, int D, double *yr, double *yi)
{
    const int K = N >= L ? (N - L) / D + 1 : 0;
    for (int k = 0; k < K; k++) {
        double sr = 0, si = 0;
        for (int j = 0; j < L; j++) {
            const int i = k * D + L - 1 - j;
            sr += h[j] * xr[i];
            si += h[j] * xi[i];
        }
        yr[k] = sr;
        yi[k] = si;
    }
    return K;
}

/* Position (0-based, parabolic interpolation) of the largest |y| (peak_position in setup_radar.m). */
static double peak_position(const double *yr, const double *yi, int n)
{
    int k = 0;
    double best = -1;
    for (int i = 0; i < n; i++) {
        const double a = hypot(yr[i], yi[i]);
        if (a > best) { best = a; k = i; }
    }
    double q = k;
    if (k > 0 && k < n - 1) {
        const double a0 = hypot(yr[k - 1], yi[k - 1]), a1 = best, a2 = hypot(yr[k + 1], yi[k + 1]);
        const double d = a0 - 2 * a1 + a2;
        if (d < 0) q += 0.5 * (a0 - a2) / d;
    }
    return q;
}

static void conv_c(const double *ar, const double *ai, int na, const double *br, const double *bi, int nb,
                   double *yr, double *yi)
{
    for (int i = 0; i < na + nb - 1; i++) yr[i] = yi[i] = 0;
    for (int i = 0; i < na; i++)
        for (int j = 0; j < nb; j++) {
            yr[i + j] += ar[i] * br[j] - ai[i] * bi[j];
            yi[i + j] += ar[i] * bi[j] + ai[i] * br[j];
        }
}

static void compression_filter(Wave *w, const Ddc *ddc, int taylorNbar, double taylorSll_dB)
{
    const double fsOut = ddc->fsOut;
    const int Lp = (int)lround(w->width_s * fsOut);
    double *win = xcalloc((size_t)Lp, sizeof(double));
    double *hmr = xcalloc((size_t)Lp, sizeof(double)), *hmi = xcalloc((size_t)Lp, sizeof(double));
    if (w->phaseCode) for (int k = 0; k < Lp; k++) win[k] = 1.0;
    else taylor_window(Lp, taylorNbar, taylorSll_dB, win);
    for (int k = 0; k < Lp; k++) {                        /* conj(flipud(ref .* win)) */
        double re, im;
        envelope(w, (k + 0.5) / fsOut, &re, &im);
        hmr[Lp - 1 - k] = re * win[k];
        hmi[Lp - 1 - k] = -im * win[k];
    }
    double *hr, *hi;
    int L;
    if (!w->mismatch) {
        hr = hmr;
        hi = hmi;
        L = Lp;
        w->mfShift = 0;
    } else {
        if (!w->phaseCode)
            fail("cfg.wave.%s.compression = 'mismatch' needs a phase code (type 'phasecode').", w->name);
        const double spc = fsOut / w->misRate_Hz;
        if (fabs(spc - round(spc)) > 1e-9 || spc < 0.5)
            fail("The %s-pulse mismatch filter (%g MHz) does not fit the DDC output rate (%g MHz): the rate must "
                 "be a whole multiple of the filter rate.", w->name, w->misRate_Hz / 1e6, fsOut / 1e6);
        const int r = (int)lround(spc);
        L = w->nMis * r;
        hr = xcalloc((size_t)L, sizeof(double));
        hi = xcalloc((size_t)L, sizeof(double));
        for (int i = 0; i < w->nMis; i++)
            for (int k = 0; k < r; k++) { hr[i * r + k] = w->misRe[i]; hi[i * r + k] = w->misIm[i]; }
        /* ddc_pulse.m: the pulse after the DDC filters at the DDC rate */
        const int nPad = ddc->L1 + ddc->D1 * (ddc->L2 - 1);
        const int N = (int)lround(w->width_s * ddc->fs) + 2 * nPad;
        double *xr = xcalloc((size_t)N, sizeof(double)), *xi = xcalloc((size_t)N, sizeof(double));
        for (int i = 0; i < N; i++) envelope(w, (i - nPad) / ddc->fs, &xr[i], &xi[i]);
        double *y1r = xcalloc((size_t)N, sizeof(double)), *y1i = xcalloc((size_t)N, sizeof(double));
        const int K1 = fir_decimate(xr, xi, N, ddc->h1, ddc->L1, ddc->D1, y1r, y1i);
        double *yr = xcalloc((size_t)N, sizeof(double)), *yi = xcalloc((size_t)N, sizeof(double));
        const int K2 = fir_decimate(y1r, y1i, K1, ddc->h2, ddc->L2, ddc->D2, yr, yi);
        const int nc = K2 + (L > Lp ? L : Lp);
        double *cr = xcalloc((size_t)nc, sizeof(double)), *ci = xcalloc((size_t)nc, sizeof(double));
        conv_c(yr, yi, K2, hmr, hmi, Lp, cr, ci);
        const double qMatched = peak_position(cr, ci, K2 + Lp - 1);
        conv_c(yr, yi, K2, hr, hi, L, cr, ci);
        const double q = peak_position(cr, ci, K2 + L - 1);
        w->mfShift = (int)lround((L - 1 - q) - (Lp - 1 - qMatched));
        free(xr); free(xi); free(y1r); free(y1i); free(yr); free(yi); free(cr); free(ci);
    }
    double nrm = 0;
    for (int k = 0; k < L; k++) nrm += hr[k] * hr[k] + hi[k] * hi[k];
    nrm = sqrt(nrm);
    w->nTaps = L;
    w->tapsRe = xcalloc((size_t)L, sizeof(double));
    w->tapsIm = xcalloc((size_t)L, sizeof(double));
    for (int k = 0; k < L; k++) { w->tapsRe[k] = hr[k] / nrm; w->tapsIm[k] = hi[k] / nrm; }
    if (hr != hmr) { free(hr); free(hi); }
    free(hmr); free(hmi); free(win);
}

/* ------------------------------------------------------------------------ */
/* Default input data of make_input_files.m (when a default file is stale)   */
/* ------------------------------------------------------------------------ */

typedef struct { double d; const double *w; int n; } ArrayCtx;

/* Azimuth array factor of make_input_files.m: |sum w_n exp(j 2 pi d n sind(az))|^2 / sum(w)^2 */
static double af_power(const ArrayCtx *a, double az)
{
    double re = 0, im = 0, sw = 0;
    const double s = sind(az);
    for (int n = 0; n < a->n; n++) {
        const double ph = 2 * PI * a->d * n * s;
        re += a->w[n] * cos(ph);
        im += a->w[n] * sin(ph);
        sw += a->w[n];
    }
    return (re * re + im * im) / (sw * sw);
}

static double half_bw(double d, const void *ctx)
{
    ArrayCtx a = *(const ArrayCtx *)ctx;
    a.d = d;
    return 10 * log10(af_power(&a, 4.1 / 2)) + 3.0103;
}

/* The default azimuth pattern (az_deg, azGain_dB: 18001 points). */
static void build_antenna(double **az, double **gain, int *n)
{
    double w[32];
    taylor_window(32, 4, -30, w);
    ArrayCtx a = { 0, w, 32 };
    a.d = fzero_brent(half_bw, &a, 0.3, 0.7);
    *n = 18001;
    *az = xcalloc((size_t)*n, sizeof(double));
    *gain = xcalloc((size_t)*n, sizeof(double));
    for (int i = 0; i < *n; i++) {
        (*az)[i] = -180.0 + i * 0.02;
        const double g = 10 * log10(af_power(&a, (*az)[i]));
        (*gain)[i] = fabs((*az)[i]) > 90 ? -50.0 : g > -50 ? g : -50.0;
    }
}

/* lowpass_fir of make_input_files.m, scaled to int16 with fracBits fractional bits. */
static void build_lowpass(int nTaps, double fCutNorm, double q, double *h)
{
    double sum = 0;
    for (int k = 0; k < nTaps; k++) {
        const double n = k - (nTaps - 1) / 2.0, x = 2 * fCutNorm * n;
        const double sinc = x != 0 ? sin(PI * x) / (PI * x) : 1.0;
        const double bl = 0.42 - 0.5 * cos(2 * PI * k / (nTaps - 1)) + 0.08 * cos(4 * PI * k / (nTaps - 1));
        h[k] = sinc * bl;
        sum += h[k];
    }
    for (int k = 0; k < nTaps; k++) {
        double v = round(q * (h[k] / sum));
        h[k] = v > 32767 ? 32767 : v < -32768 ? -32768 : v;
    }
}

/* The default elevation DBF matrix (21 beams) and the beam elevations. */
static void build_dbf(int nEl, double spacing, double tilt, double lambda, double *wRe, double *wIm, double *beamEl)
{
    const int nBeams = 21;
    const double esm = spacing * lambda;
    const double uFirst = sind(0.5 - tilt), uStep = 1.2 * 0.886 / (nEl * spacing);
    const double uEnd = uFirst + (nBeams - 1) * uStep, uMax = sind(60 - tilt);
    const double uLast = uEnd < uMax ? uEnd : uMax;
    double *taper = xcalloc((size_t)nEl, sizeof(double));
    taylor_window(nEl, 4, -30, taper);
    for (int b = 0; b < nBeams; b++) {
        const double u = b == nBeams - 1 ? uLast : uFirst + b * (uLast - uFirst) / (nBeams - 1);
        beamEl[b] = tilt + asin(u) * 180 / PI;
        double nrm = 0;
        for (int n = 0; n < nEl; n++) {
            const double ph = -2 * PI * esm / lambda * (u * n);
            wRe[b + n * nBeams] = taper[n] * cos(ph);
            wIm[b + n * nBeams] = taper[n] * sin(ph);
            nrm += taper[n] * taper[n];
        }
        nrm = sqrt(nrm);
        for (int n = 0; n < nEl; n++) { wRe[b + n * nBeams] /= nrm; wIm[b + n * nBeams] /= nrm; }
    }
    free(taper);
}

/* ------------------------------------------------------------------------ */
/* The export (setup_radar.m + export_receiver_config.m)                     */
/* ------------------------------------------------------------------------ */

static void export_config(void)
{
    const char *names[2] = { "short", "long" };
    char p[96], q[96];
    const double c = num("const.c");
    const double fs = num("adc.fs_Hz");
    const double lambda = c / num("fc");

    /* CFAR window of the selected method (cfar_window) */
    const char *method = str("cfar.method");
    int cfarMethod, guardR, trainR, guardD, trainD;
    double pfa, thr_dB = NAN;
    const char *mkey;
    if (same_text(method, "2d")) {
        cfarMethod = 0; mkey = "cfar.twoD";
        snprintf(p, sizeof p, "%s.guardRange", mkey); guardR = (int)num(p);
        snprintf(p, sizeof p, "%s.trainRange", mkey); trainR = (int)num(p);
        snprintf(p, sizeof p, "%s.guardDoppler", mkey); guardD = (int)num(p);
        snprintf(p, sizeof p, "%s.trainDoppler", mkey); trainD = (int)num(p);
    } else if (same_text(method, "range")) {
        cfarMethod = 1; mkey = "cfar.range";
        snprintf(p, sizeof p, "%s.guard", mkey); guardR = (int)num(p);
        snprintf(p, sizeof p, "%s.train", mkey); trainR = (int)num(p);
        guardD = trainD = 0;
    } else if (same_text(method, "doppler")) {
        cfarMethod = 2; mkey = "cfar.doppler";
        guardR = trainR = 0;
        snprintf(p, sizeof p, "%s.guard", mkey); guardD = (int)num(p);
        snprintf(p, sizeof p, "%s.train", mkey); trainD = (int)num(p);
    } else {
        fail("cfg.cfar.method must be '2d', 'range' or 'doppler'.");
    }
    int cfarType = 0;
    if (has("cfar.type") && !empty("cfar.type")) {
        static const char *const types[4] = { "ca", "go", "so", "os" };
        const char *t = str("cfar.type");
        cfarType = -1;
        for (int i = 0; i < 4; i++) if (same_text(t, types[i])) cfarType = i;
        if (cfarType < 0) fail("cfg.cfar.type must be 'ca', 'go', 'so' or 'os'.");
    }
    const double osRank = has("cfar.osRank") && !empty("cfar.osRank") ? num("cfar.osRank") : 0.9;
    if (!(osRank > 0 && osRank <= 1)) fail("cfg.cfar.osRank must lie in (0, 1].");
    snprintf(p, sizeof p, "%s.pfa", mkey); pfa = num(p);
    snprintf(p, sizeof p, "%s.threshold_dB", mkey);
    if (!empty(p)) thr_dB = num(p);
    int snrCfar = 0;
    if (has("report.snrNoise")) {
        const char *s = str("report.snrNoise");
        if (same_text(s, "cfar")) snrCfar = 1;
        else if (!same_text(s, "median")) fail("cfg.report.snrNoise must be 'median' or 'cfar'.");
    }
    const int dopplerInterp = has("report.dopplerInterp") ? num("report.dopplerInterp") != 0 : 0;
    const int elevationInterp = has("report.elevationInterp") ? num("report.elevationInterp") != 0 : 1;

    /* antenna file: 3 dB azimuth beamwidth */
    const MadeFor mfAnt[2] = { { "fc", { num("fc") }, 1 }, { "spacing", { num("ant.elementSpacing_lambda") }, 1 } };
    int stale;
    MatFile ant = load_input_ex("antenna", mfAnt, 2, 1, &stale);
    double *azv, *gv;
    int nAz;
    if (stale) {
        build_antenna(&azv, &gv, &nAz);
    } else {
        const MatVar *az = mat_var(&ant, "az_deg", "antenna"), *gain = mat_var(&ant, "azGain_dB", "antenna");
        if (az->numel != gain->numel || az->numel < 3) fail("cfg.files.antenna: az_deg and azGain_dB differ in length");
        nAz = (int)az->numel;
        azv = xcalloc((size_t)nAz, sizeof(double));
        gv = xcalloc((size_t)nAz, sizeof(double));
        memcpy(azv, az->re, sizeof(double) * (size_t)nAz);
        memcpy(gv, gain->re, sizeof(double) * (size_t)nAz);
    }
    double azBeamwidth;
    {
        const int n = nAz;
        int iMax = 0;
        for (int i = 1; i < n; i++) if (gv[i] > gv[iMax]) iMax = i;
        const double gmax = gv[iMax], lim = -3.0103;
        int iL = -1, iR = -1;
        for (int i = 0; i <= iMax; i++) if (gv[i] - gmax < lim) iL = i;
        for (int i = iMax; i < n && iR < 0; i++) if (gv[i] - gmax < lim) iR = i;
        if (iL < 0 || iR < 1) fail("cfg.files.antenna: the azimuth pattern has no 3 dB points");
        const double g0 = gv[iL] - gmax, g1 = gv[iL + 1] - gmax;
        const double azL = azv[iL] + (lim - g0) * (azv[iL + 1] - azv[iL]) / (g1 - g0);
        const double h0 = gv[iR - 1] - gmax, h1 = gv[iR] - gmax;
        const double azR = azv[iR - 1] + (lim - h0) * (azv[iR] - azv[iR - 1]) / (h1 - h0);
        azBeamwidth = azR - azL;
    }
    free(azv);
    free(gv);
    mat_free(&ant);
    const double rate = num("ant.rpm") * 360.0 / 60.0;
    const int nCh = (int)num("ant.nElements");
    const double tilt = num("ant.tilt_deg");

    /* DDC filters and rates (ddc_decimation.m) */
    const double outRate = num("ddc.outputRate_Hz");
    const char *ddcKey = outRate == 2e6 ? "ddc" : outRate == 4e6 ? "ddc4MHz" : NULL;
    if (!ddcKey) fail("cfg.ddc.outputRate_Hz must be 2e6 or 4e6 (is %g).", outRate);
    Ddc ddc;
    memset(&ddc, 0, sizeof ddc);
    ddc.fs = fs;
    ddc.D1 = (int)num("ddc.stage1Decimation");
    const double D2 = fs / ddc.D1 / outRate;
    if (D2 < 1 || fabs(D2 - round(D2)) > 1e-9)
        fail("The stage-1 rate %g MHz is not a whole multiple of the DDC output rate %g MHz.", fs / ddc.D1 / 1e6, outRate / 1e6);
    ddc.D2 = (int)lround(D2);
    const double fracBits = num("ddc.coefFracBits");
    const MadeFor mfDdc[3] = { { "fs_Hz", { fs }, 1 }, { "decimation", { ddc.D1, ddc.D2 }, 2 },
                               { "coefFracBits", { fracBits }, 1 } };
    MatFile ddcF = load_input_ex(ddcKey, mfDdc, 3, 1, &stale);
    const double qScale = pow(2.0, fracBits);
    if (stale) {                                          /* make_input_files.m: 70 and 66 taps */
        ddc.L1 = 70;
        ddc.L2 = 66;
        ddc.h1 = xcalloc((size_t)ddc.L1, sizeof(double));
        ddc.h2 = xcalloc((size_t)ddc.L2, sizeof(double));
        build_lowpass(70, 5.8e6 / fs, qScale, ddc.h1);
        build_lowpass(66, (outRate == 2e6 ? 1.0e6 : 1.5e6) / (fs / ddc.D1), qScale, ddc.h2);
    } else {
        const MatVar *h1 = mat_var(&ddcF, "h1", ddcKey), *h2 = mat_var(&ddcF, "h2", ddcKey);
        ddc.L1 = (int)h1->numel;
        ddc.L2 = (int)h2->numel;
        ddc.h1 = xcalloc((size_t)ddc.L1, sizeof(double));
        ddc.h2 = xcalloc((size_t)ddc.L2, sizeof(double));
        memcpy(ddc.h1, h1->re, sizeof(double) * (size_t)ddc.L1);
        memcpy(ddc.h2, h2->re, sizeof(double) * (size_t)ddc.L2);
    }
    for (int i = 0; i < ddc.L1; i++) ddc.h1[i] /= qScale;
    for (int i = 0; i < ddc.L2; i++) ddc.h2[i] /= qScale;
    mat_free(&ddcF);
    const double fs1 = fs / ddc.D1;
    ddc.fsOut = fs1 / ddc.D2;
    const double fsOut = ddc.fsOut;

    /* PRI sequence and MTI weights */
    const char *prfMode = str("prf.mode");
    int nStagger, staggered, nTaps;
    double *pri, *mti;
    if (same_text(prfMode, "fixed")) {
        staggered = 0;
        nStagger = 1;
        pri = xcalloc(1, sizeof(double));
        pri[0] = num("prf.pri_s");
        nTaps = 3;
        mti = xcalloc(3, sizeof(double));
        mti[0] = 1; mti[1] = -2; mti[2] = 1;
    } else if (same_text(prfMode, "staggered")) {
        staggered = 1;
        MatFile st = load_input("stagger", NULL, 0);
        const MatVar *pv = mat_var(&st, "priVector_s", "stagger"), *mv = mat_var(&st, "mtiCoeffs", "stagger");
        nStagger = (int)pv->numel;
        pri = xcalloc((size_t)nStagger, sizeof(double));
        memcpy(pri, pv->re, sizeof(double) * (size_t)nStagger);
        nTaps = (int)mv->cols;
        if ((int)mv->rows != 1 && (int)mv->rows != nStagger) fail("cfg.files.stagger: mtiCoeffs needs 1 or %d rows", nStagger);
        mti = xcalloc((size_t)(nStagger * nTaps), sizeof(double));      /* row-major [nStagger x nTaps] */
        for (int k = 0; k < nStagger; k++)
            for (int l = 0; l < nTaps; l++) mti[k * nTaps + l] = mv->re[(mv->rows == 1 ? 0 : k) + l * mv->rows];
        mat_free(&st);
    } else {
        fail("Unknown PRF mode \"%s\"", prfMode);
    }
    int mtiFir = 0;
    if (has("mti.type")) {
        const char *t = str("mti.type");
        if (same_text(t, "fir")) mtiFir = 1;
        else if (!same_text(t, "canceller")) fail("cfg.mti.type must be 'canceller' or 'fir'.");
    }
    if (mtiFir) {
        if (!has("files.mtiFir")) fail("cfg.mti.type = 'fir' needs the coefficient file cfg.files.mtiFir.");
        MatFile fm = load_input("mtiFir", NULL, 0);
        const MatVar *h = mat_get(&fm, "coeffs");
        if (!h || !h->re || h->numel == 0)
            fail("%s must hold the FIR coefficients in the variable 'coeffs'.", str("files.mtiFir"));
        if (h->im) fail("The FIR coefficients of %s must be real and finite.", str("files.mtiFir"));
        int rows = (int)h->rows, cols = (int)h->cols;
        if (rows == 1 || cols == 1) { cols = (int)h->numel; rows = 1; }
        if (rows != 1 && rows != nStagger)
            fail("%s: 'coeffs' has %d rows: 1 (one filter) or %d (one per stagger position).", str("files.mtiFir"), rows, nStagger);
        int nz = 0;
        for (size_t i = 0; i < h->numel; i++) {
            if (!isfinite(h->re[i])) fail("The FIR coefficients of %s must be real and finite.", str("files.mtiFir"));
            nz |= h->re[i] != 0;
        }
        if (!nz) fail("The FIR coefficients of %s are all zero.", str("files.mtiFir"));
        free(mti);
        nTaps = cols;
        mti = xcalloc((size_t)(nStagger * nTaps), sizeof(double));
        for (int k = 0; k < nStagger; k++)
            for (int l = 0; l < nTaps; l++) mti[k * nTaps + l] = h->re[(rows == 1 ? 0 : k) + (size_t)l * (rows == 1 ? 1 : rows)];
        mat_free(&fm);
    }
    double meanPri = 0, minPri = pri[0], sumPri = 0;
    for (int i = 0; i < nStagger; i++) { sumPri += pri[i]; minPri = pri[i] < minPri ? pri[i] : minPri; }
    meanPri = sumPri / nStagger;

    /* waveforms */
    Wave wv[2];
    memset(wv, 0, sizeof wv);
    MatFile codes;
    int haveCodes = 0;
    for (int i = 0; i < 2; i++) {
        Wave *w = &wv[i];
        w->name = names[i];
        snprintf(p, sizeof p, "wave.%s.type", names[i]);
        const char *type = str(p);
        if (same_text(type, "phasecode")) w->phaseCode = 1;
        else if (!same_text(type, "lfm")) fail("Unknown waveform type \"%s\"", type);
        snprintf(p, sizeof p, "wave.%s.width_s", names[i]);       w->width_s = num(p);
        snprintf(p, sizeof p, "wave.%s.bandwidth_Hz", names[i]);  w->bandwidth_Hz = num_or(p, 0);
        snprintf(p, sizeof p, "wave.%s.chip_s", names[i]);        w->chip_s = num_or(p, 1e-6);
        snprintf(p, sizeof p, "wave.%s.tTx_s", names[i]);         w->tTx_s = num(p);
        snprintf(p, sizeof p, "wave.%s.enable", names[i]);       w->enable = num_or(p, 1) != 0;
        snprintf(p, sizeof p, "wave.%s.fc_Hz", names[i]);        w->fc_Hz = num_or(p, num("fc"));
        if (!(w->fc_Hz > 0) || !isfinite(w->fc_Hz))
            fail("cfg.wave.%s.fc_Hz must be a positive carrier frequency [Hz] (or [] = cfg.fc).", names[i]);
        w->lambda_m = c / w->fc_Hz;
        {
            const double fIFp = num("adc.fIF_Hz") + (w->fc_Hz - num("fc"));
            const double bwHalf = 0.5 * (w->phaseCode ? 1.0 / w->chip_s : w->bandwidth_Hz);
            if (w->enable && (fIFp - bwHalf <= 0 || fIFp + bwHalf >= fs / 2))
                fail("cfg.wave.%s.fc_Hz = %.3f MHz puts the sub-pulse on the IF %.3f MHz: it must lie inside "
                     "(0, %.1f MHz) with its bandwidth (LO = cfg.fc - cfg.adc.fIF_Hz = %.3f MHz).", names[i],
                     w->fc_Hz / 1e6, fIFp / 1e6, fs / 2e6, (num("fc") - num("adc.fIF_Hz")) / 1e6);
        }
        snprintf(p, sizeof p, "wave.%s.compression", names[i]);
        const char *comp = has(p) ? str(p) : "matched";
        if (same_text(comp, "mismatch")) w->mismatch = 1;
        else if (!same_text(comp, "matched")) fail("cfg.wave.%s.compression must be 'matched' or 'mismatch'.", names[i]);
    }
    if (wv[0].phaseCode || wv[1].phaseCode) {
        MadeFor mf[4];
        for (int i = 0; i < 2; i++) {
            snprintf(p, sizeof p, "wave.%s.mismatchLength", names[i]);
            const int hasLen = has(p) && !empty(p);
            mf[2 * i] = (MadeFor){ i ? "longChips" : "shortChips", { round(wv[i].width_s / wv[i].chip_s) }, 1 };
            mf[2 * i + 1] = (MadeFor){ i ? "longMismatchLength" : "shortMismatchLength", { hasLen ? num(p) : 0 }, hasLen };
        }
        codes = load_input("phaseCode", mf, 4);
        haveCodes = 1;
        for (int i = 0; i < 2; i++) {
            Wave *w = &wv[i];
            if (!w->phaseCode) continue;
            snprintf(q, sizeof q, "%sCode", names[i]);
            const MatVar *cv = mat_var(&codes, q, "phaseCode");
            w->nCode = (int)cv->numel;
            if (w->nCode != (int)lround(w->width_s / w->chip_s))
                fail("The %s code of %s has %d chips, but width_s/chip_s = %d.", names[i], str("files.phaseCode"),
                     w->nCode, (int)lround(w->width_s / w->chip_s));
            w->codeRe = xcalloc((size_t)w->nCode, sizeof(double));
            w->codeIm = xcalloc((size_t)w->nCode, sizeof(double));
            for (int k = 0; k < w->nCode; k++) { w->codeRe[k] = cv->re[k]; w->codeIm[k] = cv->im ? cv->im[k] : 0; }
            if (w->mismatch) {
                snprintf(q, sizeof q, "%sMismatch", names[i]);
                const MatVar *mv = mat_var(&codes, q, "phaseCode");
                w->nMis = (int)mv->numel;
                w->misRe = xcalloc((size_t)w->nMis, sizeof(double));
                w->misIm = xcalloc((size_t)w->nMis, sizeof(double));
                for (int k = 0; k < w->nMis; k++) { w->misRe[k] = mv->re[k]; w->misIm[k] = mv->im ? mv->im[k] : 0; }
                w->misRate_Hz = 1.0 / w->chip_s;
                snprintf(q, sizeof q, "%sMismatchRate_Hz", names[i]);
                const MatVar *rv = mat_get(&codes, q);
                if (rv && rv->re && rv->numel >= 1) w->misRate_Hz = rv->re[0];
            }
        }
    }
    const int taylorNbar = (int)num("proc.taylorNbar");
    const double taylorSll = num("proc.taylorSll_dB");
    for (int i = 0; i < 2; i++) compression_filter(&wv[i], &ddc, taylorNbar, taylorSll);
    if (haveCodes) mat_free(&codes);
    Wave *ws = &wv[0], *wl = &wv[1];
    if (!ws->enable && !wl->enable) fail("Both sub-pulses are switched off (cfg.wave.short.enable, cfg.wave.long.enable).");
    if (ws->enable && wl->enable && wl->tTx_s < ws->tTx_s + ws->width_s)
        fail("The long pulse must start after the end of the short pulse.");

    /* PRI timeline, DDC sample times */
    const double maxRange = num("rx.maxRange_m");
    const double filterSpan = (ddc.L1 - 1) / fs + (ddc.L2 - 1) / fs1;
    const Wave *wLast = wl->enable ? wl : ws;             /* the last transmitted sub-pulse ends the record */
    const double rxEnd = wLast->tTx_s + wLast->width_s + 2 * maxRange / c + filterSpan + 2 / fsOut;
    if (rxEnd > minPri)
        fail("Receive window (%.1f us) does not fit in the shortest PRI (%.1f us): reduce cfg.rx.maxRange_m or "
             "start the long pulse earlier.", rxEnd * 1e6, minPri * 1e6);
    const long nAdc = lround((rxEnd - 0.0) * fs);
    const long n1 = (nAdc - ddc.L1) / ddc.D1 + 1;
    const int nDdc = (int)((n1 - ddc.L2) / ddc.D2 + 1);
    double *tDdc = xcalloc((size_t)nDdc, sizeof(double));
    for (int k = 0; k < nDdc; k++) {
        const double t1 = (double)((long)k * ddc.D2 * ddc.D1) / fs + (ddc.L1 - 1) / 2.0 / fs;
        tDdc[k] = t1 + (ddc.L2 - 1) / 2.0 / fs1;
    }

    /* CFAR window size check of the receive intervals */
    for (int i = 0; i < 2; i++) {
        Wave *w = &wv[i];
        if (!w->enable) continue;                         /* switched off: no receive interval, no cells */
        const double iv0 = i == 0 ? ws->tTx_s + ws->width_s : wl->tTx_s + wl->width_s;
        const double iv1 = i == 0 && wl->enable ? wl->tTx_s : INFINITY;
        int first = -1, count = 0;
        for (int k = 0; k < nDdc; k++) {
            const double tFirst = (double)((long)k * ddc.D2 * ddc.D1) / fs;
            const double tLast = (double)(((long)k * ddc.D2 + ddc.L2 - 1) * ddc.D1 + ddc.L1 - 1) / fs;
            if (tFirst >= iv0 && tLast < iv1) {
                if (first < 0) first = k;
                else if (k != first + count) fail("internal: the %s-pulse receive interval is not contiguous", w->name);
                count++;
            }
        }
        const int Lp = (int)lround(w->width_s * fsOut);
        w->rxFirst = first + 1;
        w->rxCount = count;
        w->nRange = count - Lp + 1;
        if (first < 0 || w->nRange < 2 * (guardR + trainR) + 1)
            fail("The %s-pulse receive interval is too short for its pulse and the CFAR window: move the long "
                 "pulse later or increase cfg.rx.maxRange_m.", w->name);
        w->ddcFirst1 = w->rxFirst - w->mfShift;
        w->ddcCount = w->nRange + w->nTaps - 1;
        w->delay_s = xcalloc((size_t)w->nRange, sizeof(double));
        w->range_m = xcalloc((size_t)w->nRange, sizeof(double));
        for (int k = 0; k < w->nRange; k++) {
            w->delay_s[k] = tDdc[first + k] - 0.5 / fsOut - w->tTx_s;
            w->range_m[k] = c * w->delay_s[k] / 2;
        }
    }
    /* handover */
    const int oR = guardR + trainR;
    const double handover = wl->enable ? wl->range_m[0] : INFINITY;
    if (ws->enable && wl->enable && ws->range_m[ws->nRange - 1] + fabs(ws->range_m[1] - ws->range_m[0]) < handover)
        warn("Gap in the range coverage: short pulse up to %.1f km, long pulse from %.1f km.",
             ws->range_m[ws->nRange - 1] / 1e3, handover / 1e3);
    {
        int last = 0;                                     /* 1-based, 0 = none */
        for (int k = 0; k < ws->nRange; k++) if (ws->range_m[k] < handover) last = k + 1;
        int nKeep = last + oR > 2 * oR + 1 ? last + oR : 2 * oR + 1;
        if (ws->enable && wl->enable && nKeep < ws->nRange) {
            ws->nRange = nKeep;
            ws->ddcCount = nKeep + ws->nTaps - 1;
        }
    }
    for (int i = 0; i < 2; i++) {                        /* CFAR-tested cells */
        Wave *w = &wv[i];
        if (!w->enable) continue;
        snprintf(p, sizeof p, "cfar.%s.minRange_m", names[i]);
        const double lo = num(p);
        snprintf(p, sizeof p, "cfar.%s.maxRange_m", names[i]);
        const double hi = num(p);
        int first = 0, last = 0;
        for (int k = 0; k < w->nRange; k++) {
            int ok = w->range_m[k] >= lo && w->range_m[k] <= hi;
            if (i == 0) ok = ok && w->range_m[k] < handover;
            if (ok) { if (!first) first = k + 1; last = k + 1; }
        }
        if (!first) {
            warn("No %s-pulse range cell is tested: check cfg.cfar.%s.minRange_m / maxRange_m.", names[i], names[i]);
            first = 1;
            last = 0;
        }
        w->testFirst = first;
        w->testLast = last;
    }

    /* DBF */
    const int dbf = nCh > 1;
    int nBeams = 1;
    double *wRe = NULL, *wIm = NULL, *beamU = NULL, *beamElKeep = NULL;
    int wCplx = 0;
    if (dbf) {
        const MadeFor mfDbf[4] = { { "fc", { num("fc") }, 1 }, { "nElements", { nCh }, 1 }, { "tilt_deg", { tilt }, 1 },
                                   { "spacing", { num("ant.elementSpacing_lambda") }, 1 } };
        MatFile df = load_input_ex("dbf", mfDbf, 4, 1, &stale);
        double *beamEl;
        if (stale) {
            nBeams = 21;
            wCplx = 1;
            wRe = xcalloc((size_t)(nBeams * nCh), sizeof(double));
            wIm = xcalloc((size_t)(nBeams * nCh), sizeof(double));
            beamEl = xcalloc((size_t)nBeams, sizeof(double));
            build_dbf(nCh, num("ant.elementSpacing_lambda"), tilt, lambda, wRe, wIm, beamEl);
        } else {
            const MatVar *W = mat_var(&df, "W", "dbf"), *el = mat_var(&df, "beamEl_deg", "dbf");
            if ((int)W->cols != nCh)
                fail("The DBF matrix has %d columns but cfg.ant.nElements = %d.", (int)W->cols, nCh);
            nBeams = (int)el->numel;
            if ((int)W->rows != nBeams) fail("cfg.files.dbf: W has %d rows but beamEl_deg %d values", (int)W->rows, nBeams);
            wCplx = W->im != NULL;
            wRe = xcalloc(W->numel, sizeof(double));
            wIm = xcalloc(W->numel, sizeof(double));
            memcpy(wRe, W->re, sizeof(double) * W->numel);
            if (W->im) memcpy(wIm, W->im, sizeof(double) * W->numel);
            beamEl = xcalloc((size_t)nBeams, sizeof(double));
            memcpy(beamEl, el->re, sizeof(double) * (size_t)nBeams);
        }
        beamU = xcalloc((size_t)nBeams, sizeof(double));
        for (int b = 0; b < nBeams; b++) beamU[b] = sind(beamEl[b] - tilt);
        beamElKeep = beamEl;
        mat_free(&df);
    } else {
        beamU = xcalloc(1, sizeof(double));
        beamU[0] = NAN;
    }
    for (int i = 0; i < 2; i++) {                        /* beam directions on the carrier of each sub-pulse */
        wv[i].beamU = xcalloc((size_t)nBeams, sizeof(double));
        for (int b = 0; b < nBeams; b++) {
            double u = beamU[b] * (num("fc") / wv[i].fc_Hz);
            if (fabs(u) > 1) u = u > 0 ? 1.0 : -1.0;
            wv[i].beamU[b] = u;
        }
    }

    /* CPI and Doppler processing */
    const double cpiTime = azBeamwidth / rate;
    /* cfg.proc.cpiMode: 'auto' (beam dwell) or 'manual' (cfg.proc.cpiLength); without the field
     * (older configurations) cpiLength = [] means auto */
    int cpiAuto = !has("proc.cpiLength") || empty("proc.cpiLength");
    if (has("proc.cpiMode") && !empty("proc.cpiMode")) {
        const char *m = str("proc.cpiMode");
        if (same_text(m, "auto")) cpiAuto = 1;
        else if (same_text(m, "manual")) cpiAuto = 0;
        else fail("cfg.proc.cpiMode must be 'auto' or 'manual'.");
    }
    int cpiLength;
    if (cpiAuto) cpiLength = (int)lround(cpiTime / meanPri);
    else {
        if (!has("proc.cpiLength") || empty("proc.cpiLength")) fail("cfg.proc.cpiMode = 'manual' needs cfg.proc.cpiLength: a whole number of pulses.");
        const double v = num("proc.cpiLength");
        if (v < 1 || v != floor(v)) fail("cfg.proc.cpiMode = 'manual' needs cfg.proc.cpiLength: a whole number of pulses.");
        cpiLength = (int)v;
    }
    if (staggered) cpiLength = nStagger * (cpiLength / nStagger);
    const int nOut = cpiLength - (nTaps - 1);
    int nRows, rowLen;
    double rowPri;
    if (staggered) { nRows = nStagger; rowLen = nOut >= 0 ? nOut / nRows : -1; rowPri = sumPri; }
    else { nRows = 1; rowLen = nOut; rowPri = pri[0]; }
    int nfft;
    if (!empty("proc.nfft")) nfft = (int)num("proc.nfft");
    else nfft = staggered ? 1 << nextpow2(rowLen) : 1 << nextpow2(cpiLength);
    if (rowLen > nfft) rowLen = nfft;
    if (rowLen < 2)
        fail("The CPI of %d pulses is too short for the MTI and the Doppler processing (cfg.proc.cpiLength).", cpiLength);
    double *win = xcalloc((size_t)rowLen, sizeof(double)), wn = 0;
    for (int j = 0; j < rowLen; j++) {
        win[j] = rowLen == 1 ? 1.0 : 0.54 - 0.46 * cos(2 * PI * j / (rowLen - 1));
        wn += win[j] * win[j];
    }
    wn = sqrt(wn);
    for (int j = 0; j < rowLen; j++) win[j] /= wn;
    double *fHz = xcalloc((size_t)nfft, sizeof(double));
    for (int k = 0; k < nfft; k++) fHz[k] = (k - nfft / 2) / (nfft * rowPri);

    /* clutter filter only up to cfg.mti.maxRange_m: the cells beyond get the bypass weights (the last tap) */
    double mtiMaxRange = INFINITY;
    if (has("mti.maxRange_m") && !empty("mti.maxRange_m")) mtiMaxRange = num("mti.maxRange_m");
    if (isnan(mtiMaxRange) || mtiMaxRange < 0)
        fail("cfg.mti.maxRange_m must be a non-negative scalar range in metres (Inf is allowed).");
    int bypassUsed = 0;
    for (int i = 0; i < 2; i++)
        for (int k = 0; k < wv[i].nRange; k++) bypassUsed |= wv[i].range_m[k] > mtiMaxRange;
    double *bypass = xcalloc((size_t)(nStagger * nTaps), sizeof(double));
    for (int k = 0; k < nStagger; k++) bypass[k * nTaps + nTaps - 1] = 1.0;

    /* noise gain and looks of every Doppler bin: white noise (impulses) through MTI, window, FFT
     * (pass 1: the same without the clutter filter, for the cells beyond mtiMaxRange) */
    double *noiseGain = xcalloc((size_t)(2 * nfft), sizeof(double));
    double *noiseGainB = noiseGain + nfft;
    double *Cre = xcalloc((size_t)(2 * nfft * nRows * nRows), sizeof(double));
    double *Cim = xcalloc((size_t)(2 * nfft * nRows * nRows), sizeof(double));
    for (int pass = 0; pass < 1 + bypassUsed; pass++) {
        const double *co = pass ? bypass : mti;
        double *gainP = pass ? noiseGainB : noiseGain;
        double *CreP = Cre + (size_t)pass * nfft * nRows * nRows, *CimP = Cim + (size_t)pass * nfft * nRows * nRows;
        double *y = xcalloc((size_t)(nOut > 0 ? nOut : 1), sizeof(double));
        double *Yr = xcalloc((size_t)(nRows * nfft), sizeof(double)), *Yi = xcalloc((size_t)(nRows * nfft), sizeof(double));
        double *cosT = xcalloc((size_t)nfft, sizeof(double)), *sinT = xcalloc((size_t)nfft, sizeof(double));
        for (int k = 0; k < nfft; k++) { cosT[k] = cos(2 * PI * k / nfft); sinT[k] = sin(2 * PI * k / nfft); }
        const int used = nRows * rowLen, off = nOut - used;
        for (int imp = 0; imp < cpiLength; imp++) {
            for (int m = 0; m < nOut; m++) {
                const int l = imp - m, st = m % nStagger;
                y[m] = (l >= 0 && l < nTaps) ? co[st * nTaps + l] : 0.0;
            }
            for (int r = 0; r < nRows; r++)
                for (int i = 0; i < nfft; i++) {
                    const int kf = ((i - nfft / 2) % nfft + nfft) % nfft;
                    double sr = 0, si = 0;
                    for (int j = 0; j < rowLen; j++) {
                        const int qq = staggered ? j * nRows + r : r * rowLen + j;
                        const double x = y[off + qq] * win[j];
                        if (x == 0) continue;
                        const int ph = (int)(((long)kf * j) % nfft);
                        sr += x * cosT[ph];
                        si -= x * sinT[ph];
                    }
                    Yr[r * nfft + i] = sr;
                    Yi[r * nfft + i] = si;
                }
            for (int i = 0; i < nfft; i++) {
                double pw = 0;
                for (int r = 0; r < nRows; r++) pw += Yr[r * nfft + i] * Yr[r * nfft + i] + Yi[r * nfft + i] * Yi[r * nfft + i];
                gainP[i] += pw / nRows;
                for (int r1 = 0; r1 < nRows; r1++)            /* C = Yk' * Yk */
                    for (int r2 = 0; r2 < nRows; r2++) {
                        const double ar = Yr[r1 * nfft + i], ai = -Yi[r1 * nfft + i];
                        const double br = Yr[r2 * nfft + i], bi = Yi[r2 * nfft + i];
                        CreP[(i * nRows + r1) * nRows + r2] += ar * br - ai * bi;
                        CimP[(i * nRows + r1) * nRows + r2] += ar * bi + ai * br;
                    }
            }
        }
        free(y); free(Yr); free(Yi); free(cosT); free(sinT);
    }
    if (!bypassUsed) memcpy(noiseGainB, noiseGain, sizeof(double) * (size_t)nfft);
    free(bypass);
    const double notchWidth = num("cfar.notchWidth");
    const double notchMaxRange_m = num("cfar.notchMaxRange_m");
    const int notchBins = (int)lround(notchWidth * 2 * nfft / rowLen);
    const int zeroBin = nfft / 2 + 1;
    double nLooks = INFINITY;
    int anyTested = 0;
    {
        const int m = 2 * nRows;
        double *A = xcalloc((size_t)(m * m), sizeof(double));
        for (int ii = 0; ii < nfft * (1 + bypassUsed); ii++) {
            const int i = ii % nfft;
            const double *CreP = Cre + (size_t)(ii / nfft) * nfft * nRows * nRows;
            const double *CimP = Cim + (size_t)(ii / nfft) * nfft * nRows * nRows;
            if (notchMaxRange_m >= maxRange && !(abs(i + 1 - zeroBin) > notchBins)) continue;
            anyTested = 1;
            double tr = 0;
            for (int r1 = 0; r1 < nRows; r1++) {
                tr += CreP[(i * nRows + r1) * nRows + r1];
                for (int r2 = 0; r2 < nRows; r2++) {       /* [Re -Im; Im Re]: the eigenvalues of C, twice */
                    const double re = CreP[(i * nRows + r1) * nRows + r2], im = CimP[(i * nRows + r1) * nRows + r2];
                    A[r1 * m + r2] = re;
                    A[(r1 + nRows) * m + (r2 + nRows)] = re;
                    A[r1 * m + (r2 + nRows)] = -im;
                    A[(r1 + nRows) * m + r2] = im;
                }
            }
            const double lmax = nRows == 1 ? tr : max_eig_sym(A, m);
            const double looks = tr / lmax;
            if (looks < nLooks) nLooks = looks;
        }
        free(A);
    }
    if (!anyTested)
        fail("The MTI notch (cfg.cfar.notchWidth = %g: +-%d bins) covers every Doppler bin (FFT length %d).",
             notchWidth, notchBins, nfft);
    MedianCtx mc = { nLooks };
    const double medianToMean = fzero_brent(median_f, &mc, 0.1, 2.0);

    /* CFAR threshold */
    const int nOuter = (2 * (guardR + trainR) + 1) * (2 * (guardD + trainD) + 1);
    const int nGuard = (2 * guardR + 1) * (2 * guardD + 1);
    const int nTrain = nOuter - nGuard;
    Detector det = { cfarType, nTrain, 0, 0, nLooks };
    det.nHalf = cfarMethod == 2 ? trainD
              : (guardR + trainR) * (2 * (guardD + trainD) + 1) - guardR * (2 * guardD + 1);
    {
        const double k = round(osRank * nTrain);
        det.osK = (int)(k < 1 ? 1 : k > nTrain ? nTrain : k);
    }
    if ((cfarType == 1 || cfarType == 2) && det.nHalf < 1)
        fail("cfg.cfar.type = '%s' needs training cells on both sides of the cell (cfg.cfar.%s).",
             cfarType == 1 ? "go" : "so", method);
    cfar_quadrature(nLooks);
    const double alpha = isnan(thr_dB) ? cfar_alpha(pfa, &det) : pow(10.0, thr_dB / 10.0);
    double noiseScale = 1;                               /* mean of the statistic for noise of mean 1 */
    if (cfarType != 0) {
        noiseScale = 0;
        for (int i = 0; i < g_nq; i++) noiseScale += g_qw[i] * (1 - cfar_statistic(g_qs[i], &det));
    }
    if (2 * (guardD + trainD) + 1 > nfft)
        fail("CFAR Doppler window (%d bins) is larger than the FFT length (%d).", 2 * (guardD + trainD) + 1, nfft);

    /* tracker (track_setup.m) */
    const int trackEnable = num("track.enable") != 0;
    const char *dimS = str("track.dim");
    int dim;
    if (same_text(dimS, "auto")) dim = 2 + dbf;
    else if (same_text(dimS, "2d")) dim = 2;
    else if (same_text(dimS, "3d")) {
        dim = 3;
        if (!dbf) { warn("cfg.track.dim = '3d' needs the elevation of the DBF (1 element: no DBF): tracking in 2-D."); dim = 2; }
    } else fail("cfg.track.dim must be 'auto', '2d' or '3d'.");
    const char *model = str("track.model");
    int ca;
    double qv[3];
    if (same_text(model, "cv")) {
        ca = 0;
        qv[0] = qv[1] = pow(num("track.sigmaAccel_mps2"), 2);
        qv[2] = pow(num("track.sigmaAccelZ_mps2"), 2);
    } else if (same_text(model, "ca")) {
        ca = 1;
        qv[0] = qv[1] = pow(num("track.sigmaJerk_mps3"), 2);
        qv[2] = pow(num("track.sigmaJerkZ_mps3"), 2);
    } else fail("cfg.track.model must be 'cv' or 'ca'.");
    const double sig[3] = { num("track.sigmaRange_m"), PI / 180 * num("track.sigmaAz_deg"), PI / 180 * num("track.sigmaEl_deg") };
    const double maxSpeed = num("track.maxSpeed_mps");
    const double gateProb = num("track.gateProb");
    double sigmaVr = empty("track.sigmaVr_mps") ? 0.4 * fabs((-fHz[1] * lambda / 2) - (-fHz[0] * lambda / 2))
                                                : num("track.sigmaVr_mps");
    const double Rvr = sigmaVr * sigmaVr;
    /* the plot velocity is known modulo the Doppler span on the carrier of its sub-pulse */
    const double vrSpanPulse[2] = { wv[0].lambda_m / (2 * rowPri), wv[1].lambda_m / (2 * rowPri) };
    const double vrSpan = !ws->enable ? vrSpanPulse[1] : !wl->enable ? vrSpanPulse[0]
                        : vrSpanPulse[0] < vrSpanPulse[1] ? vrSpanPulse[0] : vrSpanPulse[1];
    double minSpeed = 0;
    if (has("track.minSpeedEnable") && num("track.minSpeedEnable") != 0) {
        minSpeed = num("track.minSpeed_mps");
        if (!(minSpeed > 0 && isfinite(minSpeed)))
            fail("cfg.track.minSpeed_mps must be a positive speed [m/s] when cfg.track.minSpeedEnable is true.");
    }

    /* ---------------------------------------------------------------- */
    /* records, in the order of export_receiver_config.m                */
    /* ---------------------------------------------------------------- */
    rec_scalar("nCh", nCh);
    rec_scalar("nBeams", nBeams);
    rec_scalar("dbfEnabled", dbf);
    if (dbf) {
        Rec *r = add_rec("W", nBeams, nCh, wCplx);
        memcpy(r->re, wRe, sizeof(double) * (size_t)(nBeams * nCh));
        if (wCplx) memcpy(r->im, wIm, sizeof(double) * (size_t)(nBeams * nCh));
    }
    rec_vec("beamU", beamU, nBeams, 1);
    rec_scalar("tilt_deg", tilt);
    rec_scalar("height_m", num("ant.height_m"));
    rec_scalar("rate_dps", rate);
    rec_scalar("azBeamwidth_deg", azBeamwidth);
    rec_scalar("maxRange_m", maxRange);
    rec_scalar("nDdc", nDdc);
    rec_scalar("fsOut_Hz", fsOut);
    for (int i = 0; i < 2; i++) {
        Wave *w = &wv[i];
        char nm[32];
        const char *n = names[i];
        pfx(nm, n, "_enable");     rec_scalar(nm, w->enable);
        pfx(nm, n, "_fc_Hz");      rec_scalar(nm, w->fc_Hz);
        pfx(nm, n, "_lambda_m");   rec_scalar(nm, w->lambda_m);
        if (!w->enable) continue;
        pfx(nm, n, "_beamU");      rec_vec(nm, w->beamU, nBeams, 1);
        pfx(nm, n, "_rxFirst");    rec_scalar(nm, w->rxFirst - 1);
        pfx(nm, n, "_rxCount");    rec_scalar(nm, w->rxCount);
        pfx(nm, n, "_ddcFirst");   rec_scalar(nm, w->ddcFirst1 - 1);
        pfx(nm, n, "_ddcCount");   rec_scalar(nm, w->ddcCount);
        pfx(nm, n, "_nTaps");      rec_scalar(nm, w->nTaps);
        pfx(nm, n, "_nRange");     rec_scalar(nm, w->nRange);
        pfx(nm, n, "_taps");       rec_cvec(nm, w->tapsRe, w->tapsIm, w->nTaps);
        pfx(nm, n, "_range_m");    rec_vec(nm, w->range_m, w->nRange, 0);
        pfx(nm, n, "_delay_s");    rec_vec(nm, w->delay_s, w->nRange, 0);
        pfx(nm, n, "_testFirst");  rec_scalar(nm, w->testFirst - 1);
        pfx(nm, n, "_testCount");  rec_scalar(nm, w->testLast - w->testFirst + 1 > 0 ? w->testLast - w->testFirst + 1 : 0);
        pfx(nm, n, "_waveform");   rec_scalar(nm, w->phaseCode);
        pfx(nm, n, "_compression"); rec_scalar(nm, w->mismatch);
        pfx(nm, n, "_nChips");     rec_scalar(nm, w->phaseCode * round(w->width_s / w->chip_s));
        pfx(nm, n, "_mfShift");    rec_scalar(nm, w->mfShift);
        const int np = (int)lround(w->width_s * fsOut);
        double *pr = xcalloc((size_t)np, sizeof(double)), *pi = xcalloc((size_t)np, sizeof(double));
        for (int k = 0; k < np; k++) envelope(w, (k + 0.5) / fsOut, &pr[k], &pi[k]);
        pfx(nm, n, "_pulse");      rec_cvec(nm, pr, pi, np);
        free(pr);
        free(pi);
    }
    rec_scalar("cpiLength", cpiLength);
    rec_scalar("nStagger", nStagger);
    rec_scalar("staggered", staggered);
    {
        Rec *r = add_rec("mtiCoeffs", nStagger, nTaps, 0);    /* column-major */
        for (int k = 0; k < nStagger; k++)
            for (int l = 0; l < nTaps; l++) r->re[k + l * nStagger] = mti[k * nTaps + l];
    }
    rec_scalar("mtiType", mtiFir);
    rec_scalar("nRows", nRows);
    rec_scalar("rowLen", rowLen);
    rec_scalar("nfft", nfft);
    rec_vec("dopWindow", win, rowLen, 1);
    rec_vec("noiseGain", noiseGain, nfft, 0);
    rec_scalar("mtiMaxRange_m", mtiMaxRange);
    rec_vec("noiseGainBypass", noiseGainB, nfft, 0);
    rec_vec("dopFreq_Hz", fHz, nfft, 0);
    rec_scalar("notchBins", notchBins);
    rec_scalar("notchMaxRange_m", notchMaxRange_m);
    rec_scalar("medianToMean", medianToMean);
    rec_scalar("cfarMethod", cfarMethod);
    rec_scalar("guardRange", guardR);
    rec_scalar("trainRange", trainR);
    rec_scalar("guardDoppler", guardD);
    rec_scalar("trainDoppler", trainD);
    rec_scalar("nTrain", nTrain);
    rec_scalar("alpha", alpha);
    rec_scalar("cfarType", cfarType);
    rec_scalar("cfarNHalf", det.nHalf);
    rec_scalar("cfarOsK", det.osK);
    rec_scalar("cfarNoiseScale", noiseScale);
    rec_scalar("nLooks", nLooks);
    rec_scalar("sidelobeMargin_dB", num("prune.sidelobeMargin_dB"));
    rec_scalar("snrNoise", snrCfar);
    rec_scalar("dopplerInterp", dopplerInterp);
    rec_scalar("elevationInterp", elevationInterp);
    const double updateEvery = num("proc.updateEvery");
    const double minExtent = num("extract.minExtent_deg");
    const double hitsF = ceil(minExtent / (rate * meanPri * updateEvery));
    rec_scalar("minHits", hitsF > 1 ? hitsF : 1);
    rec_scalar("maxMiss", num("extract.maxMiss"));
    rec_scalar("updateEvery", updateEvery);
    rec_scalar("consolePlots", num("console.plots") != 0);
    rec_scalar("rdRate_Hz", num_or("display.rdRate_Hz", 10));   /* range-Doppler maps per second */
    rec_scalar("minExtent_deg", minExtent);
    rec_scalar("meanPri_s", meanPri);
    rec_scalar("trackEnable", trackEnable);
    rec_scalar("trackDim", dim);
    {
        double s3[3] = { 0, 0, 0 }, q3[3] = { 0, 0, 0 };
        for (int k = 0; k < dim; k++) { s3[k] = sqrt(sig[k] * sig[k]); q3[k] = qv[k]; }
        rec_vec("trackSigma", s3, 3, 0);
        rec_scalar("trackModel", ca);
        rec_vec("trackQ", q3, 3, 0);
    }
    rec_scalar("trackSigmaV0", maxSpeed / sqrt(3.0));
    rec_scalar("trackSigmaA0", num("track.maxAccel_mps2") / sqrt(3.0));
    rec_scalar("trackBatchGap_s", num("track.gnnSector_deg") / rate);
    rec_scalar("trackBatchMax_s", num("track.gnnMaxSector_deg") / rate);
    rec_scalar("trackGate", 2 * gammaincinv(gateProb, dim / 2.0));
    rec_scalar("trackUseVr", num("track.useRadialVelocity") != 0);
    rec_scalar("trackVrUnambiguous", maxSpeed < vrSpan / 2);
    rec_scalar("trackSigmaVr", sqrt(Rvr));
    rec_scalar("trackVrSpan", vrSpan);
    rec_vec("trackVrSpanPulse", vrSpanPulse, 2, 0);
    rec_scalar("trackGateVr", 2 * gammaincinv(gateProb, (dim + 1) / 2.0));
    rec_scalar("trackConfirmHits", num("track.confirmHits"));
    rec_scalar("trackMaxMissTentative", num("track.maxMissTentative"));
    rec_scalar("trackMaxMiss", num("track.maxMiss"));
    rec_scalar("trackMinSpeed", minSpeed);
    rec_scalar("scanTime_s", 360.0 / rate);
    int trkMethod = 0;                                   /* 0 EKF, 1 PMBM, 2 TOMHT */
    if (has("track.method") && !empty("track.method")) {
        const char *m = str("track.method");
        if (same_text(m, "pmbm")) trkMethod = 1;
        else if (same_text(m, "tomht")) trkMethod = 2;
        else if (!same_text(m, "ekf")) fail("cfg.track.method must be 'ekf', 'pmbm' or 'tomht'.");
    }
    const int pmbm = trkMethod != 0;                     /* the IMM-UKF of the PMBM (also the TOMHT) */
    const char *tk = trkMethod == 2 ? "track.tomht" : "track.pmbm";
#define TK(f) (snprintf(p, sizeof p, "%s.%s", tk, f), p)
    rec_scalar("trackMethod", trkMethod);
    rec_scalar("trackAz0_deg", num("ant.az0_deg"));
    if (pmbm) {                                          /* pmbm_setup (tomht_setup) of track_setup.m */
        const double pd = num(TK("pd")), ps = trkMethod == 2 ? 1.0 : num(TK("ps"));
        if (!(pd > 0 && pd < 1) || !(ps > 0 && ps <= 1)) fail("cfg.track.pmbm: pd must lie in (0, 1) and ps in (0, 1].");
        const double fpps = num(TK("falsePlotsPerScan")), und = num(TK("undetectedTargets"));
        if (!(fpps > 0) || !(und > 0)) fail("cfg.track.pmbm: falsePlotsPerScan and undetectedTargets must be positive.");
        double rLo = INFINITY, rHi = -INFINITY;          /* the range covered by the sub-pulses that are on */
        for (int i = 0; i < 2; i++)
            for (int k = 0; k < wv[i].nRange; k++) {
                if (wv[i].range_m[k] < rLo) rLo = wv[i].range_m[k];
                if (wv[i].range_m[k] > rHi) rHi = wv[i].range_m[k];
            }
        double vol = (rHi - rLo) * 2 * PI;
        if (dim == 3) {
            double emax = beamElKeep[0], emin = beamElKeep[0], de = 1;
            for (int b = 1; b < nBeams; b++) { emax = beamElKeep[b] > emax ? beamElKeep[b] : emax; emin = beamElKeep[b] < emin ? beamElKeep[b] : emin; }
            if (nBeams > 1) de = fabs(beamElKeep[nBeams - 1] - beamElKeep[nBeams - 2]);
            vol = vol * (PI / 180 * (emax - emin + de));
        }
        rec_scalar("pmbmPd", pd);
        rec_scalar("pmbmPs", ps);
        const double ll[3] = { log(fpps / vol), log(und / vol), -log(vrSpan) };
        rec_vec("pmbmLogLambda", ll, 3, 0);
        const double maxHyp = num(TK("maxHypotheses"));
        const double rPrune = trkMethod == 2 ? 1e-9 : num(TK("rPrune")), rConfirm = trkMethod == 2 ? 1.0 : num(TK("rConfirm"));
        const double rDelete = trkMethod == 2 ? 1e-6 : num(TK("rDelete"));
        if (!(maxHyp >= 1) || !(rPrune < rDelete && rDelete < rConfirm && rConfirm <= 1))
            fail("cfg.track.pmbm: maxHypotheses >= 1 and rPrune < rDelete < rConfirm <= 1 needed.");
        rec_scalar("pmbmMaxHyp", maxHyp);
        const double pr[4] = { log(num(TK("hypPrune"))), rPrune, rConfirm, rDelete };
        rec_vec("pmbmPrune", pr, 4, 0);
        /* IMM models */
        const MValue *mv = field(TK("imm.models"), 1);
        double models[3];
        int M = 0;
        static const char *const mnames[3] = { "cv", "ca", "ct" };
        const int nm = mv->type == MV_STR ? 1 : mv->type == MV_CELL ? mv->rows * mv->cols : 0;
        int okModels = nm >= 1 && nm <= 3;
        for (int i = 0; i < nm && okModels; i++) {
            const MValue *e = mv->type == MV_STR ? mv : mv->cells[i];
            int k = -1;
            if (e->type == MV_STR) for (int j = 0; j < 3; j++) if (same_text(e->str, mnames[j])) k = j;
            for (int j = 0; j < M; j++) if (models[j] == k) k = -1;
            if (k < 0) okModels = 0; else models[M++] = k;
        }
        if (!okModels) fail("cfg.track.pmbm.imm.models must be distinct models of 'cv', 'ca', 'ct'.");
        rec_vec("pmbmModels", models, M, 0);
        const MValue *muv = field(TK("imm.mu0"), 1);
        if (muv->type != MV_NUM || muv->rows * muv->cols != M) fail("cfg.track.pmbm.imm.mu0 needs one non-negative probability per model.");
        double mu[3], musum = 0;
        for (int i = 0; i < M; i++) {
            mu[i] = muv->re[i];
            if (mu[i] < 0) fail("cfg.track.pmbm.imm.mu0 needs one non-negative probability per model.");
            musum += mu[i];
        }
        if (!(musum > 0)) fail("cfg.track.pmbm.imm.mu0 needs one non-negative probability per model.");
        for (int i = 0; i < M; i++) mu[i] /= musum;
        rec_vec("pmbmMu0", mu, M, 0);
        double sw = M == 1 ? 0 : num(TK("imm.switchProb"));
        if (!(sw >= 0 && sw < (M - 1.0) / M)) fail("cfg.track.pmbm.imm.switchProb must lie in [0, (M-1)/M) for M models.");
        rec_scalar("pmbmSwitchLambda", 1 - sw * M / (M - 1 + (M == 1)));
        const double qzA = pow(num("track.sigmaAccelZ_mps2"), 2), qzJ = pow(num("track.sigmaJerkZ_mps3"), 2);
        const double sCv = num(TK("imm.sigmaAccelCv_mps2")), sCa = num(TK("imm.sigmaJerkCa_mps3"));
        const double sCt = num(TK("imm.sigmaAccelCt_mps2"));
        double qcv[3] = { sCv * sCv, sCv * sCv, qzA }, qca[3] = { sCa * sCa, sCa * sCa, qzJ }, qct[3] = { sCt * sCt, sCt * sCt, qzA };
        for (int k = dim; k < 3; k++) qcv[k] = qca[k] = qct[k] = 0;
        rec_vec("pmbmQCv", qcv, 3, 0);
        rec_vec("pmbmQCa", qca, 3, 0);
        rec_vec("pmbmQCt", qct, 3, 0);
        const double sigmaTurn0 = PI / 180 * num(TK("imm.maxTurn_dps")) / sqrt(3.0);
        const double maxClimb = num(TK("maxClimb_mps"));
        const double turn[3] = { pow(PI / 180 * num(TK("imm.sigmaTurn_dps")), 2), sigmaTurn0,
                                 (maxSpeed < maxClimb ? maxSpeed : maxClimb) / sqrt(3.0) };
        rec_vec("pmbmTurn", turn, 3, 0);
        const double sigmaA0 = num("track.maxAccel_mps2") / sqrt(3.0);
        const double unused[2] = { sigmaA0 * sigmaA0, sigmaTurn0 * sigmaTurn0 };
        rec_vec("pmbmUnused", unused, 2, 0);
        const int n = 3 * dim + 1;
        const double ua = num(TK("ukf.alpha")), ub = num(TK("ukf.beta")), uk = num(TK("ukf.kappa"));
        const double lam = ua * ua * (n + uk) - n;
        if (!(n + lam > 0)) fail("cfg.track.pmbm.ukf: alpha^2*(n + kappa) must be positive (n = %d).", n);
        const double wm0 = lam / (n + lam);
        const double ukf[4] = { sqrt(n + lam), wm0, wm0 + 1 - ua * ua + ub, 1 / (2 * (n + lam)) };
        rec_vec("pmbmUkf", ukf, 4, 0);
        const double timing[2] = { num("track.gnnMaxSector_deg") / rate + 2 * cpiTime, 1 / rate };
        rec_vec("pmbmTiming", timing, 2, 0);
    }
    if (trkMethod == 2) {                                /* tomht_setup of track_setup.m */
        const double v[8] = { num("track.tomht.pd"), num("track.tomht.maxHypotheses"), num("track.tomht.kBest"),
                              log(num("track.tomht.hypPrune")), num("track.tomht.nScan"), num("track.tomht.confirmScore"),
                              num("track.tomht.deleteScore"), num("track.tomht.pDelete") };
        const double hp = num("track.tomht.hypPrune");
        if (!(v[1] >= 1) || !(v[2] >= 1) || !(v[4] >= 1 && v[4] <= 16) || v[4] != floor(v[4]) || !(v[6] > 0) || !(v[7] > 0 && v[7] < 1) ||
            !(hp > 0 && hp < 1))
            fail("cfg.track.tomht: maxHypotheses >= 1, kBest >= 1, nScan a whole number 1 ... 16, deleteScore > 0, "
                 "pDelete and hypPrune in (0, 1) needed.");
        rec_vec("tomhtParams", v, 8, 0);
    }
#undef TK

    if (!g_quiet) {
        printf("%d element(s), %d beam(s), DDC %d samples at %.0f MHz; short %d / long %d range cells; CPI %d "
               "pulses, FFT %d x %d row(s); clutter filter %s (%d taps); CFAR %s %s, alpha %.2f dB (Pfa %.1e, %.2f looks)\n",
               nCh, nBeams, nDdc, fsOut / 1e6, ws->nRange, wl->nRange, cpiLength, nfft, nRows,
               mtiFir ? "FIR" : "canceller", nTaps, cfarType == 0 ? "CA" : cfarType == 1 ? "GO" : cfarType == 2 ? "SO" : "OS",
               method, 10 * log10(alpha), cfar_pfa(alpha, &det), nLooks);
    }
    (void)trackEnable;
    free(beamElKeep); free(g_qs); free(g_qw); g_qs = g_qw = NULL;
    free(pri); free(mti); free(tDdc); free(wRe); free(wIm); free(beamU); free(win); free(fHz);
    free(noiseGain); free(Cre); free(Cim); free(ddc.h1); free(ddc.h2);
    for (int i = 0; i < 2; i++) {
        free(wv[i].codeRe); free(wv[i].codeIm); free(wv[i].misRe); free(wv[i].misIm);
        free(wv[i].tapsRe); free(wv[i].tapsIm); free(wv[i].delay_s); free(wv[i].range_m);
    }
}

/* The export and the file (errors of export_config() come back here). */
static int run_export(const char *outPath)
{
    if (setjmp(g_jb) != 0) {
        fprintf(stderr, "radar_export_config: %s\n", g_err);
        return 1;
    }
    export_config();
    if (write_file(outPath) != 0) {
        fprintf(stderr, "radar_export_config: cannot write %s\n", outPath);
        return 1;
    }
    if (!g_quiet) printf("Receiver configuration (%d records, %s PRF) written to %s\n", g_nRec, str("prf.mode"), outPath);
    return 0;
}

static void usage(void)
{
    puts("usage: radar_export_config [--config-m FILE] [--out FILE] [--set STATEMENT]... [--quiet]\n"
         "  Writes the receiver configuration (as export_receiver_config.m) from radar_config.m.\n"
         "  --config-m FILE   MATLAB configuration (default " RADAR_SOURCE_DIR "/../radar_config.m)\n"
         "  --out FILE        output (default " RADAR_SOURCE_DIR "/config/radar_config.bin)\n"
         "  --set STATEMENT   MATLAB assignment run after the file, e.g. --set \"cfg.prf.mode = 'staggered';\"\n"
         "  --quiet           print only errors");
}

int main(int argc, char **argv)
{
    const char *mPath = RADAR_SOURCE_DIR "/../radar_config.m";
    const char *outPath = RADAR_SOURCE_DIR "/config/radar_config.bin";
    const char **sets = calloc((size_t)argc + 1, sizeof(char *));
    int nSets = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i], *v = i + 1 < argc ? argv[i + 1] : NULL;
        if (!strcmp(a, "--config-m") && v) { mPath = v; i++; }
        else if (!strcmp(a, "--out") && v) { outPath = v; i++; }
        else if (!strcmp(a, "--set") && v) { sets[nSets++] = v; i++; }
        else if (!strcmp(a, "--quiet")) g_quiet = 1;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); free(sets); return 0; }
        else { fprintf(stderr, "unknown option %s\n", a); usage(); free(sets); return 2; }
    }
    char err[1024];
    MValue *cfg = mfile_eval(mPath, sets, nSets, err, sizeof err);
    free(sets);
    if (!cfg) {
        fprintf(stderr, "radar_export_config: %s\n", err);
        return 1;
    }
    if (cfg->type != MV_STRUCT) {
        fprintf(stderr, "radar_export_config: %s does not return a struct\n", mPath);
        mvalue_free(cfg);
        return 1;
    }
    g_cfg = cfg;
    if (!g_quiet) {
        printf("Reading %s\n", mPath);
        fflush(stdout);
    }
    const int r = run_export(outPath);
    if (r != 0) {
        mvalue_free(cfg);
        return 1;
    }
    for (int i = 0; i < g_nRec; i++) { free(g_rec[i].re); free(g_rec[i].im); }
    free(g_rec);
    mvalue_free(cfg);
    return 0;
}

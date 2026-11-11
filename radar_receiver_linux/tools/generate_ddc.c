/*
 * generate_ddc.c - radar_generate_ddc: the DDC files of the simulator without
 * MATLAB, as generate_ddc_files.m does (a run configuration of the same name in
 * CLion).
 *
 *   radar_generate_ddc [--config-m FILE] [--set STATEMENT]... [--out FOLDER]
 *                      [--raw | --no-raw] [--threads N] [--targets FILE.csv]
 *                      [--max-pulses N] [--quiet]
 *     --config-m FILE  the configuration function (default ../radar_config.m of this project)
 *     --set STATEMENT  a MATLAB assignment run after the file (repeatable), e.g.
 *                      --set "cfg.sim.sectors_deg = [20 45];" --set "cfg.ddc.outputRate_Hz = 4e6;"
 *     --out FOLDER     the DDC files (default cfg.io.ddcFolder)
 *     --raw, --no-raw  also write / do not write the raw ADC data (default cfg.io.saveRaw,
 *                      into cfg.io.rawFolder)
 *     --threads N      PRIs simulated in parallel (OpenMP; default: all cores)
 *     --targets FILE   the targets from a CSV file instead of cfg.scenario (one target per
 *                      line: name,rcs_m2,x,y,z,vx,vy,vz[,tStart,turn_dps,accelAlong,accelUp]...;
 *                      positions [m] and velocities [m/s] at t = 0, x East, y North, z up from
 *                      the ground below the radar; maneuver rows of make_target.m in groups of 4;
 *                      export_targets_csv.m writes cfg.targets in this format)
 *     --max-pulses N   only the first N active PRIs (a quick test)
 *     --quiet          print only errors and warnings
 *
 * Every active PRI of the schedule (cfg.sim.nScans, cfg.sim.sectors_deg) goes
 * through steps 1 - 3 of run_radar_sim only, independently of the others, as
 * generate_ddc_files.m:
 *   1. targets (radar equation, Swerling RCS per antenna pass, maneuvers),
 *      ground clutter of the terrain map (land, mountains, hills, desert, sea;
 *      shadowing; CNR -> sigma0), transmitter leakage, STC, thermal noise on
 *      the IF of every sub-pulse (cfg.wave.<pulse>.enable / fc_Hz), receiver
 *      gain, 16-bit ADC with its noise (ENOB)
 *   2. raw ADC data -> <rawFolder>/pulse_<n>.mat (adcData int16 [nAdc x nElements], hdr)
 *   3. fixed-point DDC (2 or 4 MHz) -> <ddcFolder>/ddc_<n>.mat (I, Q int16 [nDdc x nElements], hdr)
 * The files are those of MATLAB (same variables, header and naming), so
 * radar_stream_send / replay_ddc_udp.m send them to the receiver.
 *
 * Random numbers: the target RCS of every scan, the random scenario and the
 * clutter map are drawn first from one stream seeded with cfg.sim.seed; every
 * PRI then has its own stream for the noise (seeded with cfg.sim.seed and its
 * pulse number), so the files do not depend on the number of threads or the
 * order of the PRIs. The random numbers are not those of MATLAB (another
 * generator): the same statistics, other values. Without noise, clutter and
 * Swerling fluctuation (cfg.rx.addThermalNoise = false, cfg.adc.addAdcNoise =
 * false, cfg.clutter.enable = false, cfg.sim.swerling = 0) a PRI is
 * deterministic and its DDC output equals MATLAB's up to the rounding of a few
 * ADC samples (MATLAB sums the signal in single precision).
 *
 * Not supported (MATLAB only): cfg.adc.quantize = false (analog ADC values).
 * Input files: as radar_export_config (made for other configuration values:
 * rebuild them with make_input_files.m in MATLAB).
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "mat_file.h"
#include "mat_write.h"
#include "mfile.h"

#include <errno.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#  include <omp.h>
#endif
#ifdef _WIN32
#  include <direct.h>
#  define MKDIR(p) _mkdir(p)
#else
#  include <sys/stat.h>
#  define MKDIR(p) mkdir(p, 0777)
#endif

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
    void *p = calloc(n ? n : 1, size ? size : 1);
    if (!p) fail("out of memory (%zu x %zu bytes)", n, size);
    return p;
}

static double now_s(void)
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}

/* ------------------------------------------------------------------------ */
/* Reading cfg (as radar_export_config)                                     */
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

static double num_or(const char *path, double def) { return has(path) && !empty(path) ? num(path) : def; }

static const char *str(const char *path)
{
    const MValue *v = field(path, 1);
    if (v->type != MV_STR) fail("cfg.%s must be a string", path);
    return v->str;
}

static const char *str_or(const char *path, const char *def) { return has(path) && !empty(path) ? str(path) : def; }

static int same_text(const char *a, const char *b)
{
    for (; *a && *b; a++, b++)
        if ((*a | 32) != (*b | 32)) return 0;
    return *a == *b;
}

/* A numeric matrix of cfg (values column-major), n = rows * cols. */
static const MValue *matrix(const char *path, int *rows, int *cols)
{
    const MValue *v = field(path, 1);
    if (v->type != MV_NUM) fail("cfg.%s must be numeric", path);
    *rows = v->rows;
    *cols = v->cols;
    return v;
}

/* ------------------------------------------------------------------------ */
/* Input files                                                              */
/* ------------------------------------------------------------------------ */

typedef struct { const char *name; double v[4]; int n; } MadeFor;

static MatFile load_input(const char *key, const MadeFor *mf, int nMf)
{
    char path[64];
    snprintf(path, sizeof path, "files.%s", key);
    const char *file = str(path);
    MatFile m;
    char err[256];
    if (mat_load(file, &m, err, sizeof err) != 0)
        fail("cfg.files.%s: %s (create the default input files in MATLAB with make_input_files(cfg), or give your "
             "own file)", key, err);
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
        if (!ok)
            fail("%s is a default input file made for other configuration values (madeFor.%s). Rebuild it in MATLAB "
                 "with make_input_files(cfg, '%s'), or replace it with your own file without 'madeFor'", file, bad, key);
    }
    return m;
}

static const MatVar *mat_var(const MatFile *m, const char *name, const char *key)
{
    const MatVar *v = mat_get(m, name);
    if (!v || !v->re) fail("cfg.files.%s: the variable '%s' is missing", key, name);
    return v;
}

static double *copy_var(const MatVar *v)
{
    double *p = xcalloc(v->numel, sizeof(double));
    memcpy(p, v->re, sizeof(double) * v->numel);
    return p;
}

/* ------------------------------------------------------------------------ */
/* Random numbers: xoshiro256** (splitmix64 seeding), Gaussian by the polar method */
/* ------------------------------------------------------------------------ */

typedef struct { uint64_t s[4]; int haveSpare; double spare; } Rng;

static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void rng_seed(Rng *r, uint64_t seed, uint64_t stream)
{
    uint64_t x = seed * 0x2545F4914F6CDD1DULL ^ (stream + 0x632BE59BD9B4E019ULL);
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&x);
    r->haveSpare = 0;
}

static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static uint64_t rng_next(Rng *r)
{
    uint64_t *s = r->s;
    const uint64_t result = rotl(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);
    return result;
}

static double rng_uniform(Rng *r) { return ((rng_next(r) >> 11) + 0.5) * (1.0 / 9007199254740992.0); }   /* (0, 1) */

static double rng_normal(Rng *r)
{
    if (r->haveSpare) { r->haveSpare = 0; return r->spare; }
    double u, v, s;
    do {
        u = 2.0 * rng_uniform(r) - 1.0;
        v = 2.0 * rng_uniform(r) - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    const double f = sqrt(-2.0 * log(s) / s);
    r->spare = v * f;
    r->haveSpare = 1;
    return u * f;
}

/* ------------------------------------------------------------------------ */
/* FFT (complex double, radix 2, in place)                                  */
/* ------------------------------------------------------------------------ */

typedef struct { int n, logn; int *rev; double *c, *s; } Fft;

static void fft_init(Fft *f, int n)
{
    f->n = n;
    f->logn = 0;
    while ((1 << f->logn) < n) f->logn++;
    f->rev = xcalloc((size_t)n, sizeof(int));
    f->c = xcalloc((size_t)n / 2 + 1, sizeof(double));
    f->s = xcalloc((size_t)n / 2 + 1, sizeof(double));
    for (int i = 0; i < n; i++) {
        int r = 0;
        for (int b = 0; b < f->logn; b++) if (i & (1 << b)) r |= 1 << (f->logn - 1 - b);
        f->rev[i] = r;
    }
    for (int k = 0; k < n / 2; k++) { f->c[k] = cos(2 * PI * k / n); f->s[k] = -sin(2 * PI * k / n); }
}

static void fft_free(Fft *f) { free(f->rev); free(f->c); free(f->s); memset(f, 0, sizeof *f); }

/* Forward (inverse = 0) or inverse (1, without the 1/n) transform of re/im. */
static void fft_run(const Fft *f, double *re, double *im, int inverse)
{
    const int n = f->n;
    for (int i = 0; i < n; i++) {
        const int r = f->rev[i];
        if (r > i) {
            double t = re[i]; re[i] = re[r]; re[r] = t;
            t = im[i]; im[i] = im[r]; im[r] = t;
        }
    }
    const double sg = inverse ? -1.0 : 1.0;
    for (int len = 2; len <= n; len <<= 1) {
        const int half = len >> 1, step = n / len;
        for (int i = 0; i < n; i += len)
            for (int k = 0; k < half; k++) {
                const double wr = f->c[k * step], wi = sg * f->s[k * step];
                const int a = i + k, b = a + half;
                const double xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
                re[b] = re[a] - xr;  im[b] = im[a] - xi;
                re[a] += xr;          im[a] += xi;
            }
    }
}

/* ------------------------------------------------------------------------ */
/* Antenna patterns                                                         */
/* ------------------------------------------------------------------------ */

typedef struct { int n; double *x, *g_dB, gMin_dB; } Pattern;

static Pattern pattern(const MatVar *x, const MatVar *g, const char *what)
{
    if (x->numel != g->numel || x->numel < 2) fail("cfg.files.antenna: %s: the angles and gains differ in length", what);
    Pattern p = { (int)x->numel, copy_var(x), copy_var(g), INFINITY };
    for (int i = 0; i < p.n; i++) if (p.g_dB[i] < p.gMin_dB) p.gMin_dB = p.g_dB[i];
    return p;
}

/* Linear power gain (pattern_gain.m: interp1 in dB, outside the table the lowest gain). */
static double pattern_gain(const Pattern *p, double q)
{
    if (!(q >= p->x[0] && q <= p->x[p->n - 1])) return pow(10.0, p->gMin_dB / 10.0);
    int lo = 0, hi = p->n - 1;
    while (hi - lo > 1) {
        const int mid = (lo + hi) / 2;
        if (p->x[mid] <= q) lo = mid; else hi = mid;
    }
    const double f = p->x[hi] > p->x[lo] ? (q - p->x[lo]) / (p->x[hi] - p->x[lo]) : 0.0;
    return pow(10.0, (p->g_dB[lo] + f * (p->g_dB[hi] - p->g_dB[lo])) / 10.0);
}

static double sind(double x) { const double i = x / 180.0; return i == floor(i) && isfinite(i) ? 0.0 : sin(i * PI); }
static double cosd(double x) { return sind(x + 90.0); }
static double wrap180(double a) { double r = fmod(a + 180.0, 360.0); if (r < 0) r += 360.0; return r - 180.0; }
static double mod360(double a) { double r = fmod(a, 360.0); if (r < 0) r += 360.0; return r; }

/* ------------------------------------------------------------------------ */
/* The simulated radar (setup_radar.m: what the signal generation needs)    */
/* ------------------------------------------------------------------------ */

typedef struct {
    const char *name;
    int    enable, phaseCode;
    double fc, lambda, fIF, radarConst;
    double width, bw, chip, tTx;
    int    nCode;
    double *codeRe, *codeIm;
    int    compMismatch, nMis;
    double *misRe, *misIm, misRate;
    int    rxFirst, rxLast;          /* 1-based receive interval (PRI header), 0 0 = off */
} Wave;

typedef struct {
    double c, fs, fc, lambda, fIF, fLO;
    double radarConst, N0, sigmaTh, voltGain, leakAmp;
    int    fullScale, nBits, quantize, addAdcNoise, addNoise, addLeak, stcOn;
    double lsbPerUnit, sigmaAdc, sigmaThAdc, N0eff, N0ddc;
    double height, tilt, rate, az0, spacing_m, azSpan;
    Pattern az, txEl, elElem;
    int    nCh;
    int    nStagger;
    double *pri;
    /* DDC */
    int    D1, D2, L1, L2, shift1, shift2, fixedPoint;
    double *h1, *h2;                 /* integer coefficients */
    double fracScale;                /* 2^coefFracBits */
    double fsOut, fs1;
    int    nAdc, nDdc;
    double *cosNco, *msinNco;        /* 16-bit NCO values per ADC sample */
    float  *stc;                     /* STC amplitude gain per ADC sample */
    unsigned char *rxSample;         /* receive interval samples (ADC overload) */
    int    nRx;
    double ddcT0;
    Wave   w[2];
    double overloadWarn;
} Radar;

static void envelope(const Wave *w, double t, double *re, double *im)
{
    *re = *im = 0.0;
    if (t < 0 || t >= w->width) return;
    if (!w->phaseCode) {
        const double k = w->bw / w->width, d = t - w->width / 2;
        const double ph = PI * k * (d * d);
        *re = cos(ph);
        *im = sin(ph);
    } else {
        double m = floor(t / w->chip);
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
        double nu = 1.0, de = 1.0;
        for (int n = 1; n <= nbar - 1; n++) {
            nu *= 1.0 - (m * (double)m / sp2) / (A * A + (n - 0.5) * (n - 0.5));
            if (n != m) de *= 1.0 - (double)m * m / ((double)n * n);
        }
        const double Fm = ((m + 1) % 2 ? -1.0 : 1.0) * nu / (2.0 * de);
        for (int k = 0; k < N; k++) {
            const double xi = (k - 0.5 * N + 0.5) / N;
            w[k] += 2.0 * Fm * cos(2.0 * PI * m * xi);
        }
    }
}

/* Valid FIR decimation (fir_decimate_valid.m), real filter, complex input. */
static int fir_dec(const double *xr, const double *xi, int N, const double *h, int L, int D, double *yr, double *yi)
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

/* CNR of uniform clutter of sigma0 = 1 at range R, elevation el (clutter_cnr.m). */
static double clutter_cnr(const Radar *s, double R, double el, int taylorNbar, double taylorSll)
{
    const Wave *w = s->w[0].enable ? &s->w[0] : &s->w[1];
    const double c = s->c, fs = s->fs;
    const double u = -sind(s->tilt) * cosd(el) + cosd(s->tilt) * sind(el);
    const double elRel = asin(u) * 180 / PI;
    const double gEl = pattern_gain(&s->txEl, elRel) * pattern_gain(&s->elElem, elRel);
    double azInt = 0;                                      /* trapz over the table points within azSpan */
    {
        int first = 1;
        double pa = 0, pg = 0;
        for (int i = 0; i < s->az.n; i++) {
            const double a = s->az.x[i];
            if (!(fabs(a) <= s->azSpan)) continue;
            const double g = pattern_gain(&s->az, a), g2 = g * g;
            if (!first) azInt += (a - pa) * PI / 180 * (g2 + pg) / 2;
            first = 0;
            pa = a;
            pg = g2;
        }
    }
    /* compression filter (unit norm): matched (Taylor for LFM) or mismatch */
    const int Lp = (int)lround(w->width * s->fsOut);
    double *win = xcalloc((size_t)Lp, sizeof(double));
    if (w->phaseCode) for (int k = 0; k < Lp; k++) win[k] = 1.0;
    else taylor_window(Lp, taylorNbar, taylorSll, win);
    int L = Lp;
    double *hr, *hi;
    if (w->compMismatch) {
        const int r = (int)lround(s->fsOut / w->misRate);
        L = w->nMis * r;
        hr = xcalloc((size_t)L, sizeof(double));
        hi = xcalloc((size_t)L, sizeof(double));
        for (int i = 0; i < w->nMis; i++)
            for (int k = 0; k < r; k++) { hr[i * r + k] = w->misRe[i]; hi[i * r + k] = w->misIm[i]; }
    } else {
        hr = xcalloc((size_t)L, sizeof(double));
        hi = xcalloc((size_t)L, sizeof(double));
        for (int k = 0; k < Lp; k++) {
            double re, im;
            envelope(w, (k + 0.5) / s->fsOut, &re, &im);
            hr[Lp - 1 - k] = re * win[k];
            hi[Lp - 1 - k] = -im * win[k];
        }
    }
    double nrm = 0;
    for (int k = 0; k < L; k++) nrm += hr[k] * hr[k] + hi[k] * hi[k];
    nrm = sqrt(nrm);
    for (int k = 0; k < L; k++) { hr[k] /= nrm; hi[k] /= nrm; }
    /* real filters of the DDC, the equivalent filter at the ADC rate */
    double *h1r = xcalloc((size_t)s->L1, sizeof(double)), *h2r = xcalloc((size_t)s->L2, sizeof(double));
    for (int i = 0; i < s->L1; i++) h1r[i] = s->h1[i] / s->fracScale;
    for (int i = 0; i < s->L2; i++) h2r[i] = s->h2[i] / s->fracScale;
    const int nUp = s->D1 * (s->L2 - 1) + 1, nEq = s->L1 + nUp - 1;
    double *hEq = xcalloc((size_t)nEq, sizeof(double));
    for (int i = 0; i < s->L1; i++)
        for (int j = 0; j < s->L2; j++) hEq[i + j * s->D1] += h1r[i] * h2r[j];
    const int D = s->D1 * s->D2;
    /* noise gain: sum |conv(hEq, tapsUp)|^2 */
    const int nTU = D * (L - 1) + 1, nCv = nEq + nTU - 1;
    double *cr = xcalloc((size_t)nCv, sizeof(double)), *ci = xcalloc((size_t)nCv, sizeof(double));
    for (int i = 0; i < nEq; i++)
        for (int j = 0; j < L; j++) { cr[i + j * D] += hEq[i] * hr[j]; ci[i + j * D] += hEq[i] * hi[j]; }
    double noiseGain = 0;
    for (int i = 0; i < nCv; i++) noiseGain += cr[i] * cr[i] + ci[i] * ci[i];
    /* S: energy of the response to a unit echo, summed over the delays of one ADC sample */
    const int nPad = nEq, N = (int)lround(w->width * fs) + 2 * nPad;
    double *xr = xcalloc((size_t)N, sizeof(double)), *xi = xcalloc((size_t)N, sizeof(double));
    double *y1r = xcalloc((size_t)N, sizeof(double)), *y1i = xcalloc((size_t)N, sizeof(double));
    double *yr = xcalloc((size_t)N, sizeof(double)), *yi = xcalloc((size_t)N, sizeof(double));
    double S = 0;
    for (int shift = 0; shift < D; shift++) {
        for (int k = 0; k < N; k++) envelope(w, (double)(k - nPad - shift) / fs, &xr[k], &xi[k]);
        const int K1 = fir_dec(xr, xi, N, h1r, s->L1, s->D1, y1r, y1i);
        const int K2 = fir_dec(y1r, y1i, K1, h2r, s->L2, s->D2, yr, yi);
        for (int n = 0; n < K2 + L - 1; n++) {
            double sr = 0, si = 0;
            for (int j = 0; j < L; j++) {
                const int m = n - j;
                if (m < 0 || m >= K2) continue;
                sr += yr[m] * hr[j] - yi[m] * hi[j];
                si += yr[m] * hi[j] + yi[m] * hr[j];
            }
            S += sr * sr + si * si;
        }
    }
    free(win); free(hr); free(hi); free(h1r); free(h2r); free(hEq); free(cr); free(ci);
    free(xr); free(xi); free(y1r); free(y1i); free(yr); free(yi);
    return w->radarConst * 1.0 * gEl / (R * R * R) * azInt * c / (2 * fs) * S / (s->N0ddc * noiseGain * fs);
}

/* STC attenuation [dB] of a curve [range_m attenuation_dB] at the equivalent range R (stc_gain.m). */
static double stc_att(const double *curve, int rows, double R)
{
    if (R > curve[rows - 1]) return 0.0;                  /* column-major: ranges, then dB */
    double r = R < curve[0] ? curve[0] : R;
    if (r > curve[rows - 1]) r = curve[rows - 1];
    for (int i = 0; i < rows - 1; i++)
        if (r <= curve[i + 1]) {
            const double f = curve[i + 1] > curve[i] ? (r - curve[i]) / (curve[i + 1] - curve[i]) : 0.0;
            return curve[rows + i] + f * (curve[rows + i + 1] - curve[rows + i]);
        }
    return curve[2 * rows - 1];
}

static void setup(Radar *s)
{
    static const char *names[2] = { "short", "long" };
    char p[96];
    memset(s, 0, sizeof *s);
    s->c = num("const.c");
    s->fs = num("adc.fs_Hz");
    s->fc = num("fc");
    s->lambda = s->c / s->fc;
    s->fIF = num("adc.fIF_Hz");
    s->fLO = s->fc - s->fIF;
    s->nCh = (int)num("ant.nElements");
    if (s->nCh < 1) fail("cfg.ant.nElements must be at least 1");
    s->height = num("ant.height_m");
    s->tilt = num("ant.tilt_deg");
    s->rate = num("ant.rpm") * 360.0 / 60.0;
    s->az0 = num("ant.az0_deg");

    /* antenna file */
    const MadeFor mfAnt[2] = { { "fc", { s->fc }, 1 }, { "spacing", { num("ant.elementSpacing_lambda") }, 1 } };
    MatFile ant = load_input("antenna", mfAnt, 2);
    s->az = pattern(mat_var(&ant, "az_deg", "antenna"), mat_var(&ant, "azGain_dB", "antenna"), "azimuth pattern");
    s->elElem = pattern(mat_var(&ant, "elRel_deg", "antenna"), mat_var(&ant, "elElemGain_dB", "antenna"), "element pattern");
    s->txEl = pattern(mat_var(&ant, "elRel_deg", "antenna"), mat_var(&ant, "txElGain_dB", "antenna"), "transmit pattern");
    s->spacing_m = mat_var(&ant, "elementSpacing_m", "antenna")->re[0];
    const double txGain = mat_var(&ant, "txGain_dBi", "antenna")->re[0], rxGain = mat_var(&ant, "rxChannelGain_dBi", "antenna")->re[0];
    mat_free(&ant);
    {                                                    /* 3 dB azimuth beamwidth (setup_radar.m beamwidth_3dB) */
        const int n = s->az.n;
        const double *av = s->az.x, *gv = s->az.g_dB;
        int iMax = 0;
        for (int i = 1; i < n; i++) if (gv[i] > gv[iMax]) iMax = i;
        const double gmax = gv[iMax], lim = -3.0103;
        int iL = -1, iR = -1;
        for (int i = 0; i <= iMax; i++) if (gv[i] - gmax < lim) iL = i;
        for (int i = iMax; i < n && iR < 0; i++) if (gv[i] - gmax < lim) iR = i;
        if (iL < 0 || iR < 1) fail("cfg.files.antenna: the azimuth pattern has no 3 dB points");
        const double g0 = gv[iL] - gmax, g1 = gv[iL + 1] - gmax;
        const double azL = av[iL] + (lim - g0) * (av[iL + 1] - av[iL]) / (g1 - g0);
        const double h0 = gv[iR - 1] - gmax, h1 = gv[iR] - gmax;
        const double azR = av[iR - 1] + (lim - h0) * (av[iR] - av[iR - 1]) / (h1 - h0);
        s->azSpan = 3.0 * (azR - azL);
    }
    const double Pt = num("tx.peakPower_W");
    s->radarConst = Pt * pow(10.0, txGain / 10) * pow(10.0, rxGain / 10) * s->lambda * s->lambda /
                    (pow(4 * PI, 3) * pow(10.0, num("rx.losses_dB") / 10));
    s->N0 = num("const.kB") * num("const.T0") * pow(10.0, num("rx.noiseFigure_dB") / 10);
    s->sigmaTh = sqrt(s->N0 * s->fs / 2);
    s->voltGain = pow(10.0, num("rx.gain_dB") / 20);
    s->leakAmp = sqrt(2 * Pt * pow(10.0, -num("rx.txIsolation_dB") / 10));
    s->addNoise = num("rx.addThermalNoise") != 0;
    s->addLeak = num("rx.addLeakage") != 0;
    s->nBits = (int)num("adc.nBits");
    if (s->nBits > 16) fail("cfg.adc.nBits = %d: the ADC samples are stored as int16 (at most 16 bits)", s->nBits);
    s->fullScale = 1 << (s->nBits - 1);
    const double fullScale_W = 1e-3 * pow(10.0, num("adc.fullScale_dBm") / 10);
    s->lsbPerUnit = s->fullScale / sqrt(2 * fullScale_W);
    const double sinad = pow(10.0, (6.02 * num("adc.enob") + 1.76) / 10);
    const double adcNoise2 = (double)s->fullScale * s->fullScale / 2 / sinad;
    s->sigmaAdc = sqrt(adcNoise2 - 2.0 / 12 > 0 ? adcNoise2 - 2.0 / 12 : 0);
    s->sigmaThAdc = s->sigmaTh * s->voltGain * s->lsbPerUnit;
    s->N0eff = s->N0 * (s->sigmaThAdc * s->sigmaThAdc + adcNoise2) / (s->sigmaThAdc * s->sigmaThAdc);
    s->quantize = num("adc.quantize") != 0;
    s->addAdcNoise = num("adc.addAdcNoise") != 0;
    if (!s->quantize) fail("cfg.adc.quantize = false (analog ADC values) is not supported by radar_generate_ddc");
    s->overloadWarn = num("adc.overloadWarn");

    /* PRIs */
    const char *prf = str("prf.mode");
    if (same_text(prf, "fixed")) {
        s->nStagger = 1;
        s->pri = xcalloc(1, sizeof(double));
        s->pri[0] = num("prf.pri_s");
    } else if (same_text(prf, "staggered")) {
        MatFile st = load_input("stagger", NULL, 0);
        const MatVar *pv = mat_var(&st, "priVector_s", "stagger");
        s->nStagger = (int)pv->numel;
        s->pri = copy_var(pv);
        mat_free(&st);
    } else {
        fail("Unknown PRF mode \"%s\"", prf);
    }
    double minPri = s->pri[0];
    for (int i = 1; i < s->nStagger; i++) if (s->pri[i] < minPri) minPri = s->pri[i];

    /* DDC filters, decimation (ddc_decimation.m) */
    const double outRate = num("ddc.outputRate_Hz");
    const char *ddcKey = outRate == 2e6 ? "ddc" : outRate == 4e6 ? "ddc4MHz" : NULL;
    if (!ddcKey) fail("cfg.ddc.outputRate_Hz must be 2e6 or 4e6 (is %g).", outRate);
    s->D1 = (int)num("ddc.stage1Decimation");
    const double D2 = s->fs / s->D1 / outRate;
    if (D2 < 1 || fabs(D2 - round(D2)) > 1e-9)
        fail("The stage-1 rate %g MHz is not a whole multiple of the DDC output rate %g MHz.", s->fs / s->D1 / 1e6, outRate / 1e6);
    s->D2 = (int)lround(D2);
    const double fracBits = num("ddc.coefFracBits");
    const MadeFor mfDdc[3] = { { "fs_Hz", { s->fs }, 1 }, { "decimation", { s->D1, s->D2 }, 2 }, { "coefFracBits", { fracBits }, 1 } };
    MatFile df = load_input(ddcKey, mfDdc, 3);
    const MatVar *h1 = mat_var(&df, "h1", ddcKey), *h2 = mat_var(&df, "h2", ddcKey);
    s->L1 = (int)h1->numel;  s->L2 = (int)h2->numel;
    s->h1 = copy_var(h1);    s->h2 = copy_var(h2);
    mat_free(&df);
    s->fracScale = pow(2.0, fracBits);
    {
        int r, c;
        const MValue *sh = matrix("ddc.shift", &r, &c);
        if (r * c != 2) fail("cfg.ddc.shift must hold the shifts of both stages");
        s->shift1 = (int)sh->re[0];  s->shift2 = (int)sh->re[1];
    }
    s->fixedPoint = num("ddc.fixedPoint") != 0;
    s->fs1 = s->fs / s->D1;
    s->fsOut = s->fs1 / s->D2;
    {                                                    /* DDC output noise (setup_radar.m) for the CNR */
        double e1 = 0, e2 = 0;
        const int nUp = s->D1 * (s->L2 - 1) + 1, nEq = s->L1 + nUp - 1;
        double *hEq = xcalloc((size_t)nEq, sizeof(double));
        for (int i = 0; i < s->L1; i++)
            for (int j = 0; j < s->L2; j++) hEq[i + j * s->D1] += s->h1[i] / s->fracScale * s->h2[j] / s->fracScale;
        for (int i = 0; i < nEq; i++) e1 += hEq[i] * hEq[i];
        for (int j = 0; j < s->L2; j++) e2 += pow(s->h2[j] / s->fracScale, 2);
        free(hEq);
        const double g1 = pow(2.0, fracBits - s->shift1), g2 = pow(2.0, fracBits - s->shift2);
        const double adcLevel = (s->sigmaThAdc * s->sigmaThAdc + adcNoise2) * e1 * pow(g1 * g2, 2) / 2;
        const double rounding = s->fixedPoint * (1 + e2 * g2 * g2) / 12;
        s->N0ddc = s->N0eff * (adcLevel + rounding) / adcLevel;
    }

    /* sub-pulses */
    int anyCode = 0;
    for (int i = 0; i < 2; i++) {
        Wave *w = &s->w[i];
        w->name = names[i];
        snprintf(p, sizeof p, "wave.%s.enable", names[i]);  w->enable = num_or(p, 1) != 0;
        snprintf(p, sizeof p, "wave.%s.fc_Hz", names[i]);   w->fc = num_or(p, s->fc);
        snprintf(p, sizeof p, "wave.%s.type", names[i]);
        const char *type = str(p);
        if (same_text(type, "phasecode")) w->phaseCode = 1;
        else if (!same_text(type, "lfm")) fail("Unknown waveform type \"%s\"", type);
        snprintf(p, sizeof p, "wave.%s.width_s", names[i]);       w->width = num(p);
        snprintf(p, sizeof p, "wave.%s.bandwidth_Hz", names[i]);  w->bw = num_or(p, 0);
        snprintf(p, sizeof p, "wave.%s.chip_s", names[i]);        w->chip = num_or(p, 1e-6);
        snprintf(p, sizeof p, "wave.%s.tTx_s", names[i]);         w->tTx = num(p);
        snprintf(p, sizeof p, "wave.%s.compression", names[i]);
        w->compMismatch = same_text(str_or(p, "matched"), "mismatch");
        if (!(w->fc > 0) || !isfinite(w->fc)) fail("cfg.wave.%s.fc_Hz must be a positive carrier frequency [Hz] (or [] = cfg.fc).", names[i]);
        w->lambda = s->c / w->fc;
        w->fIF = s->fIF + (w->fc - s->fc);
        const double bwHalf = 0.5 * (w->phaseCode ? 1.0 / w->chip : w->bw);
        if (w->enable && (w->fIF - bwHalf <= 0 || w->fIF + bwHalf >= s->fs / 2))
            fail("cfg.wave.%s.fc_Hz = %.3f MHz puts the sub-pulse on the IF %.3f MHz: it must lie inside (0, %.1f MHz) "
                 "with its bandwidth (LO = cfg.fc - cfg.adc.fIF_Hz = %.3f MHz).", names[i], w->fc / 1e6, w->fIF / 1e6,
                 s->fs / 2e6, s->fLO / 1e6);
        w->radarConst = s->radarConst * pow(w->lambda / s->lambda, 2);
        anyCode |= w->phaseCode;
    }
    if (!s->w[0].enable && !s->w[1].enable) fail("Both sub-pulses are switched off (cfg.wave.short.enable, cfg.wave.long.enable).");
    if (anyCode) {
        MadeFor mf[4];
        for (int i = 0; i < 2; i++) {
            snprintf(p, sizeof p, "wave.%s.mismatchLength", names[i]);
            const int hasLen = has(p) && !empty(p);
            mf[2 * i] = (MadeFor){ i ? "longChips" : "shortChips", { round(s->w[i].width / s->w[i].chip) }, 1 };
            mf[2 * i + 1] = (MadeFor){ i ? "longMismatchLength" : "shortMismatchLength", { hasLen ? num(p) : 0 }, hasLen };
        }
        MatFile codes = load_input("phaseCode", mf, 4);
        for (int i = 0; i < 2; i++) {
            Wave *w = &s->w[i];
            if (!w->phaseCode) continue;
            char q[64];
            snprintf(q, sizeof q, "%sCode", names[i]);
            const MatVar *cv = mat_var(&codes, q, "phaseCode");
            w->nCode = (int)cv->numel;
            if (w->nCode != (int)lround(w->width / w->chip))
                fail("The %s code has %d chips, but width_s/chip_s = %d.", names[i], w->nCode, (int)lround(w->width / w->chip));
            w->codeRe = xcalloc((size_t)w->nCode, sizeof(double));
            w->codeIm = xcalloc((size_t)w->nCode, sizeof(double));
            for (int k = 0; k < w->nCode; k++) { w->codeRe[k] = cv->re[k]; w->codeIm[k] = cv->im ? cv->im[k] : 0; }
            if (w->compMismatch) {
                snprintf(q, sizeof q, "%sMismatch", names[i]);
                const MatVar *mv = mat_var(&codes, q, "phaseCode");
                w->nMis = (int)mv->numel;
                w->misRe = xcalloc((size_t)w->nMis, sizeof(double));
                w->misIm = xcalloc((size_t)w->nMis, sizeof(double));
                for (int k = 0; k < w->nMis; k++) { w->misRe[k] = mv->re[k]; w->misIm[k] = mv->im ? mv->im[k] : 0; }
                w->misRate = 1.0 / w->chip;
                snprintf(q, sizeof q, "%sMismatchRate_Hz", names[i]);
                const MatVar *rv = mat_get(&codes, q);
                if (rv && rv->re && rv->numel >= 1) w->misRate = rv->re[0];
            }
        }
        mat_free(&codes);
    }
    Wave *ws = &s->w[0], *wl = &s->w[1];
    if (ws->enable && wl->enable && wl->tTx < ws->tTx + ws->width) fail("The long pulse must start after the end of the short pulse.");

    /* PRI timeline, NCO, STC, DDC samples and receive intervals */
    const double maxRange = num("rx.maxRange_m");
    const double filterSpan = (s->L1 - 1) / s->fs + (s->L2 - 1) / s->fs1;
    const Wave *wLast = wl->enable ? wl : ws;
    const double rxEnd = wLast->tTx + wLast->width + 2 * maxRange / s->c + filterSpan + 2 / s->fsOut;
    if (rxEnd > minPri)
        fail("Receive window (%.1f us) does not fit in the shortest PRI (%.1f us): reduce cfg.rx.maxRange_m or start the "
             "long pulse earlier.", rxEnd * 1e6, minPri * 1e6);
    s->nAdc = (int)lround(rxEnd * s->fs);
    const double shortEnd = wl->enable ? wl->tTx : INFINITY;
    s->cosNco = xcalloc((size_t)s->nAdc, sizeof(double));
    s->msinNco = xcalloc((size_t)s->nAdc, sizeof(double));
    s->stc = xcalloc((size_t)s->nAdc, sizeof(float));
    s->rxSample = xcalloc((size_t)s->nAdc, 1);
    s->stcOn = num("rx.stc.enable") != 0;
    int cr = 0, cc = 0, lr = 0, lc = 0;
    const double *curve = NULL, *curveL = NULL;
    if (s->stcOn) {
        curve = matrix("rx.stc.curve", &cr, &cc)->re;
        if (cc != 2 || cr < 1) fail("cfg.rx.stc.curve must have rows [range_m attenuation_dB]");
        curveL = curve;  lr = cr;  lc = cc;
        if (has("rx.stc.curveLong") && !empty("rx.stc.curveLong")) {
            curveL = matrix("rx.stc.curveLong", &lr, &lc)->re;
            if (lc != 2 || lr < 1) fail("cfg.rx.stc.curveLong must have rows [range_m attenuation_dB]");
        }
    }
    for (int k = 0; k < s->nAdc; k++) {
        const double t = k / s->fs;
        double f = ws->enable ? ws->fIF : wl->fIF;
        if (ws->enable && wl->enable && t >= wl->tTx) f = wl->fIF;
        s->cosNco[k] = round(cos(2 * PI * f * t) * 32768.0) / 32768.0;
        s->msinNco[k] = round(-sin(2 * PI * f * t) * 32768.0) / 32768.0;
        const int inShort = ws->enable && t >= ws->tTx + ws->width && t < shortEnd;
        const int inLong = wl->enable && t >= wl->tTx + wl->width;
        s->rxSample[k] = (unsigned char)(inShort || inLong);
        s->nRx += s->rxSample[k];
        double A = 0.0;
        if (s->stcOn) {
            A = curve[cr];                                /* while transmitting: the first row */
            if (inShort) A = stc_att(curve, cr, s->c * (t - ws->tTx - ws->width) / 2);
            if (inLong) A = stc_att(curveL, lr, s->c * (t - wl->tTx - wl->width) / 2);
        }
        s->stc[k] = (float)pow(10.0, -A / 20);
    }
    const int n1 = (s->nAdc - s->L1) / s->D1 + 1;
    s->nDdc = (n1 - s->L2) / s->D2 + 1;
    s->ddcT0 = (s->L1 - 1) / 2.0 / s->fs + (s->L2 - 1) / 2.0 / s->fs1;
    for (int i = 0; i < 2; i++) {
        Wave *w = &s->w[i];
        w->rxFirst = w->rxLast = 0;
        if (!w->enable) continue;
        const double iv0 = w->tTx + w->width, iv1 = i == 0 ? shortEnd : INFINITY;
        for (int k = 0; k < s->nDdc; k++) {
            const double tFirst = (double)((long)k * s->D2 * s->D1) / s->fs;
            const double tLast = (double)(((long)k * s->D2 + s->L2 - 1) * s->D1 + s->L1 - 1) / s->fs;
            if (tFirst >= iv0 && tLast < iv1) {
                if (!w->rxFirst) w->rxFirst = k + 1;
                w->rxLast = k + 1;
            }
        }
        if (!w->rxFirst) fail("The %s-pulse receive interval holds no DDC sample: move the long pulse later or increase "
                              "cfg.rx.maxRange_m.", w->name);
    }
}

/* ------------------------------------------------------------------------ */
/* Targets (make_target.m, target_state.m, target_geometry.m, scenarios)    */
/* ------------------------------------------------------------------------ */

typedef struct {
    char   name[16];
    double rcs, pos0[3], vel0[3];
    int    nMan;
    double *man;                 /* rows [tStart turn_dps accelAlong accelUp], row-major */
    double *rcsScan;             /* RCS of every antenna pass */
    int    nScan;
} Target;

static Target make_target(const char *name, double rcs, double range, double az, double alt, double vRad, double vCross,
                          const double *man, int nMan)
{
    Target t;
    memset(&t, 0, sizeof t);
    snprintf(t.name, sizeof t.name, "%s", name);
    t.rcs = rcs;
    t.pos0[0] = range * sind(az);  t.pos0[1] = range * cosd(az);  t.pos0[2] = alt;
    t.vel0[0] = vRad * sind(az) + vCross * cosd(az);
    t.vel0[1] = vRad * cosd(az) - vCross * sind(az);
    t.vel0[2] = 0.0;
    t.nMan = nMan;
    if (nMan) {
        t.man = xcalloc((size_t)nMan * 4, sizeof(double));
        memcpy(t.man, man, sizeof(double) * 4 * nMan);
    }
    return t;
}

/* Motion over tau with turn rate m[0] [deg/s], along-track m[1], vertical m[2] acceleration (target_state.m segment). */
static void segment(double *p, double *v, double tau, const double *m)
{
    double w = m[0] * PI / 180, aT = m[1];
    const double aZ = m[2];
    const double s0 = hypot(v[0], v[1]), th0 = atan2(v[0], v[1]);
    double th = th0, sp = s0;
    if (s0 > 0 || aT > 0) {
        double tauH = tau;
        if (aT < 0 && s0 + aT * tau < 0) tauH = -s0 / aT;
        th = th0 + w * tauH;
        double Is, Ic, Ius, Iuc;
        if (fabs(w * tauH) < 1e-4) {
            Is  = tauH * sin(th0) + w * tauH * tauH / 2 * cos(th0);
            Ic  = tauH * cos(th0) - w * tauH * tauH / 2 * sin(th0);
            Ius = tauH * tauH / 2 * sin(th0) + w * tauH * tauH * tauH / 3 * cos(th0);
            Iuc = tauH * tauH / 2 * cos(th0) - w * tauH * tauH * tauH / 3 * sin(th0);
        } else {
            Is  = (cos(th0) - cos(th)) / w;
            Ic  = (sin(th) - sin(th0)) / w;
            Ius = -tauH * cos(th) / w + (sin(th) - sin(th0)) / (w * w);
            Iuc =  tauH * sin(th) / w + (cos(th) - cos(th0)) / (w * w);
        }
        sp = s0 + aT * tauH;
        p[0] += s0 * Is + aT * Ius;
        p[1] += s0 * Ic + aT * Iuc;
        if (tauH < tau) { aT = 0; w = 0; }
    }
    p[2] += v[2] * tau + aZ * tau * tau / 2;
    v[0] = sp * sin(th);  v[1] = sp * cos(th);  v[2] = v[2] + aZ * tau;
}

static void target_state(const Target *t, double time, double *pos, double *vel)
{
    double p[3] = { t->pos0[0], t->pos0[1], t->pos0[2] }, v[3] = { t->vel0[0], t->vel0[1], t->vel0[2] };
    double tc = 0, cur[3] = { 0, 0, 0 };
    for (int k = 0; k < t->nMan; k++) {
        if (time <= t->man[4 * k]) break;
        segment(p, v, t->man[4 * k] - tc, cur);
        tc = t->man[4 * k];
        cur[0] = t->man[4 * k + 1];  cur[1] = t->man[4 * k + 2];  cur[2] = t->man[4 * k + 3];
    }
    segment(p, v, time - tc, cur);
    memcpy(pos, p, sizeof p);
    memcpy(vel, v, sizeof v);
}

static void target_geometry(const Target *t, double time, double hAnt, double *R, double *az, double *el, double *vr)
{
    double pos[3], vel[3];
    target_state(t, time, pos, vel);
    const double d[3] = { pos[0], pos[1], pos[2] - hAnt };
    *R = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    *az = mod360(atan2(d[0], d[1]) * 180 / PI);
    *el = asin(d[2] / *R) * 180 / PI;
    *vr = (vel[0] * d[0] + vel[1] * d[1] + vel[2] * d[2]) / *R;
}

static int scenario(const char *type, Rng *rng, Target **out)
{
    if (same_text(type, "six")) {
        Target *t = xcalloc(6, sizeof(Target));
        t[0] = make_target("T1", 1.0, 12e3,  30, 1000, -237,   0, NULL, 0);
        t[1] = make_target("T2", 0.5, 20e3,  95, 3000,  130,  40, NULL, 0);
        t[2] = make_target("T3", 5.0, 35e3, 160, 6000, -290,   0, NULL, 0);
        t[3] = make_target("T4", 10,  42e3, 215, 8000,   80, 150, NULL, 0);
        t[4] = make_target("T5", 2.0,  8e3, 270,  300,  -40,  20, NULL, 0);
        t[5] = make_target("T6", 3.0, 46e3, 330, 2000,  210, -50, NULL, 0);
        *out = t;
        return 6;
    }
    if (same_text(type, "challenge")) {
        static const double m11[12] = { 3, 5, 0, 0, 15, -5, 0, 0, 27, 0, 0, 0 }, m12[8] = { 0, 0, 8, 0, 30, 0, 0, 0 };
        static const double m13[8] = { 4, -3, -4, 4, 19, 0, 0, 0 };
        Target *t = xcalloc(13, sizeof(Target));
        t[0]  = make_target("C1",  0.3, 48e3,    50, 4500, -150,   50, NULL, 0);
        t[1]  = make_target("C2",  2.0,  8e3,   125,  300,  -25,    4, NULL, 0);
        t[2]  = make_target("C3",  5.0, 28e3,    90, 6000,    4, -220, NULL, 0);
        t[3]  = make_target("C4", 10,    5e3,   140,  400,  120,    0, NULL, 0);
        t[4]  = make_target("C5",  1.0, 23e3,   140, 3000,  150,    0, NULL, 0);
        t[5]  = make_target("C6",  2.0, 33e3,   185, 5000, -200,    0, NULL, 0);
        t[6]  = make_target("C7",  2.0, 33.3e3, 185, 5200, -160,    0, NULL, 0);
        t[7]  = make_target("C8",  3.0, 50e3,   235, 9000, -420,    0, NULL, 0);
        t[8]  = make_target("C9",  2.0, 12e3,   280, 9000,  150,  100, NULL, 0);
        t[9]  = make_target("C10", 5.0, 45e3,   325,  100, -250,   30, NULL, 0);
        t[10] = make_target("C11", 3.0, 30e3,    15, 4000, -250,    0, m11, 3);
        t[11] = make_target("C12", 2.0, 20e3,   265, 2000, -100,    0, m12, 2);
        t[12] = make_target("C13", 2.0, 25e3,   345,  600, -200,   50, m13, 2);
        *out = t;
        return 13;
    }
    if (same_text(type, "random")) {                     /* random_targets_scenario.m with this program's generator */
        const int n = (int)num("scenario.random.n");
        int r, c;
        const double *rcs = matrix("scenario.random.rcs_m2", &r, &c)->re, *rng_m = matrix("scenario.random.range_m", &r, &c)->re;
        const double *alt = matrix("scenario.random.alt_m", &r, &c)->re, *spd = matrix("scenario.random.speed_mps", &r, &c)->re;
        const double maxEl = num("scenario.random.maxEl_deg");
        const double rmax = 0.95 * num("rx.maxRange_m");
        const double r0 = rng_m[0] < rmax ? rng_m[0] : rmax, r1 = rng_m[1] < rmax ? rng_m[1] : rmax;
        if (n < 1 || r0 > r1 || rcs[0] <= 0 || rcs[0] > rcs[1] || alt[0] > alt[1] || spd[0] > spd[1])
            fail("cfg.scenario.random: n must be a positive integer and every [min max] ordered (RCS > 0).");
        Target *t = xcalloc((size_t)n, sizeof(Target));
        for (int k = 0; k < n; k++) {
            const double rc = exp(log(rcs[0]) + rng_uniform(rng) * (log(rcs[1]) - log(rcs[0])));
            const double range = r0 + rng_uniform(rng) * (r1 - r0), az = 360 * rng_uniform(rng);
            double altMax = range * tan(maxEl * PI / 180);
            if (alt[1] < altMax) altMax = alt[1];
            const double altMin = alt[0] < altMax ? alt[0] : altMax;
            const double a = altMin + rng_uniform(rng) * (altMax - altMin);
            const double speed = spd[0] + rng_uniform(rng) * (spd[1] - spd[0]), hdg = 360 * rng_uniform(rng);
            char nm[16];
            snprintf(nm, sizeof nm, "R%d", k + 1);
            t[k] = make_target(nm, rc, range, az, a, speed * cosd(hdg - az), speed * sind(hdg - az), NULL, 0);
        }
        *out = t;
        return n;
    }
    fail("cfg.scenario.type must be 'six', 'challenge' or 'random'.");
    return 0;
}

/* Targets from a CSV file: name,rcs,x,y,z,vx,vy,vz[,tStart,turn,along,up]... ('#' comments). */
static int read_targets_csv(const char *path, Target **out)
{
    FILE *f = fopen(path, "r");
    if (!f) fail("--targets %s: %s", path, strerror(errno));
    int n = 0, cap = 8;
    Target *t = xcalloc((size_t)cap, sizeof(Target));
    char line[8192];
    int lineNo = 0;
    while (fgets(line, sizeof line, f)) {
        lineNo++;
        char *s = line;
        while (*s == ' ' || *s == '\t') s++;
        if (*s == '#' || *s == '\n' || *s == '\r' || *s == 0) continue;
        double v[4 + 7 + 4 * 64];
        char name[16] = "";
        char *tok = strtok(s, ",\r\n");
        int nv = 0;
        if (tok) snprintf(name, sizeof name, "%s", tok);
        while ((tok = strtok(NULL, ",\r\n")) && nv < (int)(sizeof v / sizeof v[0])) v[nv++] = strtod(tok, NULL);
        if (nv < 7 || (nv - 7) % 4 != 0)
            fail("--targets %s line %d: name,rcs_m2,x,y,z,vx,vy,vz and maneuver rows of 4 values expected", path, lineNo);
        if (n == cap) {
            cap *= 2;
            Target *q = realloc(t, sizeof(Target) * (size_t)cap);
            if (!q) fail("out of memory");
            t = q;
        }
        Target *g = &t[n++];
        memset(g, 0, sizeof *g);
        snprintf(g->name, sizeof g->name, "%s", name);
        g->rcs = v[0];
        memcpy(g->pos0, v + 1, sizeof(double) * 3);
        memcpy(g->vel0, v + 4, sizeof(double) * 3);
        g->nMan = (nv - 7) / 4;
        if (g->nMan) {
            g->man = xcalloc((size_t)g->nMan * 4, sizeof(double));
            memcpy(g->man, v + 7, sizeof(double) * 4 * g->nMan);
        }
    }
    fclose(f);
    if (n == 0) fail("--targets %s holds no target", path);
    *out = t;
    return n;
}

/* ------------------------------------------------------------------------ */
/* Ground clutter (clutter_terrain.m, clutter_init.m, clutter_signal.m)     */
/* ------------------------------------------------------------------------ */

typedef struct { double rel_dB, icm, meanVel, meanDir; int shape; } TerrainType;   /* shape 0 flat, 1 ridge, 2 rolling */

typedef struct {
    int    nAz, nR, G;
    double azStep, dR;
    double *az, *range;
    int    *delay;               /* two-way delay of every range cell [ADC samples] */
    double *groupEl;             /* [G] */
    /* visible cells, sorted by azimuth row */
    int    nCells, *rowStart;    /* cells of row a: rowStart[a] ... rowStart[a+1]-1 */
    int    *cellRow, *cellRange, *cellGroup;
    double *rhoRe, *rhoIm, *fd, *dcRe, *dcIm;
    int    hasDc;
    double sigma0_dB, cnr_dB;
    /* per sub-pulse */
    int    nTx[2], nFft[2], start0[2];
    double *txRe[2], *txIm[2];   /* FFT of the IF pulse [nFft] */
    double fdScale[2], ampScale[2];
    int    sameCarrier;
    double *steer[2];            /* steering vectors [2][nCh * G] (re, im interleaved) */
} Clutter;

static TerrainType terrain_type(const char *name)
{
    TerrainType t = { 0, NAN, 0, 0, 0 };
    char p[96];
    snprintf(p, sizeof p, "clutter.types.%s", name);
    if (!has(p)) {
        if (!same_text(name, "land")) fail("Terrain type '%s' of cfg.clutter.regions is not in cfg.clutter.types.", name);
        return t;
    }
    snprintf(p, sizeof p, "clutter.types.%s.sigma0Rel_dB", name);  t.rel_dB = num_or(p, 0);
    snprintf(p, sizeof p, "clutter.types.%s.icmStd_mps", name);    t.icm = num_or(p, NAN);
    snprintf(p, sizeof p, "clutter.types.%s.meanVel_mps", name);   t.meanVel = num_or(p, 0);
    snprintf(p, sizeof p, "clutter.types.%s.meanDir_deg", name);   t.meanDir = num_or(p, 0);
    snprintf(p, sizeof p, "clutter.types.%s.shape", name);
    const char *sh = str_or(p, "flat");
    if (same_text(sh, "flat")) t.shape = 0;
    else if (same_text(sh, "ridge")) t.shape = 1;
    else if (same_text(sh, "rolling")) t.shape = 2;
    else fail("cfg.clutter.types.%s.shape must be 'flat', 'ridge' or 'rolling'.", name);
    return t;
}

static void clutter_build(Clutter *C, Radar *s, Rng *rng, int taylorNbar, double taylorSll)
{
    memset(C, 0, sizeof *C);
    const double c = s->c, fs = s->fs, hAnt = s->height, kae = 4.0 / 3.0 * 6371e3;
    const double elev0 = num("clutter.elevation_deg");
    /* sigma0 of land from the CNR */
    const double cnrRange = num_or("clutter.cnrRange_m", 5000);
    if (!(cnrRange > 0)) fail("cfg.clutter.cnrRange_m must be a positive range [m].");
    const double cnrPer = clutter_cnr(s, cnrRange, elev0, taylorNbar, taylorSll);
    C->sigma0_dB = has("clutter.cnr_dB") && !empty("clutter.cnr_dB") ? num("clutter.cnr_dB") - 10 * log10(cnrPer)
                                                                       : num("clutter.sigma0_dB");
    C->cnr_dB = C->sigma0_dB + 10 * log10(cnrPer);
    const double sigma0Land = pow(10.0, C->sigma0_dB / 10);
    /* terrain */
    C->azStep = num("clutter.azStep_deg");
    C->dR = num("clutter.rangeStepSamples") * c / (2 * fs);
    const char *mode = str_or("clutter.terrain", "uniform");
    const int uniform = same_text(mode, "uniform");
    if (!uniform && !same_text(mode, "regions")) fail("cfg.clutter.terrain must be 'uniform' or 'regions'.");
    int nReg = 0;
    const MValue *reg = NULL;
    if (!uniform && has("clutter.regions")) {
        reg = field("clutter.regions", 1);
        if (reg->type == MV_CELL) {
            if (reg->cols != 4) fail("cfg.clutter.regions must have rows {type, [az1 az2], [r1 r2], height_m}.");
            nReg = reg->rows;
        } else if (!(reg->type == MV_NUM && reg->rows * reg->cols == 0)) {
            fail("cfg.clutter.regions must have rows {type, [az1 az2], [r1 r2], height_m}.");
        }
    }
    const char **rType = xcalloc((size_t)(nReg ? nReg : 1), sizeof(char *));
    double *rAz = xcalloc((size_t)(2 * nReg + 1), sizeof(double)), *rR = xcalloc((size_t)(2 * nReg + 1), sizeof(double));
    double *rH = xcalloc((size_t)(nReg + 1), sizeof(double));
    for (int k = 0; k < nReg; k++) {
        const MValue *e[4];
        for (int j = 0; j < 4; j++) e[j] = reg->cells[k + j * reg->rows];   /* column-major cells */
        if (e[0]->type != MV_STR || e[1]->type != MV_NUM || e[1]->rows * e[1]->cols != 2 || e[2]->type != MV_NUM ||
            e[2]->rows * e[2]->cols != 2 || e[3]->type != MV_NUM || e[3]->rows * e[3]->cols != 1 ||
            !(e[2]->re[0] < e[2]->re[1]) || e[3]->re[0] < 0)
            fail("cfg.clutter.regions row %d: {type, [az1 az2] deg, [r1 r2] m with r1 < r2, height_m >= 0}.", k + 1);
        rType[k] = e[0]->str;
        rAz[2 * k] = e[1]->re[0];  rAz[2 * k + 1] = e[1]->re[1];
        rR[2 * k] = e[2]->re[0];   rR[2 * k + 1] = e[2]->re[1];
        rH[k] = e[3]->re[0];
    }
    double rMax = num_or("clutter.maxRange_m", NAN);
    if (isnan(rMax)) {
        rMax = sqrt(2 * kae * hAnt);
        for (int k = 0; k < nReg; k++) if (rH[k] > 0 && rR[2 * k + 1] > rMax) rMax = rR[2 * k + 1];
    }
    const int i0 = (int)ceil(num("clutter.minRange_m") / C->dR), i1 = (int)floor(rMax / C->dR);
    C->nR = i1 - i0 + 1;
    if (C->nR < 1) fail("No clutter cell: cfg.clutter.minRange_m beyond the clutter range");
    C->range = xcalloc((size_t)C->nR, sizeof(double));
    for (int r = 0; r < C->nR; r++) C->range[r] = (i0 + r) * C->dR;
    C->nAz = (int)floor((360 - C->azStep) / C->azStep + 1e-9) + 1;
    C->az = xcalloc((size_t)C->nAz, sizeof(double));
    for (int a = 0; a < C->nAz; a++) C->az[a] = a * C->azStep;
    /* types: land, then the region types in their first order */
    int nTypes = 1;
    const char **typeName = xcalloc((size_t)nReg + 1, sizeof(char *));
    typeName[0] = "land";
    int *regType = xcalloc((size_t)(nReg ? nReg : 1), sizeof(int));
    for (int k = 0; k < nReg; k++) {
        int found = -1;
        for (int j = 0; j < nTypes; j++) if (!strcmp(typeName[j], rType[k])) found = j;
        if (found < 0) { typeName[nTypes] = rType[k]; found = nTypes++; }
        regType[k] = found;
    }
    TerrainType *types = xcalloc((size_t)nTypes, sizeof(TerrainType));
    for (int j = 0; j < nTypes; j++) types[j] = terrain_type(typeName[j]);
    const size_t nCell = (size_t)C->nAz * C->nR;
    unsigned char *type = xcalloc(nCell, 1);
    double *height = xcalloc(nCell, sizeof(double));
    unsigned char *regionOf = xcalloc(nCell, 1);              /* 0 = plain land, else region k + 1 */
    for (int k = 0; k < nReg; k++) {
        double width = mod360(rAz[2 * k + 1] - rAz[2 * k]);
        if (width == 0) width = 360;
        const TerrainType *tt = &types[regType[k]];
        const double H = rH[k], r1 = rR[2 * k], r2 = rR[2 * k + 1];
        for (int a = 0; a < C->nAz; a++) {
            const double d = mod360(C->az[a] - rAz[2 * k]);
            if (!(d <= width)) continue;
            const double vv = d / width;
            for (int r = 0; r < C->nR; r++) {
                const double R = C->range[r];
                if (!(R >= r1 && R <= r2)) continue;
                const double u = (R - r1) / (r2 - r1 > 2.220446049250313e-16 ? r2 - r1 : 2.220446049250313e-16);
                double h;
                if (tt->shape == 0) h = H;
                else if (tt->shape == 1) h = H * sqrt(sin(PI * vv)) * sin(PI * u) * (1 + 0.15 * cos(2 * PI * 4 * vv) * sin(2 * PI * 3 * u));
                else {
                    const double crest = pow(sin(PI * ((R - r1) / 1500 + 0.5 * cos(2 * PI * 2 * vv))), 2);
                    h = H * (0.8 + 0.2 * cos(2 * PI * 3 * vv)) * u * (0.4 + 0.6 * crest);
                }
                const size_t q = (size_t)a * C->nR + r;
                type[q] = (unsigned char)regType[k];
                height[q] = h;
                regionOf[q] = (unsigned char)(k + 1);
            }
        }
    }
    /* elevation, visibility, sigma0, motion, groups */
    double *elDeg = xcalloc(nCell, sizeof(double)), *sigma0 = xcalloc(nCell, sizeof(double));
    unsigned char *visible = xcalloc(nCell, 1);
    for (int a = 0; a < C->nAz; a++) {
        double mx = -INFINITY;
        for (int r = 0; r < C->nR; r++) {
            const size_t q = (size_t)a * C->nR + r;
            const double R = C->range[r], tanEl = (height[q] - hAnt - R * R / (2 * kae)) / R;
            if (tanEl > mx) mx = tanEl;
            visible[q] = uniform || tanEl >= mx - 1e-12;
            elDeg[q] = height[q] == 0 ? elev0 : atan(tanEl) * 180 / PI;
            sigma0[q] = visible[q] ? sigma0Land * pow(10.0, types[type[q]].rel_dB / 10) : 0.0;
        }
    }
    int *group = xcalloc(nCell, sizeof(int));
    C->groupEl = xcalloc((size_t)nReg + 1, sizeof(double));
    C->groupEl[0] = elev0;
    C->G = 1;
    for (int k = 0; k < nReg; k++) {
        double sw = 0, swe = 0;
        for (size_t q = 0; q < nCell; q++)
            if (regionOf[q] == k + 1 && height[q] != 0 && visible[q]) {
                const double R = C->range[q % C->nR], w = sigma0[q] / (R * R * R);
                sw += w;
                swe += w * elDeg[q];
            }
        if (!(sw > 0)) {
            int any = 0;
            for (size_t q = 0; q < nCell && !any; q++) any = regionOf[q] == k + 1 && height[q] != 0 && visible[q];
            if (!any) continue;
        }
        C->groupEl[C->G] = swe / sw;
        for (size_t q = 0; q < nCell; q++) if (regionOf[q] == k + 1 && height[q] != 0 && visible[q]) group[q] = C->G;
        C->G++;
    }
    /* cells: power, complex reflectivity, internal motion (clutter_init.m) */
    C->delay = xcalloc((size_t)C->nR, sizeof(int));
    for (int r = 0; r < C->nR; r++) C->delay[r] = (int)lround(2 * C->range[r] / c * fs);
    const char *spec = str("clutter.spectrum");
    const int gauss = same_text(spec, "gaussian");
    if (!gauss && !same_text(spec, "exponential")) fail("cfg.clutter.spectrum must be 'gaussian' or 'exponential'.");
    double beta = 0, dcr = 0;
    if (!gauss) {                                        /* Billingsley (wind in mph, f in GHz) */
        double w = num("clutter.windSpeed_mps") / 0.44704;
        if (w < 0.1) w = 0.1;
        beta = num_or("clutter.expBeta_spm", 1 / (0.1 * log10(w) + 0.068));
        dcr = num_or("clutter.dcToAc", 489.8 * pow(w, -1.55) * pow(s->fc / 1e9, -1.21));
    }
    const double icm = num("clutter.icmStd_mps");
    const char *fixEnv = getenv("RADAR_GEN_FIXED_CLUTTER");     /* tests: the expected amplitude in every cell */
    const int fixedClutter = fixEnv && *fixEnv && strcmp(fixEnv, "0") != 0;
    int nVis = 0;
    for (size_t q = 0; q < nCell; q++) nVis += sigma0[q] > 0;
    C->nCells = nVis;
    C->rowStart = xcalloc((size_t)C->nAz + 1, sizeof(int));
    C->cellRow = xcalloc((size_t)nVis, sizeof(int));
    C->cellRange = xcalloc((size_t)nVis, sizeof(int));
    C->cellGroup = xcalloc((size_t)nVis, sizeof(int));
    C->rhoRe = xcalloc((size_t)nVis, sizeof(double));
    C->rhoIm = xcalloc((size_t)nVis, sizeof(double));
    C->fd = xcalloc((size_t)nVis, sizeof(double));
    C->dcRe = xcalloc((size_t)nVis, sizeof(double));
    C->dcIm = xcalloc((size_t)nVis, sizeof(double));
    int m = 0;
    for (int a = 0; a < C->nAz; a++) {
        C->rowStart[a] = m;
        for (int r = 0; r < C->nR; r++) {
            const size_t q = (size_t)a * C->nR + r;
            if (!(sigma0[q] > 0)) continue;
            const double R = C->range[r];
            const double elRel = asin(-sind(s->tilt) * cosd(elDeg[q]) + cosd(s->tilt) * sind(elDeg[q])) * 180 / PI;
            const double gEl = pattern_gain(&s->txEl, elRel) * pattern_gain(&s->elElem, elRel);
            const double area = R * (C->azStep * PI / 180) * C->dR;
            const double power = s->radarConst * sigma0[q] * area / pow(R, 4) * gEl;
            const double ampRms = sqrt(2 * power);
            double re = ampRms * rng_normal(rng) / sqrt(2.0), im = ampRms * rng_normal(rng) / sqrt(2.0);
            const TerrainType *tt = &types[type[q]];
            const double meanVel = tt->meanVel * cosd(C->az[a] - tt->meanDir);
            double v = 0;
            if (!isnan(tt->icm)) v = meanVel + tt->icm * rng_normal(rng);
            else if (gauss) v = meanVel + icm * rng_normal(rng);
            else {                                       /* Laplace beta/2 exp(-beta |v|): DC and AC part */
                const double uu = rng_uniform(rng) - 0.5;
                v = meanVel - (uu > 0 ? 1 : uu < 0 ? -1 : 0) * log(1 - 2 * fabs(uu)) / beta;
                C->dcRe[m] = re * sqrt(dcr / (1 + dcr));
                C->dcIm[m] = im * sqrt(dcr / (1 + dcr));
                re *= sqrt(1 / (1 + dcr));
                im *= sqrt(1 / (1 + dcr));
                C->hasDc = 1;
            }
            if (fixedClutter) { re = ampRms; im = 0; v = 0; }   /* fixed reflectivity, no motion (tests) */
            C->cellRow[m] = a;
            C->cellRange[m] = r;
            C->cellGroup[m] = group[q];
            C->rhoRe[m] = re;
            C->rhoIm[m] = im;
            C->fd[m] = 2 * v / s->lambda;
            m++;
        }
    }
    C->rowStart[C->nAz] = m;
    if (getenv("RADAR_GEN_DEBUG")) {                    /* expected power of the cells (debugging) */
        double sp = 0, sr = 0;
        for (int a = 0; a < C->nAz; a++)
            for (int r = 0; r < C->nR; r++) {
                const size_t q = (size_t)a * C->nR + r;
                if (!(sigma0[q] > 0)) continue;
                const double R = C->range[r];
                const double elRel = asin(-sind(s->tilt) * cosd(elDeg[q]) + cosd(s->tilt) * sind(elDeg[q])) * 180 / PI;
                const double gEl = pattern_gain(&s->txEl, elRel) * pattern_gain(&s->elElem, elRel);
                sp += 2 * s->radarConst * sigma0[q] * R * (C->azStep * PI / 180) * C->dR / pow(R, 4) * gEl;
            }
        for (int k = 0; k < m; k++) sr += C->rhoRe[k] * C->rhoRe[k] + C->rhoIm[k] * C->rhoIm[k];
        printf("debug: %d x %d cells, %d visible, groups", C->nAz, C->nR, m);
        for (int g = 0; g < C->G; g++) printf(" %.4f", C->groupEl[g]);
        printf(" deg; sum of the expected cell powers %.6e, drawn %.6e\n", sp, sr);
    }
    /* IF pulses and their FFTs; steering vectors of the groups on the carrier of every sub-pulse */
    C->sameCarrier = s->w[0].fc == s->w[1].fc;
    for (int i = 0; i < 2; i++) {
        const Wave *w = &s->w[i];
        C->fdScale[i] = w->fc / s->fc;
        C->ampScale[i] = w->lambda / s->lambda;
        C->nTx[i] = (int)lround(w->width * fs);
        int nFft = 1;
        while (nFft < s->nAdc + C->nTx[i] - 1) nFft <<= 1;
        C->nFft[i] = nFft;
        C->start0[i] = (int)lround(w->tTx * fs);
        C->steer[i] = xcalloc((size_t)2 * s->nCh * C->G, sizeof(double));
        for (int g = 0; g < C->G; g++) {
            const double u = -sind(s->tilt) * cosd(C->groupEl[g]) + cosd(s->tilt) * sind(C->groupEl[g]);
            for (int ch = 0; ch < s->nCh; ch++) {
                const double ph = 2 * PI * s->spacing_m / w->lambda * ch * u;
                C->steer[i][2 * (ch * C->G + g)] = cos(ph);
                C->steer[i][2 * (ch * C->G + g) + 1] = sin(ph);
            }
        }
        if (!w->enable) continue;
        C->txRe[i] = xcalloc((size_t)nFft, sizeof(double));
        C->txIm[i] = xcalloc((size_t)nFft, sizeof(double));
        for (int k = 0; k < C->nTx[i]; k++) {
            double er, ei;
            const double t = k / fs;
            envelope(w, t, &er, &ei);
            const double cr = cos(2 * PI * w->fIF * t), ci = sin(2 * PI * w->fIF * t);
            C->txRe[i][k] = er * cr - ei * ci;
            C->txIm[i][k] = er * ci + ei * cr;
        }
        Fft f;
        fft_init(&f, nFft);
        fft_run(&f, C->txRe[i], C->txIm[i], 0);
        fft_free(&f);
    }
    free(rType); free(rAz); free(rR); free(rH); free(typeName); free(regType); free(types);
    free(type); free(height); free(regionOf); free(elDeg); free(sigma0); free(visible); free(group);
}

/* ------------------------------------------------------------------------ */
/* One PRI                                                                  */
/* ------------------------------------------------------------------------ */

typedef struct {
    float  *x;                   /* ADC input [ch][nAdc] */
    int16_t *adc;                /* [ch][nAdc] = column-major [nAdc x nCh] */
    int16_t *I, *Q;              /* [ch][nDdc] */
    double *xc;                  /* clutter [page][g][nAdc] complex (re, im interleaved) */
    double *fr, *fi;             /* FFT work */
    double *reflRe, *reflIm;     /* [nR x G] */
    double *y1, *xI, *xQ;
    Fft    fft[2];
} Work;

static int16_t sat16(double v)
{
    return (int16_t)(v > 32767 ? 32767 : v < -32768 ? -32768 : v);
}

static double qround(double v, const Radar *s)   /* int16() and the word length (adc_sample.m quantize) */
{
    double q = round(v);
    q = q > 32767 ? 32767 : q < -32768 ? -32768 : q;
    if (s->nBits < 16) {
        if (q < -s->fullScale) q = -s->fullScale;
        if (q > s->fullScale - 1) q = s->fullScale - 1;
    }
    return q;
}

static void clutter_signal(const Clutter *C, const Radar *s, Work *wk, double az, double t_s)
{
    const int nPages = C->sameCarrier ? 1 : 2, G = C->G, nAdc = s->nAdc;
    memset(wk->xc, 0, sizeof(double) * 2 * (size_t)nPages * G * nAdc);
    for (int i = 0; i < 2; i++) {
        if (!s->w[i].enable) continue;
        const int pg = nPages == 1 ? 0 : i;
        memset(wk->reflRe, 0, sizeof(double) * (size_t)C->nR * G);
        memset(wk->reflIm, 0, sizeof(double) * (size_t)C->nR * G);
        int any = 0;
        for (int a = 0; a < C->nAz; a++) {
            const double azRel = wrap180(C->az[a] - az);
            if (!(fabs(azRel) <= s->azSpan)) continue;
            const double gRow = pattern_gain(&s->az, azRel);
            for (int m = C->rowStart[a]; m < C->rowStart[a + 1]; m++) {
                const double ph = 2 * PI * C->fdScale[i] * C->fd[m] * t_s, cr = cos(ph), ci = sin(ph);
                double wr = gRow * (C->rhoRe[m] * cr - C->rhoIm[m] * ci), wi = gRow * (C->rhoRe[m] * ci + C->rhoIm[m] * cr);
                if (C->hasDc) { wr += gRow * C->dcRe[m]; wi += gRow * C->dcIm[m]; }
                const size_t q = (size_t)C->cellRange[m] * G + C->cellGroup[m];
                wk->reflRe[q] += C->ampScale[i] * wr;
                wk->reflIm[q] += C->ampScale[i] * wi;
                any = 1;
            }
        }
        if (!any) continue;
        const int nTx = C->nTx[i], nFft = C->nFft[i];
        for (int g = 0; g < G; g++) {
            memset(wk->fr, 0, sizeof(double) * nFft);
            memset(wk->fi, 0, sizeof(double) * nFft);
            int used = 0;
            for (int r = 0; r < C->nR; r++) {
                const int n0 = C->start0[i] + C->delay[r];
                if (!(n0 > -nTx && n0 < nAdc)) continue;
                const size_t q = (size_t)r * G + g;
                if (wk->reflRe[q] == 0 && wk->reflIm[q] == 0) continue;
                wk->fr[n0 + nTx - 1] = wk->reflRe[q];
                wk->fi[n0 + nTx - 1] = wk->reflIm[q];
                used = 1;
            }
            if (!used) continue;
            fft_run(&wk->fft[i], wk->fr, wk->fi, 0);
            for (int k = 0; k < nFft; k++) {
                const double ar = wk->fr[k], ai = wk->fi[k], br = C->txRe[i][k], bi = C->txIm[i][k];
                wk->fr[k] = ar * br - ai * bi;
                wk->fi[k] = ar * bi + ai * br;
            }
            fft_run(&wk->fft[i], wk->fr, wk->fi, 1);
            double *xc = wk->xc + 2 * ((size_t)pg * G + g) * nAdc;
            for (int k = 0; k < nAdc; k++) {
                xc[2 * k] += wk->fr[nTx - 1 + k] / nFft;
                xc[2 * k + 1] += wk->fi[nTx - 1 + k] / nFft;
            }
        }
    }
}

typedef struct { int pulse, scan, staggerIdx; double t, az, pri; } Pri;

/* Received signal, ADC, DDC of one PRI (generate_rx_signal.m, adc_sample.m, ddc.m); returns the
 * fraction of the receive samples at full scale. */
static double simulate_pri(const Radar *s, const Clutter *C, const Target *tg, int nTg, const Pri *pr, uint64_t seed,
                           Work *wk)
{
    const int nAdc = s->nAdc, nCh = s->nCh;
    const double fs = s->fs, c = s->c;
    float *x = wk->x;
    memset(x, 0, sizeof(float) * (size_t)nAdc * nCh);
    /* 1) targets */
    for (int i = 0; i < nTg; i++) {
        double R, az, el, vr;
        target_geometry(&tg[i], pr->t, s->height, &R, &az, &el, &vr);
        const double azRel = wrap180(az - pr->az);
        if (fabs(azRel) > s->azSpan) continue;
        const double u = -sind(s->tilt) * cosd(el) * cosd(azRel) + cosd(s->tilt) * sind(el);
        const double elRel = asin(u) * 180 / PI;
        const double ga = pattern_gain(&s->az, azRel);
        const double gain = ga * ga * pattern_gain(&s->txEl, elRel) * pattern_gain(&s->elElem, elRel);
        /* RCS of this antenna pass (target_look.m) */
        const double rot = s->rate * pr->t, rel = mod360(az - s->az0);
        int look = (int)floor((rot - rel + 180) / 360) + 1;
        if (look < 1) look = 1;
        const double rcs = tg[i].rcsScan[(look < tg[i].nScan ? look : tg[i].nScan) - 1];
        for (int p = 0; p < 2; p++) {
            const Wave *w = &s->w[p];
            if (!w->enable) continue;
            const double amp = sqrt(2 * w->radarConst * rcs / pow(R, 4) * gain);
            const double tau0 = 2 * R / c;
            int first = (int)ceil((w->tTx + tau0) * fs) + 1, last = (int)floor((w->tTx + tau0 + w->width) * fs) + 1;
            if (first < 1) first = 1;
            if (last > nAdc) last = nAdc;
            if (first > last) continue;
            double *vre = wk->y1, *vim = wk->y1 + nCh;   /* element phases on its carrier */
            for (int ch = 0; ch < nCh; ch++) {
                const double ph = 2 * PI * s->spacing_m / w->lambda * ch * u;
                vre[ch] = cos(ph);
                vim[ch] = sin(ph);
            }
            for (int n = first; n <= last; n++) {
                const double t = (n - 1) / fs, tau = 2 * (R + vr * t) / c;
                double er, ei;
                envelope(w, t - w->tTx - tau0, &er, &ei);
                const double ph = 2 * PI * (w->fIF * (t - w->tTx - tau) - s->fLO * tau);
                const double cr = cos(ph), ci = sin(ph);
                const double zr = amp * (er * cr - ei * ci), zi = amp * (er * ci + ei * cr);
                for (int ch = 0; ch < nCh; ch++) x[(size_t)ch * nAdc + n - 1] += (float)(zr * vre[ch] - zi * vim[ch]);
            }
        }
    }
    /* 2) clutter */
    if (C) {
        clutter_signal(C, s, wk, pr->az, pr->t);
        const int nPages = C->sameCarrier ? 1 : 2;
        for (int pg = 0; pg < nPages; pg++) {
            const double *st = C->steer[nPages == 1 ? (s->w[0].enable ? 0 : 1) : pg];
            for (int g = 0; g < C->G; g++) {
                const double *xc = wk->xc + 2 * ((size_t)pg * C->G + g) * nAdc;
                for (int ch = 0; ch < nCh; ch++) {
                    const double sr = st[2 * (ch * C->G + g)], si = st[2 * (ch * C->G + g) + 1];
                    float *xr = x + (size_t)ch * nAdc;
                    for (int k = 0; k < nAdc; k++) xr[k] += (float)(xc[2 * k] * sr - xc[2 * k + 1] * si);
                }
            }
        }
    }
    /* 3) transmitter leakage */
    if (s->addLeak)
        for (int p = 0; p < 2; p++) {
            const Wave *w = &s->w[p];
            if (!w->enable) continue;
            for (int k = 0; k < nAdc; k++) {
                const double t = k / fs;
                if (!(t >= w->tTx && t < w->tTx + w->width)) continue;
                double er, ei;
                envelope(w, t - w->tTx, &er, &ei);
                const double ph = 2 * PI * w->fIF * (t - w->tTx);
                const float leak = (float)(s->leakAmp * (er * cos(ph) - ei * sin(ph)));
                for (int ch = 0; ch < nCh; ch++) x[(size_t)ch * nAdc + k] += leak;
            }
        }
    /* 4) STC, 5) thermal noise, 6) receiver gain, ADC */
    Rng rng;
    rng_seed(&rng, seed, (uint64_t)pr->pulse);
    const float gain = (float)s->voltGain, sigTh = (float)s->sigmaTh;
    long clipped = 0;
    for (int ch = 0; ch < nCh; ch++) {
        float *xr = x + (size_t)ch * nAdc;
        int16_t *ad = wk->adc + (size_t)ch * nAdc;
        for (int k = 0; k < nAdc; k++) {
            float v = xr[k];
            if (s->stcOn) v *= s->stc[k];
            if (s->addNoise) v += sigTh * (float)rng_normal(&rng);
            v *= gain;
            double q = qround((double)v * s->lsbPerUnit, s);
            if (s->addAdcNoise) q = qround((double)(float)q + s->sigmaAdc * rng_normal(&rng), s);
            ad[k] = (int16_t)q;
            if (s->rxSample[k] && (q >= s->fullScale - 1 || q <= -s->fullScale)) clipped++;
        }
    }
    /* DDC: mixer, two filter stages (valid convolution, decimation), shift, rounding, saturation */
    const int L1 = s->L1, L2 = s->L2, D1 = s->D1, D2 = s->D2;
    const int K1 = (nAdc - L1) / D1 + 1, K2 = s->nDdc;
    const double sc1 = ldexp(1.0, -s->shift1), sc2 = ldexp(1.0, -s->shift2);
    for (int ch = 0; ch < nCh; ch++) {
        const int16_t *ad = wk->adc + (size_t)ch * nAdc;
        for (int k = 0; k < nAdc; k++) { wk->xI[k] = ad[k] * s->cosNco[k]; wk->xQ[k] = ad[k] * s->msinNco[k]; }
        for (int line = 0; line < 2; line++) {
            const double *xin = line ? wk->xQ : wk->xI;
            double *y1 = wk->y1;
            for (int m = 0; m < K1; m++) {
                double acc = 0;
                const double *xp = xin + (size_t)m * D1 + L1 - 1;
                for (int j = 0; j < L1; j++) acc += s->h1[j] * xp[-j];
                double v = acc * sc1;
                if (s->fixedPoint) { v = round(v); v = v > 32767 ? 32767 : v < -32768 ? -32768 : v; }
                y1[m] = v;
            }
            int16_t *out = (line ? wk->Q : wk->I) + (size_t)ch * K2;
            for (int m = 0; m < K2; m++) {
                double acc = 0;
                const double *xp = y1 + (size_t)m * D2 + L2 - 1;
                for (int j = 0; j < L2; j++) acc += s->h2[j] * xp[-j];
                out[m] = sat16(round(acc * sc2));
            }
        }
    }
    return s->nRx ? (double)clipped / ((double)s->nRx * nCh) : 0.0;
}

/* ------------------------------------------------------------------------ */
/* Files                                                                    */
/* ------------------------------------------------------------------------ */

static void mkdir_p(const char *path)
{
    char buf[4096];
    snprintf(buf, sizeof buf, "%s", path);
    for (char *p = buf + 1; *p; p++)
        if (*p == '/' || *p == '\\') {
            const char c = *p;
            *p = 0;
            MKDIR(buf);
            *p = c;
        }
    MKDIR(buf);
}

/* pulse_header.m */
static void write_file(const char *folder, const char *prefix, int compressed, const Radar *s, const Pri *pr,
                       const char *varA, const int16_t *a, size_t rows, const char *varB, const int16_t *b)
{
    char path[4200];
    snprintf(path, sizeof path, "%s/%s_%07d.mat", folder, prefix, pr->pulse);
    MatWriter *w = mat_write_open(path, compressed);
    if (!w) { fprintf(stderr, "radar_generate_ddc: cannot write %s\n", path); exit(1); }
    mat_write_int16(w, varA, a, rows, (size_t)s->nCh);
    if (varB) mat_write_int16(w, varB, b, rows, (size_t)s->nCh);
    const Wave *ws = &s->w[0], *wl = &s->w[1];
    const double v_pulse = pr->pulse, v_scan = pr->scan, v_t = pr->t, v_az = pr->az, v_pri = pr->pri, v_st = pr->staggerIdx;
    const double v_fs = s->fs, v_fif = s->fIF, v_fc = s->fc, v_rx0 = 0.0;
    const double sTx[2] = { ws->tTx, ws->width }, lTx[2] = { wl->tTx, wl->width };
    const double v_dfs = s->fsOut, v_t0 = s->ddcT0;
    const double sIdx[2] = { ws->rxFirst, ws->rxLast }, lIdx[2] = { wl->rxFirst, wl->rxLast };
    const double sFc = ws->fc, lFc = wl->fc;
    static const char *const names[] = { "pulseNumber", "scanNumber", "time_s", "antennaAz_deg", "pri_s", "staggerIndex",
                                         "fs_Hz", "fIF_Hz", "fc_Hz", "rxStart_s", "shortTx_s", "longTx_s", "ddcFs_Hz",
                                         "ddcT0_s", "shortDdcIdx", "longDdcIdx", "shortFc_Hz", "longFc_Hz" };
    const double *vals[18] = { &v_pulse, &v_scan, &v_t, &v_az, &v_pri, &v_st, &v_fs, &v_fif, &v_fc, &v_rx0, sTx, lTx,
                               &v_dfs, &v_t0, sIdx, lIdx, &sFc, &lFc };
    const int n[18] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 1 };
    mat_write_struct(w, "hdr", 18, names, vals, n);
    if (mat_write_close(w) != 0) { fprintf(stderr, "radar_generate_ddc: error writing %s\n", path); exit(1); }
}

/* ------------------------------------------------------------------------ */
/* Main                                                                     */
/* ------------------------------------------------------------------------ */

typedef struct {
    const char *outDir, *targetsCsv;
    int raw, nThreads, maxPulses;
} Options;

static int run(const Options *opt)
{
    Radar s;
    Clutter clut, *C = NULL;
    Target *tg = NULL;
    int nTg = 0;
    const double t0 = now_s();
    setup(&s);
    const uint64_t seed = (uint64_t)(int64_t)num("sim.seed");
    Rng rng;
    rng_seed(&rng, seed, 0xFFFFFFFFFFFFULL);              /* setup stream: scenario, RCS, clutter map */
    /* pulse schedule (pulse_schedule.m) */
    const double tEnd = num("sim.nScans") * 360.0 / s.rate;
    double minPri = s.pri[0];
    for (int i = 1; i < s.nStagger; i++) if (s.pri[i] < minPri) minPri = s.pri[i];
    const int nMax = (int)ceil(tEnd / minPri) + 1;
    Pri *sched = xcalloc((size_t)nMax, sizeof(Pri));
    int nSched = 0, maxScan = 0;
    double tt = 0;
    for (int k = 0; k < nMax; k++) {
        const int st = k % s.nStagger;
        if (!(tt < tEnd)) break;
        Pri *p = &sched[nSched++];
        p->pulse = k + 1;
        p->staggerIdx = st + 1;
        p->pri = s.pri[st];
        p->t = tt;
        const double azU = s.az0 + s.rate * tt;
        p->az = mod360(azU);
        p->scan = (int)floor((azU - s.az0) / 360) + 1;
        if (p->scan > maxScan) maxScan = p->scan;
        tt += s.pri[st];
    }
    int nSect = 0, sc = 0;
    const double *sect = NULL;
    if (has("sim.sectors_deg") && !empty("sim.sectors_deg")) {
        const MValue *v = matrix("sim.sectors_deg", &nSect, &sc);
        if (sc != 2) fail("cfg.sim.sectors_deg must have rows [azStart azEnd]");
        sect = v->re;
    }
    int *active = xcalloc((size_t)nSched, sizeof(int)), nActive = 0;
    for (int k = 0; k < nSched; k++) {
        int a = nSect == 0;
        for (int i = 0; i < nSect && !a; i++) {
            double width = mod360(sect[nSect + i] - sect[i]);
            if (width == 0 && sect[nSect + i] != sect[i]) width = 360;
            a = mod360(sched[k].az - sect[i]) <= width;
        }
        if (a && (opt->maxPulses <= 0 || nActive < opt->maxPulses)) active[nActive++] = k;
    }
    /* targets: the scenario (or the file), RCS of every antenna pass (draw_target_rcs.m) */
    if (opt->targetsCsv) nTg = read_targets_csv(opt->targetsCsv, &tg);
    else nTg = scenario(str("scenario.type"), &rng, &tg);
    const int swerling = (int)num("sim.swerling");
    for (int i = 0; i < nTg; i++) {
        tg[i].nScan = maxScan + 1;
        tg[i].rcsScan = xcalloc((size_t)tg[i].nScan, sizeof(double));
        for (int k = 0; k < tg[i].nScan; k++) tg[i].rcsScan[k] = swerling == 1 ? tg[i].rcs * -log(rng_uniform(&rng)) : tg[i].rcs;
    }
    if (num("clutter.enable") != 0) {
        clutter_build(&clut, &s, &rng, (int)num("proc.taylorNbar"), num("proc.taylorSll_dB"));
        C = &clut;
    }
    /* folders */
    const char *ddcDir = opt->outDir ? opt->outDir : str("io.ddcFolder");
    const int saveRaw = opt->raw >= 0 ? opt->raw : num("io.saveRaw") != 0;
    const char *rawDir = str("io.rawFolder");
    const int compressed = same_text(str_or("io.matFormat", "-v6"), "-v7");
    mkdir_p(ddcDir);
    if (saveRaw) mkdir_p(rawDir);
    if (!g_quiet) {
        printf("Generating the DDC output of %d PRIs (%d element(s), DDC %g MHz, %d target(s)%s) -> %s\n", nActive, s.nCh,
               s.fsOut / 1e6, nTg, C ? ", terrain clutter" : "", ddcDir);
        for (int i = 0; i < 2; i++)
            printf("  %s pulse: %s\n", s.w[i].name, s.w[i].enable ? "on" : "off");
        if (C) printf("  clutter: %d visible cells, %d elevation group(s), CNR of land %.1f dB -> sigma0 %.1f dB\n",
                      C->nCells, C->G, C->cnr_dB, C->sigma0_dB);
        if (saveRaw) printf("  raw ADC data -> %s\n", rawDir);
        printf("  setup %.1f s\n", now_s() - t0);
        fflush(stdout);
    }
    /* the PRIs, in parallel */
    double *clip = xcalloc((size_t)(nActive ? nActive : 1), sizeof(double));
    int done = 0;
    const double tLoop = now_s();
#ifdef _OPENMP
    if (opt->nThreads > 0) omp_set_num_threads(opt->nThreads);
#endif
    #pragma omp parallel
    {
        Work wk;
        memset(&wk, 0, sizeof wk);
        const int nPages = C && !C->sameCarrier ? 2 : 1, G = C ? C->G : 1;
        int nFftMax = 1;
        wk.x = calloc((size_t)s.nAdc * s.nCh, sizeof(float));
        wk.adc = calloc((size_t)s.nAdc * s.nCh, sizeof(int16_t));
        wk.I = calloc((size_t)s.nDdc * s.nCh, sizeof(int16_t));
        wk.Q = calloc((size_t)s.nDdc * s.nCh, sizeof(int16_t));
        wk.xI = calloc((size_t)s.nAdc, sizeof(double));
        wk.xQ = calloc((size_t)s.nAdc, sizeof(double));
        wk.y1 = calloc((size_t)s.nAdc + 2 * (size_t)s.nCh, sizeof(double));
        if (C) {
            for (int i = 0; i < 2; i++) if (C->nFft[i] > nFftMax) nFftMax = C->nFft[i];
            wk.xc = calloc(2 * (size_t)nPages * G * s.nAdc, sizeof(double));
            wk.fr = calloc((size_t)nFftMax, sizeof(double));
            wk.fi = calloc((size_t)nFftMax, sizeof(double));
            wk.reflRe = calloc((size_t)C->nR * G, sizeof(double));
            wk.reflIm = calloc((size_t)C->nR * G, sizeof(double));
            for (int i = 0; i < 2; i++) if (s.w[i].enable) fft_init(&wk.fft[i], C->nFft[i]);
        }
        if (!wk.x || !wk.adc || !wk.I || !wk.Q || !wk.xI || !wk.xQ || !wk.y1 || (C && (!wk.xc || !wk.fr || !wk.fi))) {
            fprintf(stderr, "radar_generate_ddc: out of memory\n");
            exit(1);
        }
        #pragma omp for schedule(dynamic, 1)
        for (int j = 0; j < nActive; j++) {
            const Pri *pr = &sched[active[j]];
            clip[j] = simulate_pri(&s, C, tg, nTg, pr, seed, &wk);
            if (saveRaw) write_file(rawDir, "pulse", compressed, &s, pr, "adcData", wk.adc, (size_t)s.nAdc, NULL, NULL);
            write_file(ddcDir, "ddc", compressed, &s, pr, "I", wk.I, (size_t)s.nDdc, "Q", wk.Q);
            #pragma omp critical(progress)
            {                                            /* (OpenMP 2.0: also MSVC) */
                const int d = ++done;
                if (!g_quiet && (d % 500 == 0 || d == nActive)) {
                    printf("  %d of %d PRIs (%.0f s)\n", d, nActive, now_s() - tLoop);
                    fflush(stdout);
                }
            }
        }
        free(wk.x); free(wk.adc); free(wk.I); free(wk.Q); free(wk.xI); free(wk.xQ); free(wk.y1);
        free(wk.xc); free(wk.fr); free(wk.fi); free(wk.reflRe); free(wk.reflIm);
        for (int i = 0; i < 2; i++) if (wk.fft[i].n) fft_free(&wk.fft[i]);
    }
    double cm = 0, cx = 0;
    for (int j = 0; j < nActive; j++) { cm += clip[j]; if (clip[j] > cx) cx = clip[j]; }
    if (nActive) cm /= nActive;
    if (cx > s.overloadWarn)
        warn("ADC overload: up to %.1f%% of the receive samples of a PRI are at full scale; clipped clutter and targets "
             "give false plots. Reduce cfg.rx.gain_dB or %s.", 100 * cx,
             s.stcOn ? "raise the STC attenuation (cfg.rx.stc.curve, cfg.rx.stc.curveLong)" : "switch the STC on (cfg.rx.stc.enable = true)");
    const double tRun = now_s() - t0;
    if (!g_quiet) {
        int nThr = 1;
#ifdef _OPENMP
        nThr = omp_get_max_threads();
#endif
        printf("%d DDC files written in %.1f s (%.1f ms per PRI, %d thread(s)); ADC samples at full scale %.3f%% (mean), "
               "%.3f%% (worst PRI)\n", nActive, tRun, 1e3 * tRun / (nActive ? nActive : 1), nThr, 100 * cm, 100 * cx);
    }
    free(clip); free(active); free(sched);
    return 0;
}

static void usage(void)
{
    puts("usage: radar_generate_ddc [--config-m FILE] [--set STATEMENT]... [--out FOLDER] [--raw | --no-raw]\n"
         "                          [--threads N] [--targets FILE.csv] [--max-pulses N] [--quiet]\n"
         "  Simulates every active PRI of radar_config.m up to the DDC output and writes ddc_<n>.mat\n"
         "  (as generate_ddc_files.m).\n"
         "  --config-m FILE   MATLAB configuration (default " RADAR_SOURCE_DIR "/../radar_config.m)\n"
         "  --set STATEMENT   MATLAB assignment run after the file, e.g. --set \"cfg.sim.sectors_deg = [20 45];\"\n"
         "  --out FOLDER      the DDC files (default cfg.io.ddcFolder)\n"
         "  --raw, --no-raw   also / not the raw ADC data pulse_<n>.mat (default cfg.io.saveRaw, cfg.io.rawFolder)\n"
         "  --threads N       PRIs in parallel (default: all cores)\n"
         "  --targets FILE    targets from a CSV file: name,rcs_m2,x,y,z,vx,vy,vz[,tStart,turn_dps,accelAlong,accelUp]...\n"
         "  --max-pulses N    only the first N active PRIs\n"
         "  --quiet           print only errors and warnings");
}

int main(int argc, char **argv)
{
    const char *mPath = RADAR_SOURCE_DIR "/../radar_config.m";
    const char **sets = calloc((size_t)argc + 1, sizeof(char *));
    Options opt = { NULL, NULL, -1, 0, 0 };
    int nSets = 0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i], *v = i + 1 < argc ? argv[i + 1] : NULL;
        if (!strcmp(a, "--config-m") && v) { mPath = v; i++; }
        else if (!strcmp(a, "--set") && v) { sets[nSets++] = v; i++; }
        else if (!strcmp(a, "--out") && v) { opt.outDir = v; i++; }
        else if (!strcmp(a, "--targets") && v) { opt.targetsCsv = v; i++; }
        else if (!strcmp(a, "--threads") && v) { opt.nThreads = atoi(v); i++; }
        else if (!strcmp(a, "--max-pulses") && v) { opt.maxPulses = atoi(v); i++; }
        else if (!strcmp(a, "--raw")) opt.raw = 1;
        else if (!strcmp(a, "--no-raw")) opt.raw = 0;
        else if (!strcmp(a, "--quiet")) g_quiet = 1;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); free(sets); return 0; }
        else { fprintf(stderr, "unknown option %s\n", a); usage(); free(sets); return 2; }
    }
    char err[1024];
    MValue *cfg = mfile_eval(mPath, sets, nSets, err, sizeof err);
    free(sets);
    if (!cfg) { fprintf(stderr, "radar_generate_ddc: %s\n", err); return 1; }
    if (cfg->type != MV_STRUCT) {
        fprintf(stderr, "radar_generate_ddc: %s does not return a struct\n", mPath);
        mvalue_free(cfg);
        return 1;
    }
    g_cfg = cfg;
    if (!g_quiet) printf("Reading %s\n", mPath);
    if (setjmp(g_jb) != 0) {                            /* the errors of fail() */
        fprintf(stderr, "radar_generate_ddc: %s\n", g_err);
        mvalue_free(cfg);
        return 1;
    }
    const int r = run(&opt);
    mvalue_free(cfg);
    return r;
}

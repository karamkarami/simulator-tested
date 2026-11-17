/*
 * logger.c - logging of the processing chain (PPI tab "Controls", --log-dir):
 * the outputs of the blocks the user selects are written for every PRI as
 * MAT-files (-v6), each block in a folder of its own:
 *
 *   <folder>/<name>/radar_parameters.mat   every setting of the processing chain when the log
 *                                           started (struct params, params_mat.c)
 *   <folder>/<name>/log_info.txt           what the log holds; parameters changed during the log
 *   <folder>/<name>/screen_record.mp4      the PPI window with the microphone (LOG_SCREEN,
 *                                           recorded by ppi_qt.cpp)
 *   <folder>/<name>/ddc/ddc_<n>.mat        the DDC output as received: I, Q (int16 [nSamples x
 *                                           nCh]) and hdr, the format of the simulator's DDC files
 *                                           (save_ddc_file.m): the receiver can process them again
 *                                           (radar_stream_send <log>/ddc, replay_ddc_udp.m)
 *   <folder>/<name>/dbf/dbf_<n>.mat        the beamformer output of every sub-pulse: short, long
 *                                           (complex single [nBeams x nSamples]: the pulse
 *                                           compression input, zero outside the receive window),
 *                                           shortFirstSample, longFirstSample (1-based DDC sample of
 *                                           column 1), hdr
 *   <folder>/<name>/mf/mf_<n>.mat          the pulse compression output: short, long (complex single
 *                                           [nBeams x nRange]), shortRange_m, longRange_m, hdr (as
 *                                           cfg.io.saveMf of the simulator)
 *   <folder>/<name>/doppler/doppler_<n>.mat the range-Doppler maps of the CPI updates of the PRI:
 *                                           short, long (single [nBeams x nRange x nfft], after the
 *                                           noise-gain equalisation), f_Hz, shortV_mps, longV_mps,
 *                                           shortRange_m, longRange_m, hdr (as cfg.io.saveDoppler)
 *   <folder>/<name>/detection/detection_<n>.mat  the CFAR detections after the pruning (dets) and the
 *                                           plots closed in the PRI (plots), one row each, columns
 *                                           in 'columns', hdr
 *
 * <n> is the pulse number (7 digits). The processing thread only copies the
 * data of a PRI (logger_begin_pri ... logger_end_pri); a writer thread of its
 * own writes the files. The files of a PRI are queued when the next PRI begins:
 * then the PRI of the DDC header (hdr.pri_s, used to replay the DDC files with
 * their timing) is known. When the queue holds more than LOG_QUEUE_BYTES
 * (the disk cannot keep up), the PRI is dropped and counted.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200809L
#endif

#include "radar_rx.h"
#include "thread_compat.h"
#include "mat_write.h"

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#  include <direct.h>
#  define make_dir(p) _mkdir(p)
#else
#  include <sched.h>
#  include <sys/stat.h>
#  define make_dir(p) mkdir(p, 0775)
#endif

#define LOG_QUEUE_BYTES ((size_t)1024 * 1024 * 1024)    /* largest backlog of the writer */

typedef struct {
    int    have;
    int    nBeams, n, first;            /* DBF: nSamples, first; MF: nRange */
    float *re, *im;                     /* row-major [nBeams][n] as the processor holds them */
} CBlock;

typedef struct Job {
    uint32_t pulseNumber, scanNumber;
    uint16_t staggerIndex;
    double   time_s, pri_s, az_deg, ddcFs_Hz, ddcT0_s;
    uint32_t shortIdx[2], longIdx[2];
    int      mask;
    /* DDC */
    int      nSamples, nCh;
    int16_t *iq;
    /* DBF, MF of both sub-pulses */
    CBlock   dbf[N_PULSES], mf[N_PULSES];
    /* Doppler maps [nfft][nBeams][nRange] */
    int      haveDop[N_PULSES], dopNfft[N_PULSES], dopBeams[N_PULSES], dopRange[N_PULSES];
    float   *dop[N_PULSES];
    /* detections and plots: rows of LOG_COLS doubles */
    double  *dets, *plots;
    int      nDets, capDets, nPlots, capPlots;
    size_t   bytes;
    struct Job *next;
} Job;

#define LOG_COLS 13

struct Logger {
    char     dir[1024];
    int      mask;
    const RadarConfig *cfg;
    Job     *cur, *held;                /* the PRI being filled, the previous one (pri_s) */
    double   lastPri_s;
    /* queue of the writer */
    Job     *head, *tail;
    size_t   queuedBytes;
    int      stop;
    Mutex    lock;
    Cond     notEmpty;
#ifdef _WIN32
    HANDLE   thread;
#else
    pthread_t thread;
#endif
    /* statistics (under lock) */
    uint64_t pris, files, dropped, bytes, errors;
    double   writeSeconds;
    double   startTime;
    char     lastError[256];
};

static void job_free(Job *j)
{
    if (!j) return;
    free(j->iq);
    for (int i = 0; i < N_PULSES; i++) {
        free(j->dbf[i].re); free(j->dbf[i].im); free(j->mf[i].re); free(j->mf[i].im); free(j->dop[i]);
    }
    free(j->dets);
    free(j->plots);
    free(j);
}

/* ------------------------------------------------------------------------ */
/* Writer thread                                                            */
/* ------------------------------------------------------------------------ */

static void write_hdr(MatWriter *w, const Job *j)
{
    static const char *const names[10] = { "pulseNumber", "scanNumber", "time_s", "antennaAz_deg", "pri_s",
                                            "staggerIndex", "ddcFs_Hz", "ddcT0_s", "shortDdcIdx", "longDdcIdx" };
    const double pn = j->pulseNumber, sn = j->scanNumber, t = j->time_s, az = j->az_deg, pri = j->pri_s;
    const double st = j->staggerIndex, fs = j->ddcFs_Hz, t0 = j->ddcT0_s;
    const double si[2] = { j->shortIdx[0], j->shortIdx[1] }, li[2] = { j->longIdx[0], j->longIdx[1] };
    const double *const v[10] = { &pn, &sn, &t, &az, &pri, &st, &fs, &t0, si, li };
    const int n[10] = { 1, 1, 1, 1, 1, 1, 1, 1, 2, 2 };
    mat_write_struct(w, "hdr", 10, names, v, n);
}

static const char *const g_pulseName[N_PULSES] = { "short", "long" };

/* [nBeams][n] row-major -> MATLAB [nBeams x n] (column-major) */
static void to_matlab(const float *src, int nB, int n, float *dst)
{
    for (int b = 0; b < nB; b++)
        for (int k = 0; k < n; k++) dst[b + (size_t)nB * k] = src[(size_t)b * n + k];
}

static int finish(Logger *L, MatWriter *w, const char *path)
{
    if (!w) {
        mutex_lock(&L->lock);
        L->errors++;
        snprintf(L->lastError, sizeof L->lastError, "cannot create %s", path);
        mutex_unlock(&L->lock);
        return 0;
    }
    if (mat_write_close(w) != 0) {
        mutex_lock(&L->lock);
        L->errors++;
        snprintf(L->lastError, sizeof L->lastError, "cannot write %s (disk full?)", path);
        mutex_unlock(&L->lock);
        return 0;
    }
    return 1;
}

static void write_job(Logger *L, const Job *j)
{
    const RadarConfig *cfg = L->cfg;
    char path[1200];
    int files = 0;
    if ((j->mask & LOG_DDC) && j->iq) {
        snprintf(path, sizeof path, "%s/ddc/ddc_%07u.mat", L->dir, j->pulseNumber);
        MatWriter *w = mat_write_open(path, 0);
        if (w) {
            const size_t n = (size_t)j->nSamples * j->nCh;
            int16_t *I = malloc(sizeof(int16_t) * (n ? n : 1)), *Q = malloc(sizeof(int16_t) * (n ? n : 1));
            if (!I || !Q) abort();
            for (size_t k = 0; k < n; k++) { I[k] = j->iq[2 * k]; Q[k] = j->iq[2 * k + 1]; }
            mat_write_int16(w, "I", I, (size_t)j->nSamples, (size_t)j->nCh);
            mat_write_int16(w, "Q", Q, (size_t)j->nSamples, (size_t)j->nCh);
            free(I);
            free(Q);
            write_hdr(w, j);
        }
        files += finish(L, w, path);
    }
    for (int kind = 0; kind < 2; kind++) {                       /* DBF, MF */
        const int bit = kind ? LOG_MF : LOG_DBF;
        const CBlock *blk = kind ? j->mf : j->dbf;
        if (!(j->mask & bit) || (!blk[0].have && !blk[1].have)) continue;
        snprintf(path, sizeof path, "%s/%s/%s_%07u.mat", L->dir, kind ? "mf" : "dbf", kind ? "mf" : "dbf", j->pulseNumber);
        MatWriter *w = mat_write_open(path, 0);
        if (w) {
            for (int p = 0; p < N_PULSES; p++) {
                const CBlock *c = &blk[p];
                if (!c->have) continue;
                const size_t n = (size_t)c->nBeams * c->n;
                float *re = malloc(sizeof(float) * (n ? n : 1)), *im = malloc(sizeof(float) * (n ? n : 1));
                if (!re || !im) abort();
                to_matlab(c->re, c->nBeams, c->n, re);
                to_matlab(c->im, c->nBeams, c->n, im);
                const size_t dims[2] = { (size_t)c->nBeams, (size_t)c->n };
                mat_write_single(w, g_pulseName[p], re, im, dims, 2);
                free(re);
                free(im);
                char nm[40];
                if (kind == 0) {
                    const double first = c->first + 1;
                    snprintf(nm, sizeof nm, "%sFirstSample", g_pulseName[p]);
                    mat_write_double(w, nm, &first, 1, 1);
                } else if (cfg->pulse[p].range_m) {
                    snprintf(nm, sizeof nm, "%sRange_m", g_pulseName[p]);
                    mat_write_double(w, nm, cfg->pulse[p].range_m, 1, (size_t)cfg->pulse[p].nRange);
                }
            }
            write_hdr(w, j);
        }
        files += finish(L, w, path);
    }
    if ((j->mask & LOG_DOPPLER) && (j->haveDop[0] || j->haveDop[1])) {
        snprintf(path, sizeof path, "%s/doppler/doppler_%07u.mat", L->dir, j->pulseNumber);
        MatWriter *w = mat_write_open(path, 0);
        if (w) {
            int nfft = 0;
            for (int p = 0; p < N_PULSES; p++) {
                if (!j->haveDop[p]) continue;
                const int nB = j->dopBeams[p], nR = j->dopRange[p], nD = j->dopNfft[p];
                nfft = nD;
                const size_t n = (size_t)nB * nR * nD;
                float *x = malloc(sizeof(float) * (n ? n : 1));
                if (!x) abort();
                const float *src = j->dop[p];                    /* [d][b][r] -> [b + nB (r + nR d)] */
                for (int d = 0; d < nD; d++)
                    for (int b = 0; b < nB; b++)
                        for (int r = 0; r < nR; r++)
                            x[b + (size_t)nB * (r + (size_t)nR * d)] = src[((size_t)d * nB + b) * nR + r];
                const size_t dims[3] = { (size_t)nB, (size_t)nR, (size_t)nD };
                mat_write_single(w, g_pulseName[p], x, NULL, dims, 3);
                free(x);
                char nm[40];
                snprintf(nm, sizeof nm, "%sRange_m", g_pulseName[p]);
                mat_write_double(w, nm, cfg->pulse[p].range_m, 1, (size_t)cfg->pulse[p].nRange);
                double *v = malloc(sizeof(double) * (size_t)nD);
                if (!v) abort();
                for (int d = 0; d < nD; d++) v[d] = -cfg->dopFreq_Hz[d] * cfg->pulse[p].lambda_m / 2;
                snprintf(nm, sizeof nm, "%sV_mps", g_pulseName[p]);
                mat_write_double(w, nm, v, 1, (size_t)nD);
                free(v);
            }
            if (nfft == cfg->nfft) mat_write_double(w, "f_Hz", cfg->dopFreq_Hz, 1, (size_t)nfft);
            write_hdr(w, j);
        }
        files += finish(L, w, path);
    }
    if (j->mask & LOG_DETECTION) {
        snprintf(path, sizeof path, "%s/detection/detection_%07u.mat", L->dir, j->pulseNumber);
        MatWriter *w = mat_write_open(path, 0);
        if (w) {
            for (int t = 0; t < 2; t++) {
                const int n = t ? j->nPlots : j->nDets;
                const double *rows = t ? j->plots : j->dets;
                double *m = calloc((size_t)(n ? n : 1) * LOG_COLS, sizeof(double));
                if (!m) abort();
                for (int i = 0; i < n; i++)                      /* rows -> column-major [n x LOG_COLS] */
                    for (int c = 0; c < LOG_COLS; c++) m[i + (size_t)n * c] = rows[(size_t)i * LOG_COLS + c];
                mat_write_double(w, t ? "plots" : "dets", m, (size_t)n, LOG_COLS);
                free(m);
            }
            mat_write_char(w, "columns", "pulse (1 short, 2 long), pulseNumber, time_s, az_deg, range_m, el_deg, "
                                         "vr_mps, snr_dB, power, beam, rangeCell, dopplerBin, hits");
            write_hdr(w, j);
        }
        files += finish(L, w, path);
    }
    mutex_lock(&L->lock);
    L->files += (uint64_t)files;
    mutex_unlock(&L->lock);
}

static void writer(Logger *L)
{
    for (;;) {
        mutex_lock(&L->lock);
        while (!L->head && !L->stop) cond_wait_ms(&L->notEmpty, &L->lock, 200);
        Job *j = L->head;
        if (!j) { mutex_unlock(&L->lock); break; }       /* stopped and empty */
        L->head = j->next;
        if (!L->head) L->tail = NULL;
        mutex_unlock(&L->lock);
        const double t0 = now_seconds();
        write_job(L, j);
        const double dt = now_seconds() - t0;
        mutex_lock(&L->lock);
        L->queuedBytes -= j->bytes;
        L->bytes += j->bytes;
        L->pris++;
        L->writeSeconds += dt;
        mutex_unlock(&L->lock);
        job_free(j);
    }
}

#ifdef _WIN32
static DWORD WINAPI writer_main(LPVOID arg) { writer(arg); return 0; }
#else
static void *writer_main(void *arg) { writer(arg); return NULL; }
#endif

/* ------------------------------------------------------------------------ */
/* Processing thread                                                        */
/* ------------------------------------------------------------------------ */

static int make_dirs(const char *path)                   /* mkdir -p */
{
    char buf[1024];
    const size_t n = strlen(path);
    if (n == 0 || n >= sizeof buf) return -1;
    memcpy(buf, path, n + 1);
    for (size_t i = 1; i <= n; i++) {
        if (buf[i] == '/' || buf[i] == '\\' || buf[i] == '\0') {
            const char c = buf[i];
            buf[i] = '\0';
            if (!(i == 2 && buf[1] == ':') && make_dir(buf) != 0 && errno != EEXIST) return -1;
            buf[i] = c;
        }
    }
    return 0;
}

Logger *logger_start(const char *folder, const char *name, int mask, const RadarConfig *cfg, const RuntimeParams *rt,
                     char *err, size_t errLen)
{
    if (!folder || !folder[0]) { snprintf(err, errLen, "no log folder"); return NULL; }
    if (!(mask & (LOG_BLOCKS | LOG_SCREEN))) {
        snprintf(err, errLen, "nothing selected for the log");
        return NULL;
    }
    Logger *L = calloc(1, sizeof *L);
    if (!L) { snprintf(err, errLen, "out of memory"); return NULL; }
    char stamp[32];
    if (!name || !name[0]) {                             /* default name: the date and time */
        const time_t now = time(NULL);
        strftime(stamp, sizeof stamp, "log_%Y%m%d_%H%M%S", localtime(&now));
        name = stamp;
    }
    snprintf(L->dir, sizeof L->dir, "%s/%s", folder, name);
    static const struct { int bit; const char *sub; } subs[5] = {
        { LOG_DDC, "ddc" }, { LOG_DBF, "dbf" }, { LOG_MF, "mf" }, { LOG_DOPPLER, "doppler" }, { LOG_DETECTION, "detection" } };
    if (make_dirs(L->dir) != 0) {
        snprintf(err, errLen, "cannot create the folder %s", L->dir);
        free(L);
        return NULL;
    }
    for (int i = 0; i < 5; i++) {
        if (!(mask & subs[i].bit)) continue;
        char p[1200];
        snprintf(p, sizeof p, "%s/%s", L->dir, subs[i].sub);
        if (make_dirs(p) != 0) {
            snprintf(err, errLen, "cannot create the folder %s", p);
            free(L);
            return NULL;
        }
    }
    char started[64], blocks[96];
    {
        const time_t now = time(NULL);
        strftime(started, sizeof started, "%Y-%m-%d %H:%M:%S", localtime(&now));
        snprintf(blocks, sizeof blocks, "%s%s%s%s%s%s", mask & LOG_DDC ? " ddc" : "", mask & LOG_DBF ? " dbf" : "",
                 mask & LOG_MF ? " mf" : "", mask & LOG_DOPPLER ? " doppler" : "", mask & LOG_DETECTION ? " detection" : "",
                 mask & LOG_SCREEN ? " screen" : "");
    }
    {                                                    /* the settings of the processing chain now */
        char p[1200];
        RuntimeParams r;
        if (rt) r = *rt; else memset(&r, 0, sizeof r);
        r.started = started;
        r.logBlocks = blocks[0] ? blocks + 1 : "";
        r.screenRecord = (mask & LOG_SCREEN) != 0;
        snprintf(p, sizeof p, "%s/radar_parameters.mat", L->dir);
        if (radar_parameters_write(p, cfg, &r) != 0) {
            snprintf(err, errLen, "cannot write %s", p);
            free(L);
            return NULL;
        }
    }
    {                                                    /* what the log holds */
        char p[1200];
        snprintf(p, sizeof p, "%s/log_info.txt", L->dir);
        FILE *f = fopen(p, "w");
        if (f) {
            const time_t now = time(NULL);
            fprintf(f, "Radar receiver log, started %s", ctime(&now));
            fprintf(f, "blocks:%s\n", blocks);
            fprintf(f, "channels %d (processed %d), beams %d, DDC %d samples at %.0f Hz, short %d / long %d range cells, "
                       "Doppler %d bins (%s), CPI %d pulses, radar mode %d\n", cfg->nCh, cfg->nChUsed, cfg->nBeams,
                    cfg->nDdc, cfg->fsOut_Hz, cfg->pulse[0].nRange, cfg->pulse[1].nRange, cfg->nfft,
                    cfg->dopplerMethod == DOPPLER_NUFFT ? "NUFFT" : "FFT", cfg->cpiLength, cfg->radarMode);
            fprintf(f, "CFAR threshold %.2f dB at the start\n", 10.0 * log10(rt && rt->alpha > 0 ? rt->alpha : cfg->alpha));
            fprintf(f, "files: radar_parameters.mat (every setting at the start), <block>/<block>_<pulse number>.mat%s "
                       "(see logger.c / README)\n", mask & LOG_SCREEN ? ", screen_record.mp4" : "");
            fprintf(f, "changes during the log:\n");
            fclose(f);
        }
    }
    L->mask = mask;
    L->cfg = cfg;
    L->startTime = now_seconds();
    mutex_init(&L->lock);
    cond_init(&L->notEmpty);
#ifdef _WIN32
    L->thread = CreateThread(NULL, 0, writer_main, L, 0, NULL);
    const int ok = L->thread != NULL;
#else
    /* normal scheduling (not the real-time priority of the processing thread that starts it) */
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
    pthread_attr_setschedpolicy(&attr, SCHED_OTHER);
    struct sched_param sp;
    memset(&sp, 0, sizeof sp);
    pthread_attr_setschedparam(&attr, &sp);
    int ok = pthread_create(&L->thread, &attr, writer_main, L) == 0;
    if (!ok) ok = pthread_create(&L->thread, NULL, writer_main, L) == 0;
    pthread_attr_destroy(&attr);
#endif
    if (!ok) {
        snprintf(err, errLen, "cannot start the log writer thread");
        mutex_destroy(&L->lock);
        cond_destroy(&L->notEmpty);
        free(L);
        return NULL;
    }
    return L;
}

const char *logger_dir(const Logger *L) { return L ? L->dir : ""; }

void logger_note(Logger *L, const char *text)
{
    if (!L) return;
    char p[1200], stamp[32];
    snprintf(p, sizeof p, "%s/log_info.txt", L->dir);
    FILE *f = fopen(p, "a");
    if (!f) return;
    const time_t now = time(NULL);
    strftime(stamp, sizeof stamp, "%H:%M:%S", localtime(&now));
    fprintf(f, "  %s (%.1f s into the log): %s\n", stamp, now_seconds() - L->startTime, text);
    fclose(f);
}

int logger_mask(const Logger *L) { return L ? L->mask : 0; }

/* Queues a finished PRI (or drops it when the writer is too far behind). */
static void enqueue(Logger *L, Job *j)
{
    mutex_lock(&L->lock);
    if (L->queuedBytes + j->bytes > LOG_QUEUE_BYTES && L->head) {
        L->dropped++;
        mutex_unlock(&L->lock);
        job_free(j);
        return;
    }
    L->queuedBytes += j->bytes;
    if (L->tail) L->tail->next = j; else L->head = j;
    L->tail = j;
    cond_signal(&L->notEmpty);
    mutex_unlock(&L->lock);
}

/* The held PRI goes to the writer with the PRI to the next one (pri_s of its DDC header). */
static void release_held(Logger *L, const PriData *next)
{
    Job *h = L->held;
    if (!h) return;
    if (next && next->pulseNumber == h->pulseNumber + 1 && next->time_s > h->time_s) {
        h->pri_s = next->time_s - h->time_s;
        L->lastPri_s = h->pri_s;
    } else {
        h->pri_s = L->lastPri_s;                         /* a gap or the end: the previous PRI */
    }
    L->held = NULL;
    enqueue(L, h);
}

void logger_begin_pri(Logger *L, const PriData *pri)
{
    if (!L) return;
    release_held(L, pri);
    Job *j = calloc(1, sizeof *j);
    if (!j) return;
    j->pulseNumber = pri->pulseNumber;
    j->scanNumber = pri->scanNumber;
    j->staggerIndex = pri->staggerIndex;
    j->time_s = pri->time_s;
    j->az_deg = pri->az_deg;
    j->ddcFs_Hz = pri->ddcFs_Hz;
    j->ddcT0_s = pri->ddcT0_s;
    memcpy(j->shortIdx, pri->shortIdx, sizeof j->shortIdx);
    memcpy(j->longIdx, pri->longIdx, sizeof j->longIdx);
    j->mask = L->mask;
    if ((L->mask & LOG_DDC) && pri->iq) {                /* as received: every channel */
        const size_t n = 2 * (size_t)pri->nSamples * pri->nCh;
        j->iq = malloc(sizeof(int16_t) * (n ? n : 1));
        if (j->iq) {
            memcpy(j->iq, pri->iq, sizeof(int16_t) * n);
            j->nSamples = pri->nSamples;
            j->nCh = pri->nCh;
            j->bytes += sizeof(int16_t) * n;
        }
    }
    L->cur = j;
}

void logger_block(Logger *L, const BlockOutput *o)
{
    if (!L || !L->cur || o->pulse < 0 || o->pulse >= N_PULSES) return;
    Job *j = L->cur;
    if (o->kind == LOG_DOPPLER) {
        if (!(L->mask & LOG_DOPPLER)) return;
        const size_t n = (size_t)o->nfft * o->nBeams * o->nRange;
        float *x = realloc(j->dop[o->pulse], sizeof(float) * (n ? n : 1));
        if (!x) return;
        memcpy(x, o->re, sizeof(float) * n);
        j->dop[o->pulse] = x;
        j->haveDop[o->pulse] = 1;
        j->dopNfft[o->pulse] = o->nfft;
        j->dopBeams[o->pulse] = o->nBeams;
        j->dopRange[o->pulse] = o->nRange;
        j->bytes += sizeof(float) * n;
        return;
    }
    if (o->kind != LOG_DBF && o->kind != LOG_MF) return;
    if (!(L->mask & o->kind)) return;
    CBlock *c = o->kind == LOG_DBF ? &j->dbf[o->pulse] : &j->mf[o->pulse];
    const int n = o->kind == LOG_DBF ? o->nSamples : o->nRange;
    const size_t m = (size_t)o->nBeams * n;
    float *re = realloc(c->re, sizeof(float) * (m ? m : 1)), *im = re ? realloc(c->im, sizeof(float) * (m ? m : 1)) : NULL;
    if (re) c->re = re;
    if (!re || !im) return;
    c->im = im;
    memcpy(c->re, o->re, sizeof(float) * m);
    memcpy(c->im, o->im, sizeof(float) * m);
    c->have = 1;
    c->nBeams = o->nBeams;
    c->n = n;
    c->first = o->first;
    j->bytes += 2 * sizeof(float) * m;
}

static void add_row(double **rows, int *n, int *cap, const Detection *d)
{
    if (*n == *cap) {
        const int c = *cap ? 2 * *cap : 64;
        double *p = realloc(*rows, sizeof(double) * (size_t)c * LOG_COLS);
        if (!p) return;
        *rows = p;
        *cap = c;
    }
    double *r = *rows + (size_t)*n * LOG_COLS;
    r[0] = d->pulse + 1;  r[1] = d->pulseNumber;  r[2] = d->time_s;  r[3] = d->az_deg;  r[4] = d->range_m;
    r[5] = d->el_deg;  r[6] = d->vr_mps;  r[7] = d->snr_dB;  r[8] = d->power;  r[9] = d->beam + 1;
    r[10] = d->rIdx + 1;  r[11] = d->dIdx + 1;  r[12] = d->hits;
    (*n)++;
}

void logger_detections(Logger *L, const Detection *dets, int n)
{
    if (!L || !L->cur || !(L->mask & LOG_DETECTION)) return;
    for (int i = 0; i < n; i++) add_row(&L->cur->dets, &L->cur->nDets, &L->cur->capDets, &dets[i]);
    L->cur->bytes += sizeof(double) * LOG_COLS * (size_t)n;
}

void logger_plot(Logger *L, const Detection *plot)
{
    if (!L || !L->cur || !(L->mask & LOG_DETECTION)) return;
    add_row(&L->cur->plots, &L->cur->nPlots, &L->cur->capPlots, plot);
    L->cur->bytes += sizeof(double) * LOG_COLS;
}

void logger_end_pri(Logger *L)
{
    if (!L || !L->cur) return;
    L->held = L->cur;                                    /* queued when the next PRI begins */
    L->cur = NULL;
}

void logger_stats(Logger *L, LoggerStats *s)
{
    memset(s, 0, sizeof *s);
    if (!L) return;
    mutex_lock(&L->lock);
    s->pris = L->pris;
    s->files = L->files;
    s->dropped = L->dropped;
    s->bytes = L->bytes;
    s->errors = L->errors;
    s->queuedBytes = L->queuedBytes;
    s->writeSeconds = L->writeSeconds;
    snprintf(s->lastError, sizeof s->lastError, "%s", L->lastError);
    mutex_unlock(&L->lock);
    s->seconds = now_seconds() - L->startTime;
}

void logger_stop(Logger *L, LoggerStats *final)
{
    if (!L) return;
    if (L->cur) logger_end_pri(L);
    release_held(L, NULL);
    mutex_lock(&L->lock);
    L->stop = 1;
    cond_signal(&L->notEmpty);
    mutex_unlock(&L->lock);
#ifdef _WIN32
    WaitForSingleObject(L->thread, INFINITE);
    CloseHandle(L->thread);
#else
    pthread_join(L->thread, NULL);
#endif
    if (final) logger_stats(L, final);
    mutex_destroy(&L->lock);
    cond_destroy(&L->notEmpty);
    free(L);
}

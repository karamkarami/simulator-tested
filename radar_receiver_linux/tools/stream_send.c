/*
 * stream_send.c - sends recorded DDC data over UDP with the real PRI timing,
 * to test the receiver in real time without the MATLAB simulator.
 *
 *   radar_stream_send SOURCE [--host 127.0.0.1] [--port 50000] [--speed 1.0] [--loops 1]
 *                     [--payload 1400] [--spread 0] [--no-gaps] [--catch-up]
 *                     [--queue 256] [--timing-csv FILE] [--quiet]
 *
 * SOURCE is either
 *   - a folder of DDC files ddc_<n>.mat (run_radar_sim with cfg.io.saveDdc = true,
 *     saved with -v6, or -v7 when built with zlib): every PRI is split into
 *     the datagrams of the live stream with --payload bytes of samples each, or
 *   - a stream recorded with radar_receiver --record FILE: the recorded
 *     datagrams are sent as they are.
 * The data are not processed: the receiver gets the same datagrams as from
 * the simulator.
 *
 * Timing: the send time of every PRI comes from the PRI start time in its
 * header divided by --speed (0 = as fast as possible), on one monotonic clock
 * (no accumulated drift). The sender sleeps until shortly before a send time
 * and busy-waits for the rest, so a PRI leaves within microseconds of its
 * time. A reader thread reads the files ahead into a queue of --queue PRIs
 * (filled before the start), so reading never delays a PRI.
 *   --spread F   the datagrams of a PRI are sent evenly over the fraction F of
 *                the PRI (0 = all at the PRI start, as the simulator does)
 *   --no-gaps    a skipped sector (gap in the pulse numbers) lasts one PRI
 *                instead of its real duration
 *   --catch-up   after a late PRI, the next ones are sent at once until the
 *                original schedule is reached; by default they keep the PRI
 *                spacing (no burst for the receiver)
 * With --loops > 1 the pulse numbers and times of later loops are continued,
 * so the receiver sees one long stream. --timing-csv writes the scheduled
 * and real send time of every PRI.
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <winsock2.h>
#  include <ws2tcpip.h>
typedef SOCKET sock_t;
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <arpa/inet.h>
#  include <netinet/in.h>
#  include <sys/socket.h>
#  include <sys/stat.h>
#  include <time.h>
#  include <unistd.h>
typedef int sock_t;
#endif

#include "radar_rx.h"
#include "thread_compat.h"
#include "ddc_files.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RDDC_MAGIC    0x52444443u
#define MAX_DATAGRAM  65536
#define LATE_S        1e-4          /* a PRI sent later than this counts as late          */
#ifdef _WIN32
#  define SPIN_S      20e-3         /* Sleep() can wake up late by a timer tick (15.6 ms) */
#else
#  define SPIN_S      2e-3          /* busy-wait before every send time                   */
#endif

/* ------------------------------------------------------------------------ */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------ */

static uint32_t get_u32(const unsigned char *p) { return p[0] | p[1] << 8 | p[2] << 16 | (uint32_t)p[3] << 24; }
static void put_u32(unsigned char *p, uint32_t v) { p[0] = v & 255; p[1] = v >> 8 & 255; p[2] = v >> 16 & 255; p[3] = v >> 24; }

static void sleep_seconds(double dt)
{
#ifdef _WIN32
    Sleep((DWORD)(dt * 1e3));
#else
    struct timespec ts = { (time_t)dt, (long)((dt - (double)(time_t)dt) * 1e9) };
    nanosleep(&ts, NULL);
#endif
}

/* Waits until now_seconds() >= t: sleeps for the coarse part, busy-waits for the last SPIN_S. */
static void wait_until(double t)
{
    double dt;
    while ((dt = t - now_seconds()) > SPIN_S)
        sleep_seconds(dt - SPIN_S);
    while (now_seconds() < t) {
    }
}

static int is_folder(const char *path)
{
#ifdef _WIN32
    const DWORD a = GetFileAttributesA(path);
    return a != INVALID_FILE_ATTRIBUTES && (a & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat st;
    return stat(path, &st) == 0 && S_ISDIR(st.st_mode);
#endif
}

static int cmp_double(const void *a, const void *b)
{
    const double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y;
}

/* p-th percentile (0..100) of n values (sorts v). */
static double percentile(double *v, size_t n, double p)
{
    if (n == 0) return 0;
    qsort(v, n, sizeof *v, cmp_double);
    size_t i = (size_t)(p / 100.0 * (double)(n - 1) + 0.5);
    return v[i < n ? i : n - 1];
}

/* ------------------------------------------------------------------------ */
/* PRIs and the two sources                                                  */
/* ------------------------------------------------------------------------ */

typedef struct {
    unsigned char *data;        /* the datagrams back to back                   */
    size_t         cap;
    uint32_t      *len;         /* length of every datagram                     */
    int            lenCap, nPkt;
    uint32_t       pulse;
    double         time_s, pri_s;
} Pri;

static int pri_reserve(Pri *p, size_t bytes, int nPkt)
{
    if (bytes > p->cap) {
        unsigned char *d = realloc(p->data, bytes);
        if (!d) return -1;
        p->data = d;
        p->cap = bytes;
    }
    if (nPkt > p->lenCap) {
        uint32_t *l = realloc(p->len, (size_t)nPkt * sizeof *l);
        if (!l) return -1;
        p->len = l;
        p->lenCap = nPkt;
    }
    return 0;
}

typedef struct {
    int       folder;           /* 1: ddc_*.mat files, 0: recorded stream       */
    const char *path;
    int       payload;
    /* folder */
    char    **files;
    int       nFiles, next;
    DdcPri    ddc;
    /* recorded stream */
    FILE     *f;
    unsigned char pending[MAX_DATAGRAM];   /* first datagram of the next PRI    */
    uint32_t  pendingLen;
    int       havePending;
    /* loops: pulse numbers and times continue */
    int       loops, loop;
    uint32_t  pulseOffset, lastPulse;
    double    timeShift, t0, lastTime, lastPri;
    int       haveT0;
} Source;

static int source_open(Source *s, const char *path, int payload, int loops, char *err, size_t errLen)
{
    memset(s, 0, sizeof *s);
    s->path = path;
    s->payload = payload;
    s->loops = loops > 0 ? loops : 1;
    s->lastPri = 1e-3;
    s->folder = is_folder(path);
    if (s->folder) {
        s->nFiles = ddc_list_files(path, &s->files, err, errLen);
        return s->nFiles > 0 ? 0 : -1;
    }
    s->f = fopen(path, "rb");
    if (!s->f) {
        snprintf(err, errLen, "cannot open %s", path);
        return -1;
    }
    return 0;
}

static void source_close(Source *s)
{
    if (s->f) fclose(s->f);
    ddc_free_list(s->files, s->nFiles);
    ddc_free(&s->ddc);
}

/* Reads the next datagram of the recorded stream. Returns its length, or -1 at the end. */
static long read_datagram(FILE *f, unsigned char *buf)
{
    unsigned char len4[4];
    if (fread(len4, 1, 4, f) != 4) return -1;
    const uint32_t len = get_u32(len4);
    if (len > MAX_DATAGRAM || fread(buf, 1, len, f) != len) return -1;
    return (long)len;
}

static int is_rddc(const unsigned char *d, uint32_t len) { return len >= 64 && get_u32(d) == RDDC_MAGIC; }

static double header_time(const unsigned char *d) { double t; memcpy(&t, d + 32, 8); return t; }

/* Next PRI of the recorded stream: the datagrams up to the next pulse number
 * (foreign datagrams go with the PRI they are in). Returns 1, or 0 at the end. */
static int next_recorded(Source *s, Pri *p, char *err, size_t errLen)
{
    size_t used = 0;
    int have = 0;
    p->nPkt = 0;
    for (;;) {
        unsigned char *d;
        uint32_t len;
        if (s->havePending) {
            d = s->pending;
            len = s->pendingLen;
            s->havePending = 0;
        } else {
            if (pri_reserve(p, used + MAX_DATAGRAM, p->nPkt + 1)) goto oom;
            long n = read_datagram(s->f, p->data + used);
            if (n < 0) break;
            d = p->data + used;
            len = (uint32_t)n;
            if (have && is_rddc(d, len) && get_u32(d + 4) != p->pulse) {   /* next PRI begins */
                memcpy(s->pending, d, len);
                s->pendingLen = len;
                s->havePending = 1;
                p->pri_s = header_time(d) - p->time_s;
                return 1;
            }
        }
        if (pri_reserve(p, used + len, p->nPkt + 1)) goto oom;
        if (d == s->pending) memcpy(p->data + used, d, len);
        d = p->data + used;
        if (is_rddc(d, len) && !have) {
            have = 1;
            p->pulse = get_u32(d + 4);
            p->time_s = header_time(d);
        }
        p->len[p->nPkt++] = len;
        used += len;
    }
    if (p->nPkt == 0) return 0;
    if (!have) p->time_s = s->lastTime;               /* only foreign datagrams */
    p->pri_s = s->lastPri;                            /* last PRI of the file */
    return 1;
oom:
    snprintf(err, errLen, "out of memory");
    return -1;
}

/* Next PRI of the DDC files. Returns 1, 0 at the end or -1 on an error. */
static int next_file(Source *s, Pri *p, char *err, size_t errLen)
{
    size_t bytes;
    if (s->next >= s->nFiles) return 0;
    if (ddc_load(s->files[s->next++], &s->ddc, err, errLen)) return -1;
    const int n = ddc_datagram_count(&s->ddc, s->payload, &bytes);
    if (n > 65535) {                                  /* the packet index and count are 16-bit */
        snprintf(err, errLen, "%d datagrams per PRI with --payload %d (at most 65535): use a larger payload",
                 n, s->payload);
        return -1;
    }
    if (pri_reserve(p, bytes, n)) {
        snprintf(err, errLen, "out of memory");
        return -1;
    }
    p->nPkt = ddc_make_datagrams(&s->ddc, s->payload, p->data, p->len);
    p->pulse = s->ddc.pulseNumber;
    p->time_s = s->ddc.time_s;
    p->pri_s = s->ddc.pri_s;
    return 1;
}

/* Next PRI of the source (all loops). Returns 1, 0 at the end or -1 on an error. */
static int source_next(Source *s, Pri *p, char *err, size_t errLen)
{
    for (;;) {
        const int r = s->folder ? next_file(s, p, err, errLen) : next_recorded(s, p, err, errLen);
        if (r < 0) return -1;
        if (r == 1) break;
        if (++s->loop >= s->loops) return 0;
        /* next loop: start again, continue the pulse numbers and the time */
        s->pulseOffset = s->lastPulse;
        s->timeShift = s->lastTime + s->lastPri - s->t0;
        if (s->folder) {
            s->next = 0;
        } else {
            rewind(s->f);
            s->havePending = 0;
        }
    }
    if (!s->haveT0) {
        s->t0 = p->time_s;
        s->haveT0 = 1;
    }
    if (s->loop > 0) {                                 /* rewrite the headers of later loops */
        p->pulse += s->pulseOffset;
        p->time_s += s->timeShift;
        size_t off = 0;
        for (int i = 0; i < p->nPkt; off += p->len[i++]) {
            unsigned char *d = p->data + off;
            if (!is_rddc(d, p->len[i])) continue;
            put_u32(d + 4, get_u32(d + 4) + s->pulseOffset);
            const double t = header_time(d) + s->timeShift;
            memcpy(d + 32, &t, 8);
        }
    }
    s->lastPulse = p->pulse;
    s->lastTime = p->time_s;
    if (p->pri_s > 0) s->lastPri = p->pri_s;
    return 1;
}

/* ------------------------------------------------------------------------ */
/* Reader thread: fills a queue of PRIs ahead of the sender                  */
/* ------------------------------------------------------------------------ */

typedef struct {
    Source   *src;
    Pri      *slots;
    int       capacity, head, count, held;
    int       ended, failed, stop;
    char      err[512];
    Mutex     lock;
    Cond      notEmpty, notFull;
#ifdef _WIN32
    HANDLE    thread;
#else
    pthread_t thread;
#endif
} Reader;

static void reader_run(Reader *r)
{
    for (;;) {
        mutex_lock(&r->lock);
        while (!r->stop && r->count + r->held >= r->capacity)
            cond_wait_ms(&r->notFull, &r->lock, 200);
        if (r->stop) { mutex_unlock(&r->lock); return; }
        Pri *slot = &r->slots[(r->head + r->held + r->count) % r->capacity];
        mutex_unlock(&r->lock);

        char err[512] = "";
        const int res = source_next(r->src, slot, err, sizeof err);   /* without the lock */

        mutex_lock(&r->lock);
        if (res == 1) {
            r->count++;
        } else {
            r->ended = 1;
            if (res < 0) {
                r->failed = 1;
                memcpy(r->err, err, sizeof err);
            }
        }
        cond_signal(&r->notEmpty);
        mutex_unlock(&r->lock);
        if (res != 1) return;
    }
}

#ifdef _WIN32
static DWORD WINAPI reader_main(LPVOID arg) { reader_run(arg); return 0; }
#else
static void *reader_main(void *arg) { reader_run(arg); return NULL; }
#endif

/* Releases the PRI taken before and returns the next one, or NULL at the end. */
static Pri *reader_next(Reader *r, int *waited)
{
    Pri *p = NULL;
    mutex_lock(&r->lock);
    if (r->held) {
        r->held = 0;
        r->head = (r->head + 1) % r->capacity;
        cond_signal(&r->notFull);
    }
    *waited = 0;
    while (r->count == 0 && !r->ended) {
        *waited = 1;
        cond_wait_ms(&r->notEmpty, &r->lock, 100);
    }
    if (r->count > 0) {
        p = &r->slots[r->head];
        r->count--;
        r->held = 1;
    }
    mutex_unlock(&r->lock);
    return p;
}

/* Waits until the queue is full or the source has ended. */
static void reader_wait_full(Reader *r)
{
    mutex_lock(&r->lock);
    while (r->count < r->capacity && !r->ended)
        cond_wait_ms(&r->notEmpty, &r->lock, 100);
    mutex_unlock(&r->lock);
}

/* ------------------------------------------------------------------------ */
/* Main                                                                      */
/* ------------------------------------------------------------------------ */

static void usage(void)
{
    fprintf(stderr,
            "usage: radar_stream_send SOURCE [--host H] [--port N] [--speed X] [--loops N] [--payload BYTES]\n"
            "                         [--spread F] [--no-gaps] [--catch-up] [--queue N] [--timing-csv FILE] [--quiet]\n"
            "  SOURCE: a folder of ddc_<n>.mat files (cfg.io.saveDdc) or a stream recorded with radar_receiver --record\n");
}

int main(int argc, char **argv)
{
    const char *path = NULL, *host = "127.0.0.1", *timingCsv = NULL;
    int port = 50000, loops = 1, payload = 1400, queue = 256, keepGaps = 1, catchUp = 0, quiet = 0;
    double speed = 1.0, spread = 0.0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--host") && i + 1 < argc) host = argv[++i];
        else if (!strcmp(a, "--port") && i + 1 < argc) port = atoi(argv[++i]);
        else if (!strcmp(a, "--speed") && i + 1 < argc) speed = atof(argv[++i]);
        else if (!strcmp(a, "--loops") && i + 1 < argc) loops = atoi(argv[++i]);
        else if (!strcmp(a, "--payload") && i + 1 < argc) payload = atoi(argv[++i]);
        else if (!strcmp(a, "--spread") && i + 1 < argc) spread = atof(argv[++i]);
        else if (!strcmp(a, "--queue") && i + 1 < argc) queue = atoi(argv[++i]);
        else if (!strcmp(a, "--timing-csv") && i + 1 < argc) timingCsv = argv[++i];
        else if (!strcmp(a, "--no-gaps")) keepGaps = 0;
        else if (!strcmp(a, "--catch-up")) catchUp = 1;
        else if (!strcmp(a, "--quiet")) quiet = 1;
        else if (a[0] != '-' && !path) path = a;
        else { usage(); return 2; }
    }
    if (!path) { usage(); return 2; }
    if (payload < 4 || payload > MAX_DATAGRAM - 64 - 8 - 20) {
        fprintf(stderr, "--payload must be 4 ... %d bytes\n", MAX_DATAGRAM - 64 - 8 - 20);
        return 2;
    }
    if (queue < 2) queue = 2;
    if (spread < 0) spread = 0;
    if (spread > 1) spread = 1;

    char err[512];
    Source src;
    if (source_open(&src, path, payload, loops, err, sizeof err)) {
        fprintf(stderr, "%s\n", err);
        return 1;
    }

#ifdef _WIN32
    WSADATA wsa;
    WSAStartup(MAKEWORD(2, 2), &wsa);
#endif
    sock_t s = socket(AF_INET, SOCK_DGRAM, 0);
    int sndbuf = 32 << 20;
    setsockopt(s, SOL_SOCKET, SO_SNDBUF, (const char *)&sndbuf, sizeof sndbuf);
    struct sockaddr_in to;
    memset(&to, 0, sizeof to);
    to.sin_family = AF_INET;
    to.sin_port = htons((unsigned short)port);
    if (inet_pton(AF_INET, host, &to.sin_addr) != 1) {
        fprintf(stderr, "invalid host address %s (an IPv4 address is needed)\n", host);
        return 2;
    }

    Reader rd;
    memset(&rd, 0, sizeof rd);
    rd.src = &src;
    rd.capacity = queue;
    rd.slots = calloc((size_t)queue, sizeof *rd.slots);
    if (!rd.slots) { fprintf(stderr, "out of memory\n"); return 1; }
    mutex_init(&rd.lock);
    cond_init(&rd.notEmpty);
    cond_init(&rd.notFull);
#ifdef _WIN32
    rd.thread = CreateThread(NULL, 0, reader_main, &rd, 0, NULL);
#else
    pthread_create(&rd.thread, NULL, reader_main, &rd);
#endif
    if (!quiet)
        printf("Reading %s (%s) ...\n", path, src.folder ? "DDC files" : "recorded stream");
    reader_wait_full(&rd);
    if (rd.failed) {
        fprintf(stderr, "%s\n", rd.err);
        return 1;
    }
    if (!quiet)
        printf("Sending to %s:%d (speed %g, %d PRIs read ahead)\n", host, port, speed, rd.count);

    size_t logCap = 4096, nPri = 0;
    double *sched = malloc(logCap * sizeof *sched), *sent = malloc(logCap * sizeof *sent);
    uint32_t *pulses = malloc(logCap * sizeof *pulses);
    uint64_t datagrams = 0, bytes = 0, sendErrors = 0, starved = 0;
    double tPri = 0;
    uint32_t prevPulse = 0;
    double prevTime = 0, prevPri = 0;
    const double tStart = now_seconds() + 0.01;
    int waited;
    Pri *p;
    while ((p = reader_next(&rd, &waited)) != NULL) {
        if (nPri == 0) {
            tPri = tStart;
        } else if (speed > 0) {
            double dt = p->time_s - prevTime;
            if (!keepGaps && p->pulse != prevPulse + 1) dt = prevPri;   /* skipped sector -> one PRI */
            tPri += dt / speed;
        }
        if (waited && speed > 0 && now_seconds() > tPri + LATE_S) starved++;   /* the reader was behind */
        prevPulse = p->pulse;
        prevTime = p->time_s;
        prevPri = p->pri_s;
        const double dtPkt = speed > 0 && p->nPkt > 0 ? spread * p->pri_s / speed / p->nPkt : 0;

        size_t off = 0;
        double tSent = 0;
        for (int j = 0; j < p->nPkt; off += p->len[j++]) {
            if (speed > 0) wait_until(tPri + j * dtPkt);
            if (j == 0) {
                tSent = now_seconds();
                if (speed <= 0) tPri = tSent;
                if (nPri == logCap) {
                    logCap *= 2;
                    sched = realloc(sched, logCap * sizeof *sched);
                    sent = realloc(sent, logCap * sizeof *sent);
                    pulses = realloc(pulses, logCap * sizeof *pulses);
                    if (!sched || !sent || !pulses) { fprintf(stderr, "out of memory\n"); return 1; }
                }
                sched[nPri] = tPri;
                sent[nPri] = tSent;
                pulses[nPri] = p->pulse;
                if (!catchUp && tSent - tPri > LATE_S) tPri = tSent;    /* keep the PRI spacing from here on */
            }
            if (sendto(s, (const char *)p->data + off, (int)p->len[j], 0, (struct sockaddr *)&to, sizeof to) < 0)
                sendErrors++;
            datagrams++;
            bytes += p->len[j];
        }
        nPri++;
    }
    const double runTime = nPri ? now_seconds() - sent[0] : 0;

#ifdef _WIN32
    WaitForSingleObject(rd.thread, INFINITE);
#else
    pthread_join(rd.thread, NULL);
#endif
    if (rd.failed) fprintf(stderr, "%s\n", rd.err);

    /* Timing statistics */
    if (timingCsv) {
        FILE *f = fopen(timingCsv, "w");
        if (f) {
            fprintf(f, "pulseNumber,scheduled_s,sent_s,delay_us\n");
            for (size_t k = 0; k < nPri; k++)
                fprintf(f, "%u,%.9f,%.9f,%.3f\n", pulses[k], sched[k] - sent[0], sent[k] - sent[0],
                        1e6 * (sent[k] - sched[k]));
            fclose(f);
        } else {
            fprintf(stderr, "cannot write %s\n", timingCsv);
        }
    }
    size_t nLate = 0;
    double maxLate = 0;
    double *tmp = malloc((nPri ? nPri : 1) * sizeof *tmp);
    for (size_t k = 0; k < nPri; k++) {
        tmp[k] = sent[k] - sched[k];
        if (tmp[k] > LATE_S) nLate++;
        if (tmp[k] > maxLate) maxLate = tmp[k];
    }
    const double medLate = percentile(tmp, nPri, 50);
    /* error of the PRI intervals: real minus scheduled spacing of consecutive PRIs */
    size_t nInt = 0;
    for (size_t k = 1; k < nPri; k++) {
        const double e = (sent[k] - sent[k - 1]) - (sched[k] - sched[k - 1]);
        tmp[nInt++] = e < 0 ? -e : e;
    }
    const double int99 = percentile(tmp, nInt, 99), intMax = nInt ? tmp[nInt - 1] : 0;
    printf("Sent %llu datagrams, %llu PRIs (%.1f MB) in %.2f s (%.0f PRIs/s, %.2f Gbit/s)\n",
           (unsigned long long)datagrams, (unsigned long long)nPri, bytes / 1048576.0, runTime,
           nPri / (runTime > 0 ? runTime : 1), 8e-9 * bytes / (runTime > 0 ? runTime : 1));
    if (speed > 0)
        printf("Timing: %llu PRIs later than %.1f ms (max %.3f ms, median %.1f us); PRI interval error "
               "99%% %.1f us, max %.1f us\n", (unsigned long long)nLate, LATE_S * 1e3, maxLate * 1e3,
               medLate * 1e6, int99 * 1e6, intMax * 1e6);
    if (starved)
        printf("%llu PRIs were late because reading the files was slower than the stream: "
               "use a larger --queue or a lower --speed\n", (unsigned long long)starved);
    if (sendErrors)
        printf("Send errors: %llu datagrams\n", (unsigned long long)sendErrors);

    free(tmp);
    free(sched);
    free(sent);
    free(pulses);
    for (int i = 0; i < queue; i++) {
        free(rd.slots[i].data);
        free(rd.slots[i].len);
    }
    free(rd.slots);
    source_close(&src);
    return rd.failed ? 1 : 0;
}

/*
 * test_serial.c - tests of the serial output (serial_out.c).
 *
 *  1. CRC-16/CCITT-FALSE of "123456789" is 0x29B1.
 *  2. A detections frame and a plot frame are encoded and decoded field by field.
 *  3. Linux only: a pseudo-terminal pair stands in for the serial port. Frames
 *     are sent with serial_open/serial_send_* and read back from the master
 *     side: exact bytes, the split into frames of at most 200 detections,
 *     dropping when the ring is full (without blocking), the flush in
 *     serial_close and its time limit.
 *
 * Build (Linux):
 *   gcc -std=c11 -O2 -Wall -Wextra -Isrc src/serial_out.c tests/test_serial.c -o test_serial
 * Prints PASS and returns 0 when all checks pass.
 */
#if !defined(_WIN32)
#define _XOPEN_SOURCE 700     /* posix_openpt, grantpt, unlockpt, ptsname */
#define _DEFAULT_SOURCE       /* B921600 */
#endif

#include "radar_rx.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __linux__
#include <fcntl.h>
#include <poll.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>
#endif

#define N_DETS        450
#define MAX_PER_FRAME 200
#define FULL_FRAME    (16 + 24 * MAX_PER_FRAME)   /* detections frame with 200 records */

static int g_failures = 0;
static Detection g_dets[N_DETS];

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("  check failed (line %d): %s\n", __LINE__, #cond);       \
            g_failures++;                                                    \
        }                                                                    \
    } while (0)

/* ------------------------------------------------------------------------ */
/* Helpers                                                                   */
/* ------------------------------------------------------------------------ */

static unsigned get_u16(const uint8_t *p)
{
    return (unsigned)p[0] | ((unsigned)p[1] << 8);
}

static uint32_t get_u32(const uint8_t *p)
{
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static float get_f32(const uint8_t *p)
{
    uint32_t u = get_u32(p);
    float f;
    memcpy(&f, &u, sizeof f);
    return f;
}

static double get_f64(const uint8_t *p)
{
    uint64_t u = (uint64_t)get_u32(p) | ((uint64_t)get_u32(p + 4) << 32);
    double d;
    memcpy(&d, &u, sizeof d);
    return d;
}

/* True when a decoded float is the single-precision value of 'want' (NaN = NaN). */
static int same_f32(float got, double want)
{
    if (isnan(want))
        return isnan(got) != 0;
    return got == (float)want;
}

/* Detections with different values in every field; every 5th has no elevation. */
static void make_detections(Detection *d, int n)
{
    for (int k = 0; k < n; k++) {
        memset(&d[k], 0, sizeof d[k]);
        d[k].time_s = 100.0 + 0.001 * k;
        d[k].az_deg = 0.7 * k;
        d[k].range_m = 150.0 + 15.25 * k;
        d[k].el_deg = (k % 5 == 0) ? NAN : -2.0 + 0.125 * k;
        d[k].vr_mps = -30.0 + 0.3 * k;
        d[k].snr_dB = 10.0 + 0.1 * k;
        d[k].beam = k % 48;
        d[k].rIdx = (7 * k) % 3000;
        d[k].dIdx = k % 64;
        d[k].pulse = k % 2;
        d[k].hits = k % 9;
        d[k].pulseNumber = 5000u + (uint32_t)k;
    }
}

/* Checks the 24-byte record of one detection. */
static void check_record(const uint8_t *r, const Detection *d)
{
    CHECK(get_u16(r + 0) == (unsigned)d->beam);
    CHECK(get_u16(r + 2) == (unsigned)d->rIdx);
    CHECK(get_u16(r + 4) == (unsigned)d->dIdx);
    CHECK(get_u16(r + 6) == 0);
    CHECK(same_f32(get_f32(r + 8), d->range_m));
    CHECK(same_f32(get_f32(r + 12), d->vr_mps));
    CHECK(same_f32(get_f32(r + 16), d->el_deg));
    CHECK(same_f32(get_f32(r + 20), d->snr_dB));
}

/* True when the last two bytes of the frame hold the CRC of bytes 2 .. len-3. */
static int crc_ok(const uint8_t *frame, size_t len)
{
    return get_u16(frame + len - 2) == crc16_ccitt(frame + 2, len - 4);
}

/* ------------------------------------------------------------------------ */
/* 1. CRC                                                                    */
/* ------------------------------------------------------------------------ */

static void test_crc(void)
{
    const char *text = "123456789";
    printf("crc16_ccitt\n");
    CHECK(crc16_ccitt((const uint8_t *)text, strlen(text)) == 0x29B1);
    CHECK(crc16_ccitt((const uint8_t *)text, 0) == 0xFFFF);
}

/* ------------------------------------------------------------------------ */
/* 2. Encoders                                                               */
/* ------------------------------------------------------------------------ */

static void test_encode_detections(void)
{
    uint8_t buf[256];
    const int n = 6;
    size_t len;

    printf("serial_encode_detections\n");
    len = serial_encode_detections(buf, sizeof buf, PULSE_LONG, 0x12345678u, 271.25, g_dets, n);
    CHECK(len == 16 + 24 * (size_t)n);
    CHECK(buf[0] == 0xA5 && buf[1] == 0x5A);
    CHECK(buf[2] == 1);
    CHECK(buf[3] == PULSE_LONG);
    CHECK(buf[4] == 0x78 && buf[5] == 0x56 && buf[6] == 0x34 && buf[7] == 0x12);  /* little-endian */
    CHECK(get_u32(buf + 4) == 0x12345678u);
    CHECK(get_f32(buf + 8) == 271.25f);
    CHECK(get_u16(buf + 12) == (unsigned)n);
    for (int k = 0; k < n; k++)
        check_record(buf + 14 + 24 * k, &g_dets[k]);
    CHECK(isnan(get_f32(buf + 14 + 16)));   /* detection 0 has no elevation */
    CHECK(crc_ok(buf, len));

    /* a changed byte must change the CRC */
    buf[20] ^= 0x01;
    CHECK(!crc_ok(buf, len));

    /* buffer too small, empty frame, bad arguments */
    CHECK(serial_encode_detections(buf, len - 1, PULSE_LONG, 1, 0.0, g_dets, n) == 0);
    len = serial_encode_detections(buf, sizeof buf, PULSE_SHORT, 7, -1.5, NULL, 0);
    CHECK(len == 16);
    CHECK(buf[3] == PULSE_SHORT && get_u32(buf + 4) == 7 && get_u16(buf + 12) == 0);
    CHECK(get_f32(buf + 8) == -1.5f);
    CHECK(crc_ok(buf, len));
    CHECK(serial_encode_detections(buf, sizeof buf, PULSE_SHORT, 7, 0.0, g_dets, -1) == 0);
}

static void test_encode_plot(void)
{
    uint8_t buf[64];
    Detection p;
    size_t len;

    printf("serial_encode_plot\n");
    memset(&p, 0, sizeof p);
    p.time_s = 12345.678901234;
    p.az_deg = 123.456;
    p.range_m = 45678.9;
    p.el_deg = 3.25;
    p.vr_mps = -17.5;
    p.snr_dB = 21.75;
    p.beam = 31;
    p.hits = 4;
    p.pulse = PULSE_SHORT;
    p.pulseNumber = 987654321u;

    len = serial_encode_plot(buf, sizeof buf, &p);
    CHECK(len == 42);
    CHECK(buf[0] == 0xA5 && buf[1] == 0x5A);
    CHECK(buf[2] == 2);
    CHECK(buf[3] == PULSE_SHORT);
    CHECK(get_u32(buf + 4) == 987654321u);
    CHECK(get_f64(buf + 8) == p.time_s);
    CHECK(same_f32(get_f32(buf + 16), p.az_deg));
    CHECK(same_f32(get_f32(buf + 20), p.range_m));
    CHECK(same_f32(get_f32(buf + 24), p.el_deg));
    CHECK(same_f32(get_f32(buf + 28), p.vr_mps));
    CHECK(same_f32(get_f32(buf + 32), p.snr_dB));
    CHECK(get_u16(buf + 36) == 31);
    CHECK(get_u16(buf + 38) == 4);
    CHECK(crc_ok(buf, len));

    p.el_deg = NAN;
    p.pulse = PULSE_LONG;
    len = serial_encode_plot(buf, sizeof buf, &p);
    CHECK(len == 42 && buf[3] == PULSE_LONG && isnan(get_f32(buf + 24)) && crc_ok(buf, len));
    CHECK(serial_encode_plot(buf, 41, &p) == 0);
}

/* ------------------------------------------------------------------------ */
/* 3. Pseudo-terminal (Linux)                                                */
/* ------------------------------------------------------------------------ */

#ifdef __linux__

static double now_s(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Reads up to len bytes from fd within timeout_ms. While no data arrives,
 * 'pump' (may be NULL) gets an empty detections call, which writes what is
 * still in its ring. Returns the bytes read. */
static size_t read_bytes(int fd, uint8_t *buf, size_t len, int timeout_ms, SerialOut *pump)
{
    size_t got = 0;
    double end = now_s() + 1e-3 * timeout_ms;

    while (got < len && now_s() < end) {
        struct pollfd pfd = { fd, POLLIN, 0 };
        if (poll(&pfd, 1, 10) > 0) {
            ssize_t n = read(fd, buf + got, len - got);
            if (n > 0)
                got += (size_t)n;
        } else if (pump) {
            serial_send_detections(pump, PULSE_SHORT, 0, 0.0, NULL, 0);
        }
    }
    return got;
}

/* True when nothing arrives on fd within ms. */
static int no_data(int fd, int ms)
{
    struct pollfd pfd = { fd, POLLIN, 0 };
    return poll(&pfd, 1, ms) == 0;
}

/* Discards everything that is readable on fd. */
static void drain(int fd)
{
    uint8_t tmp[4096];
    while (!no_data(fd, 20) && read(fd, tmp, sizeof tmp) > 0) {
    }
}

/* Small detections frame, plot frame, empty update and the 200-per-frame split. */
static void test_pty_frames(int master, SerialOut *s)
{
    static uint8_t expect[16384], got[16384];
    const int counts[3] = { 200, 200, 50 };
    uint64_t sent, dropped, bytes;
    size_t len1, len2, total, off = 0;
    Detection plot = g_dets[3];

    printf("serial port: frames\n");
    len1 = serial_encode_detections(expect, sizeof expect, PULSE_SHORT, 77, 12.5, g_dets, 5);
    CHECK(serial_send_detections(s, PULSE_SHORT, 77, 12.5, g_dets, 5) == 1);
    CHECK(read_bytes(master, got, len1, 1000, s) == len1);
    CHECK(memcmp(got, expect, len1) == 0);

    plot.hits = 6;
    len2 = serial_encode_plot(expect, sizeof expect, &plot);
    CHECK(serial_send_plot(s, &plot) == 1);
    CHECK(read_bytes(master, got, len2, 1000, s) == len2);
    CHECK(memcmp(got, expect, len2) == 0);

    /* n == 0 sends nothing */
    CHECK(serial_send_detections(s, PULSE_LONG, 78, 13.0, g_dets, 0) == 1);
    CHECK(no_data(master, 50));

    /* 450 detections: frames of 200, 200 and 50 with the same header */
    total = 3 * 16 + 24 * (size_t)N_DETS;
    CHECK(serial_send_detections(s, PULSE_LONG, 99, 45.0, g_dets, N_DETS) == 1);
    CHECK(read_bytes(master, got, total, 2000, s) == total);
    for (int f = 0; f < 3; f++) {
        const uint8_t *fr = got + off;
        size_t flen = 16 + 24 * (size_t)counts[f];
        size_t elen = serial_encode_detections(expect, sizeof expect, PULSE_LONG, 99, 45.0,
                                               g_dets + MAX_PER_FRAME * f, counts[f]);
        CHECK(fr[0] == 0xA5 && fr[1] == 0x5A && fr[2] == 1 && fr[3] == PULSE_LONG);
        CHECK(get_u32(fr + 4) == 99 && get_f32(fr + 8) == 45.0f);
        CHECK(get_u16(fr + 12) == (unsigned)counts[f]);
        check_record(fr + 14, &g_dets[MAX_PER_FRAME * f]);
        CHECK(crc_ok(fr, flen));
        CHECK(elen == flen && memcmp(fr, expect, flen) == 0);
        off += flen;
    }
    CHECK(no_data(master, 50));

    serial_stats(s, &sent, &dropped, &bytes);
    CHECK(sent == 5 && dropped == 0);
    CHECK(bytes == len1 + len2 + total);
}

/* Runs in a child process: reads 'nFrames' full detections frames from the
 * master and checks them. Returns 0 when all are complete and valid. */
static int reader_child(int master, int nFrames)
{
    size_t len = (size_t)nFrames * FULL_FRAME;
    uint8_t *buf = (uint8_t *)malloc(len + 1);
    uint32_t last = 0;

    if (!buf)
        return 10;
    if (read_bytes(master, buf, len, 5000, NULL) != len)
        return 11;
    for (int f = 0; f < nFrames; f++) {
        const uint8_t *fr = buf + (size_t)f * FULL_FRAME;
        if (fr[0] != 0xA5 || fr[1] != 0x5A || fr[2] != 1 || get_u16(fr + 12) != MAX_PER_FRAME)
            return 12;
        if (!crc_ok(fr, FULL_FRAME))
            return 13;
        if (f > 0 && get_u32(fr + 4) <= last)
            return 14;
        last = get_u32(fr + 4);
    }
    if (!no_data(master, 100))
        return 15;
    free(buf);
    return 0;
}

/* Fills the ring without reading the port: frames are dropped without
 * blocking, then serial_close flushes the queued frames to a reader. */
static void test_pty_overflow(int master, SerialOut *s)
{
    uint64_t sent0, dropped0, bytes0, sent, dropped, bytes;
    int accepted = 0, status = -1;
    double t0, dt;
    pid_t pid;

    printf("serial port: full ring, close\n");
    serial_stats(s, &sent0, &dropped0, &bytes0);
    t0 = now_s();
    for (int i = 0; i < 100; i++)
        accepted += serial_send_detections(s, PULSE_SHORT, 1000u + (uint32_t)i, 1.0, g_dets,
                                           MAX_PER_FRAME);
    dt = now_s() - t0;
    serial_stats(s, &sent, &dropped, &bytes);
    printf("  100 frames of %d bytes in %.1f ms: %d queued, %d dropped\n", FULL_FRAME, 1e3 * dt,
           accepted, 100 - accepted);
    CHECK(dt < 0.5);                                   /* never blocks */
    CHECK(accepted > 50 && accepted < 100);            /* the 256 kB ring was full */
    CHECK(sent - sent0 == (uint64_t)accepted);
    CHECK(dropped - dropped0 == (uint64_t)(100 - accepted));
    CHECK(bytes - bytes0 < (uint64_t)accepted * FULL_FRAME);   /* part is still in the ring */

    /* a reader in a child process takes the data while serial_close flushes */
    fflush(stdout);
    pid = fork();
    if (pid == 0)
        _exit(reader_child(master, accepted));
    CHECK(pid > 0);
    t0 = now_s();
    serial_close(s);
    dt = now_s() - t0;
    if (pid > 0)
        waitpid(pid, &status, 0);
    printf("  serial_close took %.1f ms, reader status %d\n", 1e3 * dt,
           WIFEXITED(status) ? WEXITSTATUS(status) : -1);
    CHECK(WIFEXITED(status) && WEXITSTATUS(status) == 0);
    CHECK(dt < 2.5);
}

/* Nobody reads: serial_close gives up after about 2 s. */
static void test_pty_close_timeout(int master, const char *name)
{
    char err[256];
    SerialOut *s = serial_open(name, 921600, err, sizeof err);
    double t0, dt;
    int i;

    printf("serial port: close timeout\n");
    CHECK(s != NULL);
    if (!s)
        return;
    for (i = 0; i < 200; i++)
        if (!serial_send_detections(s, PULSE_SHORT, (uint32_t)i, 0.0, g_dets, MAX_PER_FRAME))
            break;
    CHECK(i < 200);
    t0 = now_s();
    serial_close(s);
    dt = now_s() - t0;
    printf("  serial_close took %.2f s\n", dt);
    CHECK(dt > 1.5 && dt < 3.5);
    drain(master);
}

static void test_pty(void)
{
    char name[256], err[256];
    struct termios tio;
    SerialOut *s;
    int keep, master = posix_openpt(O_RDWR | O_NOCTTY);

    printf("serial port: open\n");
    if (master < 0 || grantpt(master) != 0 || unlockpt(master) != 0 || !ptsname(master)) {
        printf("  cannot create a pseudo-terminal\n");
        g_failures++;
        return;
    }
    snprintf(name, sizeof name, "%s", ptsname(master));
    keep = open(name, O_RDWR | O_NOCTTY);   /* keeps the slave side open between the tests */
    CHECK(keep >= 0);

    err[0] = 0;
    CHECK(serial_open("/dev/no_such_serial_port", 921600, err, sizeof err) == NULL);
    CHECK(err[0] != 0);
    CHECK(serial_open(name, 12345, err, sizeof err) == NULL);
    CHECK(serial_open("/dev/null", 921600, err, sizeof err) == NULL);   /* not a tty */

    s = serial_open(name, 921600, err, sizeof err);
    CHECK(s != NULL);
    if (!s) {
        printf("  serial_open(%s): %s\n", name, err);
        close(keep);
        close(master);
        return;
    }
    CHECK(tcgetattr(keep, &tio) == 0);
    CHECK(cfgetospeed(&tio) == B921600 && cfgetispeed(&tio) == B921600);
    CHECK((tio.c_cflag & CSIZE) == CS8 && !(tio.c_cflag & PARENB) && !(tio.c_cflag & CSTOPB));
    CHECK(!(tio.c_oflag & OPOST) && !(tio.c_lflag & ICANON) && !(tio.c_iflag & IXON));

    test_pty_frames(master, s);
    test_pty_overflow(master, s);   /* closes s */
    test_pty_close_timeout(master, name);

    serial_stats(NULL, NULL, NULL, NULL);   /* must not crash */
    serial_close(NULL);
    close(keep);
    close(master);
}

#endif /* __linux__ */

int main(void)
{
    make_detections(g_dets, N_DETS);
    test_crc();
    test_encode_detections();
    test_encode_plot();
#ifdef __linux__
    test_pty();
#else
    printf("serial port test skipped (Linux only)\n");
#endif
    if (g_failures) {
        printf("FAIL (%d checks failed)\n", g_failures);
        return 1;
    }
    printf("PASS\n");
    return 0;
}

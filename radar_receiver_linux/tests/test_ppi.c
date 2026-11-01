/*
 * test_ppi.c - test of the Qt PPI display (ppi_qt.cpp).
 *
 * A worker posts antenna positions over 1.2 scans (1000 updates per second
 * of display time, faster than real time) and a few plots of both pulse
 * types, one of them without elevation (NaN), with the CFAR detections that
 * led to them, and reports a CFAR threshold, then returns. ppi_run renders
 * the PPI to a PNG file; the test checks that the file is a PNG and not tiny.
 * The display shows the plots (default mode), the detections are stored.
 *
 * Expected picture (antenna start azimuth 30 deg, clockwise, ends at 102 deg
 * in scan 2): the plot at 45 deg of scan 1 has been swept over again and is
 * gone; all other plots are shown; the antenna line points to 102 deg.
 *
 * Usage: test_ppi [snapshot.png]   (uses QT_QPA_PLATFORM=offscreen if unset)
 */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#define _POSIX_C_SOURCE 200112L
#endif

#include "radar_rx.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <time.h>
#endif

#define MAX_RANGE_M   60e3
#define AZ0_DEG       30.0      /* antenna azimuth at t = 0 (scan boundary)   */
#define RATE_DPS      36.0      /* display-time rotation rate (10 s per scan) */
#define STEP_DEG      0.036     /* 1000 antenna updates per second            */
#define N_SCANS       1.2
#define CPI_DELAY_DEG 1.0       /* plot posted when the antenna is this far past it */

typedef struct {
    int    scan;      /* scan in which the target is seen (1-based) */
    double az_deg, range_m, el_deg, vr_mps, snr_dB;
    int    pulse;
} TestPlot;

static const TestPlot kPlots[] = {
    { 1,  45.0,  8.0e3,  2.0,  -40.0, 18.0, PULSE_SHORT },   /* swept again in scan 2: removed */
    { 1, 120.0, 15.0e3,  1.5,   25.0, 21.0, PULSE_SHORT },
    { 1, 122.0, 15.4e3,  1.5,   25.0, 19.0, PULSE_SHORT },
    { 1, 200.0, 35.0e3,  3.0,  110.0, 15.0, PULSE_LONG  },
    { 1, 250.0, 42.0e3,  NAN, -150.0, 13.0, PULSE_LONG  },   /* no elevation */
    { 1, 300.0, 52.0e3,  1.0,   80.0, 14.0, PULSE_LONG  },
    { 1, 330.0, 25.0e3,  0.5,  -60.0, 20.0, PULSE_SHORT },
    { 2,  60.0, 45.0e3,  2.0,  200.0, 16.0, PULSE_LONG  },   /* scan 2 */
    { 2,  80.0,  5.0e3,  6.0,   10.0, 25.0, PULSE_SHORT },   /* scan 2 */
};
#define N_PLOTS ((int)(sizeof kPlots / sizeof kPlots[0]))

static volatile sig_atomic_t g_stop = 0;

static void sleep_ms(int ms)
{
#ifdef _WIN32
    Sleep((DWORD)ms);
#else
    struct timespec ts;
    ts.tv_sec = ms / 1000;
    ts.tv_nsec = (long)(ms % 1000) * 1000000L;
    nanosleep(&ts, NULL);
#endif
}

static void set_offscreen_platform(void)
{
    if (getenv("QT_QPA_PLATFORM")) return;
#ifdef _WIN32
    _putenv_s("QT_QPA_PLATFORM", "offscreen");
#else
    setenv("QT_QPA_PLATFORM", "offscreen", 0);
#endif
}

/* Unwrapped antenna azimuth at which the antenna points at a test plot. */
static double plot_unwrapped_az(const TestPlot *t)
{
    double rel = fmod(t->az_deg - AZ0_DEG + 360.0, 360.0);
    return AZ0_DEG + 360.0 * (t->scan - 1) + rel;
}

static void post_plot(const TestPlot *t, double time_s, uint32_t pulseNumber)
{
    Detection d;
    memset(&d, 0, sizeof d);
    d.time_s = time_s;
    d.az_deg = t->az_deg;
    d.range_m = t->range_m;
    d.el_deg = t->el_deg;
    d.vr_mps = t->vr_mps;
    d.snr_dB = t->snr_dB;
    d.pulse = t->pulse;
    d.hits = 3;
    d.pulseNumber = pulseNumber;
    ppi_post_plot(&d);
}

static int worker(void *arg)
{
    int *posted = (int *)arg;
    int done[N_PLOTS] = { 0 };
    const int nUpdates = (int)(N_SCANS * 360.0 / STEP_DEG);

    for (int k = 0; k <= nUpdates && !g_stop; ++k) {
        const double azU = AZ0_DEG + k * STEP_DEG;
        const double t = (azU - AZ0_DEG) / RATE_DPS;
        const uint32_t scan = (uint32_t)floor((azU - AZ0_DEG) / 360.0) + 1;
        ppi_post_antenna(fmod(azU, 360.0), t, scan);
        if (k == 0) ppi_post_cfar_threshold(12.8, 1e-8, 12.8);
        double request_dB;
        if (ppi_take_cfar_request(&request_dB)) ++*posted;   /* nobody pressed Apply: must not happen */
        for (int i = 0; i < N_PLOTS; ++i) {                  /* detections while the antenna passes */
            const double d = azU - plot_unwrapped_az(&kPlots[i]);
            if (d >= -2.0 && d < 2.0) {
                Detection det;
                memset(&det, 0, sizeof det);
                det.time_s = t;
                det.az_deg = fmod(azU, 360.0);
                det.range_m = kPlots[i].range_m;
                det.el_deg = kPlots[i].el_deg;
                det.snr_dB = kPlots[i].snr_dB;
                det.pulse = kPlots[i].pulse;
                ppi_post_detections(&det, 1, STEP_DEG);
            }
        }

        for (int i = 0; i < N_PLOTS; ++i) {
            if (!done[i] && azU >= plot_unwrapped_az(&kPlots[i]) + CPI_DELAY_DEG) {
                post_plot(&kPlots[i], t, (uint32_t)k);
                done[i] = 1;
                ++*posted;
            }
        }
        if (k % 20 == 0) sleep_ms(1);   /* about 1 s for the whole run */
    }
    return 0;
}

/* Returns the file size, or -1 when the file does not exist. Checks the PNG signature. */
static long png_file_size(const char *path, int *isPng)
{
    static const unsigned char sig[8] = { 0x89, 'P', 'N', 'G', '\r', '\n', 0x1A, '\n' };
    unsigned char head[8];
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    *isPng = fread(head, 1, sizeof head, f) == sizeof head && memcmp(head, sig, sizeof sig) == 0;
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fclose(f);
    return size;
}

int main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "ppi_test.png";
    int failures = 0, posted = 0, isPng = 0;

    set_offscreen_platform();
    remove(path);

    int code = ppi_run(argc, argv, MAX_RANGE_M, worker, &posted, &g_stop, path);
    long size = png_file_size(path, &isPng);

    printf("ppi_run returned %d, plots posted %d, snapshot %s: %ld bytes\n", code, posted, path, size);
    if (code != 0)          { printf("FAIL: ppi_run returned %d\n", code); failures++; }
    if (posted != N_PLOTS)  { printf("FAIL: %d of %d plots posted\n", posted, N_PLOTS); failures++; }
    if (size < 0)           { printf("FAIL: snapshot %s not written\n", path); failures++; }
    else if (!isPng)        { printf("FAIL: snapshot is not a PNG file\n"); failures++; }
    else if (size < 10000)  { printf("FAIL: snapshot is too small (%ld bytes)\n", size); failures++; }

    /* posting without a running display must be harmless */
    ppi_post_antenna(10.0, 0.0, 1);
    post_plot(&kPlots[0], 0.0, 0);
    ppi_post_detections(NULL, 0, 0.1);
    ppi_post_cfar_threshold(13.0, NAN, 13.0);

    printf(failures ? "FAIL\n" : "PASS\n");
    return failures ? 1 : 0;
}

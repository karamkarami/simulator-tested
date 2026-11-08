/*
 * bench_processor.c - processing time of the receiver per block on synthetic
 * PRIs, without the network, the display or the console output:
 *
 *   radar_bench [--config FILE] [--threads N] [--pris N] [--thresholds DB,DB,...]
 *               [--clutter DB] [--mf auto|direct|fft] [--update-every N]
 *
 * Every PRI holds new complex Gaussian noise (16 LSB rms per component, like
 * the DDC output of the simulator), a few moving targets in both sub-pulses
 * that cross the beam while the antenna turns, and optionally stationary
 * clutter at the start of both receive windows (--clutter: its power per
 * sample above the noise [dB], 0 = none). For every CFAR threshold (default:
 * the configuration's, then lower ones down to a few dB, where most cells
 * are false alarms) the mean time per PRI of every block is printed, as the
 * timing panel of the PPI shows it, with the CFAR detections per update and
 * the plots. The first threshold starts after a full CPI. Use the
 * configuration of the radar to be measured (export_receiver_config.m), e.g.
 * one element: config/radar_config_1ch.bin.
 */
#include "radar_rx.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

static uint64_t g_rng = 0x9E3779B97F4A7C15ull;

static double uniform(void)
{
    g_rng ^= g_rng << 13;
    g_rng ^= g_rng >> 7;
    g_rng ^= g_rng << 17;
    return ((g_rng >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

static void gauss2(double *a, double *b)          /* two independent N(0, 1) */
{
    const double r = sqrt(-2.0 * log(uniform())), t = 2.0 * PI * uniform();
    *a = r * cos(t);
    *b = r * sin(t);
}

static int16_t sat16(double v)
{
    v = floor(v + 0.5);
    return (int16_t)(v > 32767 ? 32767 : v < -32768 ? -32768 : v);
}

typedef struct { uint64_t dets, updates, plots; } Counts;

static void on_detections(void *user, int pulse, const Detection *d, int n)
{
    Counts *c = user;
    (void)pulse;
    (void)d;
    c->dets += (uint64_t)n;
    c->updates++;
}

static void on_plot(void *user, const Detection *plot)
{
    (void)plot;
    ((Counts *)user)->plots++;
}

typedef struct { int pulse; double cell, vr_mps, amp, az_deg; } Target;

int main(int argc, char **argv)
{
    const char *configPath = "config/radar_config.bin", *mf = "auto";
    int threads = 0, nPris = 1000, updateEvery = 0, nThr = 0;
    double thr[16], clutter_dB = 0.0;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i], *v = i + 1 < argc ? argv[i + 1] : NULL;
        if (!v && strcmp(a, "--help") && strcmp(a, "-h")) { fprintf(stderr, "%s needs a value\n", a); return 2; }
        if (!strcmp(a, "--config")) configPath = argv[++i];
        else if (!strcmp(a, "--threads")) threads = atoi(argv[++i]);
        else if (!strcmp(a, "--pris")) nPris = atoi(argv[++i]);
        else if (!strcmp(a, "--clutter")) clutter_dB = atof(argv[++i]);
        else if (!strcmp(a, "--mf")) mf = argv[++i];
        else if (!strcmp(a, "--update-every")) updateEvery = atoi(argv[++i]);
        else if (!strcmp(a, "--thresholds")) {
            char *s = argv[++i];
            for (char *tok = strtok(s, ","); tok && nThr < 16; tok = strtok(NULL, ",")) thr[nThr++] = atof(tok);
        } else {
            puts("usage: radar_bench [--config FILE] [--threads N] [--pris N] [--thresholds DB,DB,...]\n"
                 "                   [--clutter DB] [--mf auto|direct|fft] [--update-every N]");
            return strcmp(a, "--help") && strcmp(a, "-h") ? 2 : 0;
        }
    }
    RadarConfig cfg;
    char err[256];
    if (radar_config_load(configPath, &cfg, err, sizeof err) != 0) {
        fprintf(stderr, "configuration %s: %s\n", configPath, err);
        return 1;
    }
    cfg.mfMethod = !strcmp(mf, "direct") ? MF_DIRECT : !strcmp(mf, "fft") ? MF_FFT : MF_AUTO;
    if (updateEvery > 0) cfg.updateEvery = updateEvery;
    if (nThr == 0) {                                   /* the configuration's, then lower ones */
        const double t0 = 10.0 * log10(cfg.alpha), lower[] = { 10.0, 8.0, 6.0, 4.0 };
        thr[nThr++] = t0;
        for (int k = 0; k < 4; k++)
            if (lower[k] < t0 - 0.5) thr[nThr++] = lower[k];
    }
    Processor *proc = processor_create(&cfg, threads);
    if (!proc) { fprintf(stderr, "cannot create the processor\n"); return 1; }

    /* targets: in both sub-pulses, crossing the beam one after another */
    const double pri_s = cfg.meanPri_s > 0 ? cfg.meanPri_s : 1e-3, dwell_deg = cfg.azBeamwidth_deg;
    const double dv = fabs(cfg.dopFreq_Hz[1] - cfg.dopFreq_Hz[0]) * cfg.pulse[PULSE_LONG].lambda_m / 2;
    Target tg[] = {
        { PULSE_SHORT, 0.45, -(cfg.nfft / 5 + 0.3) * dv, 40.0, 0.5 * dwell_deg },
        { PULSE_LONG,  0.25, (cfg.nfft / 7 + 0.2) * dv, 10.0, 1.0 * dwell_deg },
        { PULSE_LONG,  0.60, -(cfg.nfft / 3 + 0.4) * dv, 3.0, 1.5 * dwell_deg },
        { PULSE_LONG,  0.85, (cfg.nfft / 4 + 0.1) * dv, 300.0, 2.0 * dwell_deg },
    };
    const int nT = sizeof tg / sizeof tg[0];
    const double sigma = 16.0;                          /* LSB rms per component */
    float *clutRe = NULL, *clutIm = NULL;
    if (clutter_dB > 0.0) {                             /* the first 10 % of each receive window */
        clutRe = calloc((size_t)cfg.nCh * cfg.nDdc, sizeof(float));
        clutIm = calloc((size_t)cfg.nCh * cfg.nDdc, sizeof(float));
        for (int i = 0; i < N_PULSES; i++) {
            const PulseConfig *pc = &cfg.pulse[i];
            const int n = pc->rxCount / 10;
            for (int k = 0; k < n; k++) {
                const double a = sigma * pow(10.0, clutter_dB / 20.0) * pow(1.0 + 10.0 * k / n, -1.5);
                for (int ch = 0; ch < cfg.nCh; ch++) {
                    double g1, g2;
                    gauss2(&g1, &g2);
                    clutRe[(size_t)ch * cfg.nDdc + pc->rxFirst + k] = (float)(a * g1);
                    clutIm[(size_t)ch * cfg.nDdc + pc->rxFirst + k] = (float)(a * g2);
                }
            }
        }
    }

    PriData pri;
    memset(&pri, 0, sizeof pri);
    pri.nCh = cfg.nCh;
    pri.nSamples = cfg.nDdc;
    pri.ddcFs_Hz = (float)cfg.fsOut_Hz;
    pri.shortIdx[0] = (uint32_t)cfg.pulse[PULSE_SHORT].rxFirst + 1;
    pri.shortIdx[1] = (uint32_t)(cfg.pulse[PULSE_SHORT].rxFirst + cfg.pulse[PULSE_SHORT].rxCount);
    pri.longIdx[0] = (uint32_t)cfg.pulse[PULSE_LONG].rxFirst + 1;
    pri.longIdx[1] = (uint32_t)(cfg.pulse[PULSE_LONG].rxFirst + cfg.pulse[PULSE_LONG].rxCount);
    double *re = malloc(sizeof(double) * cfg.nDdc), *im = malloc(sizeof(double) * cfg.nDdc);
    pri.iq = malloc(sizeof(int16_t) * 2 * (size_t)cfg.nCh * cfg.nDdc);
    if (!re || !im || !pri.iq) { fprintf(stderr, "out of memory\n"); return 1; }

    printf("%s: %d element(s), %d beam(s), short %d / long %d range cells, CPI %d, FFT %d, CFAR %s %s, "
           "update every %d pulse(s), %d PRIs per threshold, clutter %s\n", configPath, cfg.nCh, cfg.nBeams,
           cfg.pulse[0].nRange, cfg.pulse[1].nRange, cfg.cpiLength, cfg.nfft, cfar_type_name(cfg.cfarType),
           cfar_method_name(cfg.cfarMethod),
           cfg.updateEvery, nPris, clutter_dB > 0 ? "on" : "off");
    printf("mean time per PRI [us] (CFAR+pr: the CFAR with the pruning and the reports, as in the PPI)\n"
           "%9s %9s | %6s %6s %6s %7s %7s %6s %7s %7s | %7s | %8s %7s\n", "thr[dB]", "Pfa", "DBF",
           "MF", "MTI", "Doppler", "CFAR+pr", "interp", "extract", "track", "total", "det/upd", "plots");
    Counts cnt;
    memset(&cnt, 0, sizeof cnt);
    ProcessorCallbacks cb = { on_detections, on_plot, NULL, &cnt, NULL, NULL, NULL };
    uint32_t pulseNumber = 1;
    for (int t = 0; t < nThr; t++) {
        processor_set_cfar_alpha(proc, pow(10.0, thr[t] / 10.0));
        const int warm = t == 0 ? cfg.cpiLength : 0;
        ProcessorTiming t0;
        Counts c0 = cnt;
        processor_timing(proc, &t0);
        for (int k = 0; k < warm + nPris; k++, pulseNumber++) {
            if (k == warm) { processor_timing(proc, &t0); c0 = cnt; }
            const double time_s = pulseNumber * pri_s, az = fmod(cfg.rate_dps * time_s, 360.0);
            for (int ch = 0; ch < cfg.nCh; ch++) {
                for (int n = 0; n < cfg.nDdc; n++) {
                    gauss2(&re[n], &im[n]);
                    re[n] *= sigma;
                    im[n] *= sigma;
                }
                if (clutRe)
                    for (int n = 0; n < cfg.nDdc; n++) {
                        re[n] += clutRe[(size_t)ch * cfg.nDdc + n];
                        im[n] += clutIm[(size_t)ch * cfg.nDdc + n];
                    }
                for (int j = 0; j < nT; j++) {                 /* echoes, two-way Gaussian beam */
                    const PulseConfig *pc = &cfg.pulse[tg[j].pulse];
                    const double cycle = 4.0 * dwell_deg;          /* a target every 4 beamwidths */
                    const double daz = fmod(az, cycle) - tg[j].az_deg;
                    const double g = exp(-2.77 * daz * daz / (dwell_deg * dwell_deg));
                    if (g < 1e-3 || !pc->nPulse) continue;
                    const double ph = -4.0 * PI * tg[j].vr_mps * time_s / pc->lambda_m + 0.4 * ch;
                    const double ar = sigma * tg[j].amp * g * cos(ph), ai = sigma * tg[j].amp * g * sin(ph);
                    const int first = pc->ddcFirst + (int)(tg[j].cell * pc->nRange) + pc->mfShift;
                    for (int m = 0; m < pc->nPulse; m++) {
                        const int n = first + m;
                        if (n < 0 || n >= cfg.nDdc) continue;
                        re[n] += ar * pc->pulseRe[m] - ai * pc->pulseIm[m];
                        im[n] += ar * pc->pulseIm[m] + ai * pc->pulseRe[m];
                    }
                }
                int16_t *dst = pri.iq + 2 * (size_t)ch * cfg.nDdc;
                for (int n = 0; n < cfg.nDdc; n++) {
                    dst[2 * n] = sat16(re[n]);
                    dst[2 * n + 1] = sat16(im[n]);
                }
            }
            pri.pulseNumber = pulseNumber;
            pri.scanNumber = 1;
            pri.az_deg = (float)az;
            pri.staggerIndex = (uint16_t)(1 + (pulseNumber - 1) % (uint32_t)cfg.nStagger);
            pri.time_s = time_s;
            processor_process(proc, &pri, &cb);
        }
        ProcessorTiming t1;
        processor_timing(proc, &t1);
        const double n = (double)(t1.pris - t0.pris), us = 1e3 / (n > 0 ? n : 1);
        const double upd = (double)(cnt.updates - c0.updates);
        printf("%9.2f %9.1e | %6.1f %6.1f %6.1f %7.1f %7.1f %6.1f %7.1f %7.1f | %7.1f | %8.0f %7llu\n", thr[t],
               cfar_pfa(&cfg, pow(10.0, thr[t] / 10.0)), (t1.dbf_ms - t0.dbf_ms) * us,
               (t1.mf_ms - t0.mf_ms) * us, (t1.mti_ms - t0.mti_ms) * us, (t1.doppler_ms - t0.doppler_ms) * us,
               (t1.cfar_ms - t0.cfar_ms) * us, (t1.interp_ms - t0.interp_ms) * us, (t1.extract_ms - t0.extract_ms) * us, (t1.track_ms - t0.track_ms) * us,
               (t1.total_ms - t0.total_ms) * us, upd > 0 ? (double)(cnt.dets - c0.dets) / upd : 0.0,
               (unsigned long long)(cnt.plots - c0.plots));
        fflush(stdout);
    }
    processor_destroy(proc);
    free(re);
    free(im);
    free(pri.iq);
    free(clutRe);
    free(clutIm);
    radar_config_free(&cfg);
    return 0;
}

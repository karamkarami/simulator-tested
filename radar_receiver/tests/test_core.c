/*
 * test_core.c - unit tests of the processing core:
 *   1) batched FFT against a direct DFT,
 *   2) plot extractor (azimuth grouping, minimum hits, flush),
 *   3) the exported configuration (argv[1]) and a synthetic moving target
 *      processed end to end (DBF, matched filter, MTI, Doppler, CFAR,
 *      pruning, extraction): exactly one plot at the right range, radial
 *      velocity and elevation beam.
 */
#include "radar_rx.h"
#include "extractor.h"
#include "fft.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

static int failures = 0;
#define CHECK(cond, ...) do { if (!(cond)) { failures++; printf("FAIL %s:%d: ", __FILE__, __LINE__); \
                                              printf(__VA_ARGS__); printf("\n"); } } while (0)

/* ------------------------------------------------------------------------ */

static void test_fft(void)
{
    const int sizes[] = { 2, 8, 32, 128, 512 };
    const int nSig = 5;
    for (size_t t = 0; t < sizeof sizes / sizeof sizes[0]; t++) {
        const int n = sizes[t];
        FftPlan plan;
        CHECK(fft_plan_init(&plan, n) == 0, "plan %d", n);
        double *xr = malloc(sizeof(double) * n * nSig), *xi = malloc(sizeof(double) * n * nSig);
        float *re = malloc(sizeof(float) * n * nSig), *im = malloc(sizeof(float) * n * nSig);
        srand(1234 + n);
        for (int c = 0; c < nSig; c++)
            for (int m = 0; m < n; m++) {
                xr[c * n + m] = rand() / (double)RAND_MAX - 0.5;
                xi[c * n + m] = rand() / (double)RAND_MAX - 0.5;
                re[plan.bitrev[m] * nSig + c] = (float)xr[c * n + m];
                im[plan.bitrev[m] * nSig + c] = (float)xi[c * n + m];
            }
        fft_rows(&plan, re, im, nSig, nSig);
        double maxErr = 0;
        for (int c = 0; c < nSig; c++)
            for (int k = 0; k < n; k++) {
                double sr = 0, si = 0;
                for (int m = 0; m < n; m++) {
                    const double a = -2 * PI * k * m / n;
                    sr += xr[c * n + m] * cos(a) - xi[c * n + m] * sin(a);
                    si += xr[c * n + m] * sin(a) + xi[c * n + m] * cos(a);
                }
                const double e = hypot(re[k * nSig + c] - sr, im[k * nSig + c] - si);
                if (e > maxErr) maxErr = e;
            }
        CHECK(maxErr < 1e-4 * n, "FFT %d: max error %g", n, maxErr);
        free(xr); free(xi); free(re); free(im);
        fft_plan_free(&plan);
    }
}

/* ------------------------------------------------------------------------ */

typedef struct { Detection plots[16]; int n; } PlotList;

static void collect_plot(void *user, const Detection *plot)
{
    PlotList *l = user;
    if (l->n < 16) l->plots[l->n++] = *plot;
}

static Detection det_at(int b, int r, int d, double power, double az)
{
    Detection x;
    memset(&x, 0, sizeof x);
    x.beam = b; x.rIdx = r; x.dIdx = d; x.power = power; x.az_deg = az; x.range_m = 1000.0 * r;
    return x;
}

static void test_extractor(void)
{
    Extractor e;
    PlotList l = { .n = 0 };
    extractor_init(&e, 32, 3, 1);
    /* target A: 5 updates, peak in the third, moves one range cell and wraps in Doppler */
    for (int k = 0; k < 5; k++) {
        Detection d[2];
        int n = 0;
        d[n++] = det_at(4, 100 + (k == 3), k < 3 ? 31 : 0, k == 2 ? 50.0 : 10.0, 10.0 + k);
        if (k == 1) d[n++] = det_at(10, 300, 5, 99.0, 11.0);    /* single hit: no plot */
        extractor_update(&e, d, n, 0, collect_plot, &l);
    }
    CHECK(l.n == 0, "no plot may close while the target is still detected (%d)", l.n);
    extractor_update(&e, NULL, 0, 0, collect_plot, &l);      /* first miss: still open */
    CHECK(l.n == 0, "plot closed after the first miss (maxMiss = 1)");
    extractor_update(&e, NULL, 0, 0, collect_plot, &l);      /* second miss: closes */
    CHECK(l.n == 1, "expected 1 plot, got %d", l.n);
    if (l.n == 1) {
        CHECK(l.plots[0].hits == 5, "hits %d", l.plots[0].hits);
        CHECK(fabs(l.plots[0].az_deg - 12.0) < 1e-9, "plot azimuth %g (strongest detection)", l.plots[0].az_deg);
    }
    /* flush closes an open group with enough hits */
    l.n = 0;
    for (int k = 0; k < 3; k++) {
        Detection d = det_at(0, 50, 7, 5.0, 100.0 + k);
        extractor_update(&e, &d, 1, 0, collect_plot, &l);
    }
    extractor_update(&e, NULL, 0, 1, collect_plot, &l);
    CHECK(l.n == 1, "flush: expected 1 plot, got %d", l.n);
    extractor_free(&e);
}

/* ------------------------------------------------------------------------ */

typedef struct {
    PlotList plots;
    int updates[N_PULSES];
    int dets[N_PULSES];
} Sink;

static void on_detections(void *user, int pulse, const Detection *d, int n)
{
    Sink *s = user;
    (void)d;
    s->updates[pulse]++;
    s->dets[pulse] += n;
}

static void on_plot(void *user, const Detection *plot)
{
    collect_plot(&((Sink *)user)->plots, plot);
}

static double gauss(void)
{
    double u1 = (rand() + 1.0) / (RAND_MAX + 2.0), u2 = rand() / (RAND_MAX + 1.0);
    return sqrt(-2 * log(u1)) * cos(2 * PI * u2);
}

static int16_t sat16(double v)
{
    v = floor(v + 0.5);
    return (int16_t)(v > 32767 ? 32767 : v < -32768 ? -32768 : v);
}

/* A target of the long pulse: echo = transmitted waveform (the pulse of the
 * configuration, or else from the matched filter taps; scaled to the energy
 * of the unit-norm taps), steered to one elevation beam, moving with vr. With
 * a mismatch filter the echo starts mfShift samples after the start of the
 * target cell (its range axis). Noise 16 LSB rms per component like the
 * simulator's DDC output. */
static int test_end_to_end(const char *configPath, int mfDirect, double threshold_dB, Detection *plotOut)
{
    RadarConfig cfg;
    char err[256];
    if (radar_config_load(configPath, &cfg, err, sizeof err) != 0) {
        CHECK(0, "configuration %s: %s", configPath, err);
        return 0;
    }
    cfg.mfMethod = mfDirect ? MF_DIRECT : MF_AUTO;
    CHECK(cfg.nCh >= 1 && cfg.nDdc > 0 && cfg.cpiLength > 0, "configuration sizes");
    CHECK(cfg.pulse[PULSE_SHORT].ddcFirst < cfg.pulse[PULSE_LONG].ddcFirst, "short pulse before the long one");
    printf("config: %d channels, %d beams, %d DDC samples, CPI %d, nfft %d, %s PRF, %s matched filter\n",
           cfg.nCh, cfg.nBeams, cfg.nDdc, cfg.cpiLength, cfg.nfft, cfg.staggered ? "staggered" : "fixed",
           mfDirect ? "direct" : "automatic");
    for (int i = 0; i < N_PULSES; i++)
        printf("  %s pulse: %s, %s filter, %d taps, %d range cells\n", i == PULSE_SHORT ? "short" : "long",
               cfg.pulse[i].waveform == WAVE_PHASE_CODE ? "phase code" : "LFM",
               cfg.pulse[i].compression == COMPRESSION_MISMATCH ? "mismatch" : "matched", cfg.pulse[i].nTaps,
               cfg.pulse[i].nRange);

    const PulseConfig *pc = &cfg.pulse[PULSE_LONG];
    const int beam = cfg.nBeams > 4 ? cfg.nBeams / 3 : 0;
    const int r0 = pc->nRange / 2;                      /* range cell of the target */
    const int nEcho = pc->nPulse ? pc->nPulse : pc->nTaps;
    double pulseNorm = 0.0;
    for (int m = 0; m < pc->nPulse; m++) pulseNorm += pc->pulseRe[m] * pc->pulseRe[m] + pc->pulseIm[m] * pc->pulseIm[m];
    pulseNorm = pulseNorm > 0.0 ? sqrt(pulseNorm) : 1.0;
    const double vr = 7.3;                              /* m/s, inside the unambiguous interval */
    const double pri_s = 1e-3, amp = 30.0;              /* amplitude per element (LSB)          */
    const int nPulses = cfg.cpiLength + 60, nWith = cfg.cpiLength + 30;

    PriData pri;
    memset(&pri, 0, sizeof pri);
    pri.nCh = cfg.nCh;
    pri.nSamples = cfg.nDdc;
    pri.ddcFs_Hz = (float)cfg.fsOut_Hz;
    pri.shortIdx[0] = (uint32_t)cfg.pulse[PULSE_SHORT].ddcFirst + 1;
    pri.shortIdx[1] = (uint32_t)(cfg.pulse[PULSE_SHORT].ddcFirst + cfg.pulse[PULSE_SHORT].ddcCount);
    pri.longIdx[0] = (uint32_t)cfg.pulse[PULSE_LONG].ddcFirst + 1;
    pri.longIdx[1] = (uint32_t)(cfg.pulse[PULSE_LONG].ddcFirst + cfg.pulse[PULSE_LONG].ddcCount);
    pri.iq = malloc(sizeof(int16_t) * 2 * (size_t)cfg.nCh * cfg.nDdc);
    double *steerRe = malloc(sizeof(double) * cfg.nCh), *steerIm = malloc(sizeof(double) * cfg.nCh);
    for (int ch = 0; ch < cfg.nCh; ch++) {               /* matched to the beam weights */
        if (!cfg.dbfEnabled) { steerRe[ch] = 1; steerIm[ch] = 0; continue; }
        const double wr = cfg.wRe[beam * cfg.nCh + ch], wi = cfg.wIm[beam * cfg.nCh + ch];
        const double m = hypot(wr, wi);
        steerRe[ch] = m > 0 ? wr / m : 1;
        steerIm[ch] = m > 0 ? -wi / m : 0;
    }

    Processor *proc = processor_create(&cfg, 0);
    CHECK(proc != NULL, "processor_create");
    CHECK(processor_cfar_alpha(proc) == cfg.alpha, "the CFAR threshold starts at the configuration value");
    if (!isnan(threshold_dB)) processor_set_cfar_alpha(proc, pow(10.0, threshold_dB / 10.0));
    Sink sink;
    memset(&sink, 0, sizeof sink);
    ProcessorCallbacks cb = { on_detections, on_plot, NULL, &sink, NULL };
    srand(42);
    double tProc = 0;
    for (int k = 0; k < nPulses; k++) {
        pri.pulseNumber = 1000 + k;
        pri.scanNumber = 1;
        pri.az_deg = (float)fmod(k * cfg.rate_dps * pri_s, 360.0);
        pri.time_s = k * pri_s;
        pri.staggerIndex = 1;
        for (size_t i = 0; i < 2 * (size_t)cfg.nCh * cfg.nDdc; i++) pri.iq[i] = sat16(16.0 * gauss());
        if (k < nWith) {
            const double ph = -4 * PI * vr * pri.time_s / pc->lambda_m;   /* baseband phase -4*pi*R/lambda */
            const int first = pc->ddcFirst + r0 + pc->mfShift;
            for (int ch = 0; ch < cfg.nCh; ch++) {
                const double cr = amp * (steerRe[ch] * cos(ph) - steerIm[ch] * sin(ph));
                const double ci = amp * (steerRe[ch] * sin(ph) + steerIm[ch] * cos(ph));
                for (int m = 0; m < nEcho; m++) {
                    /* waveform s[m]: the pulse, or conj(h[nTaps-1-m]) */
                    const double sr = pc->nPulse ? pc->pulseRe[m] / pulseNorm : pc->tapsRe[pc->nTaps - 1 - m];
                    const double si = pc->nPulse ? pc->pulseIm[m] / pulseNorm : -pc->tapsIm[pc->nTaps - 1 - m];
                    int16_t *v = pri.iq + 2 * ((size_t)ch * cfg.nDdc + first + m);
                    v[0] = sat16(v[0] + cr * sr - ci * si);
                    v[1] = sat16(v[1] + cr * si + ci * sr);
                }
            }
        }
        const double t0 = now_seconds();
        processor_process(proc, &pri, &cb);
        tProc += now_seconds() - t0;
    }
    processor_flush(proc, &cb);
    ProcessorTiming tm;
    processor_timing(proc, &tm);
    CHECK(tm.rejected == 0 && tm.pris == (uint64_t)nPulses, "%llu PRIs rejected", (unsigned long long)tm.rejected);
    printf("per PRI [ms]: DBF %.2f, MF %.2f, MTI %.2f, Doppler %.2f, CFAR %.2f, extraction %.2f, total %.2f\n",
           tm.dbf_ms / tm.pris, tm.mf_ms / tm.pris, tm.mti_ms / tm.pris, tm.doppler_ms / tm.pris,
           tm.cfar_ms / tm.pris, tm.extract_ms / tm.pris, tm.total_ms / tm.pris);

    printf("end to end: %d PRIs in %.1f ms (%.2f ms per PRI), %d/%d CPI updates, %d/%d detections, %d plots\n",
           nPulses, tProc * 1e3, tProc * 1e3 / nPulses, sink.updates[0], sink.updates[1],
           sink.dets[0], sink.dets[1], sink.plots.n);
    CHECK(sink.updates[PULSE_LONG] == nPulses - cfg.cpiLength + 1 || cfg.updateEvery > 1,
          "CPI updates %d", sink.updates[PULSE_LONG]);
    if (!isnan(threshold_dB)) {          /* a threshold above the target peak: nothing is detected */
        CHECK(sink.dets[0] + sink.dets[1] == 0 && sink.plots.n == 0,
              "CFAR threshold %.1f dB: %d detections, %d plots (expected none)", threshold_dB,
              sink.dets[0] + sink.dets[1], sink.plots.n);
        processor_destroy(proc);
        free(steerRe);
        free(steerIm);
        free(pri.iq);
        radar_config_free(&cfg);
        return sink.plots.n == 0;
    }
    int found = 0;
    for (int i = 0; i < sink.plots.n; i++) {
        const Detection *p = &sink.plots.plots[i];
        printf("  plot: %s, range %.1f m (cell %d), vr %.2f m/s, beam %d, el %.2f deg, SNR %.1f dB, hits %d\n",
               p->pulse == PULSE_LONG ? "long" : "short", p->range_m, p->rIdx, p->vr_mps, p->beam,
               p->el_deg, p->snr_dB, p->hits);
        if (p->pulse == PULSE_LONG && abs(p->rIdx - r0) <= 1 && fabs(p->vr_mps - vr) < 0.5 &&
            abs(p->beam - beam) <= 1) {
            found++;
            *plotOut = *p;
        }
    }
    CHECK(found == 1, "expected the target as one long-pulse plot at cell %d, %.1f m/s, beam %d", r0, vr, beam);
    CHECK(sink.plots.n == 1, "expected exactly one plot, got %d", sink.plots.n);

    processor_destroy(proc);
    free(steerRe);
    free(steerIm);
    free(pri.iq);
    radar_config_free(&cfg);
    return found == 1;
}

/* False alarm probability of the CFAR threshold (the inverse of cfar_alpha.m). */
static void test_cfar_pfa(const char *configPath)
{
    /* one look: closed form (1 + alpha/N)^-N */
    for (int n = 8; n <= 300; n *= 3) {
        const double alpha = n * (pow(1e-6, -1.0 / n) - 1.0);
        const double pfa = cfar_pfa(alpha, n, 1.0);
        CHECK(fabs(pfa / 1e-6 - 1.0) < 1e-9, "cfar_pfa one look, N = %d: %.6e (expected 1e-6)", n, pfa);
    }
    CHECK(fabs(cfar_pfa(0.0, 100, 1.5) - 1.0) < 1e-12, "cfar_pfa(alpha = 0) = 1");
    CHECK(isnan(cfar_pfa(10.0, 100, 0.0)), "cfar_pfa without nLooks is NaN");
    /* the configuration: alpha was designed for cfg.cfar.pfa = 1e-8 with nLooks looks */
    RadarConfig cfg;
    char err[256];
    if (radar_config_load(configPath, &cfg, err, sizeof err) != 0) return;
    if (cfg.nLooks > 0) {
        const double pfa = cfar_pfa(cfg.alpha, cfg.nTrain, cfg.nLooks);
        printf("CFAR: alpha %.3f (%.2f dB), %d training cells, %.2f looks -> Pfa %.3e\n", cfg.alpha,
               10 * log10(cfg.alpha), cfg.nTrain, cfg.nLooks, pfa);
        CHECK(fabs(pfa / 1e-8 - 1.0) < 1e-3, "Pfa of the configured threshold %.4e (expected 1e-8)", pfa);
        CHECK(cfar_pfa(cfg.alpha * 2.0, cfg.nTrain, cfg.nLooks) < pfa, "a higher threshold has a lower Pfa");
    } else {
        printf("configuration without nLooks: Pfa check skipped\n");
    }
    radar_config_free(&cfg);
}

int main(int argc, char **argv)
{
    test_fft();
    test_extractor();
    if (argc > 1) {
        /* fast convolution (where cheaper) and direct matched filter: same plot */
        Detection a, b;
        if (test_end_to_end(argv[1], 0, NAN, &a) && test_end_to_end(argv[1], 1, NAN, &b)) {
            CHECK(a.rIdx == b.rIdx && a.dIdx == b.dIdx && a.beam == b.beam && a.hits == b.hits,
                  "fast and direct matched filter give different plots");
            CHECK(fabs(a.snr_dB - b.snr_dB) < 0.01 && fabs(a.power / b.power - 1) < 1e-4,
                  "fast and direct matched filter: SNR %.4f / %.4f dB, power ratio %.6f",
                  a.snr_dB, b.snr_dB, a.power / b.power);
        }
        /* the same target with the CFAR threshold raised far above its peak (run-time change) */
        test_end_to_end(argv[1], 0, 45.0, &a);
        test_cfar_pfa(argv[1]);
    } else {
        printf("no configuration file given: end-to-end test skipped\n");
    }
    if (failures) printf("%d FAILURE(S)\n", failures);
    else printf("all core tests passed\n");
    return failures ? 1 : 0;
}

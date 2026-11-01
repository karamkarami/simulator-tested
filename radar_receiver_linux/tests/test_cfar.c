/*
 * test_cfar.c - the CA-CFAR of processor.c against a direct computation.
 *
 * White-box test: it includes processor.c to call its static cfar(). On maps
 * of exponential noise (one per method), every decision of cfar() must equal
 * that of the direct sums over the training cells combined with the direct
 * local peak test (cfar() keeps only the detections that are the maximum of
 * their 3 x 3 x 3 neighbourhood), for the 2-D, the 1-D range and the 1-D
 * Doppler window and for tested range intervals in the middle, at both range
 * edges (shifted windows, also when the whole interval lies within the first
 * or last guard+training cells) and over all cells; with the beams of the
 * configuration (argv[1], which gives the window sizes) and with one beam
 * split into range chunks over several threads; the noise estimate kept for
 * the SNR of every detection (cfg.report.snrNoise = 'cfar') equals the
 * direct one. Then the exact noise floor
 * median (check_median) against a sorted copy, and the pruning of thousands
 * of detections (check_pruning) against the direct pairwise rules.
 */
#include "../src/processor.c"

#include <stdio.h>

static int failures = 0;

/* Noise estimate of cell (b, r, d) from direct sums: the outer window centred
 * on the nearest cell whose window lies inside the data (shifted inwards at
 * the range edges), minus the guard cells inside the data; Doppler circular. */
static double direct_noise(const RadarConfig *cfg, const PulseState *s, int b, int r, int d)
{
    const int nR = s->nR, nfft = cfg->nfft;
    const int gR = cfg->guardRange, oR = gR + cfg->trainRange, gD = cfg->guardDoppler, oD = gD + cfg->trainDoppler;
    const int rc = r < oR ? oR : r > nR - 1 - oR ? nR - 1 - oR : r;
    double sum = 0.0;
    int n = 0;
    for (int dd = -oD; dd <= oD; dd++)
        for (int rr = rc - oR; rr <= rc + oR; rr++) {
            const int inGuard = abs(rr - r) <= gR && abs(dd) <= gD;
            if (inGuard) continue;
            sum += P_at(s, ((d + dd) % nfft + nfft) % nfft, b, rr);
            n++;
        }
    return sum / n;
}

/* Local peak from direct comparisons with the 26 neighbours (Doppler circular). */
static int direct_peak(const RadarConfig *cfg, const PulseState *s, int b, int r, int d)
{
    const float v = P_at(s, d, b, r);
    for (int dd = -1; dd <= 1; dd++)
        for (int bb = b - 1; bb <= b + 1; bb++)
            for (int rr = r - 1; rr <= r + 1; rr++)
                if (bb >= 0 && bb < cfg->nBeams && rr >= 0 && rr < s->nR &&
                    P_at(s, (d + dd + cfg->nfft) % cfg->nfft, bb, rr) > v) return 0;
    return 1;
}

/* One CFAR window (guard and training cells per side in range and Doppler),
 * with nBeams beams (0: those of the configuration) and nThreads threads (0:
 * the default). */
static void check_window(RadarConfig cfg, const char *name, int gR, int tR, int gD, int tD, int nBeams, int nThreads)
{
    char label[64];
    cfg.guardRange = gR;  cfg.trainRange = tR;  cfg.guardDoppler = gD;  cfg.trainDoppler = tD;
    cfg.nTrain = (2 * (gR + tR) + 1) * (2 * (gD + tD) + 1) - (2 * gR + 1) * (2 * gD + 1);
    if (nBeams > 0) cfg.nBeams = nBeams;                /* the CFAR only needs the maps */
    cfg.snrNoise = SNR_NOISE_CFAR;                      /* the CFAR keeps the noise of its detections */
    Processor *p = processor_create(&cfg, nThreads);
    PulseState *s = &p->ps[PULSE_LONG];
    const int nR = s->nR, oR = gR + tR, nfft = cfg.nfft, failuresBefore = failures;
    snprintf(label, sizeof label, "%s, %d beam(s), %d thread(s)", name, cfg.nBeams, p->nThreads);
    srand(11);
    for (size_t i = 0; i < (size_t)nfft * s->nCells; i++)
        s->P[i] = (float)-log((rand() + 1.0) / (RAND_MAX + 2.0));     /* exponential, mean 1 */
    processor_set_cfar_alpha(p, 4.0);                                  /* many detections to compare */
    const int intervals[][2] = { { nR / 2 - 4, 9 }, { 0, 3 }, { 0, oR }, { 2, oR }, { nR - 3, 3 },
                                 { nR - oR, oR }, { nR - oR - 2, oR }, { 0, 1 }, { nR - 1, 1 }, { 0, nR } };
    for (size_t k = 0; k < sizeof intervals / sizeof intervals[0]; k++) {
        PulseConfig pc = cfg.pulse[PULSE_LONG];
        pc.testFirst = intervals[k][0];
        pc.testCount = intervals[k][1];
        cfar(p, s, &pc);
        int wrong = 0, detected = 0;
        for (int d = 0; d < nfft; d++)
            for (int b = 0; b < cfg.nBeams; b++)
                for (int r = 0; r < nR; r++) {
                    const int tested = r >= pc.testFirst && r < pc.testFirst + pc.testCount &&
                                       abs(d - nfft / 2) > cfg.notchBins;
                    const unsigned char got = s->det[(size_t)d * s->nCells + (size_t)b * nR + r];
                    detected += got;
                    if (!tested) { wrong += got; continue; }
                    if (b != 0 && b != cfg.nBeams - 1) continue;        /* direct sums: first and last beam */
                    const double thr = p->alpha * direct_noise(&cfg, s, b, r, d), v = P_at(s, d, b, r);
                    if (fabs(v - thr) > 1e-5 * thr && got != (v > thr && direct_peak(&cfg, s, b, r, d)))
                        wrong++;                                        /* (ties: float sums) */
                    if (got) {                                          /* its noise for the SNR */
                        const double nz = s->noise[(size_t)d * s->nCells + (size_t)b * nR + r];
                        if (fabs(nz * p->alpha - thr) > 1e-5 * thr) wrong++;
                    }
                }
        if (wrong) {
            failures++;
            printf("FAIL %s CFAR, tested cells %d ... %d of %d: %d wrong decisions (%d detections)\n", label,
                   pc.testFirst, pc.testFirst + pc.testCount - 1, nR, wrong, detected);
        }
    }
    printf("%-32s window (range %d+%d, Doppler %d+%d): %s\n", label, gR, tR, gD, tD,
           failures > failuresBefore ? "FAIL" : "every decision equals the direct sums and peak test");
    processor_destroy(p);
}

/* The noise floor median against a sorted copy: maps of several sizes (odd,
 * even, one value), with ties, zeros and a structured half, twice each (the
 * second time with the bracket from the first median), and once with a
 * wrong guess. */
static int cmp_float(const void *a, const void *b)
{
    const float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

static void check_median(void)
{
    const int failuresBefore = failures;
    float *map = malloc(sizeof(float) * 200000), *work = malloc(sizeof(float) * (200000 + 8));
    float *sorted = malloc(sizeof(float) * 200000);
    srand(3);
    for (int t = 0; t < 400; t++) {
        const int rows = 1 + rand() % 128, len = 1 + rand() % (t % 4 == 0 ? 8 : 800), stride = len + rand() % 4;
        int n = 0;
        for (int d = 0; d < rows; d++)
            for (int r = 0; r < stride; r++) {
                float v = (float)-log((rand() + 1.0) / (RAND_MAX + 2.0));
                if (t % 5 == 1) v = (float)(rand() % 4);                 /* many ties */
                if (t % 5 == 2 && rand() % 3 == 0) v = 0.0f;              /* zeros */
                if (t % 5 == 3 && d < rows / 2) v *= 1e5f;                /* structured */
                map[d * stride + r] = v;
                if (r < len) sorted[n++] = v;
            }
        qsort(sorted, n, sizeof(float), cmp_float);
        const double want = (n % 2) ? sorted[n / 2] : 0.5 * ((double)sorted[n / 2 - 1] + sorted[n / 2]);
        float guess = 0.0f;
        for (int k = 0; k < 3; k++) {
            if (k == 2) guess = 1e3f;                                    /* far from the median */
            const double got = map_median(map, rows, len, (size_t)stride, work, &guess);
            if (got != want) {
                failures++;
                if (failures - failuresBefore <= 5)
                    printf("FAIL median of %d x %d values (call %d): %.9g, sorted %.9g\n", rows, len, k, got, want);
            }
        }
    }
    printf("noise floor median: %s\n", failures > failuresBefore ? "FAIL" : "equals the sorted copy on 400 maps");
    free(map);
    free(work);
    free(sorted);
}

/* Pruning under load: a long-pulse map with a strong target (range
 * sidelobes along its line) and a low threshold gives thousands of CFAR
 * detections (local peaks); prune_and_report() must keep exactly the cells
 * of the direct rules (local peak, sidelobe blanking pair by pair). nStrong
 * more strong cells (like clutter residue) make the pairwise comparisons
 * expensive: the pruning then takes the map maxima (blank_dense), with the
 * same result and a time bounded by the map size. */
static void check_pruning(RadarConfig cfg, const char *name, int nStrong, double alpha)
{
    Processor *p = processor_create(&cfg, 0);
    PulseState *s = &p->ps[PULSE_LONG];
    const PulseConfig *pc = &cfg.pulse[PULSE_LONG];
    const int nR = s->nR, nB = cfg.nBeams, nfft = cfg.nfft, tb = nB / 2, tr = nR / 3, td = nfft / 4;
    srand(5);
    for (size_t i = 0; i < (size_t)nfft * s->nCells; i++)
        s->P[i] = (float)-log((rand() + 1.0) / (RAND_MAX + 2.0));
    for (int r = tr - 80; r <= tr + 80; r++)                              /* sidelobes of the target */
        if (r >= 0 && r < nR) s->P[(size_t)td * s->nCells + (size_t)tb * nR + r] += 1e4f * (float)(rand() + 1) / RAND_MAX;
    s->P[(size_t)td * s->nCells + (size_t)tb * nR + tr] = 1e7f;
    for (int k = 0; k < nStrong; k++)                                     /* strong scattered cells */
        s->P[(size_t)(rand() % nfft) * s->nCells + (size_t)(rand() % nB) * nR + rand() % nR] =
            1e3f + 1e5f * (float)rand() / RAND_MAX;
    processor_set_cfar_alpha(p, alpha);
    PulseConfig all = *pc;
    all.testFirst = 0;  all.testCount = nR;
    cfar(p, s, &all);
    const double t0 = now_seconds();
    prune_and_report(p, PULSE_LONG, s);
    const double dt = now_seconds() - t0;

    /* the direct rules */
    Cell *c = malloc(sizeof(Cell) * (size_t)nfft * s->nCells);
    int n = 0, nDet = 0;
    for (int d = 0; d < nfft; d++)
        for (int r = 0; r < nR; r++)
            for (int b = 0; b < nB; b++) {
                if (!s->det[(size_t)d * s->nCells + (size_t)b * nR + r]) continue;
                nDet++;
                const float v = P_at(s, d, b, r);
                int peak = 1;
                for (int dd = -1; dd <= 1; dd++)
                    for (int bb = b - 1; bb <= b + 1; bb++)
                        for (int rr = r - 1; rr <= r + 1; rr++)
                            if (bb >= 0 && bb < nB && rr >= 0 && rr < nR && P_at(s, (d + dd + nfft) % nfft, bb, rr) > v) peak = 0;
                if (peak) c[n++] = (Cell){ b, r, d, v };
            }
    const double ratio = pow(10.0, cfg.sidelobeMargin_dB / 10.0);
    int nKeep = 0, wrong = 0, k = 0;
    for (int i = 0; i < n; i++) {
        int keep = 1;
        for (int j = 0; j < n && keep; j++) {
            if (!((double)c[j].pw > ratio * c[i].pw)) continue;
            int dD = abs(c[i].d - c[j].d);  if (nfft - dD < dD) dD = nfft - dD;
            const int dR = abs(c[i].r - c[j].r);
            if ((abs(c[i].b - c[j].b) <= 1 && dD <= 1 && dR <= pc->nTaps) || (dR <= 1 && dD <= 1) ||
                (dR <= 1 && abs(c[i].b - c[j].b) <= 1)) keep = 0;
        }
        if (!keep) continue;
        nKeep++;
        if (k >= s->nReports || s->reports[k].beam != c[i].b || s->reports[k].rIdx != c[i].r || s->reports[k].dIdx != c[i].d)
            wrong++;
        k++;
    }
    if (wrong || k != s->nReports || dt > 0.5) {
        failures++;
        printf("FAIL pruning (%s): %d reports, %d expected, %d different, %.3f s\n", name, s->nReports, nKeep, wrong, dt);
    }
    printf("pruning (%s): %d CFAR detections (local peaks: %d), %d kept (the direct rules: %d) in %.1f ms: %s\n", name,
           nDet, n, s->nReports, nKeep, dt * 1e3, wrong || k != s->nReports || dt > 0.5 ? "FAIL" : "ok");
    free(c);
    processor_destroy(p);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("usage: test_cfar CONFIG\n");
        return 2;
    }
    RadarConfig cfg;
    char err[256];
    if (radar_config_load(argv[1], &cfg, err, sizeof err) != 0) {
        printf("configuration %s: %s\n", argv[1], err);
        return 1;
    }
    for (int one = 0; one < 2; one++) {                 /* the configured beams; one beam in 4 threads */
        const int nB = one ? 1 : 0, nThr = one ? 4 : 0;
        check_window(cfg, "2-D", 3, 8, 3, 4, nB, nThr);
        check_window(cfg, "1-D range", 3, 16, 0, 0, nB, nThr);
        check_window(cfg, "1-D Doppler", 0, 0, 2, 8, nB, nThr);
        check_window(cfg, "configured", cfg.guardRange, cfg.trainRange, cfg.guardDoppler, cfg.trainDoppler, nB, nThr);
    }
    check_median();
    check_pruning(cfg, "strong target, low threshold", 0, 2.5);
    check_pruning(cfg, "3000 strong cells, very low threshold", 3000, 1.2);
    radar_config_free(&cfg);
    if (failures) printf("%d FAILURE(S)\n", failures);
    else printf("all CFAR and pruning tests passed\n");
    return failures ? 1 : 0;
}

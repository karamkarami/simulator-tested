/*
 * test_cfar.c - the CA-CFAR of processor.c against a direct computation.
 *
 * White-box test: it includes processor.c to call its static cfar(). On maps
 * of exponential noise (one per method), every decision of cfar() must equal
 * that of the direct sums over the training cells, for the 2-D, the 1-D range
 * and the 1-D Doppler window and for tested range intervals in the middle, at
 * both range edges (shifted windows, also when the whole interval lies within
 * the first or last guard+training cells) and over all cells. The window of
 * the configuration (argv[1]) gives the sizes.
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

/* One CFAR window (guard and training cells per side in range and Doppler). */
static void check_window(RadarConfig cfg, const char *name, int gR, int tR, int gD, int tD)
{
    cfg.guardRange = gR;  cfg.trainRange = tR;  cfg.guardDoppler = gD;  cfg.trainDoppler = tD;
    cfg.nTrain = (2 * (gR + tR) + 1) * (2 * (gD + tD) + 1) - (2 * gR + 1) * (2 * gD + 1);
    Processor *p = processor_create(&cfg, 0);
    PulseState *s = &p->ps[PULSE_LONG];
    const int nR = s->nR, oR = gR + tR, nfft = cfg.nfft, failuresBefore = failures;
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
                    if (fabs(v - thr) > 1e-5 * thr && got != (v > thr)) wrong++;   /* (ties: float sums) */
                }
        if (wrong) {
            failures++;
            printf("FAIL %s CFAR, tested cells %d ... %d of %d: %d wrong decisions (%d detections)\n", name,
                   pc.testFirst, pc.testFirst + pc.testCount - 1, nR, wrong, detected);
        }
    }
    printf("%-12s window (range %d+%d, Doppler %d+%d): %s\n", name, gR, tR, gD, tD,
           failures > failuresBefore ? "FAIL" : "every decision equals the direct sums");
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
    check_window(cfg, "2-D", 3, 8, 3, 4);
    check_window(cfg, "1-D range", 3, 16, 0, 0);
    check_window(cfg, "1-D Doppler", 0, 0, 2, 8);
    check_window(cfg, "configured", cfg.guardRange, cfg.trainRange, cfg.guardDoppler, cfg.trainDoppler);
    radar_config_free(&cfg);
    if (failures) printf("%d FAILURE(S)\n", failures);
    else printf("all CFAR tests passed\n");
    return failures ? 1 : 0;
}

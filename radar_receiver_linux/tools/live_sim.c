/*
 * live_sim.c - radar_live_sim: the simulator of radar_generate_ddc in real time (a run
 * configuration of the same name in CLion). Included at the end of generate_ddc.c (compiled with
 * RADAR_LIVE_SIM): the same simulation of every PRI (targets, terrain clutter, leakage, STC,
 * noise, ADC, fixed-point DDC), but instead of files every PRI goes to the receiver over UDP at
 * its PRI time, as the board would send it, until Ctrl+C.
 *
 *   radar_live_sim [--config-m FILE] [--set STATEMENT]... [--host H] [--port N]
 *                  [--format rddc|vbdr] [--payload BYTES] [--threads N] [--targets FILE.csv]
 *                  [--fidelity full|fast] [--board-port N] [--board-reply HOST:PORT] [--no-board]
 *                  [--max-pulses N] [--duration S] [--speed X] [--quiet]
 *     --config-m, --set, --targets, --threads   as radar_generate_ddc
 *     --host, --port     the receiver (default cfg.udp.host, cfg.udp.port)
 *     --format           'rddc' or 'vbdr' (default cfg.udp.format)
 *     --payload BYTES    bytes per datagram (default cfg.udp.payloadBytes, VBDR: cfg.udp.vbdr.payloadBytes
 *                        of the mode)
 *     --fidelity         'full' (as radar_generate_ddc) or 'fast': the signal at the DDC rate (targets,
 *                        leakage, STC and noise, no terrain clutter, no ADC clipping; the gains of the
 *                        ADC and DDC are measured once): much less computation (default cfg.sim.liveFidelity)
 *     --board-port N     the board commands of the receiver (default cfg.board.port), answered to
 *     --board-reply H:P  (default cfg.board.pcHost:cfg.board.replyPort); --no-board: not
 *     --max-pulses N     stop after N PRIs (also those outside the sectors); --duration S: after S
 *                        seconds of simulated time (default: until Ctrl+C)
 *     --speed X          the PRIs go out X times as fast as the real time (default 1; e.g. 0.2 when the
 *                        receiver runs on the same, slow computer)
 *
 * Threads: the producers (--threads, default all cores) simulate the next PRIs ahead into a ring
 * of PRIs; the main thread sends every PRI at its time (one clock: PRI start times of the
 * schedule), answers the board commands in the idle time and prints the progress. When the
 * producers cannot keep up, the PRIs go out later and the simulated time runs slower than the
 * real time (the PRTs in the stream stay those of the radar; a warning says how much slower).
 *
 * Schedule: endless (cfg.sim.nScans is ignored); the antenna turns at the rotation speed of the
 * mode; cfg.sim.sectors_deg: only the PRIs in the sectors are simulated and sent, the time and
 * the azimuth advance outside them too. The targets move on (cfg.scenario or --targets).
 *
 * Board commands (as the board, README "VBDR stream and board control"): every command is printed
 * and answered; the registers start with the settings of the simulation. A write of the radar mode,
 * the phase code, F1, F2, the stagger pattern or (cfg.board.sendRpm) the rotation speed changes the
 * simulation (when no further command came for 0.3 s): the configuration of the new settings is
 * set up (radar_config.m with cfg.modes.active and the settings), the PRIs simulated ahead are
 * dropped, and the schedule goes on from the time and azimuth reached (pulse numbers continue).
 * STC, gain and the operation mode are only stored and answered.
 */
#if defined(_WIN32)
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <winsock2.h>
#  include <ws2tcpip.h>
typedef SOCKET lsock_t;
#else
#  include <arpa/inet.h>
#  include <netinet/in.h>
#  include <sys/socket.h>
typedef int lsock_t;
#endif
#include <signal.h>

#include "radar_rx.h"
#include "thread_compat.h"
#include "ddc_files.h"

#define LIVE_RING      64        /* PRIs simulated ahead                                   */
#define LIVE_DEBOUNCE  0.3       /* s without a board command before the settings change    */

static volatile sig_atomic_t g_liveStop = 0;
static void live_on_signal(int sig) { (void)sig; g_liveStop = 1; }

/* Waits until now_seconds() = t: sleeps for the coarse part, busy-waits the last 2 ms. */
static void live_wait_until(double t)
{
    for (;;) {
        const double dt = t - now_seconds();
        if (dt <= 0) return;
        if (dt > 3e-3) {
#ifdef _WIN32
            Sleep((DWORD)((dt - 2e-3) * 1e3));
#else
            const double sl = dt - 2e-3;
            struct timespec ts = { (time_t)sl, (long)((sl - (double)(time_t)sl) * 1e9) };
            nanosleep(&ts, NULL);
#endif
        }
    }
}

/* ------------------------------------------------------------------------ */
/* Fidelity 'fast': the signal at the DDC rate                               */
/* ------------------------------------------------------------------------ */

typedef struct {
    double gain;        /* DDC output per LSB of the IF amplitude                               */
    double sign;        /* +1: the output phase follows the IF phase, -1: its conjugate          */
    double theta;       /* output phase of an IF cosine of phase 0                               */
    double noise;       /* std of I (and Q) per LSB rms of white noise at the ADC                */
} FastCal;

/* Measures the DDC with a tone at the IF and with white noise (the run_ddc of this setup). */
static void fast_calibrate(const Radar *s, Work *wk, FastCal *cal)
{
    const int nAdc = s->nAdc, nCh = s->nCh, m = s->nDdc / 2;
    const double a = 1000.0;
    double re[2], im[2];
    for (int k2 = 0; k2 < 2; k2++) {
        const double phi = k2 ? PI / 2 : 0;
        for (int ch = 0; ch < nCh; ch++)
            for (int k = 0; k < nAdc; k++)
                wk->adc[(size_t)ch * nAdc + k] = (int16_t)lround(a * cos(2 * PI * s->fIF * (s->m0 + k) / s->fs + phi));
        run_ddc(s, wk);
        re[k2] = wk->I[m];
        im[k2] = wk->Q[m];
    }
    cal->gain = hypot(re[0], im[0]) / a;
    cal->theta = atan2(im[0], re[0]);
    const double d = atan2(im[1], re[1]) - cal->theta;
    cal->sign = sin(d) >= 0 ? 1.0 : -1.0;
    Rng r;
    rng_seed(&r, 12345, 1);
    for (int ch = 0; ch < nCh; ch++)
        for (int k = 0; k < nAdc; k++) wk->adc[(size_t)ch * nAdc + k] = (int16_t)lround(a * rng_normal(&r));
    run_ddc(s, wk);
    double acc = 0;
    int n = 0;
    for (int k = s->nDdc / 8; k < s->nDdc; k++, n++) acc += (double)wk->I[k] * wk->I[k];
    cal->noise = n ? sqrt(acc / n) / a : 0;
}

/* One PRI at the DDC rate: targets, leakage (with the STC), thermal and ADC noise. */
static void simulate_fast(const Radar *s, const Target *tg, int nTg, const Pri *pr, uint64_t seed, Work *wk,
                          const FastCal *cal)
{
    const int nDdc = s->nDdc, nCh = s->nCh;
    const double c = s->c, toLsb = s->voltGain * s->lsbPerUnit * cal->gain;
    double *zr = wk->xI, *zi = wk->xQ;              /* [ch][nDdc] (nAdc >= nCh * nDdc / ... : sized below) */
    memset(zr, 0, sizeof(double) * (size_t)nDdc * nCh);
    memset(zi, 0, sizeof(double) * (size_t)nDdc * nCh);
    /* the ADC sample of every DDC sample (STC, LO correction) */
    for (int k = 0; k < nDdc; k++) {
        const double t = s->ddcT0 + k / s->fsOut;
        int n = (int)lround(t * s->fs) - s->m0;
        if (n < 0) n = 0;
        if (n >= s->nAdc) n = s->nAdc - 1;
        wk->y1[k] = n;
    }
    for (int i = 0; i < nTg; i++) {
        double R, az, el, vr;
        target_geometry(&tg[i], pr->t, s->height, &R, &az, &el, &vr);
        const double azRel = wrap180(az - pr->az);
        if (fabs(azRel) > s->azSpan) continue;
        const double u = -sind(s->tilt) * cosd(el) * cosd(azRel) + cosd(s->tilt) * sind(el);
        const double elRel = asin(u) * 180 / PI;
        const double ga = pattern_gain(&s->az, azRel);
        const double gain = ga * ga * pattern_gain(&s->txEl, elRel) * pattern_gain(&s->elElem, elRel);
        const double rot = s->rate * pr->t, rel = mod360(az - s->az0);
        int look = (int)floor((rot - rel + 180) / 360) + 1;
        if (look < 1) look = 1;
        const double rcs = tg[i].rcsScan[(look < tg[i].nScan ? look : tg[i].nScan) - 1];
        for (int p = 0; p < 2; p++) {
            const Wave *w = &s->w[p];
            if (!w->enable) continue;
            const double amp = sqrt(2 * w->radarConst * rcs / pow(R, 4) * gain) * toLsb;
            const double tau0 = 2 * R / c;
            for (int k = 0; k < nDdc; k++) {
                const double t = s->ddcT0 + k / s->fsOut;
                if (!(t >= w->tTx + tau0 && t < w->tTx + tau0 + w->width)) continue;
                const double tau = 2 * (R + vr * t) / c;
                double er, ei;
                envelope(w, t - w->tTx - tau0, &er, &ei);
                const double psi = atan2(ei, er) - 2 * PI * (w->fIF * (w->tTx + tau) + w->fLO * tau);
                double a = amp * hypot(er, ei);
                const int n = (int)wk->y1[k];
                double ph = cal->sign * psi + cal->theta;
                if (s->loShift[p]) ph += cal->sign * atan2(s->loShift[p][2 * n + 1], s->loShift[p][2 * n]);
                if (s->stcOn) a *= s->stc[n];
                const double br = a * cos(ph), bi = a * sin(ph);
                for (int ch = 0; ch < nCh; ch++) {
                    const double e = cal->sign * 2 * PI * s->spacing_m / w->lambda * ch * u;
                    const double cr = cos(e), ci = sin(e);
                    zr[(size_t)ch * nDdc + k] += br * cr - bi * ci;
                    zi[(size_t)ch * nDdc + k] += br * ci + bi * cr;
                }
            }
        }
    }
    if (s->addLeak)
        for (int p = 0; p < 2; p++) {
            const Wave *w = &s->w[p];
            if (!w->enable) continue;
            for (int k = 0; k < nDdc; k++) {
                const double t = s->ddcT0 + k / s->fsOut;
                if (!(t >= w->tTx && t < w->tTx + w->width)) continue;
                double er, ei;
                envelope(w, t - w->tTx, &er, &ei);
                const int n = (int)wk->y1[k];
                double a = s->leakAmp * hypot(er, ei) * toLsb;
                if (s->stcOn) a *= s->stc[n];
                const double ph = cal->sign * (atan2(ei, er) - 2 * PI * w->fIF * w->tTx) + cal->theta;
                for (int ch = 0; ch < nCh; ch++) {
                    zr[(size_t)ch * nDdc + k] += a * cos(ph);
                    zi[(size_t)ch * nDdc + k] += a * sin(ph);
                }
            }
        }
    Rng rng;
    rng_seed(&rng, seed, (uint64_t)pr->pulse);
    const double sigIn2 = (s->addNoise ? pow(s->sigmaTh * s->voltGain * s->lsbPerUnit, 2) : 0) +
                          (s->addAdcNoise ? s->sigmaAdc * s->sigmaAdc : 0) + 1.0 / 12;
    const double sig = cal->noise * sqrt(sigIn2);
    for (int ch = 0; ch < nCh; ch++)
        for (int k = 0; k < nDdc; k++) {
            const size_t q = (size_t)ch * nDdc + k;
            wk->I[q] = sat16(round(zr[q] + sig * rng_normal(&rng)));
            wk->Q[q] = sat16(round(zi[q] + sig * rng_normal(&rng)));
        }
}

/* ------------------------------------------------------------------------ */
/* Work buffers, setup and freeing                                           */
/* ------------------------------------------------------------------------ */

static void work_free(Work *wk)
{
    free(wk->x); free(wk->adc); free(wk->I); free(wk->Q); free(wk->xI); free(wk->xQ); free(wk->y1);
    free(wk->xc); free(wk->fr); free(wk->fi); free(wk->reflRe); free(wk->reflIm);
    for (int i = 0; i < 2; i++) if (wk->fft[i].n) fft_free(&wk->fft[i]);
    memset(wk, 0, sizeof *wk);
}

/* The buffers of run() for this setup (also those of the fidelity 'fast'). 0 on success. */
static int work_alloc(Work *wk, const Radar *s, const Clutter *C)
{
    memset(wk, 0, sizeof *wk);
    const int nPages = C && !C->sameCarrier ? 2 : 1, G = C ? C->G : 1;
    const size_t nZ = (size_t)s->nDdc * s->nCh > (size_t)s->nAdc ? (size_t)s->nDdc * s->nCh : (size_t)s->nAdc;
    wk->x = calloc((size_t)s->nAdc * s->nCh, sizeof(float));
    wk->adc = calloc((size_t)s->nAdc * s->nCh, sizeof(int16_t));
    wk->I = calloc((size_t)s->nDdc * s->nCh, sizeof(int16_t));
    wk->Q = calloc((size_t)s->nDdc * s->nCh, sizeof(int16_t));
    wk->xI = calloc(nZ, sizeof(double));
    wk->xQ = calloc(nZ, sizeof(double));
    wk->y1 = calloc((size_t)s->nAdc + 2 * (size_t)s->nCh + (size_t)s->nDdc, sizeof(double));
    if (C) {
        int nFftMax = 1;
        for (int i = 0; i < 2; i++) if (C->nFft[i] > nFftMax) nFftMax = C->nFft[i];
        wk->xc = calloc(2 * (size_t)nPages * G * s->nAdc, sizeof(double));
        wk->fr = calloc((size_t)nFftMax, sizeof(double));
        wk->fi = calloc((size_t)nFftMax, sizeof(double));
        wk->reflRe = calloc((size_t)C->nR * G, sizeof(double));
        wk->reflIm = calloc((size_t)C->nR * G, sizeof(double));
        for (int i = 0; i < 2; i++) if (s->w[i].enable) fft_init(&wk->fft[i], C->nFft[i]);
        if (!wk->xc || !wk->fr || !wk->fi || !wk->reflRe || !wk->reflIm) return -1;
    }
    return wk->x && wk->adc && wk->I && wk->Q && wk->xI && wk->xQ && wk->y1 ? 0 : -1;
}

static void pattern_free(Pattern *p) { free(p->x); free(p->g_dB); }

static void radar_free(Radar *s)
{
    free(s->cosNco); free(s->msinNco); free(s->stc); free(s->rxSample); free(s->pri); free(s->scanPri);
    free(s->h1); free(s->h2);
    for (int p = 0; p < 2; p++) {
        free(s->loShift[p]);
        free(s->w[p].codeRe); free(s->w[p].codeIm); free(s->w[p].misRe); free(s->w[p].misIm);
    }
    pattern_free(&s->az); pattern_free(&s->txEl); pattern_free(&s->elElem);
    memset(s, 0, sizeof *s);
}

static void clutter_free(Clutter *C)
{
    free(C->az); free(C->range); free(C->delay); free(C->groupEl); free(C->rowStart); free(C->cellRow);
    free(C->cellRange); free(C->cellGroup); free(C->rhoRe); free(C->rhoIm); free(C->fd); free(C->dcRe); free(C->dcIm);
    for (int i = 0; i < 2; i++) { free(C->steer[i]); free(C->txRe[i]); free(C->txIm[i]); }
    memset(C, 0, sizeof *C);
}

/* ------------------------------------------------------------------------ */
/* The simulator                                                             */
/* ------------------------------------------------------------------------ */

typedef struct {
    Pri      pr;
    double   prevPri;
    double   azU;         /* unwrapped azimuth of the PRI start                   */
    int      k;           /* ring index of the PRI (sequence number)              */
    int      active;      /* inside the sectors: simulated and sent              */
    int      ready;
    uint64_t gen;
    int      nDdc, nCh;
    int16_t *iq;          /* I, Q pairs, channel after channel                     */
    size_t   cap;
} Slot;

typedef struct {
    /* configuration */
    const char  *mPath;
    const char **baseSets;
    int          nBase;
    BoardInfo    bi;
    int          haveBi;
    int          mode;
    ModeSettings set;
    MValue      *cfg;
    Radar        s;
    Clutter      clut;
    int          haveClut;
    Target      *tg;
    int          nTg;
    uint64_t     seed;
    int          fast;
    FastCal      cal;
    double      *sect;
    int          nSect;
    int          nThreads;
    const char  *targetsCsv;
    /* schedule */
    int          nextPulse, stPos;
    double       nextT, nextAzU, lastPri;
    /* ring */
    Slot         ring[LIVE_RING];
    Mutex        mu;
    Cond         cv;
    uint64_t     gen;
    int          nextK, sendK, busy, paused, stop;
} Live;

static int in_sectors(const Live *L, double az)
{
    if (L->nSect == 0) return 1;
    for (int i = 0; i < L->nSect; i++) {
        double width = mod360(L->sect[L->nSect + i] - L->sect[i]);
        if (width == 0 && L->sect[L->nSect + i] != L->sect[i]) width = 360;
        if (mod360(az - L->sect[i]) <= width) return 1;
    }
    return 0;
}

/* The next PRI of the schedule (under the lock). */
static void sched_next(Live *L, Pri *pr, double *prevPri, double *azUOut)
{
    const Radar *s = &L->s;
    const double t = L->nextT, azU = L->nextAzU;
    pr->pulse = L->nextPulse++;
    pr->scan = (int)floor((azU - s->az0) / 360) + 1;
    if (s->nScanPri > 0) {
        const int j = ((pr->scan - 1) % s->nScanPri + s->nScanPri) % s->nScanPri;
        pr->pri = s->scanPri[j];
        pr->staggerIdx = j + 1;
    } else {
        const int st = L->stPos++ % s->nStagger;
        pr->pri = s->pri[st];
        pr->staggerIdx = st + 1;
    }
    pr->t = t;
    pr->az = mod360(azU);
    *azUOut = azU;
    *prevPri = L->lastPri > 0 ? L->lastPri : pr->pri;
    L->lastPri = pr->pri;
    L->nextT = t + pr->pri;
    L->nextAzU = azU + s->rate * pr->pri;
}

/* Evaluates radar_config.m with the settings (mode, set) and sets everything up. On an error the
 * message is printed and the old setup stays (returns -1). */
static int live_build(Live *L, int mode, const ModeSettings *set, int first)
{
    char st[8][160], err[1024];
    const char *sets[48];
    int n = 0;
    for (int i = 0; i < L->nBase && n < 32; i++) sets[n++] = L->baseSets[i];
    if (L->haveBi && set) {
        const int ns = board_config_statements(&L->bi, mode, set, st, 8);
        for (int i = 0; i < ns; i++) sets[n++] = st[i];
    }
    MValue *cfg = mfile_eval(L->mPath, sets, n, err, sizeof err);
    if (!cfg || cfg->type != MV_STRUCT) {
        fprintf(stderr, "radar_live_sim: %s\n", cfg ? "the configuration is not a struct" : err);
        if (cfg) mvalue_free(cfg);
        return -1;
    }
    const MValue *oldCfg = g_cfg;
    g_cfg = cfg;
    Radar s;
    Clutter clut;
    memset(&clut, 0, sizeof clut);
    int haveClut = 0;
    if (setjmp(g_jb) != 0) {
        fprintf(stderr, "radar_live_sim: %s%s\n", g_err, first ? "" : " (the simulation keeps its settings)");
        g_cfg = oldCfg;
        mvalue_free(cfg);
        return -1;
    }
    setup(&s);
    Rng rng;
    rng_seed(&rng, L->seed, 0xFFFFFFFFFFFFULL);
    if (first) {                                     /* the targets stay the same for every mode */
        L->seed = (uint64_t)(int64_t)num("sim.seed");
        rng_seed(&rng, L->seed, 0xFFFFFFFFFFFFULL);
        L->nTg = L->targetsCsv ? read_targets_csv(L->targetsCsv, &L->tg) : scenario(str("scenario.type"), &rng, &L->tg);
        const int swerling = (int)num("sim.swerling");
        for (int i = 0; i < L->nTg; i++) {
            L->tg[i].nScan = 10000;
            L->tg[i].rcsScan = xcalloc((size_t)L->tg[i].nScan, sizeof(double));
            for (int k = 0; k < L->tg[i].nScan; k++)
                L->tg[i].rcsScan[k] = swerling == 1 ? L->tg[i].rcs * -log(rng_uniform(&rng)) : L->tg[i].rcs;
        }
    }
    if (!L->fast && num("clutter.enable") != 0) {
        const WinSpec rangeWin = window_spec("proc.", "rangeWindow", "taylor");
        clutter_build(&clut, &s, &rng, &rangeWin);
        haveClut = 1;
    }
    double *sect = NULL;
    int nSect = 0, sc = 0;
    if (has("sim.sectors_deg") && !empty("sim.sectors_deg")) {
        const MValue *v = matrix("sim.sectors_deg", &nSect, &sc);
        if (sc != 2) fail("cfg.sim.sectors_deg must have rows [azStart azEnd]");
        sect = xcalloc((size_t)2 * nSect, sizeof(double));
        memcpy(sect, v->re, sizeof(double) * 2 * (size_t)nSect);
    }
    FastCal cal;
    memset(&cal, 0, sizeof cal);
    if (L->fast) {
        Work wk;
        if (work_alloc(&wk, &s, NULL) != 0) fail("out of memory");
        fast_calibrate(&s, &wk, &cal);
        work_free(&wk);
    }
    /* the new setup replaces the old one */
    if (!first) {
        radar_free(&L->s);
        if (L->haveClut) clutter_free(&L->clut);
        mvalue_free(L->cfg);
        free(L->sect);
    }
    L->cfg = cfg;
    L->s = s;
    L->clut = clut;
    L->haveClut = haveClut;
    L->sect = sect;
    L->nSect = nSect;
    L->cal = cal;
    L->mode = mode;
    if (set) L->set = *set;
    if (!g_quiet) {
        printf("Mode %d: %d element(s), DDC %d samples at %g MHz, %s PRF (%d PRI(s)), %g rpm, carriers %.0f / %.0f MHz%s%s\n",
               mode, s.nCh, s.nDdc, s.fsOut / 1e6, s.nScanPri ? "scan-to-scan" : s.nStagger > 1 ? "staggered" : "fixed",
               s.nScanPri ? s.nScanPri : s.nStagger, s.rate * 60 / 360, s.w[0].fc / 1e6, s.w[1].fc / 1e6,
               haveClut ? ", terrain clutter" : "", L->fast ? ", fidelity fast" : "");
        fflush(stdout);
    }
    return 0;
}

#ifdef _WIN32
static DWORD WINAPI producer_main(LPVOID arg)
#else
static void *producer_main(void *arg)
#endif
{
    Live *L = arg;
    Work wk;
    memset(&wk, 0, sizeof wk);
    uint64_t wkGen = (uint64_t)-1;
    mutex_lock(&L->mu);
    for (;;) {
        while (!L->stop && (L->paused || L->nextK - L->sendK >= LIVE_RING)) cond_wait_ms(&L->cv, &L->mu, 50);
        if (L->stop) break;
        const int k = L->nextK++;
        Slot *sl = &L->ring[k % LIVE_RING];
        Pri pr;
        double prevPri, azU;
        sched_next(L, &pr, &prevPri, &azU);
        const uint64_t gen = L->gen;
        sl->ready = 0;
        L->busy++;
        mutex_unlock(&L->mu);

        const Radar *s = &L->s;
        const int active = in_sectors(L, pr.az);
        if (active) {
            if (wkGen != gen) {                       /* the buffers of the setup in use */
                work_free(&wk);
                if (work_alloc(&wk, s, L->haveClut ? &L->clut : NULL) != 0) {
                    fprintf(stderr, "radar_live_sim: out of memory\n");
                    exit(1);
                }
                wkGen = gen;
            }
            if (L->fast) simulate_fast(s, L->tg, L->nTg, &pr, L->seed, &wk, &L->cal);
            else simulate_pri(s, L->haveClut ? &L->clut : NULL, L->tg, L->nTg, &pr, L->seed, &wk);
            const size_t need = 2 * (size_t)s->nDdc * s->nCh;
            if (need > sl->cap) {
                int16_t *q = realloc(sl->iq, need * sizeof *q);
                if (!q) { fprintf(stderr, "radar_live_sim: out of memory\n"); exit(1); }
                sl->iq = q;
                sl->cap = need;
            }
            for (size_t q = 0; q < (size_t)s->nDdc * s->nCh; q++) { sl->iq[2 * q] = wk.I[q]; sl->iq[2 * q + 1] = wk.Q[q]; }
            sl->nDdc = s->nDdc;
            sl->nCh = s->nCh;
        }
        mutex_lock(&L->mu);
        L->busy--;
        if (gen == L->gen) {
            sl->pr = pr;
            sl->prevPri = prevPri;
            sl->azU = azU;
            sl->k = k;
            sl->active = active;
            sl->gen = gen;
            sl->ready = 1;
        }
        cond_broadcast(&L->cv);
    }
    mutex_unlock(&L->mu);
    work_free(&wk);
    return 0;
}

/* Stops the producers, sets up the new settings and goes on from the PRI after the last one sent. */
static void live_reconfigure(Live *L, int mode, const ModeSettings *set, double contT, double contAzU, int contPulse)
{
    mutex_lock(&L->mu);
    L->paused = 1;
    L->gen++;
    while (L->busy > 0) cond_wait_ms(&L->cv, &L->mu, 20);
    mutex_unlock(&L->mu);
    const double t0 = now_seconds();
    const int ok = live_build(L, mode, set, 0) == 0;
    mutex_lock(&L->mu);
    L->nextPulse = contPulse;
    L->nextT = contT;
    L->nextAzU = contAzU;
    L->stPos = 0;
    L->lastPri = 0;
    L->nextK = L->sendK;
    for (int i = 0; i < LIVE_RING; i++) L->ring[i].ready = 0;
    L->paused = 0;
    cond_broadcast(&L->cv);
    mutex_unlock(&L->mu);
    if (ok && !g_quiet) printf("  new settings in %.2f s, from pulse %d (t = %.3f s, az %.2f deg)\n", now_seconds() - t0,
                               contPulse, contT, mod360(contAzU));
}

/* The settings the board registers ask for (invalid ones: those of the mode file). */
static void settings_of_registers(const Live *L, const BoardSim *b, int *mode, ModeSettings *set)
{
    *mode = board_sim_register(b, 0x00);
    const ModeInfo *m = board_mode_info(&L->bi, *mode);
    if (!m) {
        printf("  radar mode %d is not available: mode %d stays\n", *mode, L->mode);
        *mode = L->mode;
        m = board_mode_info(&L->bi, *mode);
    }
    ModeSettings s = m->defaults, t;
    char err[256];
    s.rpm = L->set.rpm;
    t = s;  t.rpm = L->set.rpm;
    if (!board_settings_valid(&L->bi, *mode, &t, err, sizeof err)) s.rpm = m->defaults.rpm;
    struct { int addr; int *field; int base; } it[4] = {
        { 0x01, &s.code, 0 }, { 0x02, &s.f1_MHz, L->bi.fcMin_MHz }, { 0x03, &s.f2_MHz, L->bi.fcMin_MHz }, { 0x05, &s.stagger, 0 } };
    for (int i = 0; i < 4; i++) {
        const int old = *it[i].field;
        if (old == 0) continue;                         /* not used by the mode */
        *it[i].field = board_sim_register(b, it[i].addr) + it[i].base;
        if (!board_settings_valid(&L->bi, *mode, &s, err, sizeof err)) {
            printf("  %s: %s (the value of the mode file is used)\n", board_register_name(it[i].addr), err);
            *it[i].field = old;
        }
    }
    if (L->bi.sendRpm) {
        t = s;
        t.rpm = board_sim_register(b, 0x32);
        if (board_settings_valid(&L->bi, *mode, &t, err, sizeof err)) s.rpm = t.rpm;
        else printf("  RPM: %s (%g rpm stays)\n", err, s.rpm);
    }
    *set = s;
}

static void live_registers(const Live *L, uint8_t *regs)
{
    memset(regs, 0, 256);
    const ModeInfo *m = board_mode_info(&L->bi, L->mode);
    regs[0x00] = (uint8_t)L->mode;
    regs[0x01] = (uint8_t)L->set.code;
    regs[0x02] = (uint8_t)(m && m->enable[0] ? L->set.f1_MHz - L->bi.fcMin_MHz : 0);
    regs[0x03] = (uint8_t)(m && m->enable[1] ? L->set.f2_MHz - L->bi.fcMin_MHz : 0);
    regs[0x05] = (uint8_t)L->set.stagger;
    regs[0x28] = (uint8_t)L->bi.operationMode;
    for (int i = 0; i < 4; i++) { regs[0x29 + i] = (uint8_t)L->set.stc[i]; regs[0x2D + i] = (uint8_t)L->set.gain[i]; }
    regs[0x32] = (uint8_t)lround(L->set.rpm);
}

typedef struct {
    const char *host, *format, *fidelity, *boardReply;
    int port, payload, boardPort, replyPort;
    char hostBuf[128], replyHost[128];
    int vbdr, fast;
} LiveOpts;

/* The defaults of the options from the configuration (cfg.udp, cfg.board, cfg.sim.liveFidelity). */
static int live_defaults(const char *mPath, const char *const *sets, int nSets, LiveOpts *o)
{
    char err[1024];
    MValue *c0 = mfile_eval(mPath, sets, nSets, err, sizeof err);
    if (!c0 || c0->type != MV_STRUCT) {
        fprintf(stderr, "radar_live_sim: %s\n", c0 ? "the configuration is not a struct" : err);
        if (c0) mvalue_free(c0);
        return -1;
    }
    g_cfg = c0;
    if (setjmp(g_jb) != 0) {
        fprintf(stderr, "radar_live_sim: %s\n", g_err);
        mvalue_free(c0);
        g_cfg = NULL;
        return -1;
    }
    const char *format = o->format ? o->format : str_or("udp.format", "rddc");
    const char *fidelity = o->fidelity ? o->fidelity : str_or("sim.liveFidelity", "full");
    if (!same_text(format, "rddc") && !same_text(format, "vbdr")) fail("--format must be rddc or vbdr");
    if (!same_text(fidelity, "full") && !same_text(fidelity, "fast")) fail("--fidelity must be full or fast");
    o->vbdr = same_text(format, "vbdr");
    o->fast = same_text(fidelity, "fast");
    snprintf(o->hostBuf, sizeof o->hostBuf, "%s", o->host ? o->host : str_or("udp.host", "127.0.0.1"));
    if (!o->port) o->port = (int)num_or("udp.port", 50000);
    if (!o->payload) o->payload = (int)(o->vbdr ? num_or("udp.vbdr.payloadBytes", 1400) : num_or("udp.payloadBytes", 1400));
    if (!o->boardPort) o->boardPort = (int)num_or("board.port", 5000);
    o->replyPort = (int)num_or("board.replyPort", 5100);
    snprintf(o->replyHost, sizeof o->replyHost, "%s", str_or("board.pcHost", "127.0.0.1"));
    if (o->boardReply) {
        snprintf(o->replyHost, sizeof o->replyHost, "%s", o->boardReply);
        char *colon = strrchr(o->replyHost, ':');
        if (colon) { *colon = 0; o->replyPort = atoi(colon + 1); }
    }
    mvalue_free(c0);
    g_cfg = NULL;
    return 0;
}

static void live_usage(void)
{
    puts("usage: radar_live_sim [--config-m FILE] [--set STATEMENT]... [--host H] [--port N] [--format rddc|vbdr]\n"
         "                      [--payload BYTES] [--threads N] [--targets FILE.csv] [--fidelity full|fast]\n"
         "                      [--board-port N] [--board-reply HOST:PORT] [--no-board] [--max-pulses N]\n"
         "                      [--duration S] [--speed X] [--quiet]\n"
         "  Simulates the PRIs of radar_config.m (as radar_generate_ddc) and sends every one to the receiver at\n"
         "  its PRI time, until Ctrl+C; answers the board commands of the receiver and changes the mode and its\n"
         "  settings with them (README \"radar_live_sim\").");
}

int main(int argc, char **argv)
{
    const char *mPath = RADAR_SOURCE_DIR "/../radar_config.m";
    const char **sets = calloc((size_t)argc + 1, sizeof(char *));
    const char *host = NULL, *format = NULL, *fidelity = NULL, *boardReply = NULL, *targetsCsv = NULL;
    int port = 0, payload = 0, nThreads = 0, boardPort = 0, noBoard = 0, nSets = 0;
    long maxPulses = 0;
    double duration = 0, speed = 1;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i], *v = i + 1 < argc ? argv[i + 1] : NULL;
        if (!strcmp(a, "--config-m") && v) { mPath = v; i++; }
        else if (!strcmp(a, "--set") && v) { sets[nSets++] = v; i++; }
        else if (!strcmp(a, "--host") && v) { host = v; i++; }
        else if (!strcmp(a, "--port") && v) { port = atoi(v); i++; }
        else if (!strcmp(a, "--format") && v) { format = v; i++; }
        else if (!strcmp(a, "--payload") && v) { payload = atoi(v); i++; }
        else if (!strcmp(a, "--threads") && v) { nThreads = atoi(v); i++; }
        else if (!strcmp(a, "--targets") && v) { targetsCsv = v; i++; }
        else if (!strcmp(a, "--fidelity") && v) { fidelity = v; i++; }
        else if (!strcmp(a, "--board-port") && v) { boardPort = atoi(v); i++; }
        else if (!strcmp(a, "--board-reply") && v) { boardReply = v; i++; }
        else if (!strcmp(a, "--no-board")) noBoard = 1;
        else if (!strcmp(a, "--max-pulses") && v) { maxPulses = atol(v); i++; }
        else if (!strcmp(a, "--duration") && v) { duration = atof(v); i++; }
        else if (!strcmp(a, "--speed") && v) { speed = atof(v); i++; }
        else if (!strcmp(a, "--quiet")) g_quiet = 1;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { live_usage(); free(sets); return 0; }
        else { fprintf(stderr, "unknown option %s\n", a); live_usage(); free(sets); return 2; }
    }
    static Live L;
    L.mPath = mPath;
    L.baseSets = sets;
    L.nBase = nSets;
    L.targetsCsv = targetsCsv;
    char err[1024];
    /* the settings of the modes (radar_settings.c) and the common ones of the configuration */
    if (board_info_load(mPath, sets, nSets, &L.bi, err, sizeof err) != 0) {
        fprintf(stderr, "radar_live_sim: %s\n", err);
        return 1;
    }
    L.haveBi = 1;
    LiveOpts o = { host, format, fidelity, boardReply, port, payload, boardPort, 0, "", "", 0, 0 };
    if (live_defaults(mPath, sets, nSets, &o) != 0) return 2;
    if (!(speed > 0)) { fprintf(stderr, "--speed must be positive\n"); return 2; }
    const int vbdr = o.vbdr;
    L.fast = o.fast;
    port = o.port;
    payload = o.payload;
    boardPort = o.boardPort;
    const int replyPort = o.replyPort;
    char hostBuf[128], replyHost[128];
    snprintf(hostBuf, sizeof hostBuf, "%s", o.hostBuf);
    snprintf(replyHost, sizeof replyHost, "%s", o.replyHost);
    (void)&warn;
    (void)&now_s;

    /* the settings of the active mode, the setup */
    L.mode = L.bi.activeMode;
    const ModeInfo *m0 = board_mode_info(&L.bi, L.mode);
    if (!m0) { fprintf(stderr, "radar_live_sim: mode %d is not in cfg.modes.available\n", L.mode); return 1; }
    L.set = m0->defaults;
    if (!g_quiet) printf("Reading %s\n", mPath);
    if (live_build(&L, L.mode, NULL, 1) != 0) return 1;
    L.nextPulse = 1;
    L.nextAzU = L.s.az0;
#ifdef _OPENMP
    L.nThreads = nThreads > 0 ? nThreads : omp_get_num_procs();
#else
    L.nThreads = nThreads > 0 ? nThreads : 4;
#endif
    if (L.nThreads < 1) L.nThreads = 1;

#ifdef _WIN32
    WSADATA wsa;
    WSAStartup(MAKEWORD(2, 2), &wsa);
#endif
    lsock_t sock = socket(AF_INET, SOCK_DGRAM, 0);
    int sndbuf = 32 << 20;
    setsockopt(sock, SOL_SOCKET, SO_SNDBUF, (const char *)&sndbuf, sizeof sndbuf);
    struct sockaddr_in to;
    memset(&to, 0, sizeof to);
    to.sin_family = AF_INET;
    to.sin_port = htons((unsigned short)port);
    if (inet_pton(AF_INET, hostBuf, &to.sin_addr) != 1) { fprintf(stderr, "invalid host address %s\n", hostBuf); return 2; }
    BoardSim *board = NULL;
    if (!noBoard) {
        uint8_t regs[256];
        live_registers(&L, regs);
        board = board_sim_open(boardPort, replyHost, replyPort, regs, 1, err, sizeof err);
        if (!board) fprintf(stderr, "warning: the board commands are not answered: %s\n", err);
    }
    if (!g_quiet) {
        printf("Sending to %s:%d (%s, %d bytes per packet), %d producer thread(s), until %s\n", hostBuf, port,
               vbdr ? "VBDR" : "RDDC", payload, L.nThreads,
               maxPulses > 0 ? "the PRI limit" : duration > 0 ? "the duration" : "Ctrl+C");
        fflush(stdout);
    }

    mutex_init(&L.mu);
    cond_init(&L.cv);
#ifdef _WIN32
    HANDLE *thr = calloc((size_t)L.nThreads, sizeof(HANDLE));
    for (int i = 0; i < L.nThreads; i++) thr[i] = CreateThread(NULL, 0, producer_main, &L, 0, NULL);
#else
    pthread_t *thr = calloc((size_t)L.nThreads, sizeof(pthread_t));
    for (int i = 0; i < L.nThreads; i++) pthread_create(&thr[i], NULL, producer_main, &L);
#endif
    signal(SIGINT, live_on_signal);
    signal(SIGTERM, live_on_signal);

    /* pre-fill: half of the ring before the clock starts */
    mutex_lock(&L.mu);
    const double tFill = now_seconds();
    while (!g_liveStop && now_seconds() - tFill < 30) {
        int n = 0;
        for (int i = 0; i < LIVE_RING; i++) n += L.ring[i].ready;
        if (n >= LIVE_RING / 2) break;
        cond_wait_ms(&L.cv, &L.mu, 50);
    }
    mutex_unlock(&L.mu);

    size_t bufCap = 0;
    unsigned char *buf = NULL;
    uint32_t lenBuf[65536];
    uint64_t nSent = 0, nPkt = 0, sendErrors = 0, nLate = 0, nPri = 0;
    double wall0 = now_seconds() + 0.01, simT0 = -1, lateSum = 0, lastPrint = now_seconds(), lastWrite = -1;
    double lastT = 0, lastAzU = L.s.az0, lastPri = 0;
    int lastPulse = 0, dirty = 0;
    uint32_t written[8];
    while (!g_liveStop) {
        /* the board commands; a change of the simulation when no further command comes */
        if (board) {
            memset(written, 0, sizeof written);
            if (board_sim_poll(board, 0, written) > 0) {
                const int addrs[6] = { 0x00, 0x01, 0x02, 0x03, 0x05, 0x32 };
                for (int i = 0; i < 6; i++)
                    if ((written[addrs[i] / 32] >> (addrs[i] % 32) & 1) && (addrs[i] != 0x32 || L.bi.sendRpm)) {
                        dirty = 1;
                        lastWrite = now_seconds();
                    }
            }
            if (dirty && now_seconds() - lastWrite > LIVE_DEBOUNCE) {
                int mode;
                ModeSettings set;
                settings_of_registers(&L, board, &mode, &set);
                dirty = 0;
                if (mode != L.mode || set.code != L.set.code || set.stagger != L.set.stagger || set.rpm != L.set.rpm ||
                    set.f1_MHz != L.set.f1_MHz || set.f2_MHz != L.set.f2_MHz) {
                    printf("Board: mode %d, code %d, stagger %d, %g rpm, F1 %d MHz, F2 %d MHz: new settings of the simulation\n",
                           mode, set.code, set.stagger, set.rpm, set.f1_MHz, set.f2_MHz);
                    fflush(stdout);
                    const double azNext = lastPulse ? lastAzU + L.s.rate * lastPri : L.s.az0;
                    const double tNext = lastPulse ? lastT + lastPri : 0;
                    live_reconfigure(&L, mode, &set, tNext, azNext, lastPulse + 1);
                }
            }
        }
        /* the next PRI */
        mutex_lock(&L.mu);
        Slot *sl = &L.ring[L.sendK % LIVE_RING];
        if (!(sl->ready && sl->k == L.sendK && sl->gen == L.gen)) {
            cond_wait_ms(&L.cv, &L.mu, board ? 5 : 50);
            mutex_unlock(&L.mu);
            continue;
        }
        mutex_unlock(&L.mu);
        const Pri pr = sl->pr;
        if (simT0 < 0) simT0 = pr.t;
        if (duration > 0 && pr.t - simT0 >= duration) break;
        if (sl->active) {
            double tSend = wall0 + (pr.t - simT0) / speed;
            /* the idle time: the board commands */
            double left;
            while (board && (left = tSend - now_seconds()) > 5e-3 && !g_liveStop)
                board_sim_poll(board, (int)((left - 4e-3) * 1e3 < 20 ? (left - 4e-3) * 1e3 : 20), NULL);
            live_wait_until(tSend);
            const double late = now_seconds() - tSend;
            if (late > 1e-3) {                       /* the simulated time runs slower from here on */
                wall0 += late;
                lateSum += late;
                nLate++;
            }
            DdcPri d;
            memset(&d, 0, sizeof d);
            d.pulseNumber = (uint32_t)pr.pulse;
            d.scanNumber = (uint32_t)pr.scan;
            d.az_deg = (float)pr.az;
            d.staggerIndex = (uint16_t)pr.staggerIdx;
            d.time_s = pr.t;
            d.pri_s = pr.pri;
            d.prevPri_s = sl->prevPri;
            d.ddcFs_Hz = (float)L.s.fsOut;
            d.ddcT0_s = (float)L.s.ddcT0;
            d.shortIdx[0] = (uint32_t)L.s.w[0].rxFirst;  d.shortIdx[1] = (uint32_t)L.s.w[0].rxLast;
            d.longIdx[0] = (uint32_t)L.s.w[1].rxFirst;   d.longIdx[1] = (uint32_t)L.s.w[1].rxLast;
            d.nSamples = sl->nDdc;
            d.nCh = sl->nCh;
            d.iq = sl->iq;
            size_t bytes;
            const int n = vbdr ? ddc_vbdr_count(d.nSamples, payload, &bytes) : ddc_datagram_count(&d, payload, &bytes);
            if (n > 65535) { fprintf(stderr, "radar_live_sim: %d packets per PRI: use a larger --payload\n", n); break; }
            if (bytes > bufCap) {
                unsigned char *q = realloc(buf, bytes);
                if (!q) { fprintf(stderr, "radar_live_sim: out of memory\n"); break; }
                buf = q;
                bufCap = bytes;
            }
            const int np = vbdr ? ddc_make_vbdr(&d, payload, (uint32_t)sendErrors, buf, lenBuf)
                                : ddc_make_datagrams(&d, payload, buf, lenBuf);
            size_t off = 0;
            for (int j = 0; j < np; off += lenBuf[j++]) {
                if (sendto(sock, (const char *)buf + off, (int)lenBuf[j], 0, (struct sockaddr *)&to, sizeof to) < 0) sendErrors++;
                nPkt++;
            }
            nSent++;
        }
        lastT = pr.t;
        lastAzU = sl->azU;
        lastPri = pr.pri;
        lastPulse = pr.pulse;
        nPri++;
        mutex_lock(&L.mu);
        sl->ready = 0;
        L.sendK++;
        cond_broadcast(&L.cv);
        mutex_unlock(&L.mu);
        if (maxPulses > 0 && (long)nPri >= maxPulses) break;
        if (!g_quiet && now_seconds() - lastPrint > 5) {
            lastPrint = now_seconds();
            printf("t %.1f s, az %6.2f deg, pulse %d: %llu PRIs sent (%llu packets)%s", pr.t, pr.az, pr.pulse,
                   (unsigned long long)nSent, (unsigned long long)nPkt, sendErrors ? ", send errors" : "");
            if (nLate) printf("; slower than real time: %llu PRIs late, %.2f s behind", (unsigned long long)nLate, lateSum);
            printf("\n");
            fflush(stdout);
        }
    }
    mutex_lock(&L.mu);
    L.stop = 1;
    cond_broadcast(&L.cv);
    mutex_unlock(&L.mu);
    for (int i = 0; i < L.nThreads; i++) {
#ifdef _WIN32
        WaitForSingleObject(thr[i], INFINITE);
        CloseHandle(thr[i]);
#else
        pthread_join(thr[i], NULL);
#endif
    }
    printf("radar_live_sim: %llu PRIs sent (%llu packets, %llu send errors) of %llu PRIs; %llu late by more than 1 ms "
           "(%.2f s behind the real time)\n", (unsigned long long)nSent, (unsigned long long)nPkt,
           (unsigned long long)sendErrors, (unsigned long long)nPri, (unsigned long long)nLate, lateSum);
    board_sim_close(board);
    free(buf);
    free(thr);
    for (int i = 0; i < LIVE_RING; i++) free(L.ring[i].iq);
    free(sets);
    return 0;
}

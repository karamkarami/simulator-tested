/*
 * main.c - radar receiver: DDC output over UDP -> processing -> plots.
 *
 *   radar_receiver [options]
 *     --config FILE          receiver configuration from export_receiver_config.m
 *                            (default config/radar_config.bin)
 *     --port N               UDP port of the DDC stream (default 50000)
 *     --record FILE          also store the received datagrams in FILE
 *     --replay FILE          read the datagrams from FILE instead of the network
 *     --ppi                  show the PPI (Qt build only)
 *     --ppi-snapshot FILE    save the PPI as PNG at the end and close it
 *     --ppi-show MODE        PPI shows 'plots' (default, after the plot extractor),
 *                            'detections' (CFAR detections before it) or 'both'
 *     --serial DEVICE        send the CFAR detections to a serial port (e.g. /dev/ttyUSB0, COM3)
 *     --baud N               serial baud rate (default 921600)
 *     --serial-content X     'detections' (default) or 'plots'
 *     --threads N            processing threads (default: all cores but one)
 *     --mf auto|direct|fft   matched filter: direct or FFT fast convolution (default auto: the faster)
 *     --update-every N       run Doppler/CFAR every N pulses (default from the configuration)
 *     --plots-csv FILE       plot list (default plots.csv, '' = none)
 *     --detections-csv FILE  list of all pruned CFAR detections (default: none)
 *     --tracks-csv FILE      every track update and deletion (default: none; tracking is
 *                            switched on in the configuration, cfg.track.enable)
 *     --idle-timeout S       stop when no data came for S seconds after the first PRI (default 0 = never)
 *     --queue N              PRIs buffered between the input thread and the processing (default 64)
 *     --no-skip              never skip Doppler/CFAR updates when the processing falls behind
 *     --cfar-threshold DB    CA-CFAR threshold factor in dB (default from the configuration,
 *                            i.e. from cfg.cfar.pfa); can also be changed in the PPI
 *     --quiet                do not print every plot
 */
#include "radar_rx.h"

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *configPath, *recordPath, *replayPath, *snapshotPath;
    const char *serialDevice, *serialContent, *plotsCsv, *detectionsCsv, *tracksCsv, *mf, *ppiShow;
    int port, baud, threads, updateEvery, ppi, quiet, queue, noSkip;
    double idleTimeout_s;
    double cfarThreshold_dB;      /* NAN: from the configuration */
} Options;

typedef struct {
    Options opt;
    RadarConfig cfg;
    SerialOut *serial;
    FILE *plotsFile, *detFile, *trackFile;
    uint64_t nPlots, nDetections, nTrackUpdates, nConfirmed;
    double detArc_deg;            /* azimuth step between CPI updates (PPI detection arcs) */
    double console_ms, display_ms;   /* time in the console output and in the PPI calls */
    int exitCode;
} App;

static volatile sig_atomic_t g_stop = 0;

static void on_signal(int sig)
{
    (void)sig;
    g_stop = 1;
}

static void usage(void)
{
    puts("usage: radar_receiver [--config FILE] [--port N] [--record FILE] [--replay FILE]\n"
         "                      [--ppi] [--ppi-snapshot FILE] [--ppi-show plots|detections|both]\n"
         "                      [--serial DEVICE] [--baud N]\n"
         "                      [--serial-content detections|plots] [--threads N] [--mf auto|direct|fft]\n"
         "                      [--update-every N]\n"
         "                      [--plots-csv FILE] [--detections-csv FILE] [--tracks-csv FILE]\n"
         "                      [--idle-timeout S] [--quiet]\n"
         "                      [--queue N] [--no-skip] [--cfar-threshold DB]");
}

static int parse_options(int argc, char **argv, Options *o)
{
    memset(o, 0, sizeof(*o));
    o->configPath = "config/radar_config.bin";
    o->port = 50000;
    o->baud = 921600;
    o->serialContent = "detections";
    o->plotsCsv = "plots.csv";
    o->mf = "auto";
    o->ppiShow = "plots";
    o->queue = 64;
    o->cfarThreshold_dB = NAN;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
#define TAKE(field) do { if (!v) { fprintf(stderr, "%s needs a value\n", a); return -1; } field; i++; } while (0)
        if      (!strcmp(a, "--config"))         TAKE(o->configPath = v);
        else if (!strcmp(a, "--port"))           TAKE(o->port = atoi(v));
        else if (!strcmp(a, "--record"))         TAKE(o->recordPath = v);
        else if (!strcmp(a, "--replay"))         TAKE(o->replayPath = v);
        else if (!strcmp(a, "--ppi"))            o->ppi = 1;
        else if (!strcmp(a, "--ppi-snapshot"))   TAKE((o->snapshotPath = v, o->ppi = 1));
        else if (!strcmp(a, "--ppi-show"))       TAKE(o->ppiShow = v);
        else if (!strcmp(a, "--serial"))         TAKE(o->serialDevice = v);
        else if (!strcmp(a, "--baud"))           TAKE(o->baud = atoi(v));
        else if (!strcmp(a, "--serial-content")) TAKE(o->serialContent = v);
        else if (!strcmp(a, "--threads"))        TAKE(o->threads = atoi(v));
        else if (!strcmp(a, "--mf"))             TAKE(o->mf = v);
        else if (!strcmp(a, "--update-every"))   TAKE(o->updateEvery = atoi(v));
        else if (!strcmp(a, "--plots-csv"))      TAKE(o->plotsCsv = v);
        else if (!strcmp(a, "--detections-csv")) TAKE(o->detectionsCsv = v);
        else if (!strcmp(a, "--tracks-csv"))     TAKE(o->tracksCsv = v);
        else if (!strcmp(a, "--idle-timeout"))   TAKE(o->idleTimeout_s = atof(v));
        else if (!strcmp(a, "--quiet"))          o->quiet = 1;
        else if (!strcmp(a, "--queue"))          TAKE(o->queue = atoi(v));
        else if (!strcmp(a, "--no-skip"))        o->noSkip = 1;
        else if (!strcmp(a, "--cfar-threshold")) TAKE(o->cfarThreshold_dB = atof(v));
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); return 1; }
        else { fprintf(stderr, "unknown option %s\n", a); usage(); return -1; }
#undef TAKE
    }
    if (strcmp(o->mf, "auto") && strcmp(o->mf, "direct") && strcmp(o->mf, "fft")) {
        fprintf(stderr, "--mf must be 'auto', 'direct' or 'fft'\n");
        return -1;
    }
    if (!isnan(o->cfarThreshold_dB) && !(o->cfarThreshold_dB >= 0.0 && o->cfarThreshold_dB <= 50.0)) {
        fprintf(stderr, "--cfar-threshold must be 0 ... 50 dB\n");
        return -1;
    }
    if (strcmp(o->ppiShow, "plots") && strcmp(o->ppiShow, "detections") && strcmp(o->ppiShow, "both")) {
        fprintf(stderr, "--ppi-show must be 'plots', 'detections' or 'both'\n");
        return -1;
    }
    if (strcmp(o->serialContent, "detections") && strcmp(o->serialContent, "plots")) {
        fprintf(stderr, "--serial-content must be 'detections' or 'plots'\n");
        return -1;
    }
    return 0;
}

/* ------------------------------------------------------------------------ */
/* Callbacks of the processor                                               */
/* ------------------------------------------------------------------------ */

static const char *pulse_name(int pulse) { return pulse == PULSE_SHORT ? "short" : "long"; }

static double ms_since(double t0) { return (now_seconds() - t0) * 1e3; }

static void on_detections(void *user, int pulse, const Detection *dets, int n)
{
    App *app = user;
    app->nDetections += (uint64_t)n;
    /* called on every CPI update, also without detections: no frame is sent
     * then, but the bytes still queued are written to the port */
    if (app->serial) {
        if (n > 0 && !strcmp(app->opt.serialContent, "detections"))
            serial_send_detections(app->serial, pulse, dets[0].pulseNumber, dets[0].az_deg, dets, n);
        else
            serial_send_detections(app->serial, pulse, 0, 0.0, NULL, 0);
    }
    if (app->detFile) {
        for (int i = 0; i < n; i++) {
            const Detection *d = &dets[i];
            fprintf(app->detFile, "%u,%s,%.6f,%.3f,%.1f,%.3f,%.3f,%.2f,%d,%d,%d\n", d->pulseNumber, pulse_name(pulse),
                    d->time_s, d->az_deg, d->range_m, d->el_deg, d->vr_mps, d->snr_dB, d->beam, d->rIdx, d->dIdx);
        }
    }
#ifdef RADAR_WITH_QT
    if (app->opt.ppi && n > 0) {
        const double t0 = now_seconds();
        ppi_post_detections(dets, n, app->detArc_deg);
        app->display_ms += ms_since(t0);
    }
#endif
}

static void on_plot(void *user, const Detection *p)
{
    App *app = user;
    app->nPlots++;
    if (!app->opt.quiet) {
        const double t0 = now_seconds();
        printf("plot  t %8.3f s  az %7.2f deg  R %8.3f km  el %6.2f deg  vr %7.2f m/s  SNR %5.1f dB  %-5s  hits %d\n",
               p->time_s, p->az_deg, p->range_m / 1e3, p->el_deg, p->vr_mps, p->snr_dB, pulse_name(p->pulse), p->hits);
        app->console_ms += ms_since(t0);
    }
    if (app->plotsFile) {
        fprintf(app->plotsFile, "%.6f,%.4f,%.2f,%.4f,%.4f,%.2f,%s,%d,%d,%d,%d,%u\n", p->time_s, p->az_deg, p->range_m,
                p->el_deg, p->vr_mps, p->snr_dB, pulse_name(p->pulse), p->hits, p->beam, p->rIdx, p->dIdx, p->pulseNumber);
        fflush(app->plotsFile);
    }
    if (app->serial && !strcmp(app->opt.serialContent, "plots")) serial_send_plot(app->serial, p);
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        const double t0 = now_seconds();
        ppi_post_plot(p);
        app->display_ms += ms_since(t0);
    }
#endif
}

static void on_track(void *user, const Track *t)
{
    App *app = user;
    app->nTrackUpdates += !t->deleted;
    if (t->confirmed && !t->deleted && t->hits == app->cfg.trackConfirmHits) {
        app->nConfirmed++;
        if (!app->opt.quiet) {
            const double t0 = now_seconds();
            printf("track %u confirmed  t %8.3f s  R %8.3f km  az %7.2f deg  speed %6.1f m/s\n", t->id, t->time_s,
                   hypot(hypot(t->pos_m[0], t->pos_m[1]), t->pos_m[2]) / 1e3,
                   fmod(atan2(t->pos_m[0], t->pos_m[1]) * 180.0 / 3.14159265358979323846 + 360.0, 360.0),
                   hypot(t->vel_mps[0], t->vel_mps[1]));
            app->console_ms += ms_since(t0);
        }
    }
    if (app->trackFile)
        fprintf(app->trackFile, "%u,%d,%d,%.6f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%d,%d,%.4f,%.4f,%.4f\n", t->id, t->confirmed,
                t->deleted, t->time_s, t->pos_m[0], t->pos_m[1], t->pos_m[2], t->vel_mps[0], t->vel_mps[1], t->vel_mps[2],
                t->hits, t->misses, t->acc_mps2[0], t->acc_mps2[1], t->acc_mps2[2]);
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        const double t0 = now_seconds();
        ppi_post_track(t);
        app->display_ms += ms_since(t0);
    }
#endif
}

static void on_antenna(void *user, double az_deg, double time_s, uint32_t scan)
{
    App *app = user;
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        const double t0 = now_seconds();
        ppi_post_antenna(az_deg, time_s, scan);
        app->display_ms += ms_since(t0);
    }
#else
    (void)app; (void)az_deg; (void)time_s; (void)scan;
#endif
}

/* ------------------------------------------------------------------------ */
/* Worker: input -> processing                                              */
/* ------------------------------------------------------------------------ */

/* Block times of the processor and the application's console and display times. */
static void app_timing(const App *app, const Processor *proc, ProcessorTiming *t)
{
    processor_timing(proc, t);
    t->console_ms = app->console_ms;
    t->display_ms = app->display_ms;
}

static void print_timing(const App *app, const Processor *proc, const UdpInput *in, uint64_t queueDrops, double wall_s)
{
    ProcessorTiming t;
    UdpStats u;
    app_timing(app, proc, &t);
    udp_input_stats(in, &u);
    const double n = t.pris ? (double)t.pris : 1.0;
    fprintf(stderr, "\nPRIs processed: %llu (%llu CPI updates) in %.1f s; datagrams %llu, foreign %llu, "
            "incomplete PRIs %llu\n", (unsigned long long)t.pris, (unsigned long long)t.updates, wall_s,
            (unsigned long long)u.datagrams, (unsigned long long)u.foreign, (unsigned long long)u.prisIncomplete);
    fprintf(stderr, "Processing time per PRI: %.3f ms (DBF %.3f, MF %.3f, MTI %.3f, Doppler %.3f, "
            "CFAR+pruning %.3f, extraction %.3f, tracking %.3f)\n", t.total_ms / n, t.dbf_ms / n, t.mf_ms / n,
            t.mti_ms / n, t.doppler_ms / n, t.cfar_ms / n, t.extract_ms / n, t.track_ms / n);
    fprintf(stderr, "Output per PRI: console %.3f ms, display (PPI calls) %.3f ms\n", t.console_ms / n, t.display_ms / n);
    if (t.skipped || queueDrops)
        fprintf(stderr, "Overload: %llu CPI updates skipped, %llu PRIs dropped because the queue was full\n",
                (unsigned long long)t.skipped, (unsigned long long)queueDrops);
    if (t.rejected)
        fprintf(stderr, "PRIs ignored: %llu (their DDC layout differs from the configuration: export it again "
                "with export_receiver_config.m)\n", (unsigned long long)t.rejected);
}

/* Prints the CFAR threshold in use and shows it in the PPI. */
static void report_cfar(App *app, const Processor *proc, double default_dB, int changed)
{
    const double alpha = processor_cfar_alpha(proc), dB = 10.0 * log10(alpha);
    const double pfa = cfar_pfa(alpha, app->cfg.nTrain, app->cfg.nLooks);
    if (changed || fabs(dB - default_dB) > 1e-9) {
        const double t0 = now_seconds();
        if (isnan(pfa)) fprintf(stderr, "CFAR threshold %.2f dB\n", dB);
        else fprintf(stderr, "CFAR threshold %.2f dB (Pfa %.1e)\n", dB, pfa);
        app->console_ms += ms_since(t0);
    }
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        const double t0 = now_seconds();
        ppi_post_cfar_threshold(dB, pfa, default_dB);
        app->display_ms += ms_since(t0);
    }
#endif
}

static int worker(void *arg)
{
    App *app = arg;
    char err[256] = "";
    UdpInput *in = udp_input_open(app->opt.port, app->opt.recordPath, app->opt.replayPath, err, sizeof err);
    if (!in) {
        fprintf(stderr, "input: %s\n", err);
        app->exitCode = 1;
        return 1;
    }
    Processor *proc = processor_create(&app->cfg, app->opt.threads);
    if (!proc) {
        fprintf(stderr, "cannot create the processor\n");
        udp_input_close(in);
        app->exitCode = 1;
        return 1;
    }
    ProcessorCallbacks cb = { on_detections, on_plot, on_antenna, app, on_track };
#ifdef RADAR_WITH_QT
    char detName[32], trackName[32];
    snprintf(detName, sizeof detName, "CFAR %s", cfar_method_name(app->cfg.cfarMethod));
    if (app->cfg.trackEnable)
        snprintf(trackName, sizeof trackName, "Tracking %d-D %s", app->cfg.trackDim, app->cfg.trackModel == TRACK_CA ? "CA" : "CV");
    else snprintf(trackName, sizeof trackName, "Tracking (off)");
    double lastTiming = 0.0;
#endif
    const double default_dB = 10.0 * log10(app->cfg.alpha);
    if (!isnan(app->opt.cfarThreshold_dB))
        processor_set_cfar_alpha(proc, pow(10.0, app->opt.cfarThreshold_dB / 10.0));
    report_cfar(app, proc, default_dB, 0);
    /* live stream: a full queue drops PRIs; replay: the reader waits */
    const int live = app->opt.replayPath == NULL;
    InputThread *input = input_thread_start(in, app->opt.queue, live);
    if (!input) {
        fprintf(stderr, "cannot start the input thread\n");
        processor_destroy(proc);
        udp_input_close(in);
        app->exitCode = 1;
        return 1;
    }
    if (live)
        fprintf(stderr, "Waiting for the DDC stream on UDP port %d ...\n", app->opt.port);
    /* While the processing is behind a live stream (more than a quarter of the
     * queue waiting), the Doppler/CFAR updates are skipped; the fast-time
     * processing keeps the CPI buffers complete. */
    const int skipAbove = app->opt.queue / 4 > 2 ? app->opt.queue / 4 : 2;
    const double t0 = now_seconds();
    double lastData = 0;                                        /* 0: no PRI yet */
    while (!g_stop) {
        PriData *pri = NULL;
        int r = input_thread_next(input, &pri, 200);
        if (r < 0) break;                                      /* end of the replay file */
        if (r == 0) {
            if (app->opt.idleTimeout_s > 0 && lastData > 0 && now_seconds() - lastData > app->opt.idleTimeout_s) break;
            continue;
        }
        lastData = now_seconds();
#ifdef RADAR_WITH_QT
        double request_dB;
        const double tReq = now_seconds();
        const int request = app->opt.ppi && ppi_take_cfar_request(&request_dB);
        if (app->opt.ppi) app->display_ms += ms_since(tReq);
        if (request) {                                         /* Apply in the PPI */
            processor_set_cfar_alpha(proc, pow(10.0, request_dB / 10.0));
            report_cfar(app, proc, default_dB, 1);
        }
#endif
        if (live && !app->opt.noSkip)
            processor_set_slow_time(proc, input_thread_backlog(input) <= skipAbove);
        processor_process(proc, pri, &cb);
#ifdef RADAR_WITH_QT
        if (app->opt.ppi && lastData - lastTiming > 0.25) {   /* block times in the PPI panel */
            const double tPost = now_seconds();
            ProcessorTiming tm;
            app_timing(app, proc, &tm);
            ppi_post_timing(&tm, detName, trackName);
            app->display_ms += ms_since(tPost);
            lastTiming = lastData;
        }
#endif
    }
    const uint64_t queueDrops = input_thread_dropped(input);
    input_thread_stop(input);
    processor_flush(proc, &cb);
    print_timing(app, proc, in, queueDrops, now_seconds() - t0);
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        ProcessorTiming tm;
        app_timing(app, proc, &tm);
        ppi_post_timing(&tm, detName, trackName);
    }
#endif
    fprintf(stderr, "Plots: %llu, CFAR detections after pruning: %llu\n",
            (unsigned long long)app->nPlots, (unsigned long long)app->nDetections);
    if (app->cfg.trackEnable)
        fprintf(stderr, "Tracks (EKF %d-D): %llu updates, %llu confirmed\n", app->cfg.trackDim,
                (unsigned long long)app->nTrackUpdates, (unsigned long long)app->nConfirmed);
    processor_destroy(proc);
    udp_input_close(in);
    return 0;
}

/* Creates a CSV file with its header line; NULL (and a warning) if it cannot
 * be created, NULL without a warning when no file is wanted ("" or NULL). */
static FILE *open_csv(const char *path, const char *header)
{
    if (!path || !path[0]) return NULL;
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "warning: cannot create %s, it is not written\n", path);
        return NULL;
    }
    fprintf(f, "%s\n", header);
    return f;
}

int main(int argc, char **argv)
{
    App app;
    memset(&app, 0, sizeof app);
    int r = parse_options(argc, argv, &app.opt);
    if (r != 0) return r < 0 ? 2 : 0;

    char err[256] = "";
    if (radar_config_load(app.opt.configPath, &app.cfg, err, sizeof err) != 0) {
        fprintf(stderr, "configuration: %s\n", err);
        return 1;
    }
    app.cfg.mfMethod = !strcmp(app.opt.mf, "direct") ? MF_DIRECT : !strcmp(app.opt.mf, "fft") ? MF_FFT : MF_AUTO;
    if (app.opt.updateEvery > 0 && app.opt.updateEvery != app.cfg.updateEvery) {
        /* the plot extractor needs the same azimuth extent with the new update interval */
        app.cfg.updateEvery = app.opt.updateEvery;
        const double step_deg = app.cfg.rate_dps * app.cfg.meanPri_s * app.cfg.updateEvery;
        const int minHits = (int)ceil(app.cfg.minExtent_deg / step_deg);     /* as run_radar_sim.m */
        app.cfg.minHits = minHits > 1 ? minHits : 1;
    }
    app.detArc_deg = app.cfg.rate_dps * app.cfg.meanPri_s * app.cfg.updateEvery;
    fprintf(stderr, "Configuration: %d element(s), %d beam(s), DDC %d samples, short %d / long %d range cells, "
            "CPI %d pulses, FFT %d x %d row(s), update every %d pulse(s), plots need %d update(s)\n", app.cfg.nCh,
            app.cfg.nBeams, app.cfg.nDdc, app.cfg.pulse[0].nRange, app.cfg.pulse[1].nRange, app.cfg.cpiLength,
            app.cfg.nfft, app.cfg.nRows, app.cfg.updateEvery, app.cfg.minHits);
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &app.cfg.pulse[i];
        if (pc->waveform == WAVE_PHASE_CODE)
            fprintf(stderr, "%s pulse: phase code of %d chips, %s filter of %d taps (%.2f ... %.2f km)\n", pulse_name(i),
                    pc->nChips, pc->compression == COMPRESSION_MISMATCH ? "mismatch" : "matched", pc->nTaps,
                    pc->range_m[0] / 1e3, pc->range_m[pc->nRange - 1] / 1e3);
        else
            fprintf(stderr, "%s pulse: LFM, matched filter of %d taps (%.2f ... %.2f km)\n", pulse_name(i), pc->nTaps,
                    pc->range_m[0] / 1e3, pc->range_m[pc->nRange - 1] / 1e3);
    }
    fprintf(stderr, "CFAR %s (window: range %d+%d, Doppler %d+%d guard+training cells per side), ", cfar_method_name(app.cfg.cfarMethod),
            app.cfg.guardRange, app.cfg.trainRange, app.cfg.guardDoppler, app.cfg.trainDoppler);
    if (app.cfg.trackEnable)
        fprintf(stderr, "tracking: EKF %d-D, %s model, GNN association\n", app.cfg.trackDim,
                app.cfg.trackModel == TRACK_CA ? "constant-acceleration" : "constant-velocity");
    else fprintf(stderr, "tracking off\n");

    app.plotsFile = open_csv(app.opt.plotsCsv,
                             "time_s,az_deg,range_m,el_deg,vr_mps,snr_dB,pulse,hits,beam,rangeCell,dopplerBin,pulseNumber");
    app.trackFile = open_csv(app.opt.tracksCsv, "id,confirmed,deleted,time_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,"
                                                "hits,misses,ax_mps2,ay_mps2,az_mps2");
    app.detFile = open_csv(app.opt.detectionsCsv,
                           "pulseNumber,pulse,time_s,az_deg,range_m,el_deg,vr_mps,snr_dB,beam,rangeCell,dopplerBin");
    if (app.opt.serialDevice) {
        app.serial = serial_open(app.opt.serialDevice, app.opt.baud, err, sizeof err);
        if (!app.serial) {
            fprintf(stderr, "serial: %s\n", err);
            return 1;
        }
    }
    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    if (app.opt.ppi) {
#ifdef RADAR_WITH_QT
        ppi_set_display_mode(!strcmp(app.opt.ppiShow, "detections") ? PPI_SHOW_DETECTIONS
                             : !strcmp(app.opt.ppiShow, "both") ? PPI_SHOW_BOTH : PPI_SHOW_PLOTS);
        app.exitCode = ppi_run(argc, argv, app.cfg.maxRange_m, worker, &app, &g_stop, app.opt.snapshotPath);
#else
        fprintf(stderr, "this build has no Qt: --ppi is not available\n");
        app.exitCode = 1;
#endif
    } else {
        worker(&app);
    }

    if (app.serial) {
        uint64_t sent, dropped, bytes;
        serial_stats(app.serial, &sent, &dropped, &bytes);
        fprintf(stderr, "Serial: %llu frames sent, %llu dropped, %llu bytes\n", (unsigned long long)sent,
                (unsigned long long)dropped, (unsigned long long)bytes);
        serial_close(app.serial);
    }
    if (app.plotsFile) fclose(app.plotsFile);
    if (app.detFile) fclose(app.detFile);
    if (app.trackFile) fclose(app.trackFile);
    radar_config_free(&app.cfg);
    return app.exitCode;
}

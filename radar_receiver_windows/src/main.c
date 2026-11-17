/*
 * main.c - radar receiver: DDC output over UDP -> processing -> plots.
 *
 *   radar_receiver [options]
 *     --config FILE          receiver configuration from export_receiver_config.m
 *                            (default config/radar_config.bin)
 *     --port N               UDP port of the DDC stream (default 50000)
 *     --format rddc|vbdr     packet format of the DDC stream (default: of the configuration,
 *                            cfg.udp.format): 'rddc' or 'vbdr' (the VBDR board, README)
 *   Radar Settings of the PPI (tab Controls; the VBDR board, README "VBDR stream and board control"):
 *     --config-m FILE        radar_config.m of the modes (default ../radar_config.m of this project):
 *                            the settings of every mode and the configuration of a new mode
 *     --config-set STMT      a MATLAB assignment run after it (repeatable, as radar_export_config --set;
 *                            give the ones the configuration file was exported with)
 *     --no-board             no Radar Settings (no commands to the board)
 *     --record FILE          also store the received datagrams in FILE
 *     --replay FILE          read the datagrams from FILE instead of the network
 *     --ppi                  show the PPI (Qt build only)
 *     --ppi-snapshot FILE    save the PPI as PNG at the end and close it
 *     --ppi-show MODE        PPI shows 'plots' (default, after the plot extractor),
 *                            'detections' (CFAR detections before it), 'both' or
 *                            'video' (DDC output averaged over the elements), 'dbf'
 *                            (DBF output of the selected beam) or 'mf' (pulse
 *                            compression output of the selected beam)
 *     --ppi-beam N           beam shown by dbf / mf and the range-Doppler tab (1 ... beams, default 1)
 *     --ppi-tab TAB          tab at the start: 'ppi' (default) or 'rd' (range-Doppler map + small PPI)
 *     --ppi-rd-rate HZ       range-Doppler maps per second of the tab (default cfg.display.rdRate_Hz)
 *     --ppi-rd-plot P        plot below the range-Doppler map: 'spectrum' (default) or 'range' (profile)
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
 *     --log-dir DIR          log the processing chain from the start to DIR/<name> (also: tab
 *                            "Controls" of the PPI): one MAT-file per PRI and block (logger.c)
 *     --log-name NAME        name of the log folder (default log_<date>_<time>)
 *     --log-blocks LIST      blocks of the log, e.g. ddc,dbf,mf,doppler,detection (default ddc,detection);
 *                            'screen': with --ppi also a video of the window (screen_record.mp4, ffmpeg)
 *     --ffmpeg PATH          ffmpeg program of the screen recording (default: next to the program, else PATH)
 *     --record-audio SPEC    sound of the screen recording: 'auto' (default: the default microphone),
 *                            'none', or an ffmpeg input 'FORMAT:DEVICE', e.g. pulse:default, alsa:hw:1,
 *                            dshow:audio=Microphone (Realtek Audio)
 *     --record-fps N         frames per second of the screen recording (default 10)
 *     --cfar-threshold DB    CFAR threshold factor in dB (default from the configuration,
 *                            i.e. from cfg.cfar.pfa); can also be changed in the PPI
 *     --quiet                do not print every plot and confirmed track
 *     --print-plots          print them (the default comes from the configuration,
 *                            cfg.console.plots; printing can take longer than the
 *                            processing when a low CFAR threshold gives many plots)
 *   real time (Linux, e.g. with a PREEMPT_RT kernel; see README "Linux and real time"):
 *     --rt-priority N        SCHED_FIFO priority N (1 ... 98) of the processing threads
 *                            (with their OpenMP workers), N + 1 for the UDP input thread
 *     --cpus LIST            CPUs of the processing threads, e.g. 2-7 (default thread
 *                            count: their number)
 *     --input-cpus LIST      CPUs of the UDP input thread (default: those of --cpus)
 *     --mlock                lock the memory of the process (no page faults)
 */
#include "radar_rx.h"

#ifndef RADAR_SOURCE_DIR
#  define RADAR_SOURCE_DIR "."
#endif

#include <math.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *configPath, *recordPath, *replayPath, *snapshotPath;
    const char *format;           /* --format (NULL: of the configuration) */
    const char *configM;          /* --config-m */
    const char *configSets[16];   /* --config-set */
    int nConfigSets, noBoard;
    const char *serialDevice, *serialContent, *plotsCsv, *detectionsCsv, *tracksCsv, *mf, *ppiShow, *ppiTab, *ppiRdPlot;
    int ppiBeam;                  /* 1-based */
    int port, baud, threads, updateEvery, ppi, quiet, queue, noSkip;
    double rdRate;                   /* --ppi-rd-rate (0 = from the configuration) */
    int rtPriority, mlock;        /* real time (Linux) */
    const char *cpus, *inputCpus;
    double idleTimeout_s;
    double cfarThreshold_dB;      /* NAN: from the configuration */
    const char *logDir, *logName, *logBlocks;   /* logging from the start (--log-dir) */
    const char *ffmpeg, *recordAudio;             /* screen recording (PPI) */
    double recordFps;
    int logMask;
} Options;

typedef struct {
    Options opt;
    RadarConfig cfg;
    SerialOut *serial;
    FILE *plotsFile, *detFile, *trackFile;
    uint64_t nPlots, nDetections, nTrackUpdates, nConfirmed;
    uint32_t *confirmedIds;       /* confirmed track numbers already reported */
    size_t nConfirmedIds, capConfirmedIds;
    double detArc_deg;            /* azimuth step between CPI updates (PPI detection arcs) */
    double console_ms, display_ms;   /* time in the console output and in the PPI calls */
    double log_ms;                /* time of the logging in the processing thread (copies) */
    Logger *log;                  /* the running log (NULL: none) */
    BoardControl *board;          /* Radar Settings of the PPI (NULL: off) */
    int havePri;                  /* the last PRI processed (radar_parameters.mat) */
    uint32_t lastPulse;
    double lastTime_s, lastAz_deg;
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
    puts("usage: radar_receiver [--config FILE] [--port N] [--format rddc|vbdr] [--record FILE] [--replay FILE]\n"
         "                      [--config-m FILE] [--config-set STMT]... [--no-board]\n"
         "                      [--ppi] [--ppi-snapshot FILE] [--ppi-show plots|detections|both|video|dbf|mf]\n"
         "                      [--ppi-beam N] [--ppi-tab ppi|rd] [--ppi-rd-rate HZ]\n"
         "                      [--ppi-rd-plot spectrum|range]\n"
         "                      [--serial DEVICE] [--baud N]\n"
         "                      [--serial-content detections|plots] [--threads N] [--mf auto|direct|fft]\n"
         "                      [--update-every N]\n"
         "                      [--plots-csv FILE] [--detections-csv FILE] [--tracks-csv FILE]\n"
         "                      [--idle-timeout S] [--quiet | --print-plots]\n"
         "                      [--queue N] [--no-skip] [--cfar-threshold DB]\n"
         "                      [--log-dir DIR] [--log-name NAME] [--log-blocks ddc,dbf,mf,doppler,detection,screen]\n"
         "                      [--ffmpeg PATH] [--record-audio SPEC] [--record-fps N]\n"
         "                      [--rt-priority N] [--cpus LIST] [--input-cpus LIST] [--mlock]");
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
    o->ppiTab = "ppi";
    o->ppiRdPlot = "spectrum";
    o->ppiBeam = 1;
    o->queue = 64;
    o->cfarThreshold_dB = NAN;
    o->quiet = -1;                                   /* from the configuration (consolePlots) */
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
#define TAKE(field) do { if (!v) { fprintf(stderr, "%s needs a value\n", a); return -1; } field; i++; } while (0)
        if      (!strcmp(a, "--config"))         TAKE(o->configPath = v);
        else if (!strcmp(a, "--port"))           TAKE(o->port = atoi(v));
        else if (!strcmp(a, "--format"))         TAKE(o->format = v);
        else if (!strcmp(a, "--config-m"))       TAKE(o->configM = v);
        else if (!strcmp(a, "--config-set"))     TAKE((o->nConfigSets < 16 ? (void)(o->configSets[o->nConfigSets++] = v) : (void)0));
        else if (!strcmp(a, "--no-board"))       o->noBoard = 1;
        else if (!strcmp(a, "--record"))         TAKE(o->recordPath = v);
        else if (!strcmp(a, "--replay"))         TAKE(o->replayPath = v);
        else if (!strcmp(a, "--ppi"))            o->ppi = 1;
        else if (!strcmp(a, "--ppi-snapshot"))   TAKE((o->snapshotPath = v, o->ppi = 1));
        else if (!strcmp(a, "--ppi-show"))       TAKE(o->ppiShow = v);
        else if (!strcmp(a, "--ppi-beam"))       TAKE(o->ppiBeam = atoi(v));
        else if (!strcmp(a, "--ppi-tab"))        TAKE(o->ppiTab = v);
        else if (!strcmp(a, "--ppi-rd-plot"))    TAKE(o->ppiRdPlot = v);
        else if (!strcmp(a, "--ppi-rd-rate"))    TAKE(o->rdRate = atof(v));
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
        else if (!strcmp(a, "--print-plots"))    o->quiet = 0;
        else if (!strcmp(a, "--queue"))          TAKE(o->queue = atoi(v));
        else if (!strcmp(a, "--no-skip"))        o->noSkip = 1;
        else if (!strcmp(a, "--cfar-threshold")) TAKE(o->cfarThreshold_dB = atof(v));
        else if (!strcmp(a, "--log-dir"))        TAKE(o->logDir = v);
        else if (!strcmp(a, "--log-name"))       TAKE(o->logName = v);
        else if (!strcmp(a, "--log-blocks"))     TAKE(o->logBlocks = v);
        else if (!strcmp(a, "--ffmpeg"))         TAKE(o->ffmpeg = v);
        else if (!strcmp(a, "--record-audio"))   TAKE(o->recordAudio = v);
        else if (!strcmp(a, "--record-fps"))     TAKE(o->recordFps = atof(v));
        else if (!strcmp(a, "--rt-priority"))    TAKE(o->rtPriority = atoi(v));
        else if (!strcmp(a, "--cpus"))           TAKE(o->cpus = v);
        else if (!strcmp(a, "--input-cpus"))     TAKE(o->inputCpus = v);
        else if (!strcmp(a, "--mlock"))          o->mlock = 1;
        else if (!strcmp(a, "--help") || !strcmp(a, "-h")) { usage(); return 1; }
        else { fprintf(stderr, "unknown option %s\n", a); usage(); return -1; }
#undef TAKE
    }
    if (o->format && strcmp(o->format, "rddc") && strcmp(o->format, "vbdr")) {
        fprintf(stderr, "--format must be 'rddc' or 'vbdr'\n");
        return -1;
    }
    if (strcmp(o->mf, "auto") && strcmp(o->mf, "direct") && strcmp(o->mf, "fft")) {
        fprintf(stderr, "--mf must be 'auto', 'direct' or 'fft'\n");
        return -1;
    }
    if (!isnan(o->cfarThreshold_dB) && !(o->cfarThreshold_dB >= 0.0 && o->cfarThreshold_dB <= 50.0)) {
        fprintf(stderr, "--cfar-threshold must be 0 ... 50 dB\n");
        return -1;
    }
    if (o->rtPriority < 0 || o->rtPriority > 98) {
        fprintf(stderr, "--rt-priority must be 1 ... 98 (0 = normal scheduling)\n");
        return -1;
    }
    if (strcmp(o->ppiShow, "plots") && strcmp(o->ppiShow, "detections") && strcmp(o->ppiShow, "both") &&
        strcmp(o->ppiShow, "video") && strcmp(o->ppiShow, "dbf") && strcmp(o->ppiShow, "mf")) {
        fprintf(stderr, "--ppi-show must be 'plots', 'detections', 'both', 'video', 'dbf' or 'mf'\n");
        return -1;
    }
    if (strcmp(o->ppiTab, "ppi") && strcmp(o->ppiTab, "rd")) {
        fprintf(stderr, "--ppi-tab must be 'ppi' or 'rd'\n");
        return -1;
    }
    if (strcmp(o->ppiRdPlot, "spectrum") && strcmp(o->ppiRdPlot, "range")) {
        fprintf(stderr, "--ppi-rd-plot must be 'spectrum' or 'range'\n");
        return -1;
    }
    if (o->ppiBeam < 1) {
        fprintf(stderr, "--ppi-beam must be 1 or more\n");
        return -1;
    }
    if (strcmp(o->serialContent, "detections") && strcmp(o->serialContent, "plots")) {
        fprintf(stderr, "--serial-content must be 'detections' or 'plots'\n");
        return -1;
    }
    o->logMask = LOG_DDC | LOG_DETECTION;
    if (o->logBlocks) {                                  /* ddc,dbf,mf,doppler,detection */
        static const struct { const char *name; int bit; } blocks[6] = {
            { "ddc", LOG_DDC }, { "dbf", LOG_DBF }, { "mf", LOG_MF }, { "doppler", LOG_DOPPLER }, { "detection", LOG_DETECTION },
            { "screen", LOG_SCREEN } };
        o->logMask = 0;
        const char *s = o->logBlocks;
        while (*s) {
            size_t n = strcspn(s, ",");
            int found = 0;
            for (int i = 0; i < 6; i++)
                if (n == strlen(blocks[i].name) && !strncmp(s, blocks[i].name, n)) { o->logMask |= blocks[i].bit; found = 1; }
            if (!found && n > 0) {
                fprintf(stderr, "--log-blocks: unknown block '%.*s' (ddc, dbf, mf, doppler, detection, screen)\n", (int)n, s);
                return -1;
            }
            s += n + (s[n] == ',');
        }
        if (!o->logMask) {
            fprintf(stderr, "--log-blocks needs at least one block\n");
            return -1;
        }
        if ((o->logMask & LOG_SCREEN) && !o->ppi)
            fprintf(stderr, "--log-blocks screen: the screen recording needs --ppi (ignored)\n");
    }
    return 0;
}

/* The options that change a loaded configuration (also of a new radar mode). */
static void apply_options(App *app, RadarConfig *cfg)
{
    if (app->opt.format) cfg->udpFormat = !strcmp(app->opt.format, "vbdr") ? UDP_FORMAT_VBDR : UDP_FORMAT_RDDC;
    cfg->mfMethod = !strcmp(app->opt.mf, "direct") ? MF_DIRECT : !strcmp(app->opt.mf, "fft") ? MF_FFT : MF_AUTO;
    if (app->opt.rdRate > 0) cfg->rdRate_Hz = app->opt.rdRate;
    if (app->opt.updateEvery > 0 && app->opt.updateEvery != cfg->updateEvery) {
        /* the plot extractor needs the same azimuth extent with the new update interval */
        cfg->updateEvery = app->opt.updateEvery;
        const double step_deg = cfg->rate_dps * cfg->meanPri_s * cfg->updateEvery;
        const int minHits = (int)ceil(cfg->minExtent_deg / step_deg);     /* as run_radar_sim.m */
        cfg->minHits = minHits > 1 ? minHits : 1;
        for (int j = 0; j < cfg->nPrf; j++) {                              /* every PRF of the scans */
            const int h = (int)ceil(cfg->minExtent_deg / (cfg->rate_dps * cfg->prf[j].pri_s * cfg->updateEvery));
            cfg->prf[j].minHits = h > 1 ? h : 1;
        }
    }
    app->detArc_deg = cfg->rate_dps * cfg->meanPri_s * cfg->updateEvery;
}

/* Doppler axis of the range-Doppler tab [m/s] (the long pulse, else the short one). */
static void rd_axis(const RadarConfig *cfg, double *v0, double *v1)
{
    const int vp = cfg->pulse[PULSE_LONG].enable ? PULSE_LONG : PULSE_SHORT;
    *v0 = INFINITY;
    *v1 = -INFINITY;
    for (int d = 0; d < cfg->nfft; d++) {
        const double v = -cfg->dopFreq_Hz[d] * cfg->pulse[vp].lambda_m / 2;
        *v0 = v < *v0 ? v : *v0;
        *v1 = v > *v1 ? v : *v1;
    }
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
    if (app->log && n > 0) {
        const double t0 = now_seconds();
        logger_detections(app->log, dets, n);
        app->log_ms += ms_since(t0);
    }
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
    if (app->log) {
        const double t0 = now_seconds();
        logger_plot(app->log, p);
        app->log_ms += ms_since(t0);
    }
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
    int first = 0;                               /* first report of a confirmed track number */
    if (t->confirmed && !t->deleted) {
        first = 1;
        for (size_t k = 0; k < app->nConfirmedIds && first; k++) first = app->confirmedIds[k] != t->id;
        if (first) {
            if (app->nConfirmedIds == app->capConfirmedIds) {
                app->capConfirmedIds = app->capConfirmedIds ? 2 * app->capConfirmedIds : 64;
                uint32_t *p = realloc(app->confirmedIds, sizeof(uint32_t) * app->capConfirmedIds);
                if (!p) first = 0; else app->confirmedIds = p;
            }
            if (first) app->confirmedIds[app->nConfirmedIds++] = t->id;
        }
    }
    if (first) {
        app->nConfirmed++;
        if (!app->opt.quiet) {
            const double t0 = now_seconds();
            const double R = hypot(hypot(t->pos_m[0], t->pos_m[1]), t->pos_m[2]);
            const double vr = (t->pos_m[0] * t->vel_mps[0] + t->pos_m[1] * t->vel_mps[1] + t->pos_m[2] * t->vel_mps[2]) / R;
            printf("track %u confirmed  t %8.3f s  R %8.3f km  az %7.2f deg  speed %6.1f m/s  radial velocity %7.1f m/s\n",
                   t->id, t->time_s, R / 1e3,
                   fmod(atan2(t->pos_m[0], t->pos_m[1]) * 180.0 / 3.14159265358979323846 + 360.0, 360.0),
                   hypot(t->vel_mps[0], t->vel_mps[1]), vr);
            app->console_ms += ms_since(t0);
        }
    }
    if (app->trackFile) {
        double pm[3] = { NAN, NAN, NAN };        /* PMBM: probabilities of the models CV, CA, CT */
        for (int j = 0; j < t->nModes; j++) pm[t->modes[j]] = t->modeProb[j];
        fprintf(app->trackFile, "%u,%d,%d,%.6f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                t->id, t->confirmed, t->deleted, t->time_s, t->pos_m[0], t->pos_m[1], t->pos_m[2], t->vel_mps[0],
                t->vel_mps[1], t->vel_mps[2], t->hits, t->misses, t->acc_mps2[0], t->acc_mps2[1], t->acc_mps2[2],
                t->existence, pm[0], pm[1], pm[2]);
    }
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        const double t0 = now_seconds();
        ppi_post_track(t);
        app->display_ms += ms_since(t0);
    }
#endif
}

/* DBF / pulse compression output of the selected beam (PPI display modes dbf, mf). */
static void on_profile(void *user, int kind, const float *row, int nBins, double dR_m, double az_deg, double time_s)
{
    App *app = user;
#ifdef RADAR_WITH_QT
    const double t0 = now_seconds();
    ppi_post_profile(kind, row, nBins, dR_m, az_deg, time_s);
    app->display_ms += ms_since(t0);
#else
    (void)app; (void)kind; (void)row; (void)nBins; (void)dR_m; (void)az_deg; (void)time_s;
#endif
}

/* Range-Doppler map of the selected beam (range-Doppler tab of the PPI). */
static void on_rdmap(void *user, const RdMap *map)
{
    App *app = user;
#ifdef RADAR_WITH_QT
    const double t0 = now_seconds();
    ppi_post_rdmap(map);
    app->display_ms += ms_since(t0);
#else
    (void)app; (void)map;
#endif
}

/* Block outputs of the processor (DBF, MF, Doppler) for the log. */
static void on_block(void *user, const BlockOutput *out)
{
    App *app = user;
    const double t0 = now_seconds();
    logger_block(app->log, out);
    app->log_ms += ms_since(t0);
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
    t->log_ms = app->log_ms;
}

/* Starts a log (folder/name, mask of LOG_*); reports it on the console and in the PPI. */
static void start_log(App *app, Processor *proc, const char *folder, const char *name, int mask)
{
    char err[512] = "";
    if (app->log) return;
    RuntimeParams rt;                                  /* the state now, for radar_parameters.mat */
    memset(&rt, 0, sizeof rt);
    rt.configPath = app->opt.configPath;
    rt.alpha = processor_cfar_alpha(proc);
    rt.threads = processor_threads(proc);
    rt.slowTimeSkipping = app->opt.replayPath == NULL && !app->opt.noSkip;   /* live: Doppler/CFAR skipped when behind */
    rt.havePri = app->havePri;
    rt.pulseNumber = app->lastPulse;
    rt.time_s = app->lastTime_s;
    rt.az_deg = app->lastAz_deg;
    app->log = logger_start(folder, name, mask, &app->cfg, &rt, err, sizeof err);
    processor_set_log(proc, app->log ? mask & LOG_BLOCKS : 0);
    if (app->log) fprintf(stderr, "Logging to %s\n", logger_dir(app->log));
    else fprintf(stderr, "logging: %s\n", err);
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        LoggerStats st;
        logger_stats(app->log, &st);
        ppi_post_log_status(app->log != NULL, app->log ? logger_dir(app->log) : "", mask, &st, app->log ? NULL : err);
    }
#endif
}

/* Ends the log: the queued files are written first. */
static void stop_log(App *app, Processor *proc)
{
    if (!app->log) return;
    LoggerStats st;
    char dir[1024];
    snprintf(dir, sizeof dir, "%s", logger_dir(app->log));
    if (proc) processor_set_log(proc, 0);
    Logger *L = app->log;
    app->log = NULL;
    logger_stop(L, &st);                           /* writes the queued files; st: the final state */
    fprintf(stderr, "Log %s ended\n", dir);
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) ppi_post_log_status(0, dir, 0, &st, NULL);
#else
    (void)st;
#endif
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
            "CFAR+pruning %.3f, peak interpolation %.3f, extraction %.3f, tracking %.3f)\n", t.total_ms / n,
            t.dbf_ms / n, t.mf_ms / n, t.mti_ms / n, t.doppler_ms / n, t.cfar_ms / n, t.interp_ms / n,
            t.extract_ms / n, t.track_ms / n);
    fprintf(stderr, "Output per PRI: console %.3f ms, display (PPI calls) %.3f ms, logging (copies) %.3f ms\n", t.console_ms / n,
            t.display_ms / n, t.log_ms / n);
    if (t.skipped || queueDrops)
        fprintf(stderr, "Overload: %llu CPI updates skipped, %llu PRIs dropped because the queue was full\n",
                (unsigned long long)t.skipped, (unsigned long long)queueDrops);
    if (t.rejected)
        fprintf(stderr, "PRIs ignored: %llu (their DDC layout differs from the configuration: export it again "
                "with export_receiver_config.m)\n", (unsigned long long)t.rejected);
}

/* Sets the CFAR threshold [dB], not below the threshold of cfg.cfar.maxPfa
 * (cfg.alphaMin). */
static void set_cfar_threshold(App *app, Processor *proc, double dB)
{
    double alpha = pow(10.0, dB / 10.0);
    if (alpha < app->cfg.alphaMin) {
        alpha = app->cfg.alphaMin;
        fprintf(stderr, "CFAR threshold %.2f dB is below the lowest threshold %.2f dB (cfg.cfar.maxPfa): limited to it\n",
                dB, 10.0 * log10(alpha));
    }
    processor_set_cfar_alpha(proc, alpha);
}

/* Prints the CFAR threshold in use and shows it in the PPI. */
static void report_cfar(App *app, const Processor *proc, double default_dB, int changed)
{
    const double alpha = processor_cfar_alpha(proc), dB = 10.0 * log10(alpha);
    const double pfa = cfar_pfa(&app->cfg, alpha);
    if (changed || fabs(dB - default_dB) > 1e-9) {
        const double t0 = now_seconds();
        if (isnan(pfa)) fprintf(stderr, "CFAR threshold %.2f dB\n", dB);
        else fprintf(stderr, "CFAR threshold %.2f dB (Pfa %.1e)\n", dB, pfa);
        app->console_ms += ms_since(t0);
    }
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        const double t0 = now_seconds();
        ppi_post_cfar_threshold(dB, pfa, default_dB,
                                app->cfg.alphaMin > 0.0 ? 10.0 * log10(app->cfg.alphaMin) : 0.0);
        app->display_ms += ms_since(t0);
    }
#endif
}

/* Labels of the CFAR and tracking rows of the timing panel (15 characters). */
static void block_names(const RadarConfig *cfg, char *det, size_t detLen, char *trk, size_t trkLen)
{
    snprintf(det, detLen, "CFAR %s %s", cfar_type_name(cfg->cfarType), cfar_method_name(cfg->cfarMethod));
    if (!cfg->trackEnable) snprintf(trk, trkLen, "Track (off)");
    else if (cfg->trackMethod == TRACK_METHOD_PMBM) snprintf(trk, trkLen, "Track PMBM %dD", cfg->trackDim);
    else if (cfg->trackMethod == TRACK_METHOD_TOMHT) snprintf(trk, trkLen, "Track MHT %dD", cfg->trackDim);
    else snprintf(trk, trkLen, "Track EKF %dD %s", cfg->trackDim, cfg->trackModel == TRACK_CA ? "CA" : "CV");
}

/* Radar Settings of the PPI: the settings of the modes from radar_config.m and the board control
 * thread, which first reads the settings of the board (no command is written at the start). */
static void start_board(App *app)
{
    if (app->opt.noBoard) return;
    char err[1024];
    BoardInfo *bi = calloc(1, sizeof *bi);
    if (!bi) return;
    const char *mPath = app->opt.configM ? app->opt.configM : RADAR_SOURCE_DIR "/../radar_config.m";
    if (board_info_load(mPath, app->opt.configSets, app->opt.nConfigSets, bi, err, sizeof err) != 0) {
        fprintf(stderr, "Radar Settings off (the modes of %s: %s)\n", mPath, err);
        free(bi);
        return;
    }
    ModeSettings cur;
    board_settings_of_config(bi, &app->cfg, &cur);
    app->board = board_control_start(bi, app->cfg.radarMode, &cur, err, sizeof err);
    if (!app->board) fprintf(stderr, "Radar Settings off: %s\n", err);
    else {
        fprintf(stderr, "Radar Settings: board %s:%d, answers on UDP port %d, %d mode(s) of %s\n", bi->host, bi->port,
                bi->replyPort, bi->nModes, mPath);
        board_control_request(app->board, BOARD_REQ_CHECK, 0, 0, NULL);
    }
    free(bi);
}

typedef struct { App *app; RadarConfig *next; } SwapArg;

static void swap_config(void *arg)                 /* (the logger's writer is idle) */
{
    SwapArg *a = arg;
    radar_config_free(&a->app->cfg);
    a->app->cfg = *a->next;
}

/* A new radar mode (or code, stagger pair, rotation speed, carrier) from the Radar Settings: the
 * open plots are closed, the processing starts again with the new configuration. Returns the new
 * processor (the old one when the new configuration cannot be used). */
static Processor *change_config(App *app, Processor *proc, RadarConfig *next, int mode, const ProcessorCallbacks *cb,
                                int threads)
{
    apply_options(app, next);
    if (next->udpFormat == UDP_FORMAT_VBDR && next->nChUsed > 4) {
        fprintf(stderr, "Radar Settings: mode %d processes %d channels, the VBDR stream carries 4: not used\n", mode,
                next->nChUsed);
        radar_config_free(next);
        return proc;
    }
    processor_flush(proc, cb);
    const int logMask = app->log ? logger_mask(app->log) & LOG_BLOCKS : 0;
    processor_destroy(proc);
    SwapArg a = { app, next };
    logger_when_idle(app->log, swap_config, &a);      /* the logger reads app->cfg */
    proc = processor_create(&app->cfg, threads);
    if (!proc) return NULL;
    processor_set_log(proc, logMask);
    app->detArc_deg = app->cfg.rate_dps * app->cfg.meanPri_s * app->cfg.updateEvery;
    fprintf(stderr, "Radar mode %d in use: DDC %d samples, short %d / long %d range cells, CPI %d pulses, %.3g rpm\n",
            app->cfg.radarMode, app->cfg.nDdc, app->cfg.pulse[0].nRange, app->cfg.pulse[1].nRange, app->cfg.cpiLength,
            app->cfg.rate_dps * 60.0 / 360.0);
    return proc;
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
    udp_input_set_format(in, app->cfg.udpFormat);
    VbdrTiming vbdrTiming;                                 /* VBDR: time, scan and stagger of the PRIs */
    memset(&vbdrTiming, 0, sizeof vbdrTiming);
    /* real time (Linux): the processing thread gets its priority and CPUs before
     * the OpenMP workers exist, so they inherit them */
    char rtErr[512];
    int nCpus = 0, threads = app->opt.threads;
    const int rt = app->opt.rtPriority > 0 || (app->opt.cpus && app->opt.cpus[0]);
    if (rt) {
        const int r = rt_thread_setup(app->opt.rtPriority, app->opt.cpus, &nCpus, rtErr, sizeof rtErr);
        if (r == -2) {
            fprintf(stderr, "--cpus: %s\n", rtErr);
            udp_input_close(in);
            app->exitCode = 2;
            return 2;
        }
        if (r != 0) fprintf(stderr, "warning: real-time settings of the processing threads: %s\n", rtErr);
        if (threads <= 0 && nCpus > 0) threads = nCpus;     /* one thread per CPU of --cpus */
    }
    Processor *proc = processor_create(&app->cfg, threads);
    if (!proc) {
        fprintf(stderr, "cannot create the processor\n");
        udp_input_close(in);
        app->exitCode = 1;
        return 1;
    }
    ProcessorCallbacks cb = { on_detections, on_plot, on_antenna, app, on_track, on_profile, on_rdmap, on_block };
    if (app->opt.logDir) start_log(app, proc, app->opt.logDir, app->opt.logName, app->opt.logMask);
#ifdef RADAR_WITH_QT
    char detName[48], trackName[64];
    block_names(&app->cfg, detName, sizeof detName, trackName, sizeof trackName);
    double lastTiming = 0.0;
    VideoLayout video;                       /* raw video of the PPI (made when first shown) */
    memset(&video, 0, sizeof video);
    float *videoRow = NULL;
#endif
    double default_dB = 10.0 * log10(app->cfg.alpha);
    if (!isnan(app->opt.cfarThreshold_dB))
        set_cfar_threshold(app, proc, app->opt.cfarThreshold_dB);
    report_cfar(app, proc, default_dB, 0);
    /* live stream: a full queue drops PRIs; replay: the reader waits */
    const int live = app->opt.replayPath == NULL;
    /* the input thread inherits priority + 1 and its CPUs from this thread */
    const int rtInput = app->opt.rtPriority > 0 || (app->opt.inputCpus && app->opt.inputCpus[0]);
    if (rtInput) {
        const int r = rt_thread_setup(app->opt.rtPriority > 0 ? app->opt.rtPriority + 1 : 0, app->opt.inputCpus,
                                      NULL, rtErr, sizeof rtErr);
        if (r == -2) {
            fprintf(stderr, "--input-cpus: %s\n", rtErr);
            processor_destroy(proc);
            udp_input_close(in);
            app->exitCode = 2;
            return 2;
        }
        if (r != 0) fprintf(stderr, "warning: real-time settings of the input thread: %s\n", rtErr);
    }
    InputThread *input = input_thread_start(in, app->opt.queue, live);
    if (rtInput && rt_thread_setup(app->opt.rtPriority, app->opt.cpus && app->opt.cpus[0] ? app->opt.cpus : NULL,
                                   NULL, rtErr, sizeof rtErr) != 0 && app->opt.cpus)
        fprintf(stderr, "warning: real-time settings of the processing threads: %s\n", rtErr);
    if (app->opt.mlock && rt_lock_memory(rtErr, sizeof rtErr) != 0)
        fprintf(stderr, "warning: %s\n", rtErr);
    if (rt || rtInput || app->opt.mlock) {
        char prio[32] = "normal scheduling", prioIn[32] = "normal scheduling";
        if (app->opt.rtPriority > 0) {
            snprintf(prio, sizeof prio, "SCHED_FIFO %d", app->opt.rtPriority);
            snprintf(prioIn, sizeof prioIn, "SCHED_FIFO %d", app->opt.rtPriority + 1);
        }
        fprintf(stderr, "Real time: processing %s, %d thread(s), CPUs %s; input %s, CPUs %s; memory %s\n", prio,
                processor_threads(proc), app->opt.cpus ? app->opt.cpus : "all", prioIn,
                app->opt.inputCpus ? app->opt.inputCpus : app->opt.cpus ? app->opt.cpus : "all",
                app->opt.mlock ? "locked" : "not locked");
    }
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
    double lastLogStatus = 0;
    while (!g_stop) {
#ifdef RADAR_WITH_QT
        if (app->opt.ppi) {                                     /* Start / Stop logging in the tab Controls */
            int start = 0, mask = 0;
            char folder[1024], name[256];
            if (ppi_take_log_request(&start, folder, sizeof folder, name, sizeof name, &mask)) {
                if (start) start_log(app, proc, folder, name, mask);
                else stop_log(app, proc);
            }
            if (app->log && now_seconds() - lastLogStatus > 0.5) {    /* the state of the log */
                LoggerStats st;
                logger_stats(app->log, &st);
                ppi_post_log_status(1, logger_dir(app->log), logger_mask(app->log), &st, NULL);
                lastLogStatus = now_seconds();
            }
        }
#endif
        if (app->board) {                                       /* Radar Settings of the PPI */
            char text[300];
            while (board_control_take_note(app->board, text, sizeof text)) {
                fprintf(stderr, "%s\n", text);
                if (app->log) logger_note(app->log, text);
            }
            RadarConfig next;
            int mode = 0;
            if (board_control_take_config(app->board, &next, &mode, NULL)) {
                proc = change_config(app, proc, &next, mode, &cb, threads);
                if (!proc) {
                    fprintf(stderr, "cannot create the processor of radar mode %d\n", mode);
                    app->exitCode = 1;
                    break;
                }
                default_dB = 10.0 * log10(app->cfg.alpha);
                report_cfar(app, proc, default_dB, 1);
#ifdef RADAR_WITH_QT
                block_names(&app->cfg, detName, sizeof detName, trackName, sizeof trackName);
                free(videoRow);                                 /* the video layout of the new mode */
                videoRow = NULL;
                video_layout_free(&video);
                memset(&video, 0, sizeof video);
                if (app->opt.ppi) {
                    double v0, v1;
                    rd_axis(&app->cfg, &v0, &v1);
                    ppi_post_geometry(app->cfg.maxRange_m, v0, v1);
                }
#endif
            }
        }
        PriData *pri = NULL;
        int r = input_thread_next(input, &pri, 200);
        if (r < 0) break;                                      /* end of the replay file */
        if (r == 0) {
            if (app->opt.idleTimeout_s > 0 && lastData > 0 && now_seconds() - lastData > app->opt.idleTimeout_s) break;
            continue;
        }
        lastData = now_seconds();
        if (pri->format == UDP_FORMAT_VBDR) vbdr_complete_pri(&vbdrTiming, &app->cfg, pri);
#ifdef RADAR_WITH_QT
        double request_dB;
        const double tReq = now_seconds();
        const int request = app->opt.ppi && ppi_take_cfar_request(&request_dB);
        if (app->opt.ppi) app->display_ms += ms_since(tReq);
        if (request) {                                         /* Apply in the PPI */
            const double before_dB = 10.0 * log10(processor_cfar_alpha(proc));
            set_cfar_threshold(app, proc, request_dB);
            report_cfar(app, proc, default_dB, 1);
            if (app->log) {                                    /* in log_info.txt of the running log */
                char note[160];
                const double a = processor_cfar_alpha(proc), pfa = cfar_pfa(&app->cfg, a);
                snprintf(note, sizeof note, "CFAR threshold %.2f dB -> %.2f dB (Pfa %.1e) at pulse %u", before_dB,
                         10.0 * log10(a), pfa, (unsigned)pri->pulseNumber);
                logger_note(app->log, note);
            }
        }
#endif
#ifdef RADAR_WITH_QT
        if (app->opt.ppi && ppi_video_wanted() && pri->nCh >= app->cfg.nChUsed && pri->nSamples >= app->cfg.nDdc) {
            const double tVideo = now_seconds();             /* the video of the PPI (display time) */
            if (!videoRow && video_layout_init(&video, &app->cfg, app->cfg.maxRange_m) == 0)
                videoRow = malloc(sizeof(float) * (size_t)video.nBins);
            if (videoRow) {
                video_row(&video, pri, videoRow);
                ppi_post_video(videoRow, video.nBins, app->cfg.nChUsed, video.dR_m, pri->az_deg, pri->time_s);
            }
            app->display_ms += ms_since(tVideo);
        }
#endif
#ifdef RADAR_WITH_QT
        if (app->opt.ppi) {                                  /* what the PPI shows: profile, range-Doppler map */
            int beam = 0;
            const int kind = ppi_profile_wanted(&beam);
            processor_set_profile(proc, kind, beam);
            const int rd = ppi_rdmap_wanted(&beam);
            processor_set_rdmap(proc, rd, beam, 1.0 / app->cfg.rdRate_Hz);   /* cfg.display.rdRate_Hz */
        }
#endif
        if (live && !app->opt.noSkip)
            processor_set_slow_time(proc, input_thread_backlog(input) <= skipAbove);
        if (app->log) {                                        /* the DDC output of the PRI as received */
            const double tLog = now_seconds();
            logger_begin_pri(app->log, pri);
            app->log_ms += ms_since(tLog);
        }
        processor_process(proc, pri, &cb);
        if (app->log) logger_end_pri(app->log);
        app->havePri = 1;
        app->lastPulse = pri->pulseNumber;
        app->lastTime_s = pri->time_s;
        app->lastAz_deg = pri->az_deg;
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
    if (!proc) {
        stop_log(app, NULL);
        udp_input_close(in);
        return 1;
    }
    processor_flush(proc, &cb);
    stop_log(app, proc);
    print_timing(app, proc, in, queueDrops, now_seconds() - t0);
    if (app->cfg.udpFormat == UDP_FORMAT_VBDR) {
        UdpStats us;
        udp_input_stats(in, &us);
        fprintf(stderr, "VBDR stream: %llu stream(s) never arrived or were dropped, board \"Loss=\" %u\n",
                (unsigned long long)vbdrTiming.missing, (unsigned)us.boardLoss);
    }
#ifdef RADAR_WITH_QT
    if (app->opt.ppi) {
        ProcessorTiming tm;
        app_timing(app, proc, &tm);
        ppi_post_timing(&tm, detName, trackName);
    }
#endif
#ifdef RADAR_WITH_QT
    free(videoRow);
    video_layout_free(&video);
#endif
    fprintf(stderr, "Plots: %llu, CFAR detections after pruning: %llu\n",
            (unsigned long long)app->nPlots, (unsigned long long)app->nDetections);
    if (app->cfg.trackEnable) {
        char tn[40];
        tracker_name(&app->cfg, tn, sizeof tn);
        fprintf(stderr, "Tracks (%s): %llu updates, %llu confirmed\n", tn, (unsigned long long)app->nTrackUpdates,
                (unsigned long long)app->nConfirmed);
    }
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
    if (app.opt.format) app.cfg.udpFormat = !strcmp(app.opt.format, "vbdr") ? UDP_FORMAT_VBDR : UDP_FORMAT_RDDC;
    if (app.cfg.udpFormat == UDP_FORMAT_VBDR && app.cfg.nChUsed > 4) {
        fprintf(stderr, "configuration: the VBDR stream carries 4 channels, the configuration processes %d "
                "(cfg.rx.nChannelsUsed <= 4)\n", app.cfg.nChUsed);
        return 1;
    }
    apply_options(&app, &app.cfg);
    if (app.opt.quiet < 0) app.opt.quiet = !app.cfg.consolePlots;
    fprintf(stderr, "Configuration: %d element(s) (%d processed), %d beam(s), DDC %d samples, short %d / long %d range cells, "
            "CPI %d pulses, %s %d x %d row(s), update every %d pulse(s), plots need %d update(s)\n", app.cfg.nCh,
            app.cfg.nChUsed, app.cfg.nBeams, app.cfg.nDdc, app.cfg.pulse[0].nRange, app.cfg.pulse[1].nRange,
            app.cfg.cpiLength, app.cfg.dopplerMethod == DOPPLER_NUFFT ? "NUFFT" : "FFT", app.cfg.nfft, app.cfg.nRows,
            app.cfg.updateEvery, app.cfg.minHits);
    fprintf(stderr, "DDC stream: %s format, UDP port %d\n", app.cfg.udpFormat == UDP_FORMAT_VBDR ? "VBDR" : "RDDC",
            app.opt.port);
    fprintf(stderr, "Radar mode %d; DDC record %s: %d samples, the first at %.3f us", app.cfg.radarMode,
            app.cfg.ddcManual ? "manual" : "automatic", app.cfg.nDdc,
            app.cfg.ddcT0_s * 1e6);
    for (int i = 0; i < N_PULSES; i++)
        if (app.cfg.codeIndex[i] > 0) fprintf(stderr, "; %s code %d", i == PULSE_SHORT ? "short" : "long", app.cfg.codeIndex[i]);
    if (app.cfg.staggerSet > 0) fprintf(stderr, "; stagger pair %d", app.cfg.staggerSet);
    fprintf(stderr, "\n");
    if (app.cfg.scan2scan) {
        fprintf(stderr, "Scan-to-scan PRF (selected by the stagger index of the PRI header):");
        for (int j = 0; j < app.cfg.nPrf; j++)
            fprintf(stderr, " %.1f Hz (CPI %d, FFT %d)", 1.0 / app.cfg.prf[j].pri_s, app.cfg.prf[j].cpiLength, app.cfg.prf[j].nfft);
        fprintf(stderr, "\n");
    }
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &app.cfg.pulse[i];
        if (!pc->enable) {
            fprintf(stderr, "%s pulse: off\n", pulse_name(i));
            continue;
        }
        if (pc->fc_Hz > 0)
            fprintf(stderr, "%s pulse: carrier %.3f MHz, Doppler span %.1f m/s\n", pulse_name(i), pc->fc_Hz / 1e6,
                    pc->vrSpan_mps);
        if (pc->waveform == WAVE_PHASE_CODE)
            fprintf(stderr, "%s pulse: phase code of %d chips, %s filter of %d taps (%.2f ... %.2f km)\n", pulse_name(i),
                    pc->nChips, pc->compression == COMPRESSION_MISMATCH ? "mismatch" : "matched", pc->nTaps,
                    pc->range_m[0] / 1e3, pc->range_m[pc->nRange - 1] / 1e3);
        else
            fprintf(stderr, "%s pulse: LFM, matched filter of %d taps (%.2f ... %.2f km)\n", pulse_name(i), pc->nTaps,
                    pc->range_m[0] / 1e3, pc->range_m[pc->nRange - 1] / 1e3);
    }
    fprintf(stderr, "CFAR %s %s (window: range %d+%d, Doppler %d+%d guard+training cells per side), SNR noise: %s, ",
            cfar_type_name(app.cfg.cfarType), cfar_method_name(app.cfg.cfarMethod), app.cfg.guardRange, app.cfg.trainRange, app.cfg.guardDoppler,
            app.cfg.trainDoppler, app.cfg.snrNoise == SNR_NOISE_CFAR ? "CFAR window" : "median of the map");
    fprintf(stderr, "Clutter filter: %s (%d taps); peak interpolation: Doppler %s, elevation %s\n",
            app.cfg.mtiType == MTI_FIR ? "FIR" : "canceller", app.cfg.nMtiTaps, app.cfg.dopplerInterp ? "on" : "off",
            app.cfg.elevationInterp && app.cfg.dbfEnabled ? "on" : "off");
    if (app.cfg.trackEnable && app.cfg.trackMethod == TRACK_METHOD_PMBM)
        fprintf(stderr, "tracking: PMBM %d-D, IMM-UKF of %d models, %d global hypotheses\n", app.cfg.trackDim,
                app.cfg.pmbmNModels, app.cfg.pmbmMaxHyp);
    if (app.cfg.trackEnable && app.cfg.trackMethod == TRACK_METHOD_TOMHT)
        fprintf(stderr, "tracking: TOMHT %d-D, IMM-UKF of %d models, %d global hypotheses (k-best %d), %d-scan pruning\n",
                app.cfg.trackDim, app.cfg.pmbmNModels, app.cfg.tomhtMaxHyp, app.cfg.tomhtKBest, app.cfg.tomhtNScan);
    else if (app.cfg.trackEnable)
        fprintf(stderr, "tracking: EKF %d-D, %s model, GNN association\n", app.cfg.trackDim,
                app.cfg.trackModel == TRACK_CA ? "constant-acceleration" : "constant-velocity");
    else fprintf(stderr, "tracking off\n");

    app.plotsFile = open_csv(app.opt.plotsCsv,
                             "time_s,az_deg,range_m,el_deg,vr_mps,snr_dB,pulse,hits,beam,rangeCell,dopplerBin,pulseNumber");
    app.trackFile = open_csv(app.opt.tracksCsv, "id,confirmed,deleted,time_s,x_m,y_m,z_m,vx_mps,vy_mps,vz_mps,"
                                                "hits,misses,ax_mps2,ay_mps2,az_mps2,existence,pCV,pCA,pCT");
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
                             : !strcmp(app.opt.ppiShow, "both") ? PPI_SHOW_BOTH
                             : !strcmp(app.opt.ppiShow, "video") ? PPI_SHOW_VIDEO
                             : !strcmp(app.opt.ppiShow, "dbf") ? PPI_SHOW_DBF
                             : !strcmp(app.opt.ppiShow, "mf") ? PPI_SHOW_MF : PPI_SHOW_PLOTS);
        {
            double el[256];
            const int nB = app.cfg.nBeams < 256 ? app.cfg.nBeams : 256;
            for (int b = 0; b < nB; b++)
                el[b] = app.cfg.dbfEnabled ? app.cfg.tilt_deg + asin(app.cfg.beamU[b]) * 180.0 / 3.14159265358979323846 : NAN;
            ppi_set_beams(nB, el);
            if (app.opt.ppiBeam > nB) fprintf(stderr, "--ppi-beam %d: only %d beam(s), beam 1 is shown\n", app.opt.ppiBeam, nB);
            ppi_set_beam(app.opt.ppiBeam <= nB ? app.opt.ppiBeam - 1 : 0);
            ppi_set_tab(!strcmp(app.opt.ppiTab, "rd"));
            ppi_set_rd_plot(!strcmp(app.opt.ppiRdPlot, "range"));
            /* the empty range-Doppler tab before the first map: the Doppler axis of the configuration */
            double v0, v1;
            rd_axis(&app.cfg, &v0, &v1);
            ppi_set_rd_axes(v0, v1);
            start_board(&app);
            ppi_set_board_control(app.board);
            ppi_set_log_defaults(app.opt.logDir ? app.opt.logDir : "", app.opt.logName ? app.opt.logName : "", app.opt.logMask);
            ppi_set_screen_record(app.opt.ffmpeg, app.opt.recordAudio, app.opt.recordFps);
        }
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
    free(app.confirmedIds);
    if (app.trackFile) fclose(app.trackFile);
    board_control_stop(app.board);
    radar_config_free(&app.cfg);
    return app.exitCode;
}

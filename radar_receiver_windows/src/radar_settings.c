/*
 * radar_settings.c - the settings of the radar modes that the PPI can change (tab Controls,
 * Radar Settings) and the receiver configuration of a mode with such settings.
 *
 * radar_config.m is evaluated as by radar_export_config (mfile.c): once for the common settings
 * (cfg.board, cfg.fcMin_Hz / fcMax_Hz, cfg.modes) and once per mode of cfg.modes.available
 * ("cfg.modes.active = k; cfg = apply_radar_mode(cfg);"). A new mode / phase code / stagger pair /
 * rotation speed / carrier gives a new receiver configuration: the same export with these
 * settings as statements (radar_export_run of export_config.c) into a temporary file, loaded
 * with radar_config_load.
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  include <process.h>
#  define getpid _getpid
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <unistd.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mfile.h"
#include "radar_rx.h"

/* export_config.c compiled with RADAR_EXPORT_LIBRARY */
int radar_export_run(const char *mPath, const char *const *sets, int nSets, const char *outPath, char *err, size_t errLen);

static const char *g_pulseName[N_PULSES] = { "short", "long" };

static int same_text(const char *a, const char *b)
{
    for (; *a && *b; a++, b++)
        if ((*a | 32) != (*b | 32)) return 0;
    return *a == *b;
}

static const MValue *get(const MValue *cfg, const char *path)
{
    const MValue *v = mvalue_get(cfg, path);
    if (v && v->type == MV_NUM && v->rows * v->cols == 0) return NULL;   /* [] */
    return v;
}

static double num_or(const MValue *cfg, const char *path, double def)
{
    const MValue *v = get(cfg, path);
    return v && v->type == MV_NUM && v->rows * v->cols >= 1 ? v->re[0] : def;
}

/* The numbers of cfg.<path> (at most max); 0 when missing. */
static int list_of(const MValue *cfg, const char *path, double *out, int max)
{
    const MValue *v = get(cfg, path);
    if (!v || v->type != MV_NUM) return 0;
    int n = v->rows * v->cols;
    if (n > max) n = max;
    for (int i = 0; i < n; i++) out[i] = v->re[i];
    return n;
}

static int text_is(const MValue *cfg, const char *path, const char *text)
{
    const MValue *v = mvalue_get(cfg, path);
    return v && v->type == MV_STR && v->str && same_text(v->str, text);
}

static MValue *eval(const BoardInfo *bi, const char *const *extra, int nExtra, char *err, size_t errLen)
{
    const char *sets[32];
    int n = 0;
    for (int i = 0; i < bi->nSets && n < 24; i++) sets[n++] = bi->sets[i];
    for (int i = 0; i < nExtra && n < 32; i++) sets[n++] = extra[i];
    char e[1024];
    MValue *cfg = mfile_eval(bi->mPath, sets, n, e, sizeof e);
    if (!cfg) { snprintf(err, errLen, "%s", e); return NULL; }
    if (cfg->type != MV_STRUCT) {
        snprintf(err, errLen, "%s: %s", bi->mPath, cfg->type == MV_UNEVAL && cfg->str ? cfg->str : "does not return a struct");
        mvalue_free(cfg);
        return NULL;
    }
    return cfg;
}

static int carrier_MHz(const MValue *cfg, int i)
{
    char p[64];
    snprintf(p, sizeof p, "wave.%s.fc_Hz", g_pulseName[i]);
    return (int)lround(num_or(cfg, p, num_or(cfg, "fc", 0)) / 1e6);
}

/* The settings of one mode from its evaluated configuration. */
static void mode_info(const MValue *cfg, int mode, ModeInfo *m)
{
    memset(m, 0, sizeof *m);
    m->mode = mode;
    const MValue *nm = mvalue_get(cfg, "modes.name");
    snprintf(m->name, sizeof m->name, "%s", nm && nm->type == MV_STR && nm->str ? nm->str : "");
    double buf[64];
    int firstCode = 0;
    for (int i = 0; i < N_PULSES; i++) {
        char p[64];
        snprintf(p, sizeof p, "wave.%s.enable", g_pulseName[i]);
        m->enable[i] = num_or(cfg, p, 1) != 0;
        snprintf(p, sizeof p, "wave.%s.type", g_pulseName[i]);
        m->phaseCode[i] = m->enable[i] && text_is(cfg, p, "phasecode");
        if (!m->phaseCode[i]) continue;
        snprintf(p, sizeof p, "wave.%s.availableCodes", g_pulseName[i]);
        int n = list_of(cfg, p, buf, 64);
        snprintf(p, sizeof p, "wave.%s.codeIndex", g_pulseName[i]);
        const int code = (int)num_or(cfg, p, n > 0 ? buf[0] : 1);
        if (n == 0) { buf[0] = code; n = 1; }
        if (!firstCode) {                                              /* the first phase-coded sub-pulse */
            firstCode = code;
            for (int k = 0; k < n; k++) m->codes[k] = (int)buf[k];
            m->nCodes = n;
        } else {                                                       /* both: the codes of both lists */
            int k2 = 0;
            for (int k = 0; k < m->nCodes; k++) {
                int found = 0;
                for (int j = 0; j < n; j++) found |= (int)buf[j] == m->codes[k];
                if (found) m->codes[k2++] = m->codes[k];
            }
            m->nCodes = k2;
        }
    }
    m->defaults.code = firstCode;
    m->defaults.f1_MHz = m->enable[PULSE_SHORT] ? carrier_MHz(cfg, PULSE_SHORT) : 0;
    m->defaults.f2_MHz = m->enable[PULSE_LONG] ? carrier_MHz(cfg, PULSE_LONG) : 0;
    m->staggered = text_is(cfg, "prf.mode", "staggered");
    if (m->staggered) {
        int n = list_of(cfg, "prf.availableStaggerSets", buf, 64);
        if (n == 0) { buf[0] = 1; n = 1; }
        for (int k = 0; k < n; k++) m->staggerSets[k] = (int)buf[k];
        m->nStagger = n;
        m->defaults.stagger = (int)num_or(cfg, "prf.staggerSet", buf[0]);
    }
    int n = list_of(cfg, "ant.availableRpm", m->rpms, 64);
    if (n == 0) { m->rpms[0] = num_or(cfg, "ant.rpm", 6); n = 1; }
    m->nRpm = n;
    m->defaults.rpm = num_or(cfg, "ant.rpm", m->rpms[0]);
    for (int k = 0; k < 4; k++) m->stcOptions[k] = m->gainOptions[k] = 4;
    for (int k = 0; k < 4; k++) m->defaults.stc[k] = m->defaults.gain[k] = 1;
    if (list_of(cfg, "board.stcOptions", buf, 4) == 4) for (int k = 0; k < 4; k++) m->stcOptions[k] = (int)buf[k];
    if (list_of(cfg, "board.gainOptions", buf, 4) == 4) for (int k = 0; k < 4; k++) m->gainOptions[k] = (int)buf[k];
    if (list_of(cfg, "board.stc", buf, 4) == 4) for (int k = 0; k < 4; k++) m->defaults.stc[k] = (int)buf[k];
    if (list_of(cfg, "board.gain", buf, 4) == 4) for (int k = 0; k < 4; k++) m->defaults.gain[k] = (int)buf[k];
}

int board_info_load(const char *mPath, const char *const *sets, int nSets, BoardInfo *bi, char *err, size_t errLen)
{
    memset(bi, 0, sizeof *bi);
    snprintf(bi->mPath, sizeof bi->mPath, "%s", mPath);
    if (nSets > 16) nSets = 16;
    for (int i = 0; i < nSets; i++) snprintf(bi->sets[i], sizeof bi->sets[i], "%s", sets[i]);
    bi->nSets = nSets;
    MValue *cfg = eval(bi, NULL, 0, err, errLen);
    if (!cfg) return -1;
    double av[8];
    int nAv = list_of(cfg, "modes.available", av, 8);
    bi->activeMode = (int)num_or(cfg, "modes.active", 0);
    if (nAv == 0) { av[0] = bi->activeMode; nAv = 1; }
    bi->fcMin_MHz = (int)lround(num_or(cfg, "fcMin_Hz", 150e6) / 1e6);
    bi->fcMax_MHz = (int)lround(num_or(cfg, "fcMax_Hz", 405e6) / 1e6);
    bi->nOperationModes = (int)num_or(cfg, "board.nOperationModes", 1);
    bi->operationMode = (int)num_or(cfg, "board.operationMode", 1);
    bi->sendRpm = num_or(cfg, "board.sendRpm", 0) != 0;
    const MValue *h = mvalue_get(cfg, "board.host");
    snprintf(bi->host, sizeof bi->host, "%s", h && h->type == MV_STR && h->str ? h->str : "127.0.0.1");
    bi->port = (int)num_or(cfg, "board.port", 5000);
    bi->replyPort = (int)num_or(cfg, "board.replyPort", 5100);
    bi->timeout_s = num_or(cfg, "board.timeout_s", 0.5);
    bi->retries = (int)num_or(cfg, "board.retries", 1);
    mvalue_free(cfg);
    for (int k = 0; k < nAv; k++) {
        char stmt[128];
        snprintf(stmt, sizeof stmt, "cfg.modes.active = %d; cfg = apply_radar_mode(cfg);", (int)av[k]);
        const char *extra[1] = { stmt };
        MValue *c = eval(bi, extra, 1, err, errLen);
        if (!c) return -1;
        mode_info(c, (int)av[k], &bi->modes[bi->nModes++]);
        mvalue_free(c);
    }
    return 0;
}

const ModeInfo *board_mode_info(const BoardInfo *bi, int mode)
{
    for (int i = 0; bi && i < bi->nModes; i++)
        if (bi->modes[i].mode == mode) return &bi->modes[i];
    return NULL;
}

int board_settings_valid(const BoardInfo *bi, int mode, const ModeSettings *s, char *err, size_t errLen)
{
    const ModeInfo *m = board_mode_info(bi, mode);
    if (!m) { snprintf(err, errLen, "radar mode %d is not available", mode); return 0; }
    int ok = m->nCodes == 0 ? s->code == 0 : 0;
    for (int k = 0; k < m->nCodes; k++) ok |= s->code == m->codes[k];
    if (!ok) { snprintf(err, errLen, "mode %d: phase code %d is not available", mode, s->code); return 0; }
    ok = !m->staggered ? s->stagger == 0 : 0;
    for (int k = 0; k < m->nStagger; k++) ok |= s->stagger == m->staggerSets[k];
    if (!ok) { snprintf(err, errLen, "mode %d: stagger pattern %d is not available", mode, s->stagger); return 0; }
    ok = 0;
    for (int k = 0; k < m->nRpm; k++) ok |= fabs(s->rpm - m->rpms[k]) < 1e-9;
    if (!ok) { snprintf(err, errLen, "mode %d: %g rpm is not available", mode, s->rpm); return 0; }
    const int f[2] = { s->f1_MHz, s->f2_MHz };
    for (int i = 0; i < N_PULSES; i++) {
        if (!m->enable[i] ? f[i] != 0 : f[i] < bi->fcMin_MHz || f[i] > bi->fcMax_MHz) {
            snprintf(err, errLen, "mode %d: F%d = %d MHz (%s; %d ... %d MHz)", mode, i + 1, f[i],
                     m->enable[i] ? "outside the carrier range" : "the sub-pulse is off: 0", bi->fcMin_MHz, bi->fcMax_MHz);
            return 0;
        }
    }
    for (int k = 0; k < 4; k++) {
        if (s->stc[k] < 1 || s->stc[k] > m->stcOptions[k] || s->gain[k] < 1 || s->gain[k] > m->gainOptions[k]) {
            snprintf(err, errLen, "mode %d: STC / gain of channel %d outside 1 ... %d / %d", mode, k + 1, m->stcOptions[k],
                     m->gainOptions[k]);
            return 0;
        }
    }
    return 1;
}

int board_config_statements(const BoardInfo *bi, int mode, const ModeSettings *s, char st[][160], int max)
{
    const ModeInfo *m = board_mode_info(bi, mode);
    if (!m || max < 8) return 0;
    int n = 0;
    snprintf(st[n++], 160, "cfg.modes.active = %d; cfg = apply_radar_mode(cfg);", mode);
    for (int i = 0; i < N_PULSES; i++)
        if (m->phaseCode[i] && s->code > 0) snprintf(st[n++], 160, "cfg.wave.%s.codeIndex = %d;", g_pulseName[i], s->code);
    if (m->staggered && s->stagger > 0) snprintf(st[n++], 160, "cfg.prf.staggerSet = %d;", s->stagger);
    snprintf(st[n++], 160, "cfg.ant.rpm = %.17g;", s->rpm);
    /* a carrier of the mode file that is not a whole MHz stays as it is unless it was changed */
    const int f[2] = { s->f1_MHz, s->f2_MHz }, f0[2] = { m->defaults.f1_MHz, m->defaults.f2_MHz };
    for (int i = 0; i < N_PULSES; i++)
        if (m->enable[i] && f[i] > 0 && f[i] != f0[i]) snprintf(st[n++], 160, "cfg.wave.%s.fc_Hz = %de6;", g_pulseName[i], f[i]);
    return n;
}

int board_build_config(const BoardInfo *bi, int mode, const ModeSettings *s, RadarConfig *cfg, char *err, size_t errLen)
{
    const ModeInfo *m = board_mode_info(bi, mode);
    if (!m) { snprintf(err, errLen, "radar mode %d is not available", mode); return -1; }
    char st[8][160];
    const char *extra[8];
    const int n = board_config_statements(bi, mode, s, st, 8);
    for (int i = 0; i < n; i++) extra[i] = st[i];

    const char *sets[32];
    int nSets = 0;
    for (int i = 0; i < bi->nSets; i++) sets[nSets++] = bi->sets[i];
    for (int i = 0; i < n; i++) sets[nSets++] = extra[i];
    const char *tmp = getenv("TMPDIR");
    if (!tmp || !tmp[0]) tmp = getenv("TEMP");
    if (!tmp || !tmp[0]) tmp = getenv("TMP");
#ifdef _WIN32
    if (!tmp || !tmp[0]) tmp = ".";
#else
    if (!tmp || !tmp[0]) tmp = "/tmp";
#endif
    char path[1100];
    snprintf(path, sizeof path, "%s/radar_receiver_%d_mode.bin", tmp, (int)getpid());
    if (radar_export_run(bi->mPath, sets, nSets, path, err, errLen) != 0) return -1;
    const int r = radar_config_load(path, cfg, err, errLen);
    remove(path);
    return r;
}

void board_settings_of_config(const BoardInfo *bi, const RadarConfig *cfg, ModeSettings *s)
{
    const ModeInfo *m = board_mode_info(bi, cfg->radarMode);
    memset(s, 0, sizeof *s);
    if (m) *s = m->defaults;
    s->code = 0;
    for (int i = 0; i < N_PULSES && !s->code; i++) s->code = cfg->codeIndex[i];
    if (m && m->nCodes == 0) s->code = 0;
    s->stagger = cfg->staggerSet;
    if (cfg->rate_dps > 0) s->rpm = cfg->rate_dps * 60.0 / 360.0;
    if (m) {                                   /* the rpm of the list (the exact number) */
        for (int k = 0; k < m->nRpm; k++)
            if (fabs(m->rpms[k] - s->rpm) < 1e-6) s->rpm = m->rpms[k];
    }
    const int f[2] = { cfg->pulse[0].enable && cfg->pulse[0].fc_Hz > 0 ? (int)lround(cfg->pulse[0].fc_Hz / 1e6) : 0,
                       cfg->pulse[1].enable && cfg->pulse[1].fc_Hz > 0 ? (int)lround(cfg->pulse[1].fc_Hz / 1e6) : 0 };
    if (f[0] > 0) s->f1_MHz = f[0];
    if (f[1] > 0) s->f2_MHz = f[1];
}

/*
 * board_control.c - the Radar Settings of the PPI (tab Controls): a thread of its own sends the
 * settings to the VBDR board and reads them back (board_link.c), keeps the confirmed settings of
 * every mode while the program runs, and builds the receiver configuration of a new mode, phase
 * code, stagger pair, rotation speed or carrier (radar_settings.c) for the worker.
 *
 * Every setting is written ("CNTL W"), then read ("CNTL R"): it is confirmed when the read gives
 * the written value. Only confirmed settings are kept for a mode; a mode chosen again gets its
 * kept settings (or the defaults of its mode file) sent again. The rotation speed goes to the
 * board only with cfg.board.sendRpm (address 0x32, the simulator); else only the receiver uses it.
 * README "VBDR stream and board control".
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#endif

#include "thread_compat.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mat_file.h"
#include "mat_write.h"
#include "radar_rx.h"

#define QUEUE 32
#define NOTES 64

typedef struct {
    int    kind, item;
    double value;
    char   path[1024];
} Request;

struct BoardControl {
    BoardInfo    bi;
    BoardLink   *link;
    Mutex        mu;
    Cond         cv;
#ifdef _WIN32
    HANDLE       thread;
#else
    pthread_t    thread;
#endif
    int          stop;
    Request      q[QUEUE];
    int          qHead, qCount;
    BoardState   st;            /* the state for the GUI (mu)                 */
    int          have[8];       /* confirmed settings of mode k (control thread) */
    ModeSettings mem[8];
    int          haveCfg;       /* a configuration for the worker (mu)         */
    RadarConfig  cfg;
    int          cfgMode;
    ModeSettings cfgS;
    char         notes[NOTES][256];
    int          nHead, nCount;
    char         status[2048];  /* status being built (control thread)         */
    int          level;
};

/* ------------------------------------------------------------------------ */
/* Helpers                                                                   */
/* ------------------------------------------------------------------------ */

static void note(BoardControl *bc, const char *fmt, ...)
{
    char t[256];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(t, sizeof t, fmt, ap);
    va_end(ap);
    mutex_lock(&bc->mu);
    if (bc->nCount == NOTES) { bc->nHead = (bc->nHead + 1) % NOTES; bc->nCount--; }   /* drop the oldest */
    snprintf(bc->notes[(bc->nHead + bc->nCount) % NOTES], 256, "%s", t);
    bc->nCount++;
    mutex_unlock(&bc->mu);
}

static void say(BoardControl *bc, int level, const char *fmt, ...)
{
    char t[512];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(t, sizeof t, fmt, ap);
    va_end(ap);
    const size_t n = strlen(bc->status);
    snprintf(bc->status + n, sizeof bc->status - n, "%s%s", n ? "\n" : "", t);
    if (level > bc->level) bc->level = level;
}

/* Publishes the state for the GUI. */
static void publish(BoardControl *bc, int busy, int mode, const ModeSettings *cur, int opMode)
{
    mutex_lock(&bc->mu);
    bc->st.serial++;
    bc->st.busy = busy;
    if (mode >= 0) bc->st.mode = mode;
    if (cur) bc->st.cur = *cur;
    if (opMode > 0) bc->st.opMode = opMode;
    snprintf(bc->st.status, sizeof bc->st.status, "%s", bc->status);
    bc->st.level = bc->level;
    mutex_unlock(&bc->mu);
}

static int item_applies(const ModeInfo *m, int sendRpm, int item)
{
    switch (item) {
    case BOARD_CODE:    return m->nCodes > 0;
    case BOARD_F1:      return m->enable[PULSE_SHORT];
    case BOARD_F2:      return m->enable[PULSE_LONG];
    case BOARD_STAGGER: return m->staggered;
    case BOARD_RPM:     return sendRpm;
    case BOARD_EXT_TRIG: return 0;
    default:            return 1;
    }
}

/* The register value of an item of the settings. */
static int item_value(const BoardControl *bc, int mode, const ModeSettings *s, int opMode, int item)
{
    switch (item) {
    case BOARD_MODE:    return mode;
    case BOARD_CODE:    return s->code;
    case BOARD_F1:      return s->f1_MHz - bc->bi.fcMin_MHz;
    case BOARD_F2:      return s->f2_MHz - bc->bi.fcMin_MHz;
    case BOARD_STAGGER: return s->stagger;
    case BOARD_OPMODE:  return opMode;
    case BOARD_RPM:     return (int)lround(s->rpm);
    default:
        if (item >= BOARD_STC_MAIN && item <= BOARD_STC_J3) return s->stc[item - BOARD_STC_MAIN];
        if (item >= BOARD_GAIN_MAIN && item <= BOARD_GAIN_J3) return s->gain[item - BOARD_GAIN_MAIN];
        return 0;
    }
}

static void set_item(ModeSettings *s, int item, double v)
{
    switch (item) {
    case BOARD_CODE:    s->code = (int)v; break;
    case BOARD_F1:      s->f1_MHz = (int)v; break;
    case BOARD_F2:      s->f2_MHz = (int)v; break;
    case BOARD_STAGGER: s->stagger = (int)v; break;
    case BOARD_RPM:     s->rpm = v; break;
    default:
        if (item >= BOARD_STC_MAIN && item <= BOARD_STC_J3) s->stc[item - BOARD_STC_MAIN] = (int)v;
        if (item >= BOARD_GAIN_MAIN && item <= BOARD_GAIN_J3) s->gain[item - BOARD_GAIN_MAIN] = (int)v;
    }
}

static int processing_item(int item)
{
    return item == BOARD_MODE || item == BOARD_CODE || item == BOARD_F1 || item == BOARD_F2 || item == BOARD_STAGGER ||
           item == BOARD_RPM;
}

/* Text of the value of an item (GUI units). */
static void value_text(const ModeSettings *s, int mode, int opMode, int item, char *buf, size_t len)
{
    switch (item) {
    case BOARD_MODE:   snprintf(buf, len, "%d", mode); break;
    case BOARD_F1:     snprintf(buf, len, "%d MHz", s->f1_MHz); break;
    case BOARD_F2:     snprintf(buf, len, "%d MHz", s->f2_MHz); break;
    case BOARD_OPMODE: snprintf(buf, len, "%d", opMode); break;
    case BOARD_RPM:    snprintf(buf, len, "%g rpm", s->rpm); break;
    default: {
        ModeSettings t = *s;
        snprintf(buf, len, "%d", item_value(NULL, mode, &t, opMode, item));
    }
    }
}

/* Write, then read back. 1 confirmed, 0 no answer, 2 the board holds another value (*got), -1 socket error. */
static int send_item(BoardControl *bc, int item, int value, int *got)
{
    if (item == BOARD_RPM && !bc->bi.sendRpm) return 1;           /* only the receiver uses it */
    const int addr = board_item_address(item);
    uint32_t r = 0;
    int k = board_link_exchange(bc->link, 1, addr, value, bc->bi.timeout_s, bc->bi.retries, &r);
    if (k != 1) return k;
    k = board_link_exchange(bc->link, 0, addr, 0, bc->bi.timeout_s, bc->bi.retries, &r);
    if (k != 1) return k;
    *got = (int)(r & 0xFF);
    return *got == (value & 0xFF) ? 1 : 2;
}

static void report_failure(BoardControl *bc, int item, int result, int got)
{
    if (result == 0)
        say(bc, 2, "%s not set: the board did not respond (%s:%d)", board_item_name(item), bc->bi.host, bc->bi.port);
    else if (result == 2)
        say(bc, 2, "%s not set: the board reads back %d", board_item_name(item), got);
    else
        say(bc, 2, "%s not set: network error", board_item_name(item));
}

/* Builds the configuration of (mode, s) and hands it to the worker. Returns 1 on success. */
static int post_config(BoardControl *bc, int mode, const ModeSettings *s)
{
    RadarConfig cfg;
    char err[1024];
    memset(&cfg, 0, sizeof cfg);
    if (board_build_config(&bc->bi, mode, s, &cfg, err, sizeof err) != 0) {
        say(bc, 2, "the receiver keeps its processing: configuration of mode %d: %s", mode, err);
        return 0;
    }
    mutex_lock(&bc->mu);
    if (bc->haveCfg) radar_config_free(&bc->cfg);
    bc->cfg = cfg;
    bc->cfgMode = mode;
    bc->cfgS = *s;
    bc->haveCfg = 1;
    mutex_unlock(&bc->mu);
    return 1;
}

static const ModeSettings *settings_for(BoardControl *bc, int mode)
{
    const ModeInfo *m = board_mode_info(&bc->bi, mode);
    return bc->have[mode & 7] ? &bc->mem[mode & 7] : &m->defaults;
}

/* ------------------------------------------------------------------------ */
/* Requests                                                                  */
/* ------------------------------------------------------------------------ */

/* Sends all settings s of the mode (with the mode register when withMode); the confirmed ones are
 * kept, the others reported (and kept as the defaults). Then the receiver uses (mode, s). */
static void send_mode(BoardControl *bc, int mode, ModeSettings s, int withMode, const char *why)
{
    const ModeInfo *m = board_mode_info(&bc->bi, mode);
    const int oldMode = bc->st.mode;
    const ModeSettings old = bc->st.cur;
    int got = 0, r;
    if (withMode) {
        r = send_item(bc, BOARD_MODE, mode, &got);
        if (r != 1) {
            report_failure(bc, BOARD_MODE, r, got);
            say(bc, 2, "the radar stays in mode %d", oldMode);
            return;
        }
    }
    ModeSettings kept = s;
    int nFail = 0;
    for (int item = BOARD_CODE; item < BOARD_EXT_TRIG; item++) {
        if (item == BOARD_OPMODE || !item_applies(m, 1, item)) continue;
        r = send_item(bc, item, item_value(bc, mode, &s, 0, item), &got);
        if (r == 1) continue;
        report_failure(bc, item, r, got);
        set_item(&kept, item, item == BOARD_RPM ? m->defaults.rpm : item_value(bc, mode, &m->defaults, 0, item) +
                                                        (item == BOARD_F1 || item == BOARD_F2 ? bc->bi.fcMin_MHz : 0));
        nFail++;
        if (r == 0 && nFail >= 2) { say(bc, 2, "the other settings were not sent"); break; }
    }
    bc->mem[mode & 7] = kept;
    bc->have[mode & 7] = 1;
    if (!post_config(bc, mode, &s)) return;
    if (mode != oldMode) note(bc, "Radar Settings: radar mode %d -> %d (%s)", oldMode, mode, why);
    else note(bc, "Radar Settings: mode %d %s", mode, why);
    char t[256];
    snprintf(t, sizeof t, "mode %d: code %d, stagger %d, %g rpm, F1 %d MHz, F2 %d MHz, STC %d %d %d %d, gain %d %d %d %d", mode,
             s.code, s.stagger, s.rpm, s.f1_MHz, s.f2_MHz, s.stc[0], s.stc[1], s.stc[2], s.stc[3], s.gain[0], s.gain[1],
             s.gain[2], s.gain[3]);
    note(bc, "Radar Settings: %s", t);
    if (!nFail) say(bc, 0, "%s: radar mode %d %s set", why, mode, mode != oldMode ? "and its settings" : "settings");
    (void)old;
    publish(bc, 1, mode, &s, -1);
}

static void do_set(BoardControl *bc, int item, double value)
{
    const int mode = bc->st.mode;
    const int opMode = bc->st.opMode;
    if (item == BOARD_MODE) {
        const int k = (int)value;
        if (!board_mode_info(&bc->bi, k)) { say(bc, 2, "radar mode %d is not available", k); return; }
        send_mode(bc, k, *settings_for(bc, k), 1, "Radar Mode");
        return;
    }
    if (item == BOARD_EXT_TRIG) {
        int got = 0;
        const int r = send_item(bc, BOARD_EXT_TRIG, (int)value, &got);
        if (r == 1) say(bc, 0, "Ext Trig = %d set", (int)value);
        else report_failure(bc, item, r, got);
        return;
    }
    ModeSettings s = bc->st.cur;
    char before[48], after[48];
    value_text(&s, mode, opMode, item, before, sizeof before);
    int newOp = opMode;
    if (item == BOARD_OPMODE) newOp = (int)value;
    else set_item(&s, item, value);
    char err[256];
    if (item != BOARD_OPMODE && !board_settings_valid(&bc->bi, mode, &s, err, sizeof err)) {
        say(bc, 2, "%s", err);
        return;
    }
    if (item == BOARD_OPMODE && (newOp < 1 || newOp > bc->bi.nOperationModes)) {
        say(bc, 2, "Operation Mode %d is outside 1 ... %d", newOp, bc->bi.nOperationModes);
        return;
    }
    value_text(&s, mode, newOp, item, after, sizeof after);
    int got = 0;
    const int r = send_item(bc, item, item_value(bc, mode, &s, newOp, item), &got);
    if (r != 1) {
        report_failure(bc, item, r, got);
        return;
    }
    if (processing_item(item) && !post_config(bc, mode, &s)) return;
    if (item != BOARD_OPMODE) {
        bc->mem[mode & 7] = s;
        bc->have[mode & 7] = 1;
    }
    say(bc, 0, "%s = %s set%s", board_item_name(item), after,
        item == BOARD_RPM && !bc->bi.sendRpm ? " (receiver only)" : "");
    note(bc, "Radar Settings: %s %s -> %s (mode %d)", board_item_name(item), before, after, mode);
    publish(bc, 1, mode, &s, newOp);
}

static void do_check(BoardControl *bc)
{
    const int mode = bc->st.mode;
    const ModeInfo *m = board_mode_info(&bc->bi, mode);
    const ModeSettings s = bc->st.cur;
    int nDiff = 0, nRead = 0;
    char diff[1024] = "";
    for (int item = BOARD_MODE; item < BOARD_EXT_TRIG; item++) {
        if (!m || !item_applies(m, bc->bi.sendRpm, item)) continue;
        uint32_t r = 0;
        const int k = board_link_exchange(bc->link, 0, board_item_address(item), 0, bc->bi.timeout_s, bc->bi.retries, &r);
        if (k != 1) {
            if (nRead == 0) {
                say(bc, 2, "the board did not respond (%s:%d): the settings could not be read", bc->bi.host, bc->bi.port);
                return;
            }
            say(bc, 1, "%s: the board did not respond", board_item_name(item));
            continue;
        }
        nRead++;
        const int want = item_value(bc, mode, &s, bc->st.opMode, item) & 0xFF;
        if ((int)(r & 0xFF) != want) {
            const size_t n = strlen(diff);
            snprintf(diff + n, sizeof diff - n, "%s%s: board %u, receiver %d", nDiff ? "; " : "", board_item_name(item),
                     (unsigned)(r & 0xFF), want);
            nDiff++;
        }
    }
    if (nDiff) say(bc, 1, "the board differs from the receiver (nothing was sent): %s", diff);
    else say(bc, 0, "the board holds the settings of the receiver (%d registers read)", nRead);
}

static const char *g_columns = "mode code stagger rpm F1_MHz F2_MHz stcMain stcJ1 stcJ2 stcJ3 gainMain gainJ1 gainJ2 gainJ3";

static void do_export(BoardControl *bc, const char *path)
{
    double rows[8][14];
    int n = 0;
    for (int k = 0; k < 8; k++) {
        const ModeSettings *s = bc->have[k] ? &bc->mem[k] : k == bc->st.mode ? &bc->st.cur : NULL;
        if (!s) continue;
        const double r[14] = { k, s->code, s->stagger, s->rpm, s->f1_MHz, s->f2_MHz, s->stc[0], s->stc[1], s->stc[2], s->stc[3],
                               s->gain[0], s->gain[1], s->gain[2], s->gain[3] };
        memcpy(rows[n++], r, sizeof r);
    }
    double *cm = calloc((size_t)(n ? n : 1) * 14, sizeof(double));      /* column-major [n x 14] */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < 14; j++) cm[i + (size_t)j * n] = rows[i][j];
    MatWriter *w = mat_write_open(path, 0);
    if (!w) {
        free(cm);
        say(bc, 2, "cannot write %s", path);
        return;
    }
    const double active = bc->st.mode, op = bc->st.opMode;
    mat_write_double(w, "activeMode", &active, 1, 1);
    mat_write_double(w, "operationMode", &op, 1, 1);
    mat_write_double(w, "settings", cm, (size_t)n, 14);
    mat_write_char(w, "columns", g_columns);
    const int r = mat_write_close(w);
    free(cm);
    if (r != 0) say(bc, 2, "cannot write %s", path);
    else say(bc, 0, "settings of %d mode(s) exported to %s", n, path);
}

static void do_import(BoardControl *bc, const char *path)
{
    MatFile mf;
    char err[512];
    if (mat_load(path, &mf, err, sizeof err) != 0) {
        say(bc, 2, "import: %s", err);
        return;
    }
    const MatVar *set = mat_get(&mf, "settings"), *op = mat_get(&mf, "operationMode");
    int ok = set && set->re && (set->cols == 14 || set->numel == 0) && op && op->re && op->numel == 1;
    if (!ok) {
        mat_free(&mf);
        say(bc, 2, "import: %s is not a settings file (variables settings [n x 14], operationMode)", path);
        return;
    }
    const int nOp = (int)op->re[0];
    if (op->re[0] != nOp || nOp < 1 || nOp > bc->bi.nOperationModes) {
        mat_free(&mf);
        say(bc, 2, "import rejected: Operation Mode %g is outside 1 ... %d", op->re[0], bc->bi.nOperationModes);
        return;
    }
    int have[8] = { 0 };
    ModeSettings mem[8];
    const size_t n = set->numel ? set->rows : 0;
    for (size_t i = 0; i < n; i++) {
        double r[14];
        for (int j = 0; j < 14; j++) r[j] = set->re[i + (size_t)j * n];
        const int mode = (int)r[0];
        ModeSettings s;
        s.code = (int)r[1]; s.stagger = (int)r[2]; s.rpm = r[3]; s.f1_MHz = (int)r[4]; s.f2_MHz = (int)r[5];
        for (int k = 0; k < 4; k++) { s.stc[k] = (int)r[6 + k]; s.gain[k] = (int)r[10 + k]; }
        int whole = r[0] == mode && mode >= 0 && mode <= 7;
        for (int j = 1; j < 14; j++) whole &= j == 3 || r[j] == floor(r[j]);
        if (!whole || !board_settings_valid(&bc->bi, mode, &s, err, sizeof err)) {
            mat_free(&mf);
            say(bc, 2, "import rejected (row %d): %s", (int)i + 1, whole ? err : "not whole numbers / mode 0 ... 7");
            return;
        }
        have[mode] = 1;
        mem[mode] = s;
    }
    mat_free(&mf);
    memcpy(bc->have, have, sizeof have);
    for (int k = 0; k < 8; k++) if (have[k]) bc->mem[k] = mem[k];
    say(bc, 0, "imported %s: %d mode(s)", path, (int)n);
    note(bc, "Radar Settings: imported %s", path);
    int got = 0;
    const int r = send_item(bc, BOARD_OPMODE, nOp, &got);
    if (r == 1) {
        if (nOp != bc->st.opMode) note(bc, "Radar Settings: Operation Mode %d -> %d", bc->st.opMode, nOp);
        publish(bc, 1, -1, NULL, nOp);
    } else {
        report_failure(bc, BOARD_OPMODE, r, got);
    }
    const int mode = bc->st.mode;
    send_mode(bc, mode, *settings_for(bc, mode), 0, "Import");
}

static void run_request(BoardControl *bc, const Request *q)
{
    bc->status[0] = 0;
    bc->level = 0;
    publish(bc, 1, -1, NULL, -1);
    const int mode = bc->st.mode;
    switch (q->kind) {
    case BOARD_REQ_SET:
        do_set(bc, q->item, q->value);
        break;
    case BOARD_REQ_RESET_MODE:
        bc->have[mode & 7] = 0;
        send_mode(bc, mode, board_mode_info(&bc->bi, mode)->defaults, 0, "Reset defaults of this mode");
        break;
    case BOARD_REQ_RESET_ALL:
        memset(bc->have, 0, sizeof bc->have);
        send_mode(bc, mode, board_mode_info(&bc->bi, mode)->defaults, 0, "Reset defaults of all modes");
        break;
    case BOARD_REQ_EXPORT:
        do_export(bc, q->path);
        break;
    case BOARD_REQ_IMPORT:
        do_import(bc, q->path);
        break;
    case BOARD_REQ_CHECK:
        do_check(bc);
        break;
    }
    publish(bc, 0, -1, NULL, -1);
}

#ifdef _WIN32
static DWORD WINAPI control_main(LPVOID arg)
#else
static void *control_main(void *arg)
#endif
{
    BoardControl *bc = arg;
    for (;;) {
        Request q;
        mutex_lock(&bc->mu);
        while (!bc->stop && bc->qCount == 0) cond_wait_ms(&bc->cv, &bc->mu, 200);
        if (bc->stop) { mutex_unlock(&bc->mu); break; }
        q = bc->q[bc->qHead];
        bc->qHead = (bc->qHead + 1) % QUEUE;
        bc->qCount--;
        mutex_unlock(&bc->mu);
        run_request(bc, &q);
    }
    return 0;
}

/* ------------------------------------------------------------------------ */
/* Public functions                                                          */
/* ------------------------------------------------------------------------ */

BoardControl *board_control_start(const BoardInfo *bi, int mode, const ModeSettings *cur, char *err, size_t errLen)
{
    if (!board_mode_info(bi, mode)) {
        snprintf(err, errLen, "the radar mode %d of the configuration is not in cfg.modes.available", mode);
        return NULL;
    }
    BoardControl *bc = calloc(1, sizeof *bc);
    if (!bc) { snprintf(err, errLen, "out of memory"); return NULL; }
    bc->bi = *bi;
    bc->link = board_link_open(bi->host, bi->port, bi->replyPort, err, errLen);
    if (!bc->link) { free(bc); return NULL; }
    mutex_init(&bc->mu);
    cond_init(&bc->cv);
    bc->st.mode = mode;
    bc->st.cur = *cur;
    bc->st.opMode = bi->operationMode;
    bc->st.serial = 1;
    snprintf(bc->st.status, sizeof bc->st.status, "reading the settings of the board ...");
#ifdef _WIN32
    bc->thread = CreateThread(NULL, 0, control_main, bc, 0, NULL);
    const int ok = bc->thread != NULL;
#else
    const int ok = pthread_create(&bc->thread, NULL, control_main, bc) == 0;
#endif
    if (!ok) {
        board_link_close(bc->link);
        free(bc);
        snprintf(err, errLen, "cannot start the board control thread");
        return NULL;
    }
    return bc;
}

const BoardInfo *board_control_info(const BoardControl *bc)
{
    return bc ? &bc->bi : NULL;
}

int board_control_request(BoardControl *bc, int kind, int item, double value, const char *path)
{
    if (!bc) return -1;
    mutex_lock(&bc->mu);
    if (bc->qCount == QUEUE) { mutex_unlock(&bc->mu); return -1; }
    Request *q = &bc->q[(bc->qHead + bc->qCount) % QUEUE];
    q->kind = kind;
    q->item = item;
    q->value = value;
    snprintf(q->path, sizeof q->path, "%s", path ? path : "");
    bc->qCount++;
    bc->st.busy = 1;
    bc->st.serial++;
    cond_signal(&bc->cv);
    mutex_unlock(&bc->mu);
    return 0;
}

void board_control_state(BoardControl *bc, BoardState *st)
{
    mutex_lock(&bc->mu);
    *st = bc->st;
    if (bc->qCount) st->busy = 1;
    mutex_unlock(&bc->mu);
}

int board_control_take_config(BoardControl *bc, RadarConfig *cfg, int *mode, ModeSettings *s)
{
    if (!bc) return 0;
    mutex_lock(&bc->mu);
    const int have = bc->haveCfg;
    if (have) {
        *cfg = bc->cfg;
        if (mode) *mode = bc->cfgMode;
        if (s) *s = bc->cfgS;
        bc->haveCfg = 0;
        memset(&bc->cfg, 0, sizeof bc->cfg);
    }
    mutex_unlock(&bc->mu);
    return have;
}

int board_control_take_note(BoardControl *bc, char *text, size_t len)
{
    if (!bc) return 0;
    mutex_lock(&bc->mu);
    const int have = bc->nCount > 0;
    if (have) {
        snprintf(text, len, "%s", bc->notes[bc->nHead]);
        bc->nHead = (bc->nHead + 1) % NOTES;
        bc->nCount--;
    }
    mutex_unlock(&bc->mu);
    return have;
}

void board_control_stop(BoardControl *bc)
{
    if (!bc) return;
    mutex_lock(&bc->mu);
    bc->stop = 1;
    cond_signal(&bc->cv);
    mutex_unlock(&bc->mu);
#ifdef _WIN32
    WaitForSingleObject(bc->thread, INFINITE);
    CloseHandle(bc->thread);
#else
    pthread_join(bc->thread, NULL);
#endif
    if (bc->haveCfg) radar_config_free(&bc->cfg);
    board_link_close(bc->link);
    mutex_destroy(&bc->mu);
    cond_destroy(&bc->cv);
    free(bc);
}

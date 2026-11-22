/*
 * test_board.c - the VBDR board commands and the Radar Settings (board_link.c, radar_settings.c,
 * board_control.c):
 *   1) the texts of the commands and answers (the slides) and their parsing
 *   2) the simulator as the board (board_sim_*) and the client (board_link_*) over localhost
 *   3) the settings of the modes of radar_config.m (argv[1])
 *   4) the control thread with the simulated board: check at the start, a phase code, a mode
 *      change (new receiver configuration), memory of a mode, export / import, reset, and a
 *      board that does not answer; grouped auxiliary STC/gain with per-register feedback
 * Prints PASS or FAIL and returns 0 on success.
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  include <windows.h>
static void sleep_ms(int ms) { Sleep((DWORD)ms); }
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <time.h>
static void sleep_ms(int ms)
{
    struct timespec ts = { ms / 1000, (long)(ms % 1000) * 1000000L };
    nanosleep(&ts, NULL);
}
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "radar_rx.h"

static int g_failures = 0;
#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);               \
            g_failures++;                                                        \
        }                                                                        \
    } while (0)

/* Observe the control thread's exchanges and inject a fault at a specific
 * register. Successful requests still run through the UDP client and board.
 * The test sets and reads this trace only while the control thread is idle. */
static int g_writes[256], g_reads[256], g_commands, g_sequence[12];
static int g_timeoutAddr = -1, g_heldAddr = -1;
static int test_exchange(BoardLink *link, int write, int address, int value,
                         double timeout, int retries, uint32_t *returned)
{
    if (g_commands < 12) g_sequence[g_commands] = write ? address : -address;
    g_commands++;
    if (write) {
        g_writes[address]++;
        if (address == g_timeoutAddr) return 0;
        if (address == g_heldAddr) value = 1;   /* this register stays at 1 */
    } else {
        g_reads[address]++;
    }
    return board_link_exchange(link, write, address, value, timeout, retries, returned);
}
#define board_link_exchange test_exchange
#include "../src/board_control.c"
#undef board_link_exchange

static BoardSim *g_sim;
static int g_simOn = 1;

/* Waits until the control thread is idle (the simulated board answering meanwhile). */
static void wait_idle(BoardControl *bc, BoardState *st)
{
    for (int i = 0; i < 2000; i++) {
        if (g_simOn) board_sim_poll(g_sim, 5, NULL);
        else sleep_ms(5);
        board_control_state(bc, st);
        if (!st->busy) return;
    }
    printf("timeout of the board control\n");
    g_failures++;
}

static void set_aux(BoardControl *bc, BoardState *st, int item, double value, int expectedWrites)
{
    memset(g_writes, 0, sizeof g_writes);
    memset(g_reads, 0, sizeof g_reads);
    g_commands = 0;
    CHECK(board_control_request(bc, BOARD_REQ_SET_AUX, item, value, NULL) == 0);
    wait_idle(bc, st);
    int expectedCommands = 0;
    for (int k = 0; k < 3; k++) {
        const int addr = board_item_address(item + k);
        const int expectedReads = expectedWrites && addr != g_timeoutAddr;
        CHECK(g_writes[addr] == expectedWrites && g_reads[addr] == expectedReads);
        expectedCommands += expectedWrites + expectedReads;
    }
    CHECK(g_commands == expectedCommands);             /* no command to Main or the other group */
}

static void check_channels(const int state[4], int mainItem, int j1, int j2, int j3)
{
    const int expected[] = {j1, j2, j3};
    for (int k = 1; k <= 3; k++) {
        CHECK(state[k] == expected[k - 1]);
        CHECK(board_sim_register(g_sim, board_item_address(mainItem + k)) == expected[k - 1]);
    }
}

int main(int argc, char **argv)
{
    const char *mPath = argc > 1 ? argv[1] : "../radar_config.m";
    char buf[160], err[1024];
    int w, a, v;
    uint32_t r;

    /* 1) texts */
    board_format_command(buf, sizeof buf, 1, 0x02, 0x1C);
    CHECK(!strcmp(buf, "CNTL W 0x02 0x1C"));
    board_format_command(buf, sizeof buf, 0, 0x02, 0x55);
    CHECK(!strcmp(buf, "CNTL R 0x02 0x00"));
    board_format_reply(buf, sizeof buf, 1, 0x02, 0x1C, 0);
    CHECK(!strcmp(buf, "CNTL W 0x02 0x1C Returned Value = 0x00000000"));
    board_format_reply(buf, sizeof buf, 0, 0x02, 0, 0x1C);
    CHECK(!strcmp(buf, "CNTL R 0x02 0x00 Returned Value = 0x0000001C"));
    CHECK(board_parse_reply(buf, 0, 0x02, 0, &r) && r == 0x1C);
    CHECK(!board_parse_reply(buf, 0, 0x03, 0, &r));
    CHECK(board_parse_command("CNTL W 0x2D 0x03", &w, &a, &v) && w == 1 && a == 0x2D && v == 3);
    CHECK(board_parse_command(" CNTL R 0x00 0x00\r\n", &w, &a, &v) && w == 0 && a == 0);
    CHECK(!board_parse_command("CNTL X 0x00 0x00", &w, &a, &v));
    CHECK(!board_parse_command("CNTL W 0x00", &w, &a, &v));
    CHECK(!board_parse_command("CNTL W 0x00 0x01 junk", &w, &a, &v));
    CHECK(board_item_address(BOARD_STAGGER) == 0x05 && board_item_address(BOARD_OPMODE) == 0x28 &&
          board_item_address(BOARD_GAIN_J3) == 0x30 && board_item_address(BOARD_RPM) == 0x32);
    CHECK(!strcmp(board_item_name(BOARD_STAGGER), "Stagger Pattern") && !strcmp(board_item_name(BOARD_STC_MAIN), "STC Main"));

    /* 2) localhost */
    const int boardPort = 50211, replyPort = 50212;
    uint8_t regs[256] = { 0 };
    regs[0x28] = 1;
    g_sim = board_sim_open(boardPort, "127.0.0.1", replyPort, regs, 0, err, sizeof err);
    CHECK(g_sim != NULL);
    if (!g_sim) { printf("board_sim_open: %s\nFAIL\n", err); return 1; }
    {
        BoardLink *l = board_link_open("127.0.0.1", boardPort, replyPort + 10, err, sizeof err);
        CHECK(l != NULL);
        board_sim_close(g_sim);                       /* the answers go to replyPort: reopen for this client */
        g_sim = board_sim_open(boardPort, "127.0.0.1", replyPort + 10, regs, 0, err, sizeof err);
        /* the client waits in its own call: the simulated board answers from a second attempt on */
        CHECK(board_link_exchange(l, 1, 0x2D, 3, 0.05, 0, &r) == 0);   /* nobody polled: no answer */
        board_sim_poll(g_sim, 0, NULL);                                 /* answers the stale command */
        uint32_t written[8] = { 0 };
        board_sim_poll(g_sim, 0, written);
        CHECK(board_sim_register(g_sim, 0x2D) == 3);
        board_link_close(l);
        board_sim_close(g_sim);
        g_sim = board_sim_open(boardPort, "127.0.0.1", replyPort, regs, 0, err, sizeof err);
        CHECK(g_sim != NULL);
    }

    /* 3) the modes of radar_config.m */
    BoardInfo bi;
    const char *sets[3] = { "cfg.board.host = '127.0.0.1';", "cfg.board.port = 50211;", "cfg.board.replyPort = 50212;" };
    if (board_info_load(mPath, sets, 3, &bi, err, sizeof err) != 0) {
        printf("board_info_load: %s\nFAIL\n", err);
        return 1;
    }
    CHECK(bi.nModes == 3 && bi.modes[0].mode == 0 && bi.modes[2].mode == 2 && bi.port == 50211);
    CHECK(bi.fcMin_MHz == 150 && bi.fcMax_MHz == 200 && bi.nOperationModes == 4);
    const ModeInfo *m0 = board_mode_info(&bi, 0), *m1 = board_mode_info(&bi, 1);
    CHECK(m0 && m0->nCodes == 4 && m0->defaults.code == 1 && m0->defaults.f1_MHz == 170 && m0->nRpm == 2);
    CHECK(m1 && m1->nCodes == 0 && m1->defaults.code == 0);           /* LFM: no phase code */
    CHECK(m0 && m0->stcOptions[0] == 4 && m0->defaults.gain[3] == 1);
    ModeSettings s = m0->defaults;
    CHECK(board_settings_valid(&bi, 0, &s, err, sizeof err));
    s.code = 9;
    CHECK(!board_settings_valid(&bi, 0, &s, err, sizeof err));
    s = m0->defaults;
    s.f1_MHz = 220;
    CHECK(!board_settings_valid(&bi, 0, &s, err, sizeof err));

    /* 4) the control thread */
    RadarConfig cfg0;
    memset(&cfg0, 0, sizeof cfg0);
    CHECK(board_build_config(&bi, 0, &m0->defaults, &cfg0, err, sizeof err) == 0);
    ModeSettings cur;
    board_settings_of_config(&bi, &cfg0, &cur);
    CHECK(cur.code == 1 && cur.f1_MHz == 170 && cur.f2_MHz == 170 && fabs(cur.rpm - 6) < 1e-9);
    bi.modes[0].stcOptions[2] = 2;    /* J1 and J3 accept 3, J2 does not */
    BoardControl *bc = board_control_start(&bi, 0, &cur, err, sizeof err);
    CHECK(bc != NULL);
    if (!bc) { printf("board_control_start: %s\nFAIL\n", err); return 1; }
    BoardState st;
    RadarConfig next;
    int mode;
    char text[300];

    board_sim_set_register(g_sim, 0x01, 1);              /* the board as the receiver, but STC Main 2 */
    board_sim_set_register(g_sim, 0x02, 20);
    board_sim_set_register(g_sim, 0x03, 20);
    board_sim_set_register(g_sim, 0x29, 2);
    for (int k = 0x2A; k <= 0x30; k++) board_sim_set_register(g_sim, k, 1);
    board_control_request(bc, BOARD_REQ_CHECK, 0, 0, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 1 && strstr(st.status, "STC Main: board 2, receiver 1"));
    CHECK(!board_control_take_config(bc, &next, &mode, NULL));    /* the check sends nothing */

    board_control_request(bc, BOARD_REQ_SET, BOARD_CODE, 3, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 0 && st.cur.code == 3 && board_sim_register(g_sim, 0x01) == 3);
    CHECK(board_control_take_config(bc, &next, &mode, NULL) && mode == 0 && next.codeIndex[0] == 3);
    radar_config_free(&next);

    const double invalid[] = {3, 1.5, NAN, 256};
    for (size_t i = 0; i < sizeof invalid / sizeof invalid[0]; i++) {
        set_aux(bc, &st, BOARD_STC_J1, invalid[i], 0);
        CHECK(st.level == 2);
    }
    check_channels(st.cur.stc, BOARD_STC_MAIN, 1, 1, 1);

    g_timeoutAddr = board_item_address(BOARD_STC_J1);
    set_aux(bc, &st, BOARD_STC_J1, 2, 1);
    CHECK(st.level == 2 && strstr(st.status, "STC J1"));
    check_channels(st.cur.stc, BOARD_STC_MAIN, 1, 2, 2);
    g_timeoutAddr = -1;
    set_aux(bc, &st, BOARD_STC_J1, 2, 1);
    CHECK(st.level == 0);
    check_channels(st.cur.stc, BOARD_STC_MAIN, 2, 2, 2);
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 1, 1, 1);

    g_commands = 0;
    CHECK(board_control_request(bc, BOARD_REQ_SET_AUX, BOARD_STC_J1, 2, NULL) == 0);
    CHECK(board_control_request(bc, BOARD_REQ_SET_AUX, BOARD_GAIN_J1, 2, NULL) == 0);
    wait_idle(bc, &st);
    CHECK(st.level == 0 && g_commands == 12);
    for (int i = 0; i < 6; i++) {
        const int first = i < 3 ? BOARD_STC_J1 : BOARD_GAIN_J1;
        const int addr = board_item_address(first + i % 3);
        CHECK(g_sequence[2 * i] == addr && g_sequence[2 * i + 1] == -addr);
    }
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 2, 2, 2);
    g_heldAddr = board_item_address(BOARD_GAIN_J2);
    set_aux(bc, &st, BOARD_GAIN_J1, 3, 1);
    CHECK(st.level == 2 && strstr(st.status, "Gain J2"));
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 3, 1, 3);      /* observed J2 differs from both 2 and 3 */
    g_heldAddr = -1;
    set_aux(bc, &st, BOARD_GAIN_J1, 3, 1);
    CHECK(st.level == 0);
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 3, 3, 3);
    check_channels(st.cur.stc, BOARD_STC_MAIN, 2, 2, 2);
    CHECK(st.cur.stc[0] == 1 && board_sim_register(g_sim, 0x29) == 2);
    CHECK(st.cur.gain[0] == 1 && board_sim_register(g_sim, 0x2D) == 1);
    CHECK(!board_control_take_config(bc, &next, &mode, NULL));

    board_control_request(bc, BOARD_REQ_SET, BOARD_GAIN_J2, 4, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 0 && st.cur.gain[2] == 4 && board_sim_register(g_sim, 0x2F) == 4);
    CHECK(!board_control_take_config(bc, &next, &mode, NULL));    /* the processing does not change */

    board_control_request(bc, BOARD_REQ_SET, BOARD_F2, 175, NULL);
    wait_idle(bc, &st);
    CHECK(st.cur.f2_MHz == 175 && board_sim_register(g_sim, 0x03) == 25);
    CHECK(board_control_request(bc, BOARD_REQ_SET, BOARD_F1, 1000, NULL) == -2);   /* above fcMax: refused */
    CHECK(board_control_request(bc, BOARD_REQ_SET, BOARD_F2, 10, NULL) == -2);     /* below fcMin: refused */
    CHECK(board_control_take_config(bc, &next, &mode, NULL) && fabs(next.pulse[PULSE_LONG].fc_Hz - 175e6) < 1);
    radar_config_free(&next);

    board_control_request(bc, BOARD_REQ_SET, BOARD_MODE, 2, NULL);           /* long-range mode */
    wait_idle(bc, &st);
    CHECK(st.mode == 2 && board_sim_register(g_sim, 0x00) == 2 && st.cur.code == 0);
    CHECK(board_control_take_config(bc, &next, &mode, NULL) && mode == 2 && next.radarMode == 2 && next.nDdc == 3556);
    radar_config_free(&next);

    board_control_request(bc, BOARD_REQ_SET, BOARD_MODE, 0, NULL);           /* back: the kept settings */
    wait_idle(bc, &st);
    CHECK(st.mode == 0 && st.cur.code == 3 && st.cur.gain[2] == 4 && st.cur.f2_MHz == 175);
    CHECK(board_sim_register(g_sim, 0x01) == 3);
    check_channels(st.cur.stc, BOARD_STC_MAIN, 2, 2, 2);
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 3, 4, 3);
    CHECK(board_control_take_config(bc, &next, &mode, NULL) && mode == 0 && next.codeIndex[1] == 3);
    radar_config_free(&next);

    board_control_request(bc, BOARD_REQ_SET, BOARD_OPMODE, 3, NULL);
    wait_idle(bc, &st);
    CHECK(st.opMode == 3 && board_sim_register(g_sim, 0x28) == 3);

    const char *file = "test_board_settings.mat";
    board_control_request(bc, BOARD_REQ_EXPORT, 0, 0, file);
    wait_idle(bc, &st);
    CHECK(st.level == 0);
    board_control_request(bc, BOARD_REQ_RESET_ALL, 0, 0, NULL);
    wait_idle(bc, &st);
    CHECK(st.mode == 0 && st.cur.code == 1 && st.cur.gain[2] == 1 && board_sim_register(g_sim, 0x2F) == 1);
    while (board_control_take_config(bc, &next, &mode, NULL)) radar_config_free(&next);
    board_control_request(bc, BOARD_REQ_IMPORT, 0, 0, file);
    wait_idle(bc, &st);
    CHECK(st.level == 0 && st.cur.code == 3 && st.opMode == 3);
    check_channels(st.cur.stc, BOARD_STC_MAIN, 2, 2, 2);
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 3, 4, 3);
    while (board_control_take_config(bc, &next, &mode, NULL)) radar_config_free(&next);
    remove(file);
    board_control_request(bc, BOARD_REQ_IMPORT, 0, 0, "no_such_file.mat");
    wait_idle(bc, &st);
    CHECK(st.level == 2 && st.cur.code == 3);

    g_timeoutAddr = board_item_address(BOARD_GAIN_J2);
    board_control_request(bc, BOARD_REQ_RESET_MODE, 0, 0, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 2 && strstr(st.status, "Gain J2"));
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 1, 4, 1);
    CHECK(board_control_take_config(bc, &next, &mode, &cur));
    CHECK(cur.gain[1] == 1 && cur.gain[2] == 4 && cur.gain[3] == 1);
    radar_config_free(&next);
    g_timeoutAddr = -1;
    board_control_request(bc, BOARD_REQ_SET, BOARD_MODE, 2, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 0 && st.mode == 2);
    board_control_request(bc, BOARD_REQ_SET, BOARD_MODE, 0, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 0 && st.mode == 0);
    check_channels(st.cur.gain, BOARD_GAIN_MAIN, 1, 4, 1);       /* remembered across modes */
    while (board_control_take_config(bc, &next, &mode, NULL)) radar_config_free(&next);

    g_simOn = 0;                                         /* the board does not answer */
    board_control_request(bc, BOARD_REQ_SET, BOARD_STC_J1, 2, NULL);
    wait_idle(bc, &st);
    CHECK(st.level == 2 && strstr(st.status, "did not respond") && st.cur.stc[1] == 1);
    board_control_request(bc, BOARD_REQ_SET, BOARD_MODE, 1, NULL);
    wait_idle(bc, &st);
    CHECK(st.mode == 0 && strstr(st.status, "stays in mode 0"));
    CHECK(!board_control_take_config(bc, &next, &mode, NULL));

    int nNotes = 0;
    while (board_control_take_note(bc, text, sizeof text)) nNotes++;
    CHECK(nNotes >= 6);
    board_control_stop(bc);
    board_sim_close(g_sim);
    radar_config_free(&cfg0);
    if (g_failures) {
        printf("FAIL (%d checks failed)\n", g_failures);
        return 1;
    }
    printf("PASS\n");
    return 0;
}

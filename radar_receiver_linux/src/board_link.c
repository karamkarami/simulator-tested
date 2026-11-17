/*
 * board_link.c - the ASCII commands of the VBDR board over UDP: the client of the receiver
 * (board_link_*) and the simulator acting as the board (board_sim_*), with the text of the
 * commands and answers (README "VBDR stream and board control"):
 *
 *   write:  "CNTL W 0x02 0x1C"  ->  "CNTL W 0x02 0x1C Returned Value = 0x00000000"
 *   read:   "CNTL R 0x02 0x00"  ->  "CNTL R 0x02 0x00 Returned Value = 0x0000001C"
 *
 * The board listens on its command port (5000) and answers to the receiver PC (port 5100).
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <winsock2.h>
#  include <ws2tcpip.h>
typedef SOCKET sock_t;
#  define BAD_SOCKET   INVALID_SOCKET
#  define close_socket closesocket
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <arpa/inet.h>
#  include <errno.h>
#  include <fcntl.h>
#  include <netdb.h>
#  include <netinet/in.h>
#  include <sys/select.h>
#  include <sys/socket.h>
#  include <unistd.h>
typedef int sock_t;
#  define BAD_SOCKET   (-1)
#  define close_socket close
#endif

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "radar_rx.h"

/* ------------------------------------------------------------------------ */
/* Registers and texts                                                       */
/* ------------------------------------------------------------------------ */

static const struct { int address; const char *name, *reg; } g_items[BOARD_N_ITEMS] = {
    { 0x00, "Radar Mode",      "Radar_Mode" },
    { 0x01, "Phase Code",      "Code_index" },
    { 0x02, "F1",              "F1_Freq" },
    { 0x03, "F2",              "F2_Freq" },
    { 0x05, "Stagger Pattern", "PRF_index" },
    { 0x28, "Operation Mode",  "Operation_Mode" },
    { 0x29, "STC Main",        "STC_Main" },
    { 0x2A, "STC J1",          "STC_J1" },
    { 0x2B, "STC J2",          "STC_J2" },
    { 0x2C, "STC J3",          "STC_J3" },
    { 0x2D, "Gain Main",       "Gain_Main" },
    { 0x2E, "Gain J1",         "Gain_J1" },
    { 0x2F, "Gain J2",         "Gain_J2" },
    { 0x30, "Gain J3",         "Gain_J3" },
    { 0x32, "RPM",             "RPM (simulator)" },
    { 0x31, "Ext Trig",        "Ext_Trig" },
};

int board_item_address(int item)
{
    return item >= 0 && item < BOARD_N_ITEMS ? g_items[item].address : -1;
}

const char *board_item_name(int item)
{
    return item >= 0 && item < BOARD_N_ITEMS ? g_items[item].name : "unknown";
}

const char *board_register_name(int address)
{
    for (int i = 0; i < BOARD_N_ITEMS; i++)
        if (g_items[i].address == address) return g_items[i].reg;
    return "unknown";
}

void board_format_command(char *buf, size_t len, int write, int address, int value)
{
    snprintf(buf, len, "CNTL %c 0x%02X 0x%02X", write ? 'W' : 'R', address & 0xFF, write ? value & 0xFF : 0);
}

void board_format_reply(char *buf, size_t len, int write, int address, int value, uint32_t returned)
{
    char cmd[64];
    board_format_command(cmd, sizeof cmd, write, address, value);
    snprintf(buf, len, "%s Returned Value = 0x%08X", cmd, (unsigned)returned);
}

/* Skips spaces; reads "0x<hex>" (1 ... 8 digits). Returns the position after it, NULL if none. */
static const char *read_hex(const char *p, unsigned long *v)
{
    while (*p == ' ' || *p == '\t') p++;
    if (p[0] != '0' || (p[1] != 'x' && p[1] != 'X') || !isxdigit((unsigned char)p[2])) return NULL;
    char *end;
    *v = strtoul(p + 2, &end, 16);
    return end - (p + 2) <= 8 ? end : NULL;
}

int board_parse_command(const char *text, int *write, int *address, int *value)
{
    const char *p = text;
    while (*p == ' ' || *p == '\t') p++;
    if (strncmp(p, "CNTL", 4) != 0) return 0;
    p += 4;
    while (*p == ' ' || *p == '\t') p++;
    const char op = (char)toupper((unsigned char)*p);
    if (op != 'W' && op != 'R') return 0;
    p++;
    unsigned long a, v;
    if (!(p = read_hex(p, &a)) || !(p = read_hex(p, &v)) || a > 0xFF) return 0;
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' || *p == 0) {
        if (*p == 0) {
            *write = op == 'W';
            *address = (int)a;
            *value = (int)(v & 0xFF);
            return 1;
        }
        p++;
    }
    return 0;
}

int board_parse_reply(const char *text, int write, int address, int value, uint32_t *returned)
{
    char cmd[64];
    board_format_command(cmd, sizeof cmd, write, address, value);
    const size_t n = strlen(cmd);
    if (strncmp(text, cmd, n) != 0) return 0;
    const char *p = text + n;
    while (*p == ' ') p++;
    static const char key[] = "Returned Value";
    if (strncmp(p, key, sizeof key - 1) != 0) return 0;
    p += sizeof key - 1;
    while (*p == ' ') p++;
    if (*p != '=') return 0;
    unsigned long r;
    if (!read_hex(p + 1, &r)) return 0;
    *returned = (uint32_t)r;
    return 1;
}

/* ------------------------------------------------------------------------ */
/* Sockets                                                                   */
/* ------------------------------------------------------------------------ */

static int net_start(void)
{
#ifdef _WIN32
    WSADATA wsa;
    return WSAStartup(MAKEWORD(2, 2), &wsa) == 0;
#else
    return 1;
#endif
}

static void net_stop(void)
{
#ifdef _WIN32
    WSACleanup();
#endif
}

static int resolve(const char *host, int port, struct sockaddr_in *addr)
{
    memset(addr, 0, sizeof *addr);
    addr->sin_family = AF_INET;
    addr->sin_port = htons((uint16_t)port);
    if (inet_pton(AF_INET, host, &addr->sin_addr) == 1) return 1;
    struct addrinfo hints, *res = NULL;
    memset(&hints, 0, sizeof hints);
    hints.ai_family = AF_INET;
    hints.ai_socktype = SOCK_DGRAM;
    if (getaddrinfo(host, NULL, &hints, &res) != 0 || !res) return 0;
    addr->sin_addr = ((struct sockaddr_in *)res->ai_addr)->sin_addr;
    freeaddrinfo(res);
    return 1;
}

/* A UDP socket bound to port (all interfaces); BAD_SOCKET with err on failure. */
static sock_t bound_socket(int port, char *err, size_t errLen)
{
    sock_t s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (s == BAD_SOCKET) {
        snprintf(err, errLen, "cannot create a UDP socket");
        return BAD_SOCKET;
    }
    struct sockaddr_in a;
    memset(&a, 0, sizeof a);
    a.sin_family = AF_INET;
    a.sin_addr.s_addr = htonl(INADDR_ANY);
    a.sin_port = htons((uint16_t)port);
    if (bind(s, (struct sockaddr *)&a, (socklen_t)sizeof a) != 0) {
        snprintf(err, errLen, "cannot bind UDP port %d (in use?)", port);
        close_socket(s);
        return BAD_SOCKET;
    }
    return s;
}

/* Waits up to ms for a datagram (ms <= 0: does not wait). 1 readable, 0 timeout, -1 error. */
static int wait_readable(sock_t s, double ms)
{
    fd_set set;
    FD_ZERO(&set);
    FD_SET(s, &set);
    struct timeval tv;
    if (ms < 0) ms = 0;
    tv.tv_sec = (long)(ms / 1000.0);
    tv.tv_usec = (long)((ms - 1000.0 * (double)tv.tv_sec) * 1000.0);
#ifdef _WIN32
    const int r = select(0, &set, NULL, NULL, &tv);
    return r == SOCKET_ERROR ? -1 : r > 0;
#else
    const int r = select(s + 1, &set, NULL, NULL, &tv);
    if (r < 0) return errno == EINTR ? 0 : -1;
    return r > 0;
#endif
}

/* Receives one datagram as a 0-terminated text (at most len - 1 bytes). Returns its length or -1. */
static int recv_text(sock_t s, char *buf, size_t len, struct sockaddr_in *from)
{
    struct sockaddr_in a;
    socklen_t al = (socklen_t)sizeof a;
    const int n = (int)recvfrom(s, buf, (int)len - 1, 0, (struct sockaddr *)&a, &al);
    if (n < 0) return -1;
    buf[n] = 0;
    if (from) *from = a;
    return n;
}

/* ------------------------------------------------------------------------ */
/* Client                                                                    */
/* ------------------------------------------------------------------------ */

struct BoardLink {
    sock_t             sock;
    struct sockaddr_in board;
};

BoardLink *board_link_open(const char *host, int port, int replyPort, char *err, size_t errLen)
{
    if (!net_start()) {
        snprintf(err, errLen, "WSAStartup failed");
        return NULL;
    }
    BoardLink *l = calloc(1, sizeof *l);
    if (!l) { net_stop(); snprintf(err, errLen, "out of memory"); return NULL; }
    if (!resolve(host, port, &l->board)) {
        snprintf(err, errLen, "unknown board address %s", host);
        free(l);
        net_stop();
        return NULL;
    }
    l->sock = bound_socket(replyPort, err, errLen);
    if (l->sock == BAD_SOCKET) {
        free(l);
        net_stop();
        return NULL;
    }
    return l;
}

int board_link_exchange(BoardLink *l, int write, int address, int value, double timeout_s, int retries,
                        uint32_t *returned)
{
    char cmd[64], buf[512];
    board_format_command(cmd, sizeof cmd, write, address, value);
    while (wait_readable(l->sock, 0) > 0 && recv_text(l->sock, buf, sizeof buf, NULL) >= 0) {}   /* old answers */
    for (int attempt = 0; attempt <= retries; attempt++) {
        if (sendto(l->sock, cmd, (int)strlen(cmd), 0, (const struct sockaddr *)&l->board, (socklen_t)sizeof l->board) < 0)
            return -1;
        const double deadline = now_seconds() + timeout_s;
        for (;;) {
            const double left_ms = (deadline - now_seconds()) * 1e3;
            if (left_ms <= 0) break;
            const int r = wait_readable(l->sock, left_ms);
            if (r < 0) return -1;
            if (r == 0) continue;
            if (recv_text(l->sock, buf, sizeof buf, NULL) < 0) continue;   /* e.g. ICMP unreachable (Windows) */
            if (board_parse_reply(buf, write, address, value, returned)) return 1;
        }
    }
    return 0;
}

void board_link_close(BoardLink *l)
{
    if (!l) return;
    close_socket(l->sock);
    free(l);
    net_stop();
}

/* ------------------------------------------------------------------------ */
/* The simulator as the board                                                */
/* ------------------------------------------------------------------------ */

struct BoardSim {
    sock_t             sock;
    struct sockaddr_in reply;
    uint8_t            regs[256];
    int                verbose;
};

BoardSim *board_sim_open(int port, const char *replyHost, int replyPort, const uint8_t *regs, int verbose,
                         char *err, size_t errLen)
{
    if (!net_start()) {
        snprintf(err, errLen, "WSAStartup failed");
        return NULL;
    }
    BoardSim *b = calloc(1, sizeof *b);
    if (!b) { net_stop(); snprintf(err, errLen, "out of memory"); return NULL; }
    if (!resolve(replyHost, replyPort, &b->reply)) {
        snprintf(err, errLen, "unknown address %s of the receiver PC", replyHost);
        free(b);
        net_stop();
        return NULL;
    }
    b->sock = bound_socket(port, err, errLen);
    if (b->sock == BAD_SOCKET) {
        free(b);
        net_stop();
        return NULL;
    }
    if (regs) memcpy(b->regs, regs, sizeof b->regs);
    b->verbose = verbose;
    if (verbose)
        fprintf(stderr, "Board commands: listening on UDP port %d, answers to %s:%d\n", port, replyHost, replyPort);
    return b;
}

int board_sim_poll(BoardSim *b, int wait_ms, uint32_t *written)
{
    if (!b) return 0;
    int n = 0;
    char buf[512], reply[160];
    for (int guard = 0; guard < 256; guard++) {
        if (wait_readable(b->sock, n == 0 ? wait_ms : 0) <= 0) break;
        if (recv_text(b->sock, buf, sizeof buf, NULL) < 0) continue;
        n++;
        int write, address, value;
        if (!board_parse_command(buf, &write, &address, &value)) {
            if (b->verbose) fprintf(stderr, "board command: unknown \"%.80s\"\n", buf);
            continue;
        }
        if (write) {
            b->regs[address] = (uint8_t)value;
            if (written) written[address / 32] |= 1u << (address % 32);
            board_format_reply(reply, sizeof reply, 1, address, value, 0);
            if (b->verbose) fprintf(stderr, "board command: write %s (0x%02X) = %d\n", board_register_name(address), address, value);
        } else {
            board_format_reply(reply, sizeof reply, 0, address, 0, b->regs[address]);
            if (b->verbose) fprintf(stderr, "board command: read  %s (0x%02X) -> %d\n", board_register_name(address), address,
                                    b->regs[address]);
        }
        sendto(b->sock, reply, (int)strlen(reply), 0, (const struct sockaddr *)&b->reply, (socklen_t)sizeof b->reply);
    }
    return n;
}

uint8_t board_sim_register(const BoardSim *b, int address)
{
    return b && address >= 0 && address < 256 ? b->regs[address] : 0;
}

void board_sim_set_register(BoardSim *b, int address, uint8_t value)
{
    if (b && address >= 0 && address < 256) b->regs[address] = value;
}

void board_sim_close(BoardSim *b)
{
    if (!b) return;
    close_socket(b->sock);
    free(b);
    net_stop();
}

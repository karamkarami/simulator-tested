/*
 * test_udp_input.c - test of the UDP input (udp_input.c) and of now_seconds().
 *
 * Network part: the input is opened on a free localhost port with a record
 * file and datagrams are sent from a second UDP socket of this program:
 *   - PRI A complete (packets shuffled, one duplicated, one late duplicate),
 *   - PRI B with a missing packet (dropped when PRI C starts),
 *   - PRI C complete (shuffled, one duplicated), with another layout,
 *   - foreign datagrams (bad size, magic, layout, older pulse number),
 *   - a restarted sender (pulse numbers start again from 1).
 * Replay part: the record file is replayed and must give the same PRIs and
 * statistics, then -1. Fuzz part: random datagrams are replayed (no crash).
 * VBDR part: streams of the VBDR board are written to a record file (packets
 * shuffled and duplicated, a stream with a missing packet, a lost stream,
 * a late packet) and replayed with the format UDP_FORMAT_VBDR: the 4 channels,
 * the system parameters and vbdr_complete_pri() (time from the PRTs with a lost
 * stream, scan from the azimuth wrap, stagger position from the PRT, layout of
 * the configuration).
 *
 * Prints PASS or FAIL and returns 0 on success.
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
#  include <netinet/in.h>
#  include <sys/socket.h>
#  include <unistd.h>
typedef int sock_t;
#  define BAD_SOCKET   (-1)
#  define close_socket close
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "radar_rx.h"

#define RECORD_FILE "test_udp_input.rec"
#define FUZZ_FILE   "test_udp_input_fuzz.rec"
#define VBDR_FILE   "test_udp_input_vbdr.rec"
#define MAGIC       0x52444443u
#define RESYNC_RUN  64          /* RESYNC_AFTER of udp_input.c */

static int g_failures = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);           \
            g_failures++;                                                    \
        }                                                                    \
    } while (0)

/* Description of a test PRI. */
typedef struct {
    uint32_t pulse, scan;
    float    az;
    uint16_t nCh, stagger;
    uint32_t nSamples;
    double   time;
    float    fs, t0;
    uint32_t idx[4];
    uint32_t perPacket;          /* payload bytes per packet (multiple of 4) */
} TestPri;

/* ------------------------------------------------------------------------ */
/* Datagram construction                                                    */
/* ------------------------------------------------------------------------ */

static uint32_t g_rand = 12345u;
static uint32_t rnd(void)
{
    g_rand = g_rand * 1664525u + 1013904223u;
    return g_rand >> 8;
}

/* Deterministic sample value (I: q = 0, Q: q = 1), covering the int16 range. */
static int16_t sample_value(const TestPri *p, uint32_t ch, uint32_t n, uint32_t q)
{
    uint32_t x = p->pulse * 2654435761u ^ (ch * 40503u + n * 2u + q) * 2246822519u;
    x ^= x >> 15;
    x *= 0x2C1B3C6Du;
    x ^= x >> 12;
    return (int16_t)((int32_t)(x & 0xFFFFu) - 32768);
}

static uint32_t pri_bytes(const TestPri *p) { return 4u * p->nSamples * p->nCh; }
static int packet_count(const TestPri *p) { return (int)((pri_bytes(p) + p->perPacket - 1) / p->perPacket); }

static void put16(uint8_t *b, int off, uint16_t v) { memcpy(b + off, &v, sizeof v); }
static void put32(uint8_t *b, int off, uint32_t v) { memcpy(b + off, &v, sizeof v); }
static void putf(uint8_t *b, int off, float v)     { memcpy(b + off, &v, sizeof v); }
static void putd(uint8_t *b, int off, double v)    { memcpy(b + off, &v, sizeof v); }

static void write_header(uint8_t *b, const TestPri *p, uint16_t pktIdx, uint16_t pktCount, uint32_t first)
{
    int i;
    put32(b, 0, MAGIC);
    put32(b, 4, p->pulse);
    put32(b, 8, p->scan);
    putf(b, 12, p->az);
    put16(b, 16, pktIdx);
    put16(b, 18, pktCount);
    put16(b, 20, p->nCh);
    put16(b, 22, p->stagger);
    put32(b, 24, p->nSamples);
    put32(b, 28, first);
    putd(b, 32, p->time);
    putf(b, 40, p->fs);
    putf(b, 44, p->t0);
    for (i = 0; i < 4; i++)
        put32(b, 48 + 4 * i, p->idx[i]);
}

/* Builds packet 'pktIdx' of PRI p in b. Returns the datagram length. */
static size_t make_packet(uint8_t *b, const TestPri *p, int pktIdx)
{
    uint32_t total = p->nSamples * p->nCh;          /* complex samples */
    uint32_t first = (uint32_t)pktIdx * (p->perPacket / 4);
    uint32_t count = total - first < p->perPacket / 4 ? total - first : p->perPacket / 4;
    uint32_t k;

    write_header(b, p, (uint16_t)pktIdx, (uint16_t)packet_count(p), first);
    for (k = 0; k < count; k++) {
        uint32_t flat = first + k, ch = flat / p->nSamples, n = flat % p->nSamples;
        int16_t iq[2];
        iq[0] = sample_value(p, ch, n, 0);
        iq[1] = sample_value(p, ch, n, 1);
        memcpy(b + 64 + 4 * k, iq, 4);
    }
    return 64 + 4 * (size_t)count;
}

/* ------------------------------------------------------------------------ */
/* Sender                                                                    */
/* ------------------------------------------------------------------------ */

typedef struct {
    sock_t s;
    struct sockaddr_in dest;
    int sent;                    /* datagrams sent */
} Sender;

static void send_bytes(Sender *tx, const uint8_t *b, size_t len)
{
#ifdef _WIN32
    int r = sendto(tx->s, (const char *)b, (int)len, 0,
                   (const struct sockaddr *)&tx->dest, (int)sizeof tx->dest);
#else
    int r = (int)sendto(tx->s, b, len, 0, (const struct sockaddr *)&tx->dest, (socklen_t)sizeof tx->dest);
#endif
    CHECK(r == (int)len);
    tx->sent++;
}

static void send_packet(Sender *tx, const TestPri *p, int pktIdx)
{
    static uint8_t b[65536];
    send_bytes(tx, b, make_packet(b, p, pktIdx));
}

/* Sends the packets of p in a shuffled order with one duplicate in the
 * middle; 'skip' (or -1) is a packet that is not sent. */
static void send_pri(Sender *tx, const TestPri *p, int skip)
{
    int order[1024], n = packet_count(p), i;
    for (i = 0; i < n; i++)
        order[i] = i;
    for (i = n - 1; i > 0; i--) {
        int j = (int)(rnd() % (uint32_t)(i + 1)), t = order[i];
        order[i] = order[j];
        order[j] = t;
    }
    for (i = 0; i < n; i++) {
        if (order[i] != skip)
            send_packet(tx, p, order[i]);
        if (i == n / 2 && order[0] != skip)
            send_packet(tx, p, order[0]);              /* duplicate */
    }
}

static int sender_open(Sender *tx, int port)
{
    memset(tx, 0, sizeof *tx);
    tx->s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    tx->dest.sin_family = AF_INET;
    tx->dest.sin_port = htons((uint16_t)port);
    tx->dest.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    return tx->s != BAD_SOCKET;
}

/* Returns a UDP port that is currently free. */
static int free_port(void)
{
    struct sockaddr_in a;
    socklen_t len = (socklen_t)sizeof a;
    sock_t s = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    int port = -1;

    memset(&a, 0, sizeof a);
    a.sin_family = AF_INET;
    a.sin_addr.s_addr = htonl(INADDR_ANY);
    a.sin_port = 0;
    if (s != BAD_SOCKET && bind(s, (struct sockaddr *)&a, len) == 0 &&
        getsockname(s, (struct sockaddr *)&a, &len) == 0)
        port = ntohs(a.sin_port);
    if (s != BAD_SOCKET)
        close_socket(s);
    return port;
}

/* ------------------------------------------------------------------------ */
/* Checks                                                                    */
/* ------------------------------------------------------------------------ */

static void check_pri(const PriData *d, const TestPri *p)
{
    uint32_t ch, n;
    long bad = 0;

    CHECK(d->pulseNumber == p->pulse);
    CHECK(d->scanNumber == p->scan);
    CHECK(d->az_deg == p->az);
    CHECK(d->staggerIndex == p->stagger);
    CHECK(d->time_s == p->time);
    CHECK(d->ddcFs_Hz == p->fs);
    CHECK(d->ddcT0_s == p->t0);
    CHECK(d->shortIdx[0] == p->idx[0] && d->shortIdx[1] == p->idx[1]);
    CHECK(d->longIdx[0] == p->idx[2] && d->longIdx[1] == p->idx[3]);
    CHECK(d->nSamples == (int)p->nSamples);
    CHECK(d->nCh == (int)p->nCh);
    CHECK(d->iq != NULL);
    CHECK(d->capacity >= 2 * (size_t)p->nSamples * p->nCh);
    if (!d->iq || d->nSamples != (int)p->nSamples || d->nCh != (int)p->nCh)
        return;
    for (ch = 0; ch < p->nCh; ch++)
        for (n = 0; n < p->nSamples; n++) {
            size_t k = 2 * ((size_t)ch * p->nSamples + n);
            bad += d->iq[k] != sample_value(p, ch, n, 0);
            bad += d->iq[k + 1] != sample_value(p, ch, n, 1);
        }
    if (bad)
        printf("  pulse %u: %ld wrong sample values\n", (unsigned)p->pulse, bad);
    CHECK(bad == 0);
}

static void check_stats(const UdpStats *s, uint64_t datagrams, uint64_t foreign,
                        uint64_t complete, uint64_t incomplete)
{
    printf("  stats: datagrams %llu, foreign %llu, complete %llu, incomplete %llu\n",
           (unsigned long long)s->datagrams, (unsigned long long)s->foreign,
           (unsigned long long)s->prisComplete, (unsigned long long)s->prisIncomplete);
    CHECK(s->datagrams == datagrams);
    CHECK(s->foreign == foreign);
    CHECK(s->prisComplete == complete);
    CHECK(s->prisIncomplete == incomplete);
}

/* ------------------------------------------------------------------------ */
/* Test PRIs                                                                 */
/* ------------------------------------------------------------------------ */

static TestPri g_a, g_b, g_c, g_restart[RESYNC_RUN + 1];

static void make_test_pris(void)
{
    TestPri a = { 100, 3, 12.5f, 4, 2, 1500, 1.25, 5.0e6f, 1.0e-6f, { 1, 10, 20, 300 }, 3500 };
    TestPri c = { 102, 3, 12.75f, 2, 3, 777, 1.2534567890123, 2.5e6f, 2.0e-6f, { 5, 50, 60, 700 }, 1000 };
    int i;

    g_a = a;                       /* 24000 bytes: 7 packets, the last one shorter */
    g_b = a;
    g_b.pulse = 101;
    g_b.time = 1.251;
    g_c = c;                       /* 6216 bytes: 7 packets */
    for (i = 0; i <= RESYNC_RUN; i++) {
        TestPri r = { (uint32_t)(i + 1), 0, 1.0f, 1, 1, 16, 0.001 * i, 1.0e6f, 0.0f, { 1, 2, 3, 4 }, 64 };
        g_restart[i] = r;          /* single-packet PRIs of a restarted sender */
    }
}

/* ------------------------------------------------------------------------ */
/* Network part                                                              */
/* ------------------------------------------------------------------------ */

/* Sends the foreign datagrams. Returns how many of them must be counted as foreign. */
static int send_foreign(Sender *tx)
{
    static uint8_t b[65536];
    TestPri p = g_c;
    size_t len;
    int sent0 = tx->sent;

    memset(b, 0x5A, 10);                                 /* too short */
    send_bytes(tx, b, 10);

    p.pulse = 103;
    len = make_packet(b, &p, 0);
    put32(b, 0, 0x12345678u);                            /* wrong magic */
    send_bytes(tx, b, len);

    make_packet(b, &p, 0);                               /* payload not a multiple of 4 */
    send_bytes(tx, b, 64 + 6);

    make_packet(b, &p, 0);                               /* payload does not fit */
    put32(b, 28, p.nSamples * p.nCh - 1);
    send_bytes(tx, b, 64 + 8);

    len = make_packet(b, &p, 0);                         /* packet index >= packet count */
    put16(b, 16, (uint16_t)packet_count(&p));
    send_bytes(tx, b, len);

    len = make_packet(b, &p, 0);                         /* no channel */
    put16(b, 20, 0);
    send_bytes(tx, b, len);

    len = make_packet(b, &p, 0);                         /* absurd size */
    put32(b, 24, 0x7FFFFFFFu);
    put16(b, 20, 48);
    send_bytes(tx, b, len);

    send_packet(tx, &g_b, 0);                            /* older pulse number */

    send_packet(tx, &g_c, 1);                            /* late duplicate: ignored, not foreign */
    return tx->sent - sent0 - 1;
}

static void test_network(int port, UdpStats *finalStats)
{
    char err[256];
    PriData pri;
    UdpStats st;
    Sender tx;
    UdpInput *in;
    int foreign, r, i;
    double t0, dt;

    memset(&pri, 0, sizeof pri);
    in = udp_input_open(port, RECORD_FILE, NULL, err, sizeof err);
    if (!in) {
        printf("FAIL udp_input_open: %s\n", err);
        g_failures++;
        return;
    }
    CHECK(sender_open(&tx, port));

    /* nothing sent yet: timeout */
    t0 = now_seconds();
    CHECK(udp_input_next(in, &pri, 50) == 0);
    dt = now_seconds() - t0;
    CHECK(dt >= 0.045 && dt < 2.0);

    /* PRI A */
    send_pri(&tx, &g_a, -1);
    r = udp_input_next(in, &pri, 2000);
    CHECK(r == 1);
    if (r == 1)
        check_pri(&pri, &g_a);
    send_packet(&tx, &g_a, 0);                           /* late duplicate */

    /* PRI B (incomplete), then PRI C */
    send_pri(&tx, &g_b, 3);
    send_pri(&tx, &g_c, -1);
    r = udp_input_next(in, &pri, 2000);
    CHECK(r == 1);
    if (r == 1)
        check_pri(&pri, &g_c);

    /* foreign datagrams: nothing is delivered, the call times out */
    foreign = send_foreign(&tx);
    t0 = now_seconds();
    CHECK(udp_input_next(in, &pri, 300) == 0);
    dt = now_seconds() - t0;
    printf("  timeout of 300 ms took %.1f ms\n", 1e3 * dt);
    CHECK(dt >= 0.29 && dt < 2.0);
    udp_input_stats(in, &st);
    check_stats(&st, (uint64_t)tx.sent, (uint64_t)foreign, 2, 1);

    /* restarted sender: older pulse numbers are foreign until RESYNC_RUN of them */
    for (i = 0; i <= RESYNC_RUN; i++)
        send_packet(&tx, &g_restart[i], 0);
    r = udp_input_next(in, &pri, 2000);
    CHECK(r == 1);
    if (r == 1)
        check_pri(&pri, &g_restart[RESYNC_RUN]);
    CHECK(udp_input_next(in, &pri, 100) == 0);
    udp_input_stats(in, &st);
    check_stats(&st, (uint64_t)tx.sent, (uint64_t)(foreign + RESYNC_RUN), 3, 1);
    *finalStats = st;

    close_socket(tx.s);
    udp_input_close(in);
    pri_data_free(&pri);
    CHECK(pri.iq == NULL && pri.capacity == 0);
}

/* ------------------------------------------------------------------------ */
/* Replay part                                                               */
/* ------------------------------------------------------------------------ */

static void test_replay(const UdpStats *expected)
{
    char err[256];
    PriData pri;
    UdpStats st;
    UdpInput *in;
    int r;

    memset(&pri, 0, sizeof pri);
    in = udp_input_open(0, NULL, RECORD_FILE, err, sizeof err);
    if (!in) {
        printf("FAIL udp_input_open (replay): %s\n", err);
        g_failures++;
        return;
    }
    r = udp_input_next(in, &pri, 0);
    CHECK(r == 1);
    if (r == 1)
        check_pri(&pri, &g_a);
    r = udp_input_next(in, &pri, 0);
    CHECK(r == 1);
    if (r == 1)
        check_pri(&pri, &g_c);
    r = udp_input_next(in, &pri, 0);
    CHECK(r == 1);
    if (r == 1)
        check_pri(&pri, &g_restart[RESYNC_RUN]);
    CHECK(udp_input_next(in, &pri, 0) == -1);
    CHECK(udp_input_next(in, &pri, 0) == -1);
    udp_input_stats(in, &st);
    check_stats(&st, expected->datagrams, expected->foreign,
                expected->prisComplete, expected->prisIncomplete);
    udp_input_close(in);
    pri_data_free(&pri);
}

/* ------------------------------------------------------------------------ */
/* Fuzz part: random datagrams in a replay file must never crash            */
/* ------------------------------------------------------------------------ */

static void write_record(FILE *f, const uint8_t *b, uint32_t len)
{
    uint8_t prefix[4] = { (uint8_t)len, (uint8_t)(len >> 8), (uint8_t)(len >> 16), (uint8_t)(len >> 24) };
    fwrite(prefix, 1, 4, f);
    fwrite(b, 1, len, f);
}

/* Random PRIs whose datagrams are sometimes dropped, duplicated or damaged,
 * with junk datagrams in between, written to a replay file. A PRI without
 * any lost or damaged datagram ("intact") must be delivered with its exact
 * samples; for the others the input just must not crash. */
enum { FUZZ_PRIS = 1500, FUZZ_BASE = 100000 };

typedef struct {
    TestPri pri[FUZZ_PRIS];
    int     intact[FUZZ_PRIS];
    int     delivered[FUZZ_PRIS];
    uint32_t records;
} Fuzz;

static void fuzz_layout(TestPri *p, int k)
{
    TestPri t = { (uint32_t)(FUZZ_BASE + k), 7, 0.5f * (float)k, 1, 1, 1, 0.001 * k,
                  4.0e6f, 1.0e-6f, { 1, 2, 3, 4 }, 4 };
    t.nCh = (uint16_t)(1 + rnd() % 3);
    t.nSamples = 1 + rnd() % 200;
    t.stagger = (uint16_t)(1 + k % 3);
    t.perPacket = 4 * (8 + rnd() % 64);
    *p = t;
}

/* Damages a datagram; never gives it the pulse number of another test PRI. */
static uint32_t fuzz_damage(uint8_t *b, uint32_t len, uint32_t pulse)
{
    uint32_t off;
    switch (rnd() % 4) {
    case 0:  return rnd() % len;                               /* truncated */
    case 1:  put32(b, 4, pulse - 10000 - rnd() % 1000); break; /* old pulse */
    case 2:  put32(b, 0, rnd()); break;                        /* magic     */
    default:
        do off = rnd() % 64; while (off >= 4 && off < 8);
        b[off] = (uint8_t)rnd();                               /* header    */
        break;
    }
    return len;
}

static void fuzz_record(FILE *f, Fuzz *z, const uint8_t *b, uint32_t len)
{
    write_record(f, b, len);
    z->records++;
}

static void fuzz_write_pri(FILE *f, Fuzz *z, int k)
{
    static uint8_t b[65536];
    const TestPri *p = &z->pri[k];
    int n = packet_count(p), i, order[1024];

    z->intact[k] = 1;
    for (i = 0; i < n; i++)
        order[i] = i;
    for (i = n - 1; i > 0; i--) {
        int j = (int)(rnd() % (uint32_t)(i + 1)), t = order[i];
        order[i] = order[j];
        order[j] = t;
    }
    for (i = 0; i < n; i++) {
        uint32_t len = (uint32_t)make_packet(b, p, order[i]), r = rnd() % 100, j;
        if (r < 1) {                                           /* lost      */
            z->intact[k] = 0;
            continue;
        }
        if (r < 3) {                                           /* damaged   */
            z->intact[k] = 0;
            fuzz_record(f, z, b, fuzz_damage(b, len, p->pulse));
            len = (uint32_t)make_packet(b, p, order[i]);
        }
        if (r >= 3 && r < 5)                                   /* duplicated */
            fuzz_record(f, z, b, len);
        if (r >= 5 && r < 6) {                                 /* junk      */
            uint32_t junk = rnd() % 300;
            static uint8_t jb[300];
            for (j = 0; j < junk; j++)
                jb[j] = (uint8_t)rnd();
            fuzz_record(f, z, jb, junk);
        }
        fuzz_record(f, z, b, len);
    }
}

static void fuzz_check(Fuzz *z, const PriData *pri)
{
    uint32_t k = pri->pulseNumber - FUZZ_BASE;
    CHECK(pri->iq != NULL && pri->nSamples > 0 && pri->nCh > 0);
    CHECK(pri->capacity >= 2 * (size_t)pri->nSamples * (size_t)pri->nCh);
    if (k >= FUZZ_PRIS) {
        CHECK(k < FUZZ_PRIS);
        return;
    }
    z->delivered[k]++;
    if (z->intact[k])
        check_pri(pri, &z->pri[k]);
}

static void test_fuzz(void)
{
    static Fuzz z;
    char err[256];
    PriData pri;
    UdpStats st;
    UdpInput *in;
    FILE *f = fopen(FUZZ_FILE, "wb");
    int k, r, nIntact = 0, missing = 0, twice = 0;
    const uint8_t corrupt[8] = { 0x70, 0x11, 0x01, 0x00, 1, 2, 3, 4 };

    CHECK(f != NULL);
    if (!f)
        return;
    for (k = 0; k < FUZZ_PRIS; k++) {
        fuzz_layout(&z.pri[k], k);
        fuzz_write_pri(f, &z, k);
    }
    fwrite(corrupt, 1, sizeof corrupt, f);               /* length 70000: corrupted, end */
    fclose(f);

    memset(&pri, 0, sizeof pri);
    in = udp_input_open(0, NULL, FUZZ_FILE, err, sizeof err);
    CHECK(in != NULL);
    if (!in)
        return;
    while ((r = udp_input_next(in, &pri, 0)) == 1)
        fuzz_check(&z, &pri);
    CHECK(r == -1);
    for (k = 0; k < FUZZ_PRIS; k++) {
        nIntact += z.intact[k];
        missing += z.intact[k] && z.delivered[k] != 1;
        twice += z.delivered[k] > 1;
    }
    udp_input_stats(in, &st);
    printf("  fuzz: %u datagrams, %llu foreign, %llu PRIs delivered (%d intact), %llu incomplete\n",
           (unsigned)z.records, (unsigned long long)st.foreign, (unsigned long long)st.prisComplete,
           nIntact, (unsigned long long)st.prisIncomplete);
    CHECK(st.datagrams == z.records);
    CHECK(missing == 0);
    CHECK(twice == 0);
    CHECK(nIntact > FUZZ_PRIS / 4);
    udp_input_close(in);
    pri_data_free(&pri);
    remove(FUZZ_FILE);
}

/* ------------------------------------------------------------------------ */

static void test_open_errors(void)
{
    char err[256];
    UdpInput *in = udp_input_open(0, NULL, "does/not/exist.rec", err, sizeof err);
    CHECK(in == NULL && err[0] != '\0');
    in = udp_input_open(70000, NULL, NULL, err, sizeof err);
    CHECK(in == NULL && err[0] != '\0');
    udp_input_close(NULL);
    pri_data_free(NULL);
}

/* ------------------------------------------------------------------------ */
/* VBDR stream                                                               */
/* ------------------------------------------------------------------------ */

static int16_t vbdr_value(uint32_t stream, int ch, int n, int q)
{
    return (int16_t)((int)((stream * 131u + (uint32_t)ch * 7919u + (uint32_t)n * 3u + (uint32_t)q) % 60001u) - 30000);
}

static void put_rec(FILE *f, const uint8_t *b, size_t len)
{
    uint8_t pre[4] = { (uint8_t)len, (uint8_t)(len >> 8), (uint8_t)(len >> 16), (uint8_t)(len >> 24) };
    fwrite(pre, 1, 4, f);
    fwrite(b, 1, len, f);
}

/* The packets of one stream: nS samples, payload bytes per packet, azimuth [0, 65535], PRTs [us].
 * Returns the number of packets; packet k (0-based) in pk[k] with length len[k]. */
static int vbdr_stream(uint32_t stream, int nS, int payload, uint16_t az, uint16_t prev, uint16_t cur, unsigned loss,
                       uint8_t pk[][1500], size_t *len)
{
    const int halfBytes = 64 + 8 * nS, nHalf = (halfBytes + payload - 1) / payload, count = 2 * nHalf;
    uint8_t *half = malloc((size_t)halfBytes);
    for (int h = 0; h < 2; h++) {
        uint16_t par[32];
        memset(par, 0, sizeof par);
        par[0] = (uint16_t)stream; par[1] = az; par[3] = 0x5555; par[4] = prev; par[5] = cur;
        memcpy(half, par, 64);
        for (int n = 0; n < nS; n++)
            for (int c = 0; c < 2; c++)
                for (int q = 0; q < 2; q++) {
                    int16_t v = vbdr_value(stream, 2 * h + c, n, q);
                    memcpy(half + 64 + 8 * n + 4 * c + 2 * q, &v, 2);
                }
        for (int j = 0; j < nHalf; j++) {
            const int k = h * nHalf + j, n = (j + 1) * payload <= halfBytes ? payload : halfBytes - j * payload;
            put32(pk[k], 0, stream);
            put32(pk[k], 4, (uint32_t)k + 1);
            put32(pk[k], 8, (uint32_t)count);
            memset(pk[k] + 12, 0, 16);
            snprintf((char *)pk[k] + 12, 16, "Loss=%u", loss);
            memcpy(pk[k] + 28, half + j * payload, (size_t)n);
            len[k] = 28 + (size_t)n;
        }
    }
    free(half);
    return count;
}

static void check_vbdr_samples(const PriData *d, uint32_t stream, int nS)
{
    int bad = 0;
    CHECK(d->nCh == 4 && d->nSamples == nS && d->format == UDP_FORMAT_VBDR);
    for (int ch = 0; ch < 4 && d->nCh == 4; ch++)
        for (int n = 0; n < nS; n++)
            for (int q = 0; q < 2; q++)
                bad += d->iq[2 * ((size_t)ch * nS + n) + q] != vbdr_value(stream, ch, n, q);
    CHECK(bad == 0);
}

static void test_vbdr(void)
{
    static uint8_t pk[64][1500];
    size_t len[64];
    const int nS = 300, payload = 700;            /* half = 2464 bytes: 4 packets, the last 364 */
    FILE *f = fopen(VBDR_FILE, "wb");
    CHECK(f != NULL);
    if (!f) return;
    /* stream 10: shuffled with a duplicate; staggered PRTs 1000, 1100, 1250 us (positions 1, 2, 3) */
    int n = vbdr_stream(10, nS, payload, 65000, 1250, 1000, 0, pk, len);
    CHECK(n == 8 && len[3] == 28 + 364);
    const int order1[9] = { 7, 0, 3, 3, 1, 6, 2, 5, 4 };
    for (int i = 0; i < 9; i++) put_rec(f, pk[order1[i]], len[order1[i]]);
    put_rec(f, pk[2], len[2]);                    /* late duplicate of the delivered stream */
    /* stream 11: one packet missing (dropped when stream 12 starts) */
    n = vbdr_stream(11, nS, payload, 65300, 1000, 1100, 0, pk, len);
    for (int i = 0; i < n; i++) if (i != 5) put_rec(f, pk[i], len[i]);
    /* stream 12: complete, the azimuth wraps (scan 2) */
    n = vbdr_stream(12, nS, payload, 100, 1100, 1250, 2, pk, len);
    for (int i = n - 1; i >= 0; i--) put_rec(f, pk[i], len[i]);
    /* stream 13 never sent; stream 14 (position 2: previous PRT 1000 = position 1) */
    n = vbdr_stream(14, nS, payload, 400, 1000, 1100, 2, pk, len);
    for (int i = 0; i < n; i++) put_rec(f, pk[i], len[i]);
    {   /* a packet of a foreign size and a label of a wrong count */
        uint8_t junk[40];
        memset(junk, 0, sizeof junk);
        put_rec(f, junk, 20);
        put32(junk, 0, 15); put32(junk, 4, 1); put32(junk, 8, 3);
        put_rec(f, junk, 40);
    }
    fclose(f);

    char err[256];
    UdpInput *in = udp_input_open(0, NULL, VBDR_FILE, err, sizeof err);
    CHECK(in != NULL);
    if (!in) return;
    udp_input_set_format(in, UDP_FORMAT_VBDR);
    RadarConfig cfg;
    memset(&cfg, 0, sizeof cfg);
    double seq[3] = { 1000e-6, 1100e-6, 1250e-6 };
    cfg.fsOut_Hz = 2e6; cfg.ddcT0_s = 3e-6; cfg.staggered = 1; cfg.nStagger = 3; cfg.nPriSequence = 3;
    cfg.priSequence_s = seq; cfg.meanPri_s = (seq[0] + seq[1] + seq[2]) / 3;
    cfg.pulse[PULSE_SHORT].enable = 1; cfg.pulse[PULSE_SHORT].rxFirst = 10; cfg.pulse[PULSE_SHORT].rxCount = 90;
    cfg.pulse[PULSE_LONG].enable = 1;  cfg.pulse[PULSE_LONG].rxFirst = 120; cfg.pulse[PULSE_LONG].rxCount = 180;
    VbdrTiming vt;
    memset(&vt, 0, sizeof vt);
    PriData d;
    memset(&d, 0, sizeof d);

    CHECK(udp_input_next(in, &d, 0) == 1);       /* stream 10 */
    vbdr_complete_pri(&vt, &cfg, &d);
    CHECK(d.pulseNumber == 10 && d.time_s == 0 && d.scanNumber == 1 && d.staggerIndex == 1);
    CHECK(fabs(d.az_deg - 65000 * 360.0 / 65536) < 1e-3 && d.prfNumber == 10 && d.loss == 0);
    CHECK(fabs(d.prevPrt_s - 1250e-6) < 1e-12 && fabs(d.curPrt_s - 1000e-6) < 1e-12);
    CHECK(d.ddcFs_Hz == 2e6f && d.shortIdx[0] == 11 && d.shortIdx[1] == 100 && d.longIdx[0] == 121 && d.longIdx[1] == 300);
    check_vbdr_samples(&d, 10, nS);

    CHECK(udp_input_next(in, &d, 0) == 1);       /* stream 12 (11 is incomplete) */
    vbdr_complete_pri(&vt, &cfg, &d);
    /* stream 10 lasted 1000 us, 11 lasted 1100 (previous PRT of 12) */
    CHECK(d.pulseNumber == 12 && fabs(d.time_s - 2100e-6) < 1e-12 && d.scanNumber == 2 && d.staggerIndex == 3 && d.loss == 2);
    check_vbdr_samples(&d, 12, nS);

    CHECK(udp_input_next(in, &d, 0) == 1);       /* stream 14 (13 lost) */
    vbdr_complete_pri(&vt, &cfg, &d);
    /* 12 lasted 1250, 13 lasted 1000 (previous PRT of 14) */
    CHECK(d.pulseNumber == 14 && fabs(d.time_s - 4350e-6) < 1e-12 && d.scanNumber == 2 && d.staggerIndex == 2);
    check_vbdr_samples(&d, 14, nS);
    CHECK(vt.missing == 2);

    CHECK(udp_input_next(in, &d, 0) == -1);
    UdpStats st;
    udp_input_stats(in, &st);
    CHECK(st.prisComplete == 3 && st.prisIncomplete == 1 && st.foreign == 2 && st.boardLoss == 2);
    udp_input_close(in);
    pri_data_free(&d);
    remove(VBDR_FILE);
}

int main(void)
{
    UdpStats netStats;
    double t1, t2;
    int port;

#ifdef _WIN32
    WSADATA wsa;
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        printf("FAIL WSAStartup\n");
        return 1;
    }
#endif
    t1 = now_seconds();
    t2 = now_seconds();
    CHECK(t2 >= t1);

    make_test_pris();
    memset(&netStats, 0, sizeof netStats);
    port = free_port();
    CHECK(port > 0);
    printf("network test on UDP port %d\n", port);
    if (port > 0)
        test_network(port, &netStats);
    printf("replay test\n");
    test_replay(&netStats);
    printf("fuzz test\n");
    test_fuzz();
    test_open_errors();
    remove(RECORD_FILE);
    printf("VBDR test\n");
    test_vbdr();

#ifdef _WIN32
    WSACleanup();
#endif
    if (g_failures) {
        printf("FAIL (%d checks failed)\n", g_failures);
        return 1;
    }
    printf("PASS\n");
    return 0;
}

/*
 * udp_input.c - reception of the DDC output over UDP (or from a record file)
 * and reassembly of the datagrams of every PRI.
 *
 * Datagram = 64-byte header + int16 (I, Q) pairs (little-endian):
 *    0  uint32 magic 0x52444443 ('RDDC')    4  uint32 pulse number
 *    8  uint32 scan number                  12 float  antenna azimuth [deg]
 *   16  uint16 packet index (0-based)       18 uint16 packet count
 *   20  uint16 number of channels           22 uint16 stagger index (1-based)
 *   24  uint32 samples per channel          28 uint32 first complex sample of this packet
 *   32  double PRI start time [s]           40 float  DDC sample rate [Hz]
 *   44  float  time of the first DDC sample 48 uint32 x4 short/long first/last sample
 *   64  int16  I, Q, I, Q, ... channel after channel
 *
 * VBDR stream (udp_input_set_format(in, UDP_FORMAT_VBDR)): every PRI is one stream of N packets,
 * each with a 28-byte label (uint32 stream number, uint32 packet number 1 ... N, uint32 N, 16 bytes
 * ASCII "Loss=<n>") and a part of its half of the stream: packets 1 ... N/2 = half 1 (CH1 & CH2),
 * N/2+1 ... N = half 2 (CH3 & CH4). A half = 32 uint16 system parameters (1 PRF number, 2 azimuth
 * 0 ... 65535 = 0 ... 360 deg, 4 0x5555, 5 previous PRT [us], 6 current PRT [us], ...) + int16
 * CHa_I, CHa_Q, CHb_I, CHb_Q of every DDC sample. The packets of a stream are kept until all N
 * have arrived, then the halves are joined (their sizes found from the packets themselves) and the
 * samples sorted into the channel-after-channel order of PriData (4 channels).
 *
 * Record file: every datagram as [uint32 length (little-endian)][bytes]; it
 * is created (truncated) by udp_input_open.
 *
 * udp_input_next(timeout_ms): 0 does not wait (only the queued datagrams are
 * read), a negative timeout waits forever.
 *
 * Received samples are assembled in an internal buffer. When a PRI is
 * complete, that buffer is handed to the caller's PriData and the caller's
 * previous buffer becomes the new assembly buffer (no copy).
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
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <arpa/inet.h>
#  include <errno.h>
#  include <fcntl.h>
#  include <netinet/in.h>
#  include <sys/select.h>
#  include <sys/socket.h>
#  include <sys/types.h>
#  include <unistd.h>
#endif

#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "radar_rx.h"

#ifdef _WIN32
typedef SOCKET sock_t;
#  define BAD_SOCKET   INVALID_SOCKET
#  define close_socket closesocket
#else
typedef int sock_t;
#  define BAD_SOCKET   (-1)
#  define close_socket close
#endif

#define RDDC_MAGIC      0x52444443u
#define HEADER_BYTES    64
#define MAX_DATAGRAM    65536                  /* receive buffer per datagram          */
#define MAX_PRI_BYTES   ((uint64_t)256 << 20)  /* sanity limit of the samples of a PRI */
#define RCVBUF_BYTES    (32 << 20)             /* wanted socket receive buffer         */
#define FILE_BUF_BYTES  (1 << 20)              /* stdio buffer of record/replay files  */
/* After this many consecutive datagrams of older pulses the sender is assumed
 * to have been restarted (pulse numbers start again) and they are accepted. */
#define RESYNC_AFTER    64
/* Deadline check interval (datagrams) while data keeps arriving. */
#define DEADLINE_EVERY  64

enum { RECV_ERROR = -1, RECV_EMPTY = -2 };

typedef struct {
    uint32_t magic, pulse, scan;
    float    az;
    uint16_t pktIdx, pktCount, nCh, stagger;
    uint32_t nSamples, first;
    double   time;
    float    fs, t0;
    uint32_t idx[4];
} PacketHeader;

struct UdpInput {
    sock_t       sock;           /* BAD_SOCKET in replay mode                       */
    FILE        *record;         /* may be NULL                                     */
    FILE        *replay;         /* NULL in network mode                            */
    UdpStats     stats;
#ifdef _WIN32
    int          wsaStarted;
#endif

    /* PRI being assembled */
    int          havePulse;      /* hdr.pulse holds the current pulse number        */
    int          active;         /* the current PRI is still being received         */
    PacketHeader hdr;            /* header of the first packet of the current PRI   */
    uint32_t     nReceived;      /* distinct packets received                       */
    uint64_t     bytesReceived;  /* payload bytes of the distinct packets           */
    uint64_t     priBytes;       /* expected payload bytes = 4 * nSamples * nCh     */
    int          olderRun;       /* consecutive datagrams of older pulse numbers    */
    int16_t     *buf;            /* assembly buffer                                 */
    size_t       bufCap;         /* allocated int16 values in buf                   */
    uint8_t      got[65536];     /* got[i] = 1 when packet i has been received      */
    /* VBDR stream being assembled */
    int          format;         /* UDP_FORMAT_RDDC or UDP_FORMAT_VBDR              */
    int          vHave, vActive; /* vStream holds a stream number / it is still open */
    uint32_t     vStream, vCount, vReceived, vLoss;
    uint8_t     *vBuf;           /* payloads of the packets in arrival order        */
    size_t       vLen, vCap;
    uint32_t    *vOff, *vSize;   /* [65536] offset and size of packet k (0-based)   */
    uint8_t      dgram[MAX_DATAGRAM];
};

/* ------------------------------------------------------------------------ */
/* Small helpers                                                             */
/* ------------------------------------------------------------------------ */

static void set_err(char *err, size_t errLen, const char *fmt, ...)
{
    va_list ap;
    if (!err || errLen == 0)
        return;
    va_start(ap, fmt);
    vsnprintf(err, errLen, fmt, ap);
    va_end(ap);
}

static int last_socket_error(void)
{
#ifdef _WIN32
    return WSAGetLastError();
#else
    return errno;
#endif
}

static const char *socket_error_text(int code)
{
#ifdef _WIN32
    (void)code;
    return "winsock error";
#else
    return strerror(code);
#endif
}

static uint16_t get_u16(const uint8_t *p) { uint16_t v; memcpy(&v, p, sizeof v); return v; }
static uint32_t get_u32(const uint8_t *p) { uint32_t v; memcpy(&v, p, sizeof v); return v; }
static float    get_f32(const uint8_t *p) { float v;    memcpy(&v, p, sizeof v); return v; }
static double   get_f64(const uint8_t *p) { double v;   memcpy(&v, p, sizeof v); return v; }

/* Explicit little-endian length prefix of the record file. */
static void put_le32(uint8_t *p, uint32_t v)
{
    p[0] = (uint8_t)v;
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
}

static uint32_t get_le32(const uint8_t *p)
{
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

/* ------------------------------------------------------------------------ */
/* Datagram parsing and PRI reassembly                                      */
/* ------------------------------------------------------------------------ */

static void parse_header(const uint8_t *d, PacketHeader *h)
{
    int i;
    h->magic    = get_u32(d + 0);
    h->pulse    = get_u32(d + 4);
    h->scan     = get_u32(d + 8);
    h->az       = get_f32(d + 12);
    h->pktIdx   = get_u16(d + 16);
    h->pktCount = get_u16(d + 18);
    h->nCh      = get_u16(d + 20);
    h->stagger  = get_u16(d + 22);
    h->nSamples = get_u32(d + 24);
    h->first    = get_u32(d + 28);
    h->time     = get_f64(d + 32);
    h->fs       = get_f32(d + 40);
    h->t0       = get_f32(d + 44);
    for (i = 0; i < 4; i++)
        h->idx[i] = get_u32(d + 48 + 4 * i);
}

/* Parses a datagram; returns 0 when it is foreign (bad magic, size or layout). */
static int parse_valid(const uint8_t *d, size_t len, PacketHeader *h)
{
    uint64_t payload, priBytes;
    if (len < HEADER_BYTES)
        return 0;
    parse_header(d, h);
    payload = len - HEADER_BYTES;
    if (h->magic != RDDC_MAGIC || payload % 4 != 0)
        return 0;
    if (h->pktCount == 0 || h->pktIdx >= h->pktCount || h->nCh == 0 || h->nSamples == 0)
        return 0;
    priBytes = (uint64_t)h->nSamples * h->nCh * 4;
    if (priBytes > MAX_PRI_BYTES)
        return 0;
    return (uint64_t)h->first * 4 + payload <= priBytes;
}

static int same_layout(const PacketHeader *a, const PacketHeader *b)
{
    return a->pktCount == b->pktCount && a->nSamples == b->nSamples && a->nCh == b->nCh;
}

/* Starts the reception of a new PRI. Returns 0 if its buffer cannot be allocated. */
static int start_pri(UdpInput *in, const PacketHeader *h)
{
    size_t need = (size_t)h->nSamples * h->nCh * 2;   /* int16 values */

    in->hdr = *h;
    in->havePulse = 1;
    in->active = 0;
    if (need > in->bufCap) {
        int16_t *p = (int16_t *)realloc(in->buf, need * sizeof *p);
        if (!p) {
            in->stats.prisIncomplete++;   /* this PRI is lost */
            return 0;
        }
        in->buf = p;
        in->bufCap = need;
    }
    memset(in->got, 0, h->pktCount);
    in->nReceived = 0;
    in->bytesReceived = 0;
    in->priBytes = (uint64_t)h->nSamples * h->nCh * 4;
    in->active = 1;
    return 1;
}

/* Decides from the pulse number what to do with a valid packet. Returns 1
 * when it belongs to the PRI being received (possibly just started). */
static int accept_pulse(UdpInput *in, const PacketHeader *h)
{
    int32_t age = in->havePulse ? (int32_t)(h->pulse - in->hdr.pulse) : 1;

    if (age < 0 && ++in->olderRun <= RESYNC_AFTER) {
        in->stats.foreign++;              /* late packet of an older PRI */
        return 0;
    }
    if (age == 0) {
        in->olderRun = 0;
        if (!in->active)
            return 0;                     /* late duplicate of a delivered PRI */
        if (!same_layout(&in->hdr, h)) {
            in->stats.foreign++;
            return 0;
        }
        return 1;
    }
    /* newer pulse, or a long run of older ones (the sender was restarted) */
    in->olderRun = 0;
    if (in->active)
        in->stats.prisIncomplete++;
    return start_pri(in, h);
}

/* Fills 'pri' with the completed PRI and swaps the sample buffers. */
static void deliver(UdpInput *in, PriData *pri)
{
    const PacketHeader *h = &in->hdr;
    int16_t *oldIq = pri->iq;
    size_t oldCap = pri->capacity;

    pri->pulseNumber  = h->pulse;
    pri->scanNumber   = h->scan;
    pri->az_deg       = h->az;
    pri->staggerIndex = h->stagger;
    pri->time_s       = h->time;
    pri->ddcFs_Hz     = h->fs;
    pri->ddcT0_s      = h->t0;
    pri->shortIdx[0]  = h->idx[0];
    pri->shortIdx[1]  = h->idx[1];
    pri->longIdx[0]   = h->idx[2];
    pri->longIdx[1]   = h->idx[3];
    pri->nSamples     = (int)h->nSamples;
    pri->nCh          = (int)h->nCh;
    pri->format       = UDP_FORMAT_RDDC;
    pri->iq           = in->buf;
    pri->capacity     = in->bufCap;

    in->buf = oldIq;               /* reused (and grown if needed) for the next PRI */
    in->bufCap = oldIq ? oldCap : 0;
    in->stats.prisComplete++;
}

/* Stores the payload of a packet of the current PRI. Returns 1 when the PRI
 * is complete and has been delivered in 'pri'. */
static int add_packet(UdpInput *in, const PacketHeader *h, size_t len, PriData *pri)
{
    size_t payload = len - HEADER_BYTES;

    if (in->got[h->pktIdx])
        return 0;                                  /* duplicate packet */
    in->got[h->pktIdx] = 1;
    memcpy(in->buf + 2 * (size_t)h->first, in->dgram + HEADER_BYTES, payload);
    in->nReceived++;
    in->bytesReceived += payload;
    if (in->nReceived < in->hdr.pktCount)
        return 0;

    in->active = 0;
    if (in->bytesReceived != in->priBytes) {       /* packets overlap or leave gaps */
        in->stats.prisIncomplete++;
        return 0;
    }
    deliver(in, pri);
    return 1;
}

/* ------------------------------------------------------------------------ */
/* VBDR stream                                                               */
/* ------------------------------------------------------------------------ */

#define VBDR_LABEL 28

/* Joins the packets of the complete stream and delivers it in 'pri'. Returns 1, or 0 when the
 * halves do not fit (counted as incomplete). */
static int vbdr_deliver(UdpInput *in, PriData *pri)
{
    const uint32_t nHalf = in->vCount / 2;
    uint64_t halfBytes[2] = { 0, 0 };
    for (uint32_t k = 0; k < in->vCount; k++) halfBytes[k >= nHalf] += in->vSize[k];
    if (halfBytes[0] != halfBytes[1] || halfBytes[0] < 64 || (halfBytes[0] - 64) % 8 != 0 ||
        halfBytes[0] > MAX_PRI_BYTES / 2) {
        in->stats.prisIncomplete++;
        return 0;
    }
    const size_t nS = (size_t)((halfBytes[0] - 64) / 8);
    const size_t need = nS * 4 * 2;                          /* int16 values of 4 channels */
    if (need > in->bufCap) {
        int16_t *p = (int16_t *)realloc(in->buf, (need ? need : 1) * sizeof *p);
        if (!p) {
            in->stats.prisIncomplete++;
            return 0;
        }
        in->buf = p;
        in->bufCap = need;
    }
    uint8_t par[64];
    for (int h = 0; h < 2; h++) {
        /* walk through the bytes of half h: 64 parameter bytes, then 8 bytes per sample */
        uint64_t pos = 0;                                    /* byte of the half */
        int16_t *ch0 = in->buf + 2 * (size_t)(2 * h) * nS, *ch1 = in->buf + 2 * (size_t)(2 * h + 1) * nS;
        uint8_t smp[8];
        int nb = 0;
        size_t n = 0;
        for (uint32_t k = h * nHalf; k < (uint32_t)(h + 1) * nHalf; k++) {
            const uint8_t *b = in->vBuf + in->vOff[k];
            uint32_t len = in->vSize[k], i = 0;
            if (pos < 64) {                                  /* system parameters */
                const uint32_t m = (uint32_t)(64 - pos) < len ? (uint32_t)(64 - pos) : len;
                if (h == 0) memcpy(par + pos, b, m);
                pos += m;
                i = m;
            }
            for (; i < len; i++) {                           /* samples (also across packet edges) */
                smp[nb++] = b[i];
                if (nb == 8) {
                    memcpy(&ch0[2 * n], smp, 4);             /* CHa I, Q */
                    memcpy(&ch1[2 * n], smp + 4, 4);         /* CHb I, Q */
                    n++;
                    nb = 0;
                }
            }
        }
    }
    int16_t *oldIq = pri->iq;
    size_t oldCap = pri->capacity;
    pri->scanNumber   = 0;                                   /* vbdr_complete_pri() */
    pri->time_s       = 0;
    pri->ddcFs_Hz     = 0;
    pri->ddcT0_s      = 0;
    pri->shortIdx[0]  = pri->shortIdx[1] = pri->longIdx[0] = pri->longIdx[1] = 0;
    pri->pulseNumber  = in->vStream;
    pri->az_deg       = (float)(get_u16(par + 2) * (360.0 / 65536.0));
    pri->staggerIndex = 1;
    pri->nSamples     = (int)nS;
    pri->nCh          = 4;
    pri->format       = UDP_FORMAT_VBDR;
    pri->prfNumber    = get_u16(par + 0);
    pri->prevPrt_s    = get_u16(par + 8) * 1e-6;
    pri->curPrt_s     = get_u16(par + 10) * 1e-6;
    pri->loss         = in->vLoss;
    pri->iq           = in->buf;
    pri->capacity     = in->bufCap;
    in->buf = oldIq;
    in->bufCap = oldIq ? oldCap : 0;
    in->stats.prisComplete++;
    return 1;
}

/* Handles a datagram of the VBDR stream. Returns 1 when a PRI was delivered. */
static int vbdr_datagram(UdpInput *in, size_t len, PriData *pri)
{
    if (len <= VBDR_LABEL) { in->stats.foreign++; return 0; }
    const uint32_t stream = get_u32(in->dgram), pk = get_u32(in->dgram + 4), count = get_u32(in->dgram + 8);
    if (count < 2 || count % 2 || count > 65534 || pk < 1 || pk > count) { in->stats.foreign++; return 0; }
    if (!in->vOff) {
        in->vOff = (uint32_t *)malloc(sizeof(uint32_t) * 65536);
        in->vSize = (uint32_t *)malloc(sizeof(uint32_t) * 65536);
        if (!in->vOff || !in->vSize) return 0;
    }
    const int32_t age = in->vHave ? (int32_t)(stream - in->vStream) : 1;
    if (age < 0 && ++in->olderRun <= RESYNC_AFTER) { in->stats.foreign++; return 0; }   /* late packet */
    if (age == 0) {
        in->olderRun = 0;
        if (!in->vActive) return 0;                          /* late duplicate of a delivered stream */
        if (count != in->vCount) { in->stats.foreign++; return 0; }
    } else {                                                 /* a new stream (or the sender restarted) */
        in->olderRun = 0;
        if (in->vActive) in->stats.prisIncomplete++;
        in->vHave = 1;
        in->vActive = 1;
        in->vStream = stream;
        in->vCount = count;
        in->vReceived = 0;
        in->vLen = 0;
        memset(in->got, 0, count);
        char txt[17];
        memcpy(txt, in->dgram + 12, 16);
        txt[16] = 0;
        unsigned loss = 0;
        in->vLoss = sscanf(txt, "Loss=%u", &loss) == 1 ? loss : 0;
        in->stats.boardLoss = in->vLoss;
    }
    if (in->got[pk - 1]) return 0;                           /* duplicate packet */
    const size_t payload = len - VBDR_LABEL;
    if (in->vLen + payload > in->vCap) {
        size_t cap = in->vCap ? 2 * in->vCap : (size_t)1 << 20;
        while (cap < in->vLen + payload) cap *= 2;
        uint8_t *p = (uint8_t *)realloc(in->vBuf, cap);
        if (!p) return 0;
        in->vBuf = p;
        in->vCap = cap;
    }
    memcpy(in->vBuf + in->vLen, in->dgram + VBDR_LABEL, payload);
    in->got[pk - 1] = 1;
    in->vOff[pk - 1] = (uint32_t)in->vLen;
    in->vSize[pk - 1] = (uint32_t)payload;
    in->vLen += payload;
    if (++in->vReceived < in->vCount) return 0;
    in->vActive = 0;
    return vbdr_deliver(in, pri);
}

static void record_datagram(UdpInput *in, size_t len)
{
    uint8_t prefix[4];
    if (!in->record)
        return;
    put_le32(prefix, (uint32_t)len);
    if (fwrite(prefix, 1, 4, in->record) != 4 || fwrite(in->dgram, 1, len, in->record) != len) {
        fprintf(stderr, "udp_input: write error on the record file, recording stopped\n");
        fclose(in->record);
        in->record = NULL;
    }
}

/* Handles the datagram in in->dgram. Returns 1 when a PRI was delivered. */
static int handle_datagram(UdpInput *in, size_t len, PriData *pri)
{
    PacketHeader h;

    in->stats.datagrams++;
    record_datagram(in, len);
    if (in->format == UDP_FORMAT_VBDR)
        return vbdr_datagram(in, len, pri);
    if (!parse_valid(in->dgram, len, &h)) {
        in->stats.foreign++;
        return 0;
    }
    if (!accept_pulse(in, &h))
        return 0;
    return add_packet(in, &h, len, pri);
}

/* ------------------------------------------------------------------------ */
/* Network                                                                   */
/* ------------------------------------------------------------------------ */

static void set_receive_buffer(sock_t s)
{
    int want = RCVBUF_BYTES, got = 0, forced = 0;
    socklen_t len = (socklen_t)sizeof got;

#ifdef SO_RCVBUFFORCE
    /* Linux: exceeds net.core.rmem_max when the process has CAP_NET_ADMIN */
    forced = setsockopt(s, SOL_SOCKET, SO_RCVBUFFORCE, &want, sizeof want) == 0;
#endif
    if (!forced)
        setsockopt(s, SOL_SOCKET, SO_RCVBUF, (const char *)&want, (socklen_t)sizeof want);

    if (getsockopt(s, SOL_SOCKET, SO_RCVBUF, (char *)&got, &len) != 0)
        return;
#ifdef __linux__
    got /= 2;                     /* Linux reports twice the size it was given */
#endif
    if (got < want)
        fprintf(stderr, "udp_input: warning: UDP receive buffer is %d KB instead of %d KB, "
                        "packets may be lost (Linux: sysctl -w net.core.rmem_max=%d)\n",
                got / 1024, want / 1024, want);
}

static int set_nonblocking(sock_t s)
{
#ifdef _WIN32
    u_long on = 1;
    return ioctlsocket(s, FIONBIO, &on) == 0;
#else
    int flags = fcntl(s, F_GETFL, 0);
    return flags >= 0 && fcntl(s, F_SETFL, flags | O_NONBLOCK) == 0;
#endif
}

static int open_socket(UdpInput *in, int port, char *err, size_t errLen)
{
    struct sockaddr_in addr;

#ifdef _WIN32
    WSADATA wsa;
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        set_err(err, errLen, "WSAStartup failed");
        return 0;
    }
    in->wsaStarted = 1;
#endif
    in->sock = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (in->sock == BAD_SOCKET) {
        int e = last_socket_error();
        set_err(err, errLen, "socket(): %s (%d)", socket_error_text(e), e);
        return 0;
    }
#ifndef _WIN32
    if (in->sock >= FD_SETSIZE) {
        set_err(err, errLen, "socket descriptor %d too large for select()", in->sock);
        return 0;
    }
#endif
    set_receive_buffer(in->sock);

    memset(&addr, 0, sizeof addr);
    addr.sin_family = AF_INET;
    addr.sin_addr.s_addr = htonl(INADDR_ANY);
    addr.sin_port = htons((uint16_t)port);
    if (bind(in->sock, (struct sockaddr *)&addr, (socklen_t)sizeof addr) != 0) {
        int e = last_socket_error();
        set_err(err, errLen, "bind() to UDP port %d: %s (%d)", port, socket_error_text(e), e);
        return 0;
    }
    if (!set_nonblocking(in->sock)) {
        set_err(err, errLen, "cannot make the UDP socket non-blocking");
        return 0;
    }
    return 1;
}

/* Reads one datagram into in->dgram without waiting. Returns its length,
 * RECV_EMPTY when none is queued or RECV_ERROR. */
static long recv_datagram(UdpInput *in)
{
    for (;;) {
#ifdef _WIN32
        int n = recv(in->sock, (char *)in->dgram, MAX_DATAGRAM, 0);
        int e;
        if (n >= 0)
            return n;
        e = WSAGetLastError();
        if (e == WSAEWOULDBLOCK)
            return RECV_EMPTY;
        if (e == WSAEMSGSIZE) {          /* truncated datagram: not ours */
            in->stats.datagrams++;
            in->stats.foreign++;
            continue;
        }
        if (e == WSAECONNRESET || e == WSAEINTR)
            continue;
        return RECV_ERROR;
#else
        ssize_t n = recv(in->sock, in->dgram, MAX_DATAGRAM, 0);
        if (n >= 0)
            return (long)n;
        if (errno == EAGAIN || errno == EWOULDBLOCK)
            return RECV_EMPTY;
        if (errno == EINTR || errno == ECONNREFUSED)
            continue;
        return RECV_ERROR;
#endif
    }
}

/* Waits until the socket is readable. ms < 0 waits forever.
 * Returns 1 when readable (or interrupted), 0 on timeout, -1 on error. */
static int wait_readable(sock_t s, double ms)
{
    fd_set set;
    struct timeval tv, *ptv = NULL;
    int r;

    FD_ZERO(&set);
    FD_SET(s, &set);
    if (ms >= 0) {
        long sec = (long)(ms / 1000.0);
        tv.tv_sec = sec;
        tv.tv_usec = (long)((ms - 1000.0 * (double)sec) * 1000.0);
        ptv = &tv;
    }
#ifdef _WIN32
    r = select(0, &set, NULL, NULL, ptv);
    if (r == SOCKET_ERROR)
        return WSAGetLastError() == WSAEINTR ? 1 : -1;
#else
    r = select(s + 1, &set, NULL, NULL, ptv);
    if (r < 0)
        return errno == EINTR ? 1 : -1;
#endif
    return r > 0 ? 1 : 0;
}

static int next_from_network(UdpInput *in, PriData *pri, int timeout_ms)
{
    double deadline = now_seconds() + 1e-3 * timeout_ms;
    unsigned handled = 0;

    for (;;) {
        long n = recv_datagram(in);
        if (n == RECV_ERROR)
            return -1;
        if (n == RECV_EMPTY) {
            double left_ms = timeout_ms < 0 ? -1.0 : 1e3 * (deadline - now_seconds());
            int r;
            if (timeout_ms >= 0 && left_ms <= 0)
                return 0;
            r = wait_readable(in->sock, left_ms);
            if (r < 0)
                return -1;
            /* readable, or select() timed out: the deadline is checked again
             * above, so an early wake-up (Windows timer ticks) waits on */
            continue;
        }
        if (handle_datagram(in, (size_t)n, pri))
            return 1;
        /* keep the timeout even when foreign datagrams keep arriving */
        if (timeout_ms >= 0 && ++handled % DEADLINE_EVERY == 0 && now_seconds() >= deadline)
            return 0;
    }
}

/* ------------------------------------------------------------------------ */
/* Replay                                                                    */
/* ------------------------------------------------------------------------ */

/* Reads the next datagram of the replay file into in->dgram. Returns its
 * length, or -1 at the end of the file (or when the file is corrupted). */
static long read_replay(UdpInput *in)
{
    uint8_t prefix[4];
    uint32_t len;

    if (fread(prefix, 1, 4, in->replay) != 4)
        return -1;
    len = get_le32(prefix);
    if (len > MAX_DATAGRAM)
        return -1;
    if (len > 0 && fread(in->dgram, 1, len, in->replay) != len)
        return -1;
    return (long)len;
}

static int next_from_replay(UdpInput *in, PriData *pri)
{
    for (;;) {
        long n = read_replay(in);
        if (n < 0) {
            if (in->active || in->vActive) {   /* the last PRI never completed */
                in->stats.prisIncomplete++;
                in->active = 0;
                in->vActive = 0;
            }
            return -1;
        }
        if (handle_datagram(in, (size_t)n, pri))
            return 1;
    }
}

/* ------------------------------------------------------------------------ */
/* Public functions                                                          */
/* ------------------------------------------------------------------------ */

UdpInput *udp_input_open(int port, const char *recordPath, const char *replayPath, char *err, size_t errLen)
{
    UdpInput *in;

    set_err(err, errLen, "%s", "");
    if (!replayPath && (port < 0 || port > 65535)) {
        set_err(err, errLen, "invalid UDP port %d", port);
        return NULL;
    }
    if (recordPath && replayPath && strcmp(recordPath, replayPath) == 0) {
        set_err(err, errLen, "the record and the replay file must differ");
        return NULL;
    }
    in = (UdpInput *)calloc(1, sizeof *in);
    if (!in) {
        set_err(err, errLen, "out of memory");
        return NULL;
    }
    in->sock = BAD_SOCKET;

    if (replayPath) {
        in->replay = fopen(replayPath, "rb");
        if (!in->replay) {
            set_err(err, errLen, "cannot open the replay file '%s'", replayPath);
            goto fail;
        }
        setvbuf(in->replay, NULL, _IOFBF, FILE_BUF_BYTES);
    } else if (!open_socket(in, port, err, errLen)) {
        goto fail;
    }

    if (recordPath) {
        in->record = fopen(recordPath, "wb");
        if (!in->record) {
            set_err(err, errLen, "cannot create the record file '%s'", recordPath);
            goto fail;
        }
        setvbuf(in->record, NULL, _IOFBF, FILE_BUF_BYTES);
    }
    return in;

fail:
    udp_input_close(in);
    return NULL;
}

int udp_input_next(UdpInput *in, PriData *pri, int timeout_ms)
{
    if (!in || !pri)
        return -1;
    if (in->replay)
        return next_from_replay(in, pri);
    return next_from_network(in, pri, timeout_ms);
}

void udp_input_stats(const UdpInput *in, UdpStats *stats)
{
    if (!stats)
        return;
    if (in)
        *stats = in->stats;
    else
        memset(stats, 0, sizeof *stats);
}

void udp_input_close(UdpInput *in)
{
    if (!in)
        return;
    if (in->sock != BAD_SOCKET)
        close_socket(in->sock);
#ifdef _WIN32
    if (in->wsaStarted)
        WSACleanup();
#endif
    if (in->record)
        fclose(in->record);
    if (in->replay)
        fclose(in->replay);
    free(in->buf);
    free(in->vBuf);
    free(in->vOff);
    free(in->vSize);
    free(in);
}

void udp_input_set_format(UdpInput *in, int format)
{
    if (in) in->format = format == UDP_FORMAT_VBDR ? UDP_FORMAT_VBDR : UDP_FORMAT_RDDC;
}

void pri_data_free(PriData *pri)
{
    if (!pri)
        return;
    free(pri->iq);
    memset(pri, 0, sizeof *pri);
}

/* ------------------------------------------------------------------------ */
/* VBDR: the fields of a PRI the stream does not carry                       */
/* ------------------------------------------------------------------------ */

/* Stagger position (1-based) of a PRT: the expected next position when it fits, else the closest
 * PRI of the sequence (within 1 us: the PRT is sent in whole us); 'expected' when none fits. */
static int stagger_of(const double *seq, int n, double prt, int expected)
{
    const double tol = 1.0e-6;
    if (expected >= 1 && expected <= n && fabs(seq[expected - 1] - prt) <= tol) return expected;
    int best = 0;
    double bestErr = tol;
    for (int k = 0; k < n; k++) {
        const double e = fabs(seq[k] - prt);
        if (e <= bestErr) { bestErr = e; best = k + 1; if (e == 0) break; }
    }
    return best ? best : expected;
}

void vbdr_complete_pri(VbdrTiming *st, const RadarConfig *cfg, PriData *pri)
{
    /* layout of the DDC record (the configuration of the mode) */
    pri->ddcFs_Hz = (float)cfg->fsOut_Hz;
    pri->ddcT0_s = (float)cfg->ddcT0_s;
    const PulseConfig *ps = &cfg->pulse[PULSE_SHORT], *pl = &cfg->pulse[PULSE_LONG];
    pri->shortIdx[0] = ps->enable && ps->rxCount > 0 ? (uint32_t)ps->rxFirst + 1 : 0;
    pri->shortIdx[1] = ps->enable && ps->rxCount > 0 ? (uint32_t)(ps->rxFirst + ps->rxCount) : 0;
    pri->longIdx[0]  = pl->enable && pl->rxCount > 0 ? (uint32_t)pl->rxFirst + 1 : 0;
    pri->longIdx[1]  = pl->enable && pl->rxCount > 0 ? (uint32_t)(pl->rxFirst + pl->rxCount) : 0;

    /* streams between the last PRI and this one */
    const int32_t gap = st->have ? (int32_t)(pri->pulseNumber - st->lastStream) : 1;
    const int nSeq = cfg->nPriSequence;
    const double meanPrt = cfg->meanPri_s > 0 ? cfg->meanPri_s : pri->curPrt_s;

    /* stagger index from the current PRT */
    int stagger = 1;
    if (cfg->scan2scan && cfg->nPrf > 0) {
        stagger = st->have ? st->lastStagger : 1;
        double bestErr = 1.0e-6;
        for (int j = 0; j < cfg->nPrf; j++) {
            const double e = fabs(cfg->prf[j].pri_s - pri->curPrt_s);
            if (e <= bestErr) { bestErr = e; stagger = j + 1; }
        }
    } else if (cfg->staggered && nSeq > 1) {
        const int expected = st->have && gap > 0 ? (int)((st->lastStagger - 1 + (int64_t)gap) % nSeq) + 1 : 1;
        stagger = stagger_of(cfg->priSequence_s, nSeq, pri->curPrt_s, expected);
    }
    pri->staggerIndex = (uint16_t)(stagger > 0 ? stagger : 1);

    /* PRI start time: the sum of the PRTs (the first PRI at 0) */
    double t = 0;
    if (st->have) {
        if (gap == 1) {
            t = st->lastTime_s + (pri->prevPrt_s > 0 ? pri->prevPrt_s : st->lastCurPrt_s);
        } else if (gap > 1) {
            /* the last PRI lasted lastCurPrt; the PRIs that never arrived (all but the one just before
             * this, which lasted prevPrt) are estimated: the stagger sequence or the last PRT */
            t = st->lastTime_s + st->lastCurPrt_s;
            const int64_t nUnknown = (int64_t)gap - 2;
            if (cfg->staggered && nSeq > 1 && nUnknown < 1000000) {
                for (int64_t j = 1; j <= nUnknown; j++) t += cfg->priSequence_s[(st->lastStagger - 1 + j) % nSeq];
            } else {
                t += (double)nUnknown * (cfg->staggered ? meanPrt : st->lastCurPrt_s);
            }
            t += pri->prevPrt_s > 0 ? pri->prevPrt_s : st->lastCurPrt_s;
            st->missing += (uint64_t)(gap - 1);
        } else {                                             /* older stream number: the sender restarted */
            t = st->lastTime_s + st->lastCurPrt_s;
        }
    }
    pri->time_s = t;

    /* scan number from the azimuth wrap */
    if (!st->have) st->scan = 1;
    else if (pri->az_deg < st->lastAz_deg - 180.0f) st->scan++;
    pri->scanNumber = st->scan;

    st->have = 1;
    st->lastStream = pri->pulseNumber;
    st->lastTime_s = t;
    st->lastCurPrt_s = pri->curPrt_s > 0 ? pri->curPrt_s : meanPrt;
    st->lastAz_deg = pri->az_deg;
    st->lastStagger = pri->staggerIndex;
}

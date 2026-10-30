/*
 * serial_out.c - serial output of the CFAR detections and plots.
 *
 * Frame format (all values little-endian, see README):
 *
 *   Detections frame (type 1), 16 + 24*N bytes, at most 200 detections:
 *      0: A5 5A (sync)   2: u8 type = 1   3: u8 pulse (0 = short, 1 = long)
 *      4: u32 pulse number   8: f32 azimuth [deg]   12: u16 N
 *     14: N records of 24 bytes:
 *         u16 beam, u16 range cell, u16 Doppler bin, u16 reserved (0),
 *         f32 range [m], f32 radial velocity [m/s], f32 elevation [deg], f32 SNR [dB]
 *     14+24N: u16 CRC
 *
 *   Plot frame (type 2), 42 bytes:
 *      0: A5 5A   2: u8 type = 2   3: u8 pulse   4: u32 pulse number
 *      8: f64 time [s]   16: f32 azimuth [deg]   20: f32 range [m]
 *     24: f32 elevation [deg]   28: f32 radial velocity [m/s]   32: f32 SNR [dB]
 *     36: u16 beam   38: u16 hits   40: u16 CRC
 *
 *   CRC-16/CCITT-FALSE (poly 0x1021, init 0xFFFF, no reflection, xorout 0)
 *   over the bytes after the sync up to the byte before the CRC.
 *
 * The output never blocks: frames are queued in a ring buffer (a frame is
 * queued completely or dropped and counted), then the ring is written to the
 * port as far as the port accepts data without blocking. What is left is
 * written by the next call. serial_send_detections() with n == 0 sends no
 * frame but still writes pending bytes, so calling it every CPI update keeps
 * the ring moving. A SerialOut is not thread-safe: use it from one thread.
 */
#if !defined(_WIN32) && !defined(_DEFAULT_SOURCE)
#define _DEFAULT_SOURCE   /* termios extensions (B921600, CRTSCTS) with -std=c11 */
#endif

#include "radar_rx.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#else
#include <errno.h>
#include <fcntl.h>
#include <poll.h>
#include <sys/ioctl.h>
#include <termios.h>
#include <time.h>
#include <unistd.h>
#endif

enum {
    SYNC0 = 0xA5,
    SYNC1 = 0x5A,
    TYPE_DETECTIONS = 1,
    TYPE_PLOT = 2,
    DET_HEADER_LEN = 14,          /* sync .. number of detections */
    DET_RECORD_LEN = 24,
    CRC_LEN = 2,
    PLOT_FRAME_LEN = 42,
    MAX_DETS_PER_FRAME = 200,
    MAX_DET_FRAME_LEN = DET_HEADER_LEN + MAX_DETS_PER_FRAME * DET_RECORD_LEN + CRC_LEN,
    RING_SIZE = 256 * 1024,
    WIN_IN_QUEUE = 4096,          /* Windows driver queues */
    WIN_OUT_QUEUE = 64 * 1024,
    WIN_WRITE_TIMEOUT_MS = 2      /* longest time a WriteFile call may take */
};

#define CLOSE_TIMEOUT_S 2.0       /* longest time serial_close waits for the port */

struct SerialOut {
#ifdef _WIN32
    HANDLE   h;
#else
    int      fd;
#endif
    uint8_t *ring;                /* RING_SIZE bytes waiting for the port */
    size_t   head;                /* index of the oldest byte in the ring */
    size_t   count;               /* bytes in the ring */
    uint64_t framesSent;          /* frames queued */
    uint64_t framesDropped;       /* frames that did not fit in the ring */
    uint64_t bytesSent;           /* bytes accepted by the port */
    uint8_t  frame[MAX_DET_FRAME_LEN];   /* encoding buffer */
};

static void set_err(char *err, size_t errLen, const char *fmt, ...)
{
    va_list ap;
    if (!err || errLen == 0)
        return;
    va_start(ap, fmt);
    vsnprintf(err, errLen, fmt, ap);
    va_end(ap);
}

/* ------------------------------------------------------------------------ */
/* CRC and frame encoding                                                    */
/* ------------------------------------------------------------------------ */

uint16_t crc16_ccitt(const uint8_t *data, size_t len)
{
    uint16_t crc = 0xFFFF;
    for (size_t i = 0; i < len; i++) {
        crc ^= (uint16_t)(data[i] << 8);
        for (int bit = 0; bit < 8; bit++)
            crc = (uint16_t)((crc & 0x8000u) ? ((unsigned)crc << 1) ^ 0x1021u : (unsigned)crc << 1);
    }
    return crc;
}

/* The put_* helpers write one little-endian value and return the next position. */
static uint8_t *put_u8(uint8_t *p, unsigned v)
{
    p[0] = (uint8_t)v;
    return p + 1;
}

static uint8_t *put_u16(uint8_t *p, unsigned v)
{
    p[0] = (uint8_t)v;
    p[1] = (uint8_t)(v >> 8);
    return p + 2;
}

static uint8_t *put_u32(uint8_t *p, uint32_t v)
{
    p[0] = (uint8_t)v;
    p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
    return p + 4;
}

static uint8_t *put_f32(uint8_t *p, double v)
{
    float f = (float)v;
    uint32_t u;
    memcpy(&u, &f, sizeof u);
    return put_u32(p, u);
}

static uint8_t *put_f64(uint8_t *p, double v)
{
    uint64_t u;
    memcpy(&u, &v, sizeof u);
    p = put_u32(p, (uint32_t)u);
    return put_u32(p, (uint32_t)(u >> 32));
}

static unsigned to_u16(int v)
{
    return v < 0 ? 0u : v > 0xFFFF ? 0xFFFFu : (unsigned)v;
}

/* Sync, type, pulse and pulse number (bytes 0..7 of both frame types). */
static uint8_t *put_header(uint8_t *p, unsigned type, int pulse, uint32_t pulseNumber)
{
    p = put_u8(p, SYNC0);
    p = put_u8(p, SYNC1);
    p = put_u8(p, type);
    p = put_u8(p, (unsigned)pulse & 0xFFu);
    return put_u32(p, pulseNumber);
}

static uint8_t *put_detection(uint8_t *p, const Detection *d)
{
    p = put_u16(p, to_u16(d->beam));
    p = put_u16(p, to_u16(d->rIdx));
    p = put_u16(p, to_u16(d->dIdx));
    p = put_u16(p, 0);                    /* reserved */
    p = put_f32(p, d->range_m);
    p = put_f32(p, d->vr_mps);
    p = put_f32(p, d->el_deg);
    return put_f32(p, d->snr_dB);
}

/* Writes the CRC of a complete frame of 'len' bytes into its last two bytes. */
static void put_crc(uint8_t *frame, size_t len)
{
    put_u16(frame + len - CRC_LEN, crc16_ccitt(frame + 2, len - 2 - CRC_LEN));
}

size_t serial_encode_detections(uint8_t *buf, size_t bufLen, int pulse, uint32_t pulseNumber,
                                double az_deg, const Detection *dets, int n)
{
    size_t len;
    uint8_t *p;

    if (!buf || n < 0 || n > 0xFFFF || (n > 0 && !dets))
        return 0;
    len = DET_HEADER_LEN + (size_t)n * DET_RECORD_LEN + CRC_LEN;
    if (bufLen < len)
        return 0;

    p = put_header(buf, TYPE_DETECTIONS, pulse, pulseNumber);
    p = put_f32(p, az_deg);
    p = put_u16(p, (unsigned)n);
    for (int k = 0; k < n; k++)
        p = put_detection(p, &dets[k]);
    put_crc(buf, len);
    return len;
}

size_t serial_encode_plot(uint8_t *buf, size_t bufLen, const Detection *plot)
{
    uint8_t *p;

    if (!buf || !plot || bufLen < PLOT_FRAME_LEN)
        return 0;
    p = put_header(buf, TYPE_PLOT, plot->pulse, plot->pulseNumber);
    p = put_f64(p, plot->time_s);
    p = put_f32(p, plot->az_deg);
    p = put_f32(p, plot->range_m);
    p = put_f32(p, plot->el_deg);
    p = put_f32(p, plot->vr_mps);
    p = put_f32(p, plot->snr_dB);
    p = put_u16(p, to_u16(plot->beam));
    put_u16(p, to_u16(plot->hits));
    put_crc(buf, PLOT_FRAME_LEN);
    return PLOT_FRAME_LEN;
}

/* ------------------------------------------------------------------------ */
/* Port access, Windows                                                      */
/* ------------------------------------------------------------------------ */

#ifdef _WIN32

/* QueryPerformanceCounter works on every Windows version; GetTickCount64 is
   not declared when _WIN32_WINNT < 0x0600. */
static double mono_seconds(void)
{
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    return (double)count.QuadPart / (double)freq.QuadPart;
}

static int configure_port(HANDLE h, int baud)
{
    DCB dcb;
    COMMTIMEOUTS to;

    SetupComm(h, WIN_IN_QUEUE, WIN_OUT_QUEUE);   /* a hint, some drivers ignore it */
    memset(&dcb, 0, sizeof dcb);
    dcb.DCBlength = sizeof dcb;
    if (!GetCommState(h, &dcb))
        return 0;
    dcb.BaudRate = (DWORD)baud;
    dcb.ByteSize = 8;
    dcb.Parity = NOPARITY;
    dcb.StopBits = ONESTOPBIT;
    dcb.fBinary = TRUE;
    dcb.fParity = FALSE;
    dcb.fOutxCtsFlow = FALSE;                    /* no flow control */
    dcb.fOutxDsrFlow = FALSE;
    dcb.fDsrSensitivity = FALSE;
    dcb.fOutX = FALSE;
    dcb.fInX = FALSE;
    dcb.fNull = FALSE;
    dcb.fAbortOnError = FALSE;
    dcb.fDtrControl = DTR_CONTROL_ENABLE;
    dcb.fRtsControl = RTS_CONTROL_ENABLE;
    if (!SetCommState(h, &dcb))
        return 0;

    memset(&to, 0, sizeof to);
    to.ReadIntervalTimeout = MAXDWORD;           /* reads return at once (not used) */
    to.WriteTotalTimeoutMultiplier = 0;
    to.WriteTotalTimeoutConstant = WIN_WRITE_TIMEOUT_MS;   /* WriteFile never blocks long */
    return SetCommTimeouts(h, &to) ? 1 : 0;
}

static int port_open(SerialOut *s, const char *device, int baud, char *err, size_t errLen)
{
    char path[MAX_PATH];

    if (strncmp(device, "\\\\.\\", 4) == 0)
        snprintf(path, sizeof path, "%s", device);
    else
        snprintf(path, sizeof path, "\\\\.\\%s", device);   /* needed for COM10 and up */

    s->h = CreateFileA(path, GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
    if (s->h == INVALID_HANDLE_VALUE) {
        set_err(err, errLen, "cannot open %s (error %lu)", path, (unsigned long)GetLastError());
        return -1;
    }
    if (!configure_port(s->h, baud)) {
        set_err(err, errLen, "cannot configure %s at %d baud (error %lu)", path, baud,
                (unsigned long)GetLastError());
        CloseHandle(s->h);
        return -1;
    }
    PurgeComm(s->h, PURGE_TXCLEAR | PURGE_RXCLEAR);
    return 0;
}

/* Bytes in the driver output queue (not yet transmitted). */
static size_t port_pending(const SerialOut *s)
{
    COMSTAT st;
    DWORD errors = 0;
    return ClearCommError(s->h, &errors, &st) ? (size_t)st.cbOutQue : 0;
}

/* Writes what the port takes without blocking; returns the bytes written. */
static size_t port_write(SerialOut *s, const uint8_t *data, size_t len)
{
    COMSTAT st;
    DWORD errors = 0, written = 0;

    /* Limit the write to the free space of the driver queue; the short write
       timeout covers drivers that report no queue size. */
    if (ClearCommError(s->h, &errors, &st)) {
        size_t room = st.cbOutQue < WIN_OUT_QUEUE ? WIN_OUT_QUEUE - st.cbOutQue : 0;
        if (len > room)
            len = room;
    }
    if (len == 0)
        return 0;
    if (!WriteFile(s->h, data, (DWORD)len, &written, NULL) && written == 0)
        return 0;   /* timeout or error: try again later */
    return (size_t)written;
}

static void port_wait(const SerialOut *s)
{
    (void)s;
    Sleep(1);
}

static void port_discard(SerialOut *s)
{
    PurgeComm(s->h, PURGE_TXABORT | PURGE_TXCLEAR);
}

static void port_close(SerialOut *s)
{
    CloseHandle(s->h);
}

/* ------------------------------------------------------------------------ */
/* Port access, Linux and other POSIX systems                               */
/* ------------------------------------------------------------------------ */

#else

static double mono_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* termios constant of a baud rate, B0 if not supported. */
static speed_t baud_constant(int baud)
{
    switch (baud) {
    case 1200:    return B1200;
    case 2400:    return B2400;
    case 4800:    return B4800;
    case 9600:    return B9600;
    case 19200:   return B19200;
    case 38400:   return B38400;
#ifdef B57600
    case 57600:   return B57600;
#endif
#ifdef B115200
    case 115200:  return B115200;
#endif
#ifdef B230400
    case 230400:  return B230400;
#endif
#ifdef B460800
    case 460800:  return B460800;
#endif
#ifdef B500000
    case 500000:  return B500000;
#endif
#ifdef B576000
    case 576000:  return B576000;
#endif
#ifdef B921600
    case 921600:  return B921600;
#endif
#ifdef B1000000
    case 1000000: return B1000000;
#endif
#ifdef B1500000
    case 1500000: return B1500000;
#endif
#ifdef B2000000
    case 2000000: return B2000000;
#endif
#ifdef B3000000
    case 3000000: return B3000000;
#endif
    default:      return B0;
    }
}

/* Raw 8N1 without flow control, reads return at once. */
static void make_raw(struct termios *t)
{
    t->c_iflag &= ~(tcflag_t)(IGNBRK | BRKINT | PARMRK | ISTRIP | INLCR | IGNCR | ICRNL |
                              IXON | IXOFF | IXANY);
    t->c_oflag &= ~(tcflag_t)OPOST;
    t->c_lflag &= ~(tcflag_t)(ECHO | ECHONL | ICANON | ISIG | IEXTEN);
    t->c_cflag &= ~(tcflag_t)(CSIZE | PARENB | CSTOPB);
#ifdef CRTSCTS
    t->c_cflag &= ~(tcflag_t)CRTSCTS;
#endif
    t->c_cflag |= CS8 | CLOCAL | CREAD;
    t->c_cc[VMIN] = 0;
    t->c_cc[VTIME] = 0;
}

static int port_open(SerialOut *s, const char *device, int baud, char *err, size_t errLen)
{
    struct termios tio;
    speed_t speed = baud_constant(baud);

    if (speed == B0) {
        set_err(err, errLen, "unsupported baud rate %d", baud);
        return -1;
    }
    s->fd = open(device, O_RDWR | O_NOCTTY | O_NONBLOCK);
    if (s->fd < 0) {
        set_err(err, errLen, "cannot open %s: %s", device, strerror(errno));
        return -1;
    }
    if (tcgetattr(s->fd, &tio) != 0) {
        set_err(err, errLen, "%s is not a serial port: %s", device, strerror(errno));
        close(s->fd);
        return -1;
    }
    make_raw(&tio);
    if (cfsetispeed(&tio, speed) != 0 || cfsetospeed(&tio, speed) != 0 ||
        tcsetattr(s->fd, TCSANOW, &tio) != 0) {
        set_err(err, errLen, "cannot configure %s at %d baud: %s", device, baud, strerror(errno));
        close(s->fd);
        return -1;
    }
    return 0;
}

/* Bytes accepted by the driver but not yet transmitted (0 if unknown). */
static size_t port_pending(const SerialOut *s)
{
#ifdef TIOCOUTQ
    int n = 0;
    if (ioctl(s->fd, TIOCOUTQ, &n) == 0 && n > 0)
        return (size_t)n;
#else
    (void)s;
#endif
    return 0;
}

/* Writes what the port takes without blocking; returns the bytes written. */
static size_t port_write(SerialOut *s, const uint8_t *data, size_t len)
{
    for (;;) {
        ssize_t n = write(s->fd, data, len);
        if (n >= 0)
            return (size_t)n;
        if (errno != EINTR)
            return 0;   /* EAGAIN: the port is full; other errors: try again later */
    }
}

/* Waits up to a few ms until the port can take data (or just sleeps). */
static void port_wait(const SerialOut *s)
{
    struct pollfd pfd;
    pfd.fd = s->fd;
    pfd.events = POLLOUT;
    pfd.revents = 0;
    if (s->count > 0)
        poll(&pfd, 1, 10);
    else
        poll(NULL, 0, 2);   /* the driver is still transmitting */
}

static void port_discard(SerialOut *s)
{
    tcflush(s->fd, TCOFLUSH);
}

static void port_close(SerialOut *s)
{
    close(s->fd);
}

#endif

/* ------------------------------------------------------------------------ */
/* Ring buffer                                                               */
/* ------------------------------------------------------------------------ */

static size_t ring_free(const SerialOut *s)
{
    return RING_SIZE - s->count;
}

/* Appends len bytes; the caller has checked that they fit. */
static void ring_put(SerialOut *s, const uint8_t *data, size_t len)
{
    size_t tail = (s->head + s->count) % RING_SIZE;
    size_t first = RING_SIZE - tail;
    if (first > len)
        first = len;
    memcpy(s->ring + tail, data, first);
    memcpy(s->ring, data + first, len - first);
    s->count += len;
}

/* Writes as much of the ring as the port takes without blocking. */
static void ring_flush(SerialOut *s)
{
    while (s->count > 0) {
        size_t chunk = RING_SIZE - s->head;   /* contiguous part */
        size_t n;
        if (chunk > s->count)
            chunk = s->count;
        n = port_write(s, s->ring + s->head, chunk);
        s->head = (s->head + n) % RING_SIZE;
        s->count -= n;
        s->bytesSent += n;
        if (n < chunk)
            break;                            /* the port is full */
    }
    if (s->count == 0)
        s->head = 0;                          /* keep the next frames contiguous */
}

/* Makes room for 'len' bytes if possible; returns 1 when they fit. */
static int ring_reserve(SerialOut *s, size_t len)
{
    if (len > ring_free(s))
        ring_flush(s);
    return len <= ring_free(s);
}

/* ------------------------------------------------------------------------ */
/* Public functions                                                          */
/* ------------------------------------------------------------------------ */

SerialOut *serial_open(const char *device, int baud, char *err, size_t errLen)
{
    SerialOut *s;

    if (!device || !device[0] || baud <= 0) {
        set_err(err, errLen, "invalid serial device or baud rate");
        return NULL;
    }
    s = (SerialOut *)calloc(1, sizeof *s);
    if (s)
        s->ring = (uint8_t *)malloc(RING_SIZE);
    if (!s || !s->ring) {
        if (s)
            free(s);
        set_err(err, errLen, "out of memory");
        return NULL;
    }
    if (port_open(s, device, baud, err, errLen) != 0) {
        free(s->ring);
        free(s);
        return NULL;
    }
    return s;
}

int serial_send_detections(SerialOut *s, int pulse, uint32_t pulseNumber, double az_deg,
                           const Detection *dets, int n)
{
    int nFrames, ok = 1;
    size_t total;

    if (!s || n < 0 || (n > 0 && !dets))
        return 0;

    /* All frames of one CPI update are queued, or none of them. */
    nFrames = n / MAX_DETS_PER_FRAME + (n % MAX_DETS_PER_FRAME != 0);   /* no overflow */
    total = (size_t)nFrames * (DET_HEADER_LEN + CRC_LEN) + (size_t)n * DET_RECORD_LEN;
    if (!ring_reserve(s, total)) {
        s->framesDropped += (uint64_t)nFrames;
        ok = 0;
    } else {
        for (int first = 0; first < n; first += MAX_DETS_PER_FRAME) {
            int count = n - first < MAX_DETS_PER_FRAME ? n - first : MAX_DETS_PER_FRAME;
            size_t len = serial_encode_detections(s->frame, sizeof s->frame, pulse, pulseNumber,
                                                  az_deg, dets + first, count);
            ring_put(s, s->frame, len);
            s->framesSent++;
        }
    }
    ring_flush(s);   /* also when n == 0: write what is left from earlier calls */
    return ok;
}

int serial_send_plot(SerialOut *s, const Detection *plot)
{
    int ok = 0;

    if (!s || !plot)
        return 0;
    if (ring_reserve(s, PLOT_FRAME_LEN)) {
        size_t len = serial_encode_plot(s->frame, sizeof s->frame, plot);
        ring_put(s, s->frame, len);
        s->framesSent++;
        ok = 1;
    } else {
        s->framesDropped++;
    }
    ring_flush(s);
    return ok;
}

void serial_stats(const SerialOut *s, uint64_t *framesSent, uint64_t *framesDropped, uint64_t *bytesSent)
{
    if (framesSent)
        *framesSent = s ? s->framesSent : 0;
    if (framesDropped)
        *framesDropped = s ? s->framesDropped : 0;
    if (bytesSent)
        *bytesSent = s ? s->bytesSent : 0;
}

void serial_close(SerialOut *s)
{
    double deadline;

    if (!s)
        return;
    /* Write the ring and let the driver transmit it, for at most CLOSE_TIMEOUT_S. */
    deadline = mono_seconds() + CLOSE_TIMEOUT_S;
    for (;;) {
        ring_flush(s);
        if (s->count == 0 && port_pending(s) == 0)
            break;
        if (mono_seconds() >= deadline) {
            port_discard(s);   /* so that closing the port does not wait */
            break;
        }
        port_wait(s);
    }
    port_close(s);
    free(s->ring);
    free(s);
}

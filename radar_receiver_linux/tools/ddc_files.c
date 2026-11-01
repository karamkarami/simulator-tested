/*
 * ddc_files.c - reading the DDC files of the simulator and splitting a PRI
 * into the datagrams of the live stream (see ddc_files.h).
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#else
#  include <dirent.h>
#endif

#include "ddc_files.h"
#include "mat_file.h"

#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RDDC_MAGIC   0x52444443u
#define HEADER_BYTES 64

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
/* File list                                                                 */
/* ------------------------------------------------------------------------ */

static int is_ddc_name(const char *name)
{
    const size_t n = strlen(name);
    return n > 8 && !strncmp(name, "ddc_", 4) && !strcmp(name + n - 4, ".mat");
}

static int cmp_str(const void *a, const void *b)
{
    return strcmp(*(char *const *)a, *(char *const *)b);
}

static int add_path(char ***paths, int *n, int *cap, const char *folder, const char *name)
{
    if (*n == *cap) {
        int newCap = *cap ? 2 * *cap : 256;
        char **p = realloc(*paths, (size_t)newCap * sizeof *p);
        if (!p) return -1;
        *paths = p;
        *cap = newCap;
    }
    const size_t len = strlen(folder) + strlen(name) + 2;
    char *s = malloc(len);
    if (!s) return -1;
    snprintf(s, len, "%s/%s", folder, name);
    (*paths)[(*n)++] = s;
    return 0;
}

int ddc_list_files(const char *folder, char ***paths, char *err, size_t errLen)
{
    int n = 0, cap = 0, failed = 0;
    *paths = NULL;
#ifdef _WIN32
    char pattern[4096];
    WIN32_FIND_DATAA fd;
    snprintf(pattern, sizeof pattern, "%s\\ddc_*.mat", folder);
    HANDLE h = FindFirstFileA(pattern, &fd);
    if (h == INVALID_HANDLE_VALUE) {
        set_err(err, errLen, "no ddc_*.mat files in %s", folder);
        return -1;
    }
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) && is_ddc_name(fd.cFileName))
            failed |= add_path(paths, &n, &cap, folder, fd.cFileName);
    } while (!failed && FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(folder);
    struct dirent *e;
    if (!d) {
        set_err(err, errLen, "cannot open the folder %s", folder);
        return -1;
    }
    while (!failed && (e = readdir(d)) != NULL)
        if (is_ddc_name(e->d_name))
            failed |= add_path(paths, &n, &cap, folder, e->d_name);
    closedir(d);
#endif
    if (failed) {
        ddc_free_list(*paths, n);
        *paths = NULL;
        set_err(err, errLen, "out of memory");
        return -1;
    }
    if (n == 0) {
        free(*paths);
        *paths = NULL;
        set_err(err, errLen, "no ddc_*.mat files in %s", folder);
        return -1;
    }
    qsort(*paths, (size_t)n, sizeof **paths, cmp_str);    /* zero-padded numbers: name order = pulse order */
    return n;
}

void ddc_free_list(char **paths, int n)
{
    for (int i = 0; i < n; i++)
        free(paths[i]);
    free(paths);
}

/* ------------------------------------------------------------------------ */
/* One file                                                                  */
/* ------------------------------------------------------------------------ */

/* MATLAB's conversion of a double to an integer type: round half away from
 * zero, saturate, NaN -> 0. */
static double to_int(double x, double lo, double hi)
{
    if (x != x) return 0;
    x = round(x);
    return x < lo ? lo : x > hi ? hi : x;
}

static int scalar(const MatVar *hdr, const char *name, double *v, int count, char *err, size_t errLen)
{
    const MatVar *f = mat_field(hdr, name);
    if (!f || !f->re || f->numel < (size_t)count) {
        set_err(err, errLen, "hdr.%s is missing", name);
        return -1;
    }
    for (int i = 0; i < count; i++)
        v[i] = f->re[i];
    return 0;
}

int ddc_load(const char *path, DdcPri *p, char *err, size_t errLen)
{
    MatFile mf;
    char why[256];
    double v[2];
    int r = -1;

    if (mat_load(path, &mf, why, sizeof why)) {
        set_err(err, errLen, "%s: %s", path, why);
        return -1;
    }
    const MatVar *I = mat_get(&mf, "I"), *Q = mat_get(&mf, "Q"), *hdr = mat_get(&mf, "hdr");
    if (!I || !Q || !hdr || hdr->mxClass != MAT_CLASS_STRUCT) {
        set_err(err, errLen, "%s: the variables I, Q and hdr are needed", path);
        goto done;
    }
    if (!I->re || !Q->re || I->rows != Q->rows || I->cols != Q->cols || I->rows > 0xFFFFFFFFu || I->cols > 0xFFFFu) {
        set_err(err, errLen, "%s: I and Q must be numeric matrices of the same size", path);
        goto done;
    }

#define GET(name, n) if (scalar(hdr, name, v, n, why, sizeof why)) { set_err(err, errLen, "%s: %s", path, why); goto done; }
    GET("pulseNumber", 1)   p->pulseNumber  = (uint32_t)to_int(v[0], 0, 4294967295.0);
    GET("scanNumber", 1)    p->scanNumber   = (uint32_t)to_int(v[0], 0, 4294967295.0);
    GET("antennaAz_deg", 1) p->az_deg       = (float)v[0];
    GET("staggerIndex", 1)  p->staggerIndex = (uint16_t)to_int(v[0], 0, 65535.0);
    GET("time_s", 1)        p->time_s       = v[0];
    GET("pri_s", 1)         p->pri_s        = v[0];
    GET("ddcFs_Hz", 1)      p->ddcFs_Hz     = (float)v[0];
    GET("ddcT0_s", 1)       p->ddcT0_s      = (float)v[0];
    GET("shortDdcIdx", 2)   p->shortIdx[0]  = (uint32_t)to_int(v[0], 0, 4294967295.0);
                            p->shortIdx[1]  = (uint32_t)to_int(v[1], 0, 4294967295.0);
    GET("longDdcIdx", 2)    p->longIdx[0]   = (uint32_t)to_int(v[0], 0, 4294967295.0);
                            p->longIdx[1]   = (uint32_t)to_int(v[1], 0, 4294967295.0);
#undef GET

    {
        const size_t n = I->numel;                    /* complex samples */
        if (2 * n > p->capacity) {
            int16_t *q = realloc(p->iq, 2 * n * sizeof *q);
            if (!q) {
                set_err(err, errLen, "out of memory");
                goto done;
            }
            p->iq = q;
            p->capacity = 2 * n;
        }
        for (size_t k = 0; k < n; k++) {              /* column-major: sample after sample, then channel */
            p->iq[2 * k]     = (int16_t)to_int(I->re[k], -32768.0, 32767.0);
            p->iq[2 * k + 1] = (int16_t)to_int(Q->re[k], -32768.0, 32767.0);
        }
        p->nSamples = (int)I->rows;
        p->nCh = (int)I->cols;
    }
    r = 0;
done:
    mat_free(&mf);
    return r;
}

void ddc_free(DdcPri *p)
{
    free(p->iq);
    memset(p, 0, sizeof *p);
}

/* ------------------------------------------------------------------------ */
/* Datagrams                                                                 */
/* ------------------------------------------------------------------------ */

static unsigned char *put_u16(unsigned char *p, uint16_t v) { p[0] = v & 255; p[1] = v >> 8; return p + 2; }
static unsigned char *put_u32(unsigned char *p, uint32_t v)
{
    p[0] = v & 255; p[1] = v >> 8 & 255; p[2] = v >> 16 & 255; p[3] = v >> 24;
    return p + 4;
}
static unsigned char *put_f32(unsigned char *p, float v) { uint32_t u; memcpy(&u, &v, 4); return put_u32(p, u); }
static unsigned char *put_f64(unsigned char *p, double v) { memcpy(p, &v, 8); return p + 8; }   /* little-endian host */

static size_t per_packet(int payloadBytes)
{
    size_t n = 4 * (size_t)(payloadBytes / 4);        /* whole complex samples */
    return n ? n : 4;
}

int ddc_datagram_count(const DdcPri *p, int payloadBytes, size_t *totalBytes)
{
    const size_t bytes = 4 * (size_t)p->nSamples * (size_t)p->nCh;
    const size_t per = per_packet(payloadBytes);
    const int n = (int)((bytes + per - 1) / per);
    if (totalBytes) *totalBytes = bytes + (size_t)n * HEADER_BYTES;
    return n;
}

int ddc_make_datagrams(const DdcPri *p, int payloadBytes, unsigned char *out, uint32_t *len)
{
    const size_t bytes = 4 * (size_t)p->nSamples * (size_t)p->nCh;
    const size_t per = per_packet(payloadBytes);
    const int n = ddc_datagram_count(p, payloadBytes, NULL);
    const unsigned char *samples = (const unsigned char *)p->iq;   /* int16 little-endian */

    for (int i = 0; i < n; i++) {
        const size_t first = (size_t)i * per;
        const size_t payload = first + per <= bytes ? per : bytes - first;
        unsigned char *h = out;
        h = put_u32(h, RDDC_MAGIC);
        h = put_u32(h, p->pulseNumber);
        h = put_u32(h, p->scanNumber);
        h = put_f32(h, p->az_deg);
        h = put_u16(h, (uint16_t)i);
        h = put_u16(h, (uint16_t)n);
        h = put_u16(h, (uint16_t)p->nCh);
        h = put_u16(h, p->staggerIndex);
        h = put_u32(h, (uint32_t)p->nSamples);
        h = put_u32(h, (uint32_t)(first / 4));
        h = put_f64(h, p->time_s);
        h = put_f32(h, p->ddcFs_Hz);
        h = put_f32(h, p->ddcT0_s);
        h = put_u32(h, p->shortIdx[0]);
        h = put_u32(h, p->shortIdx[1]);
        h = put_u32(h, p->longIdx[0]);
        h = put_u32(h, p->longIdx[1]);
        memcpy(h, samples + first, payload);
        len[i] = (uint32_t)(HEADER_BYTES + payload);
        out += len[i];
    }
    return n;
}

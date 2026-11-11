/*
 * mat_write.c - minimal MAT-file writer (see mat_write.h). The elements are
 * written as MATLAB writes them: 128-byte header, then per variable an
 * miMATRIX element (array flags, dimensions, name, data), 8-byte aligned;
 * with -v7 each variable inside an miCOMPRESSED element.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "mat_write.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef RADAR_HAVE_ZLIB
#  include <zlib.h>
#endif

enum { miINT8 = 1, miINT16 = 3, miINT32 = 5, miUINT32 = 6, miDOUBLE = 9, miMATRIX = 14, miCOMPRESSED = 15 };
enum { mxSTRUCT = 2, mxDOUBLE = 6, mxINT16 = 10 };

typedef struct { unsigned char *p; size_t n, cap; } Buf;

struct MatWriter {
    FILE *f;
    int   compressed, failed;
};

static void put(Buf *b, const void *d, size_t n)
{
    if (b->n + n > b->cap) {
        b->cap = 2 * (b->n + n) + 256;
        unsigned char *p = realloc(b->p, b->cap);
        if (!p) abort();
        b->p = p;
    }
    if (n) memcpy(b->p + b->n, d, n);
    b->n += n;
}
static void put32(Buf *b, uint32_t v) { put(b, &v, 4); }              /* little-endian hosts */
static void pad8(Buf *b) { static const unsigned char z[8]; put(b, z, (8 - b->n % 8) % 8); }

static void element(Buf *b, uint32_t type, const void *d, uint32_t n)
{
    if (n <= 4) {                                                      /* small data element */
        unsigned char z[4] = { 0 };
        put32(b, n << 16 | type);
        memcpy(z, d, n);
        put(b, z, 4);
    } else {
        put32(b, type);
        put32(b, n);
        put(b, d, n);
        pad8(b);
    }
}

static void matrix_head(Buf *b, int mxClass, uint32_t rows, uint32_t cols, const char *name)
{
    const uint32_t flags[2] = { (uint32_t)mxClass, 0 };
    const int32_t dims[2] = { (int32_t)rows, (int32_t)cols };
    element(b, miUINT32, flags, 8);
    element(b, miINT32, dims, 8);
    element(b, miINT8, name, (uint32_t)strlen(name));
}

static void matrix_end(Buf *out, Buf *body)
{
    put32(out, miMATRIX);
    put32(out, (uint32_t)body->n);
    put(out, body->p, body->n);
    free(body->p);
}

MatWriter *mat_write_open(const char *path, int compressed)
{
    FILE *f = fopen(path, "wb");
    if (!f) return NULL;
    MatWriter *w = calloc(1, sizeof *w);
    if (!w) { fclose(f); return NULL; }
    w->f = f;
    w->compressed = compressed;
    unsigned char head[128];
    memset(head, ' ', 116);
    const char *text = "MATLAB 5.0 MAT-file, written by radar_generate_ddc";
    memcpy(head, text, strlen(text));
    memset(head + 116, 0, 8);
    head[124] = 0x00; head[125] = 0x01; head[126] = 'I'; head[127] = 'M';
    if (fwrite(head, 1, 128, f) != 128) w->failed = 1;
    return w;
}

/* Writes one variable (a complete miMATRIX element), compressed with -v7. */
static void emit(MatWriter *w, Buf *var)
{
    if (!w->compressed) {
        if (fwrite(var->p, 1, var->n, w->f) != var->n) w->failed = 1;
        free(var->p);
        return;
    }
    Buf out = { 0 };
#ifdef RADAR_HAVE_ZLIB
    uLongf zlen = compressBound((uLong)var->n);
    unsigned char *z = malloc(zlen);
    if (!z || compress2(z, &zlen, var->p, (uLong)var->n, 1) != Z_OK) abort();
    put32(&out, miCOMPRESSED);
    put32(&out, (uint32_t)zlen);
    put(&out, z, zlen);
    free(z);
#else
    {                                                                  /* stored blocks */
        const size_t n = var->n, nBlk = (n + 65534) / 65535;
        const uint32_t zlen = (uint32_t)(2 + n + 5 * (nBlk ? nBlk : 1) + 4);
        uint32_t a = 1, b = 0;
        put32(&out, miCOMPRESSED);
        put32(&out, zlen);
        const unsigned char zh[2] = { 0x78, 0x01 };
        put(&out, zh, 2);
        size_t off = 0;
        do {
            const size_t len = n - off < 65535 ? n - off : 65535;
            const unsigned char bh[5] = { (unsigned char)(off + len >= n), (unsigned char)(len & 0xFF),
                                          (unsigned char)(len >> 8), (unsigned char)(~len & 0xFF),
                                          (unsigned char)((~len >> 8) & 0xFF) };
            put(&out, bh, 5);
            put(&out, var->p + off, len);
            off += len;
        } while (off < n);
        for (size_t i = 0; i < n; i++) { a = (a + var->p[i]) % 65521; b = (b + a) % 65521; }
        const uint32_t adler = b << 16 | a;
        const unsigned char ad[4] = { (unsigned char)(adler >> 24), (unsigned char)(adler >> 16),
                                      (unsigned char)(adler >> 8), (unsigned char)adler };
        put(&out, ad, 4);
    }
#endif
    if (fwrite(out.p, 1, out.n, w->f) != out.n) w->failed = 1;
    free(out.p);
    free(var->p);
}

void mat_write_int16(MatWriter *w, const char *name, const int16_t *data, size_t rows, size_t cols)
{
    Buf m = { 0 }, v = { 0 };
    matrix_head(&m, mxINT16, (uint32_t)rows, (uint32_t)cols, name);
    element(&m, miINT16, data, (uint32_t)(sizeof(int16_t) * rows * cols));
    matrix_end(&v, &m);
    emit(w, &v);
}

void mat_write_struct(MatWriter *w, const char *name, int nFields, const char *const *fieldNames,
                      const double *const *values, const int *n)
{
    Buf m = { 0 }, v = { 0 };
    int32_t len = 1;
    for (int k = 0; k < nFields; k++) {
        const int32_t l = (int32_t)strlen(fieldNames[k]) + 1;
        if (l > len) len = l;
    }
    char *names = calloc((size_t)nFields, (size_t)len);
    if (!names) abort();
    for (int k = 0; k < nFields; k++) memcpy(names + (size_t)k * len, fieldNames[k], strlen(fieldNames[k]));
    matrix_head(&m, mxSTRUCT, 1, 1, name);
    element(&m, miINT32, &len, 4);
    element(&m, miINT8, names, (uint32_t)(len * nFields));
    free(names);
    for (int k = 0; k < nFields; k++) {
        Buf f = { 0 };
        matrix_head(&f, mxDOUBLE, 1, (uint32_t)n[k], "");
        element(&f, miDOUBLE, values[k], (uint32_t)(sizeof(double) * n[k]));
        matrix_end(&m, &f);
    }
    matrix_end(&v, &m);
    emit(w, &v);
}

int mat_write_close(MatWriter *w)
{
    if (!w) return -1;
    const int bad = w->failed || fclose(w->f) != 0;
    free(w);
    return bad ? -1 : 0;
}

int mat_write_has_zlib(void)
{
#ifdef RADAR_HAVE_ZLIB
    return 1;
#else
    return 0;
#endif
}

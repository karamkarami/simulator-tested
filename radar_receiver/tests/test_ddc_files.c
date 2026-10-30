/*
 * test_ddc_files.c - test of the MAT-file reader and of the DDC files
 * (tools/mat_file.c, tools/ddc_files.c) used by radar_stream_send.
 *
 * MAT-files are written byte by byte as MATLAB writes them, including what
 * Octave does not do: numeric data stored in a smaller type than its class
 * (int16 data as int8, doubles as uint8/uint16/uint32), small data elements,
 * a skipped char field and, when built with zlib, compressed variables
 * (-v7). Checks:
 *   - folder listing (only ddc_*.mat, sorted by name),
 *   - header fields and I/Q values of a plain and of a compressed file,
 *   - the datagrams of a PRI, replayed through udp_input.c, give back the
 *     same PRI (the receiver accepts them),
 *   - errors: not a MAT-file, -v7.3, big-endian, truncated, missing variables.
 *
 * Prints PASS or FAIL and returns 0 on success.
 */
#if defined(_WIN32)
#  ifndef _CRT_SECURE_NO_WARNINGS
#    define _CRT_SECURE_NO_WARNINGS
#  endif
#  include <direct.h>
#  define make_dir(p) _mkdir(p)
#else
#  ifndef _POSIX_C_SOURCE
#    define _POSIX_C_SOURCE 200809L
#  endif
#  include <sys/stat.h>
#  define make_dir(p) mkdir(p, 0755)
#endif

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef RADAR_HAVE_ZLIB
#  include <zlib.h>
#endif

#include "ddc_files.h"
#include "mat_file.h"
#include "radar_rx.h"

#define DIR_NAME    "test_ddc_files_dir"
#define RECORD_FILE "test_ddc_files.rec"
#define NS          37          /* DDC samples per channel */
#define NCH         3           /* channels                */

static int g_failures = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        if (!(cond)) {                                                       \
            printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);           \
            g_failures++;                                                    \
        }                                                                    \
    } while (0)

/* ------------------------------------------------------------------------ */
/* MAT-file writer (as MATLAB writes them)                                   */
/* ------------------------------------------------------------------------ */

enum { miINT8 = 1, miUINT8 = 2, miINT16 = 3, miUINT16 = 4, miINT32 = 5, miUINT32 = 6, miDOUBLE = 9,
       miMATRIX = 14, miCOMPRESSED = 15 };
enum { mxCHAR = 4, mxSTRUCT = 2, mxDOUBLE = 6, mxINT16 = 10 };

typedef struct { unsigned char *p; size_t n, cap; } Buf;

static void put(Buf *b, const void *d, size_t n)
{
    if (b->n + n > b->cap) {
        b->cap = 2 * (b->n + n) + 256;
        b->p = realloc(b->p, b->cap);
    }
    memcpy(b->p + b->n, d, n);
    b->n += n;
}
static void put32(Buf *b, uint32_t v) { put(b, &v, 4); }   /* little-endian host */
static void pad8(Buf *b) { static const unsigned char z[8]; put(b, z, (8 - b->n % 8) % 8); }

/* A data element; small format when it fits into 4 bytes. */
static void element(Buf *b, uint32_t type, const void *d, uint32_t n)
{
    if (n <= 4) {
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

/* miMATRIX header: flags, dimensions, name. */
static void matrix_head(Buf *b, int mxClass, uint32_t rows, uint32_t cols, const char *name)
{
    uint32_t flags[2] = { (uint32_t)mxClass, 0 };
    int32_t dims[2] = { (int32_t)rows, (int32_t)cols };
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

/* A numeric matrix whose data is stored as 'type' (raw bytes given). */
static void numeric(Buf *out, int mxClass, uint32_t rows, uint32_t cols, const char *name,
                    uint32_t type, const void *data, uint32_t bytes)
{
    Buf m = { 0 };
    matrix_head(&m, mxClass, rows, cols, name);
    element(&m, type, data, bytes);
    matrix_end(out, &m);
}

static void scalar_u8(Buf *out, const char *name, uint8_t v)   { numeric(out, mxDOUBLE, 1, 1, name, miUINT8, &v, 1); }
static void scalar_u16(Buf *out, const char *name, uint16_t v) { numeric(out, mxDOUBLE, 1, 1, name, miUINT16, &v, 2); }
static void scalar_u32(Buf *out, const char *name, uint32_t v) { numeric(out, mxDOUBLE, 1, 1, name, miUINT32, &v, 4); }
static void scalar_f64(Buf *out, const char *name, double v)   { numeric(out, mxDOUBLE, 1, 1, name, miDOUBLE, &v, 8); }

/* The header struct of a DDC file. */
static void hdr_struct(Buf *out, uint32_t pulse, double time_s)
{
    static const char *names[] = { "pulseNumber", "scanNumber", "time_s", "antennaAz_deg", "pri_s",
                                   "staggerIndex", "comment", "ddcFs_Hz", "ddcT0_s", "shortDdcIdx", "longDdcIdx" };
    const int nf = (int)(sizeof names / sizeof names[0]);
    const int32_t len = 32;
    char fn[32 * 11];
    Buf m = { 0 };
    memset(fn, 0, sizeof fn);
    for (int i = 0; i < nf; i++) strcpy(fn + 32 * i, names[i]);
    matrix_head(&m, mxSTRUCT, 1, 1, "hdr");
    element(&m, miINT32, &len, 4);
    element(&m, miINT8, fn, (uint32_t)(32 * nf));
    scalar_u32(&m, "", pulse);                         /* fields have empty names */
    scalar_u16(&m, "", 2);
    scalar_f64(&m, "", time_s);
    scalar_f64(&m, "", 123.25);
    scalar_f64(&m, "", 1e-3);
    scalar_u8(&m, "", 3);
    {
        Buf c = { 0 };                                 /* char field: skipped */
        const uint16_t text[3] = { 'a', 'b', 'c' };
        matrix_head(&c, mxCHAR, 1, 3, "");
        element(&c, miUINT16, text, 6);
        matrix_end(&m, &c);
    }
    scalar_u32(&m, "", 2000000);
    scalar_f64(&m, "", 28e-6);
    {
        const uint8_t idx[2] = { 1, 12 };              /* [1 12] as uint8 */
        numeric(&m, mxDOUBLE, 1, 2, "", miUINT8, idx, 2);
    }
    {
        const uint16_t idx[2] = { 20, NS };
        numeric(&m, mxDOUBLE, 1, 2, "", miUINT16, idx, 4);
    }
    matrix_end(out, &m);
}

static int16_t I_value(int s, int ch) { return (int16_t)((s * 977 + ch * 3001) % 60001 - 30000); }
static int8_t  Q_value(int s, int ch) { return (int8_t)((s * 7 + ch * 13) % 255 - 127); }

/* Writes a DDC file; compressed: every variable in a zlib element (-v7). */
static void write_ddc_file(const char *path, uint32_t pulse, double time_s, int compressed)
{
    Buf vars[3] = { { 0 } }, file = { 0 };
    int16_t I[NS * NCH];
    int8_t Q[NS * NCH];
    for (int ch = 0; ch < NCH; ch++)
        for (int s = 0; s < NS; s++) {
            I[ch * NS + s] = I_value(s, ch);
            Q[ch * NS + s] = Q_value(s, ch);
        }
    numeric(&vars[0], mxINT16, NS, NCH, "I", miINT16, I, sizeof I);
    numeric(&vars[1], mxINT16, NS, NCH, "Q", miINT8, Q, sizeof Q);   /* int16 data stored as int8 */
    hdr_struct(&vars[2], pulse, time_s);

    unsigned char head[128];
    memset(head, ' ', 116);
    memcpy(head, "MATLAB 5.0 MAT-file, written by test_ddc_files", 46);
    memset(head + 116, 0, 8);
    head[124] = 0x00; head[125] = 0x01; head[126] = 'I'; head[127] = 'M';
    put(&file, head, 128);
    for (int v = 0; v < 3; v++) {
#ifdef RADAR_HAVE_ZLIB
        if (compressed) {
            uLongf zlen = compressBound((uLong)vars[v].n);
            unsigned char *z = malloc(zlen);
            compress(z, &zlen, vars[v].p, (uLong)vars[v].n);
            put32(&file, miCOMPRESSED);
            put32(&file, (uint32_t)zlen);
            put(&file, z, zlen);                       /* not padded */
            free(z);
            free(vars[v].p);
            continue;
        }
#else
        if (compressed) {                              /* the reader must refuse it: content does not matter */
            put32(&file, miCOMPRESSED);
            put32(&file, (uint32_t)vars[v].n);
            put(&file, vars[v].p, vars[v].n);
            free(vars[v].p);
            continue;
        }
#endif
        put(&file, vars[v].p, vars[v].n);
        free(vars[v].p);
    }
    FILE *f = fopen(path, "wb");
    fwrite(file.p, 1, file.n, f);
    fclose(f);
    free(file.p);
}

static void write_bytes(const char *path, const void *d, size_t n)
{
    FILE *f = fopen(path, "wb");
    fwrite(d, 1, n, f);
    fclose(f);
}

/* ------------------------------------------------------------------------ */
/* Tests                                                                     */
/* ------------------------------------------------------------------------ */

static void check_pri(const DdcPri *p, uint32_t pulse, double time_s)
{
    CHECK(p->pulseNumber == pulse);
    CHECK(p->scanNumber == 2);
    CHECK(p->time_s == time_s);
    CHECK(p->az_deg == 123.25f);
    CHECK(p->pri_s == 1e-3);
    CHECK(p->staggerIndex == 3);
    CHECK(p->ddcFs_Hz == 2e6f);
    CHECK(p->ddcT0_s == (float)28e-6);
    CHECK(p->shortIdx[0] == 1 && p->shortIdx[1] == 12);
    CHECK(p->longIdx[0] == 20 && p->longIdx[1] == NS);
    CHECK(p->nSamples == NS && p->nCh == NCH);
    int bad = 0;
    for (int ch = 0; ch < NCH; ch++)
        for (int s = 0; s < NS; s++) {
            bad += p->iq[2 * (ch * NS + s)] != I_value(s, ch);
            bad += p->iq[2 * (ch * NS + s) + 1] != Q_value(s, ch);
        }
    CHECK(bad == 0);
}

static void test_files(void)
{
    char err[512];
    char **paths;
    DdcPri p;
    memset(&p, 0, sizeof p);

    make_dir(DIR_NAME);
    write_ddc_file(DIR_NAME "/ddc_0000002.mat", 2, 1.001, 1);
    write_ddc_file(DIR_NAME "/ddc_0000001.mat", 1, 1.000, 0);
    write_bytes(DIR_NAME "/other.mat", "x", 1);                  /* not listed */
    write_bytes(DIR_NAME "/ddc_notes.txt", "x", 1);              /* not listed */

    const int n = ddc_list_files(DIR_NAME, &paths, err, sizeof err);
    CHECK(n == 2);
    if (n != 2) return;
    CHECK(strstr(paths[0], "ddc_0000001.mat") != NULL);
    CHECK(strstr(paths[1], "ddc_0000002.mat") != NULL);

    CHECK(ddc_load(paths[0], &p, err, sizeof err) == 0);
    check_pri(&p, 1, 1.000);
    if (mat_has_zlib()) {
        CHECK(ddc_load(paths[1], &p, err, sizeof err) == 0);
        check_pri(&p, 2, 1.001);
    } else {
        CHECK(ddc_load(paths[1], &p, err, sizeof err) != 0 && strstr(err, "-v6") != NULL);
        printf("(built without zlib: compressed file rejected as expected)\n");
    }

    /* Datagrams of PRI 1 -> record file -> udp_input: the same PRI */
    CHECK(ddc_load(paths[0], &p, err, sizeof err) == 0);
    for (int payload = 100; payload <= 1000; payload += 900) {
        size_t bytes;
        const int nPkt = ddc_datagram_count(&p, payload, &bytes);
        unsigned char *d = malloc(bytes);
        uint32_t *len = malloc((size_t)nPkt * sizeof *len);
        CHECK(ddc_make_datagrams(&p, payload, d, len) == nPkt);
        CHECK(nPkt == (4 * NS * NCH + 4 * (payload / 4) - 1) / (4 * (payload / 4)));
        FILE *f = fopen(RECORD_FILE, "wb");
        size_t off = 0;
        for (int i = nPkt - 1; i >= 0; i--) {                     /* reverse order */
            off = 0;
            for (int k = 0; k < i; k++) off += len[k];
            unsigned char l4[4] = { (unsigned char)(len[i] & 255), (unsigned char)(len[i] >> 8), 0, 0 };
            fwrite(l4, 1, 4, f);
            fwrite(d + off, 1, len[i], f);
        }
        fclose(f);
        UdpInput *in = udp_input_open(0, NULL, RECORD_FILE, err, sizeof err);
        CHECK(in != NULL);
        if (in) {
            PriData pri;
            memset(&pri, 0, sizeof pri);
            CHECK(udp_input_next(in, &pri, 0) == 1);
            CHECK(pri.pulseNumber == 1 && pri.scanNumber == 2 && pri.staggerIndex == 3);
            CHECK(pri.time_s == 1.000 && pri.az_deg == 123.25f && pri.ddcFs_Hz == 2e6f);
            CHECK(pri.shortIdx[1] == 12 && pri.longIdx[0] == 20);
            CHECK(pri.nSamples == NS && pri.nCh == NCH);
            CHECK(pri.iq && memcmp(pri.iq, p.iq, 4 * NS * NCH) == 0);
            CHECK(udp_input_next(in, &pri, 0) == -1);
            pri_data_free(&pri);
            udp_input_close(in);
        }
        free(d);
        free(len);
    }
    ddc_free(&p);
    ddc_free_list(paths, n);
    remove(DIR_NAME "/ddc_0000001.mat");
    remove(DIR_NAME "/ddc_0000002.mat");
    remove(DIR_NAME "/other.mat");
    remove(DIR_NAME "/ddc_notes.txt");
    remove(RECORD_FILE);
}

static void test_errors(void)
{
    char err[512];
    MatFile mf;
    unsigned char head[128];
    DdcPri p;
    memset(&p, 0, sizeof p);

    CHECK(mat_parse((const unsigned char *)"hello", 5, &mf, err, sizeof err) != 0);
    memset(head, ' ', sizeof head);
    memcpy(head, "MATLAB 7.3 MAT-file", 19);
    CHECK(mat_parse(head, sizeof head, &mf, err, sizeof err) != 0 && strstr(err, "7.3") != NULL);
    memcpy(head, "MATLAB 5.0 MAT-file", 19);
    head[124] = 0x01; head[125] = 0x00; head[126] = 'M'; head[127] = 'I';
    CHECK(mat_parse(head, sizeof head, &mf, err, sizeof err) != 0 && strstr(err, "big-endian") != NULL);

    /* a truncated variable and a file without hdr */
    write_ddc_file("test_ddc_full.mat", 1, 1.0, 0);
    FILE *f = fopen("test_ddc_full.mat", "rb");
    unsigned char buf[8192];
    const size_t n = fread(buf, 1, sizeof buf, f);
    fclose(f);
    CHECK(mat_parse(buf, n, &mf, err, sizeof err) == 0 && mf.nVars == 3);
    mat_free(&mf);
    CHECK(mat_parse(buf, n - 40, &mf, err, sizeof err) != 0);
    write_bytes("test_ddc_cut.mat", buf, 128 + 300);          /* I and part of Q */
    CHECK(ddc_load("test_ddc_cut.mat", &p, err, sizeof err) != 0);
    CHECK(ddc_load("does_not_exist.mat", &p, err, sizeof err) != 0);
    remove("test_ddc_full.mat");
    remove("test_ddc_cut.mat");
    ddc_free(&p);
}

int main(void)
{
    test_files();
    test_errors();
    printf("%s\n", g_failures ? "FAIL" : "PASS");
    return g_failures ? 1 : 0;
}

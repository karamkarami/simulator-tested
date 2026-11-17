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
 *   - errors: not a MAT-file, -v7.3, big-endian, truncated, missing variables,
 *   - (--log CONFIG) the logging of the receiver (logger.c): a log of every block over a
 *     CPI, the files counted and read back (the DDC files as the replay reads them).
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

#include <math.h>
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
        if (compressed) {                              /* a zlib stream of stored blocks (no compressor) */
            const size_t n = vars[v].n, nBlk = (n + 65534) / 65535;
            const uint32_t zlen = (uint32_t)(2 + n + 5 * (nBlk ? nBlk : 1) + 4);
            uint32_t a = 1, b = 0;
            put32(&file, miCOMPRESSED);
            put32(&file, zlen);
            const unsigned char zh[2] = { 0x78, 0x01 };
            put(&file, zh, 2);
            size_t off = 0;
            do {
                const size_t len = n - off < 65535 ? n - off : 65535;
                const unsigned char bh[5] = { (unsigned char)(off + len >= n), (unsigned char)(len & 0xFF),
                                              (unsigned char)(len >> 8), (unsigned char)(~len & 0xFF),
                                              (unsigned char)((~len >> 8) & 0xFF) };
                put(&file, bh, 5);
                put(&file, vars[v].p + off, len);
                off += len;
            } while (off < n);
            for (size_t i = 0; i < n; i++) { a = (a + vars[v].p[i]) % 65521; b = (b + a) % 65521; }
            const uint32_t adler = b << 16 | a;
            const unsigned char ad[4] = { (unsigned char)(adler >> 24), (unsigned char)(adler >> 16),
                                          (unsigned char)(adler >> 8), (unsigned char)adler };
            put(&file, ad, 4);
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
    CHECK(ddc_load(paths[1], &p, err, sizeof err) == 0);       /* compressed: zlib or the built-in inflate */
    check_pri(&p, 2, 1.001);
#ifdef RADAR_HAVE_ZLIB
    {   /* the built-in inflate (used without zlib) against zlib: dynamic and fixed Huffman blocks */
        const size_t n = 200000;
        unsigned char *src = malloc(n);
        uint32_t x = 12345;
        for (size_t i = 0; i < n; i++) {
            x = x * 1103515245u + 12345u;
            src[i] = (unsigned char)(i % 1000 < 600 ? (x >> 24) % 7 : (x >> 16));   /* repeats and noise */
        }
        for (int level = 1; level <= 9; level += 8) {
            uLongf zlen = compressBound((uLong)n);
            unsigned char *z = malloc(zlen);
            compress2(z, &zlen, src, (uLong)n, level);
            size_t outLen = 0;
            unsigned char *out = mat_inflate_builtin(z, zlen, &outLen);
            CHECK(out && outLen == n && memcmp(out, src, n) == 0);
            free(out);
            free(z);
        }
        free(src);
    }
#endif

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

/* The files of radar_generate_ddc (ctest generate_ddc): every file of the
 * folder loads, consecutive pulse numbers, the layout of the configuration
 * (channels, samples), the receive intervals in the header, signal present. */
static int check_generated(const char *folder, int nCh)
{
    char **paths = NULL, err[256];
    const int n = ddc_list_files(folder, &paths, err, sizeof err);
    CHECK(n > 0);
    DdcPri p;
    memset(&p, 0, sizeof p);
    uint32_t prev = 0;
    for (int i = 0; i < n; i++) {
        CHECK(ddc_load(paths[i], &p, err, sizeof err) == 0);
        CHECK(p.nCh == nCh && p.nSamples > 0);
        CHECK(i == 0 || p.pulseNumber == prev + 1);
        prev = p.pulseNumber;
        CHECK(p.shortIdx[0] >= 1 && p.shortIdx[1] > p.shortIdx[0] && p.longIdx[0] > p.shortIdx[1] &&
              p.longIdx[1] <= (uint32_t)p.nSamples);
        double e = 0.0;
        for (size_t k = 0; k < 2 * (size_t)p.nSamples * p.nCh; k++) e += (double)p.iq[k] * p.iq[k];
        CHECK(e > 0.0);
    }
    printf("generated DDC files: %d in %s, %d channels, %d samples\n", n, folder, p.nCh, p.nSamples);
    ddc_free(&p);
    ddc_free_list(paths, n > 0 ? n : 0);
    return n;
}

/* ------------------------------------------------------------------------ */
/* Logging of the receiver (logger.c): the files of every block, read back   */
/* ------------------------------------------------------------------------ */

typedef struct { Logger *log; int nDet, nPlot; } LogSink;
static void log_block(void *u, const BlockOutput *o) { logger_block(((LogSink *)u)->log, o); }
static void log_dets(void *u, int pulse, const Detection *d, int n) { (void)pulse; logger_detections(((LogSink *)u)->log, d, n); ((LogSink *)u)->nDet += n; }
static void log_plot(void *u, const Detection *p) { logger_plot(((LogSink *)u)->log, p); ((LogSink *)u)->nPlot++; }

static int count_files(const char *folder, const char *prefix)
{
    char **paths;
    char err[256];
    (void)prefix;
    const int n = ddc_list_files(folder, &paths, err, sizeof err);   /* (lists ddc_*.mat only) */
    if (n > 0) ddc_free_list(paths, n);
    return n;
}

static void test_log(const char *configPath)
{
    RadarConfig cfg;
    char err[512];
    if (radar_config_load(configPath, &cfg, err, sizeof err) != 0) { printf("configuration: %s\n", err); g_failures++; return; }
    Processor *proc = processor_create(&cfg, 1);
    CHECK(proc != NULL);
    if (!proc) return;
    make_dir("test_log_dir");
    const int mask = LOG_DDC | LOG_DBF | LOG_MF | LOG_DOPPLER | LOG_DETECTION;
    LogSink sink = { NULL, 0, 0 };
    sink.log = logger_start("test_log_dir", "run", mask, &cfg, NULL, err, sizeof err);
    CHECK(sink.log != NULL);
    if (!sink.log) { printf("logger_start: %s\n", err); processor_destroy(proc); radar_config_free(&cfg); return; }
    processor_set_log(proc, mask);
    ProcessorCallbacks cb = { log_dets, log_plot, NULL, &sink, NULL, NULL, NULL, log_block };
    PriData pri;
    memset(&pri, 0, sizeof pri);
    pri.nCh = cfg.nCh;
    pri.nSamples = cfg.nDdc;
    pri.ddcFs_Hz = (float)cfg.fsOut_Hz;
    for (int i = 0; i < N_PULSES; i++) {
        uint32_t *idx = i == PULSE_SHORT ? pri.shortIdx : pri.longIdx;
        if (cfg.pulse[i].enable) { idx[0] = (uint32_t)cfg.pulse[i].rxFirst + 1; idx[1] = (uint32_t)(cfg.pulse[i].rxFirst + cfg.pulse[i].rxCount); }
    }
    pri.iq = malloc(sizeof(int16_t) * 2 * (size_t)cfg.nCh * cfg.nDdc);
    const int nPri = cfg.cpiLength + 2;
    unsigned seed = 7;
    for (int k = 0; k < nPri; k++) {
        pri.pulseNumber = 500 + k;
        pri.scanNumber = 1;
        pri.az_deg = (float)(0.01 * k);
        pri.time_s = 2.0 + k * 1e-3;
        pri.staggerIndex = (uint16_t)(cfg.staggered ? k % cfg.nStagger + 1 : 1);
        for (size_t i = 0; i < 2 * (size_t)cfg.nCh * cfg.nDdc; i++) {
            seed = seed * 1103515245u + 12345u;
            pri.iq[i] = (int16_t)((int)(seed >> 16 & 0xFF) - 128);
        }
        logger_begin_pri(sink.log, &pri);
        processor_process(proc, &pri, &cb);
        logger_end_pri(sink.log);
    }
    processor_flush(proc, &cb);
    logger_note(sink.log, "CFAR threshold 12.00 dB -> 13.00 dB");
    LoggerStats st;
    logger_stop(sink.log, &st);
    const int updates = nPri - cfg.cpiLength + 1;
    const uint64_t expFiles = 4u * (uint64_t)nPri + (uint64_t)updates;   /* ddc, dbf, mf, detection + doppler */
    printf("log: %llu PRIs, %llu files, %.2f MB, %llu dropped, %llu errors (expected %llu files)\n",
           (unsigned long long)st.pris, (unsigned long long)st.files, st.bytes / 1e6, (unsigned long long)st.dropped,
           (unsigned long long)st.errors, (unsigned long long)expFiles);
    CHECK(st.pris == (uint64_t)nPri && st.files == expFiles && st.dropped == 0 && st.errors == 0);
    /* the DDC files: as received, readable by the replay (ddc_files.c) */
    CHECK(count_files("test_log_dir/run/ddc", "ddc") == nPri);
    DdcPri p;
    memset(&p, 0, sizeof p);
    char path[256];
    snprintf(path, sizeof path, "test_log_dir/run/ddc/ddc_%07u.mat", (unsigned)(500 + nPri - 1));
    CHECK(ddc_load(path, &p, err, sizeof err) == 0);
    CHECK(p.nCh == cfg.nCh && p.nSamples == cfg.nDdc && p.pulseNumber == (uint32_t)(500 + nPri - 1));
    CHECK(p.nSamples > 0 && memcmp(p.iq, pri.iq, sizeof(int16_t) * 2 * (size_t)cfg.nCh * cfg.nDdc) == 0);
    CHECK(fabs(p.pri_s - 1e-3) < 1e-9 && fabs(p.time_s - (2.0 + (nPri - 1) * 1e-3)) < 1e-9);
    ddc_free(&p);
    /* the other blocks: sizes of the variables */
    const int tp = cfg.pulse[PULSE_LONG].enable ? PULSE_LONG : PULSE_SHORT;
    const char *pn = tp == PULSE_LONG ? "long" : "short";
    const PulseConfig *pc = &cfg.pulse[tp];
    MatFile m;
    snprintf(path, sizeof path, "test_log_dir/run/dbf/dbf_%07u.mat", 500u);
    CHECK(mat_load(path, &m, err, sizeof err) == 0);
    const MatVar *v = mat_get(&m, pn);
    CHECK(v && v->im && (int)v->rows == cfg.nBeams && (int)v->cols == pc->ddcCount && mat_get(&m, "hdr"));
    mat_free(&m);
    snprintf(path, sizeof path, "test_log_dir/run/mf/mf_%07u.mat", 500u);
    CHECK(mat_load(path, &m, err, sizeof err) == 0);
    v = mat_get(&m, pn);
    CHECK(v && v->im && (int)v->rows == cfg.nBeams && (int)v->cols == pc->nRange);
    v = mat_get(&m, tp == PULSE_LONG ? "longRange_m" : "shortRange_m");
    CHECK(v && (int)v->numel == pc->nRange && fabs(v->re[0] - pc->range_m[0]) < 1e-9);
    mat_free(&m);
    snprintf(path, sizeof path, "test_log_dir/run/doppler/doppler_%07u.mat", (unsigned)(500 + nPri - 1));
    CHECK(mat_load(path, &m, err, sizeof err) == 0);
    v = mat_get(&m, pn);
    CHECK(v && !v->im && (int)v->rows == cfg.nBeams && (int)v->cols == pc->nRange * cfg.nfft);
    v = mat_get(&m, "f_Hz");
    CHECK(v && (int)v->numel == cfg.nfft);
    mat_free(&m);
    snprintf(path, sizeof path, "test_log_dir/run/detection/detection_%07u.mat", (unsigned)(500 + nPri - 1));
    CHECK(mat_load(path, &m, err, sizeof err) == 0);
    v = mat_get(&m, "dets");
    CHECK(v && v->cols == 13 && mat_get(&m, "plots") && mat_get(&m, "columns"));
    mat_free(&m);
    /* the settings at the start (struct params) and the change noted during the log */
    CHECK(mat_load("test_log_dir/run/radar_parameters.mat", &m, err, sizeof err) == 0);
    CHECK(mat_get(&m, "params") != NULL);
    mat_free(&m);
    {
        FILE *f = fopen("test_log_dir/run/log_info.txt", "r");
        char text[4096] = "";
        const size_t n = f ? fread(text, 1, sizeof text - 1, f) : 0;
        if (f) fclose(f);
        text[n] = 0;
        CHECK(strstr(text, "changes during the log:") && strstr(text, "CFAR threshold 12.00 dB -> 13.00 dB"));
    }
    printf("log files read back: DDC (replayable), DBF, MF, Doppler, detections (%d in the CPI updates), "
           "radar_parameters.mat, log_info.txt\n", sink.nDet);
    free(pri.iq);
    processor_destroy(proc);
    radar_config_free(&cfg);
}

int main(int argc, char **argv)
{
    if (argc > 2 && !strcmp(argv[1], "--log")) {         /* test_ddc_files --log CONFIG */
        test_log(argv[2]);
        printf("%s\n", g_failures ? "FAIL" : "PASS");
        return g_failures ? 1 : 0;
    }
    if (argc > 2) {                                      /* test_ddc_files FOLDER NCH */
        check_generated(argv[1], atoi(argv[2]));
        printf("%s\n", g_failures ? "FAIL" : "PASS");
        return g_failures ? 1 : 0;
    }
    test_files();
    test_errors();
    printf("%s\n", g_failures ? "FAIL" : "PASS");
    return g_failures ? 1 : 0;
}

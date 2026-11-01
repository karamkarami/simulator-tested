/*
 * mat_file.c - minimal reader of MATLAB level-5 MAT-files (see mat_file.h).
 *
 * File layout: a 128-byte header (text, version 0x0100, endian indicator
 * 'IM'), then data elements. Every element has an 8-byte tag (uint32 type,
 * uint32 byte count) followed by its data, padded to 8 bytes; a small
 * element (at most 4 bytes) packs the type and count into the first 4 tag
 * bytes and its data into the other 4. A variable is a miMATRIX element
 * with the sub-elements array flags, dimensions, name and data (numeric:
 * real part; struct: field name length, field names, one miMATRIX per
 * field). save -v7 wraps every variable in a zlib-compressed element.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "mat_file.h"

#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef RADAR_HAVE_ZLIB
#  include <zlib.h>
#endif

enum {
    miINT8 = 1, miUINT8 = 2, miINT16 = 3, miUINT16 = 4, miINT32 = 5, miUINT32 = 6,
    miSINGLE = 7, miDOUBLE = 9, miINT64 = 12, miUINT64 = 13, miMATRIX = 14, miCOMPRESSED = 15
};

#define MAX_DEPTH 8

static void set_err(char *err, size_t errLen, const char *fmt, ...)
{
    va_list ap;
    if (!err || errLen == 0)
        return;
    va_start(ap, fmt);
    vsnprintf(err, errLen, fmt, ap);
    va_end(ap);
}

static uint32_t get_le32(const unsigned char *p)
{
    return (uint32_t)p[0] | (uint32_t)p[1] << 8 | (uint32_t)p[2] << 16 | (uint32_t)p[3] << 24;
}

typedef struct {
    uint32_t             type, nbytes;
    const unsigned char *data;
    size_t               total;     /* bytes of the element including tag and padding */
} Element;

/* Reads the element at p (avail bytes left). Returns 0 if it does not fit. */
static int read_element(const unsigned char *p, size_t avail, Element *e)
{
    if (avail < 8)
        return 0;
    const uint32_t w0 = get_le32(p);
    if (w0 >> 16) {                                   /* small data element */
        e->type = w0 & 0xFFFFu;
        e->nbytes = w0 >> 16;
        if (e->nbytes > 4)
            return 0;
        e->data = p + 4;
        e->total = 8;
        return 1;
    }
    e->type = w0;
    e->nbytes = get_le32(p + 4);
    if (e->nbytes > avail - 8)
        return 0;
    e->data = p + 8;
    if (e->type == miCOMPRESSED) {
        e->total = 8 + (size_t)e->nbytes;             /* compressed elements are not padded */
    } else {
        const size_t padded = 8 + (((size_t)e->nbytes + 7) & ~(size_t)7);
        e->total = padded < avail ? padded : avail;
    }
    return 1;
}

static size_t type_size(uint32_t type)
{
    switch (type) {
    case miINT8: case miUINT8: return 1;
    case miINT16: case miUINT16: return 2;
    case miINT32: case miUINT32: case miSINGLE: return 4;
    case miDOUBLE: case miINT64: case miUINT64: return 8;
    default: return 0;
    }
}

/* Converts n values of storage type 'type' to double (little-endian host). */
static void convert(uint32_t type, const unsigned char *d, double *out, size_t n)
{
    size_t i;
    switch (type) {
    case miINT8:   for (i = 0; i < n; i++) out[i] = (int8_t)d[i]; break;
    case miUINT8:  for (i = 0; i < n; i++) out[i] = d[i]; break;
    case miINT16:  for (i = 0; i < n; i++) { int16_t v;  memcpy(&v, d + 2 * i, 2); out[i] = v; } break;
    case miUINT16: for (i = 0; i < n; i++) { uint16_t v; memcpy(&v, d + 2 * i, 2); out[i] = v; } break;
    case miINT32:  for (i = 0; i < n; i++) { int32_t v;  memcpy(&v, d + 4 * i, 4); out[i] = v; } break;
    case miUINT32: for (i = 0; i < n; i++) { uint32_t v; memcpy(&v, d + 4 * i, 4); out[i] = v; } break;
    case miSINGLE: for (i = 0; i < n; i++) { float v;    memcpy(&v, d + 4 * i, 4); out[i] = v; } break;
    case miDOUBLE: memcpy(out, d, 8 * n); break;
    case miINT64:  for (i = 0; i < n; i++) { int64_t v;  memcpy(&v, d + 8 * i, 8); out[i] = (double)v; } break;
    case miUINT64: for (i = 0; i < n; i++) { uint64_t v; memcpy(&v, d + 8 * i, 8); out[i] = (double)v; } break;
    default: break;
    }
}

static void free_var(MatVar *v)
{
    free(v->re);
    for (int i = 0; i < v->nFields; i++)
        free_var(&v->fields[i]);
    free(v->fields);
    memset(v, 0, sizeof *v);
}

/* Parses the content of a miMATRIX element into v. */
static int parse_matrix(const unsigned char *d, size_t n, MatVar *v, int depth, char *err, size_t errLen)
{
    Element flags, dims, name, e;
    size_t pos = 0;

    memset(v, 0, sizeof *v);
    if (n == 0)
        return 0;                                     /* empty matrix (e.g. an empty field) */
    if (depth > MAX_DEPTH) {
        set_err(err, errLen, "variables nested too deeply");
        return -1;
    }
    if (!read_element(d, n, &flags) || flags.nbytes < 8) goto bad;
    pos += flags.total;
    if (!read_element(d + pos, n - pos, &dims) || dims.type != miINT32 || dims.nbytes < 8) goto bad;
    pos += dims.total;
    if (!read_element(d + pos, n - pos, &name)) goto bad;
    pos += name.total;

    v->mxClass = (int)(get_le32(flags.data) & 0xFFu);
    v->rows = get_le32(dims.data);
    v->cols = 1;
    for (uint32_t k = 1; k < dims.nbytes / 4; k++)
        v->cols *= get_le32(dims.data + 4 * k);
    v->numel = v->rows * v->cols;
    {
        size_t len = name.nbytes < sizeof v->name - 1 ? name.nbytes : sizeof v->name - 1;
        memcpy(v->name, name.data, len);
        v->name[len] = '\0';
    }

    if (v->mxClass >= MAT_CLASS_DOUBLE && v->mxClass <= 15) {      /* numeric classes */
        if (!read_element(d + pos, n - pos, &e)) goto bad;
        const size_t sz = type_size(e.type);
        if (sz == 0 || e.nbytes / sz < v->numel) goto bad;
        if (v->numel > 0) {
            v->re = malloc(v->numel * sizeof *v->re);
            if (!v->re) {
                set_err(err, errLen, "out of memory");
                return -1;
            }
            convert(e.type, e.data, v->re, v->numel);
        }
    } else if (v->mxClass == MAT_CLASS_STRUCT) {
        Element fnLen, fnNames;
        if (!read_element(d + pos, n - pos, &fnLen) || fnLen.nbytes < 4) goto bad;
        pos += fnLen.total;
        if (!read_element(d + pos, n - pos, &fnNames)) goto bad;
        pos += fnNames.total;
        const uint32_t len = get_le32(fnLen.data);
        const int nf = len ? (int)(fnNames.nbytes / len) : 0;
        if (v->numel > 0 && nf > 0) {
            v->fields = calloc((size_t)nf, sizeof *v->fields);
            if (!v->fields) {
                set_err(err, errLen, "out of memory");
                return -1;
            }
            v->nFields = nf;
            for (int f = 0; f < nf; f++) {           /* the values of the first element */
                if (!read_element(d + pos, n - pos, &e) || e.type != miMATRIX) goto bad;
                if (parse_matrix(e.data, e.nbytes, &v->fields[f], depth + 1, err, errLen))
                    return -1;
                size_t fl = len < sizeof v->name - 1 ? len : sizeof v->name - 1;
                memcpy(v->fields[f].name, fnNames.data + (size_t)f * len, fl);
                v->fields[f].name[fl] = '\0';
                pos += e.total;
            }
        }
    } else {
        v->mxClass = 0;                               /* cell, char, sparse, object: skipped */
    }
    return 0;

bad:
    set_err(err, errLen, "corrupted variable '%s'", v->name);
    return -1;
}

#ifdef RADAR_HAVE_ZLIB
/* Inflates a miCOMPRESSED element. Returns a malloc'ed buffer or NULL. */
static unsigned char *inflate_element(const unsigned char *d, size_t n, size_t *outLen)
{
    size_t cap = 4 * n + 1024;
    unsigned char *out = malloc(cap);
    z_stream z;
    int r;

    if (!out)
        return NULL;
    memset(&z, 0, sizeof z);
    if (inflateInit(&z) != Z_OK) {
        free(out);
        return NULL;
    }
    z.next_in = (Bytef *)d;
    z.avail_in = (uInt)n;
    for (;;) {
        z.next_out = out + z.total_out;
        z.avail_out = (uInt)(cap - z.total_out);
        r = inflate(&z, Z_NO_FLUSH);
        if (r == Z_STREAM_END)
            break;
        if (r != Z_OK && r != Z_BUF_ERROR) {
            inflateEnd(&z);
            free(out);
            return NULL;
        }
        if (z.avail_out == 0) {                       /* output buffer full: grow it */
            unsigned char *p = realloc(out, 2 * cap);
            if (!p) {
                inflateEnd(&z);
                free(out);
                return NULL;
            }
            out = p;
            cap *= 2;
        } else if (z.avail_in == 0) {                 /* truncated stream */
            inflateEnd(&z);
            free(out);
            return NULL;
        }
    }
    *outLen = z.total_out;
    inflateEnd(&z);
    return out;
}
#endif

static int add_var(MatFile *mf, const unsigned char *d, size_t n, char *err, size_t errLen)
{
    MatVar *p = realloc(mf->vars, (size_t)(mf->nVars + 1) * sizeof *p);
    if (!p) {
        set_err(err, errLen, "out of memory");
        return -1;
    }
    mf->vars = p;
    if (parse_matrix(d, n, &mf->vars[mf->nVars], 0, err, errLen)) {
        free_var(&mf->vars[mf->nVars]);
        return -1;
    }
    mf->nVars++;
    return 0;
}

int mat_parse(const unsigned char *data, size_t size, MatFile *mf, char *err, size_t errLen)
{
    size_t pos = 128;

    memset(mf, 0, sizeof *mf);
    set_err(err, errLen, "%s", "");
    if (size >= 10 && memcmp(data, "MATLAB 7.3", 10) == 0) {
        set_err(err, errLen, "MAT-file -v7.3 (HDF5) is not supported: save with -v6 or -v7");
        return -1;
    }
    if (size < 128 || data[124] != 0x00 || data[125] != 0x01) {
        if (size >= 128 && data[126] == 'M' && data[127] == 'I')
            set_err(err, errLen, "big-endian MAT-file is not supported");
        else
            set_err(err, errLen, "not a MAT-file (level 5)");
        return -1;
    }
    if (data[126] != 'I' || data[127] != 'M') {
        set_err(err, errLen, "big-endian MAT-file is not supported");
        return -1;
    }
    while (pos + 8 <= size) {
        Element e;
        if (!read_element(data + pos, size - pos, &e)) {
            set_err(err, errLen, "truncated MAT-file");
            mat_free(mf);
            return -1;
        }
        if (e.type == miMATRIX) {
            if (add_var(mf, e.data, e.nbytes, err, errLen)) {
                mat_free(mf);
                return -1;
            }
        } else if (e.type == miCOMPRESSED) {
#ifdef RADAR_HAVE_ZLIB
            size_t len = 0;
            Element inner;
            unsigned char *raw = inflate_element(e.data, e.nbytes, &len);
            int r = -1;
            if (!raw)
                set_err(err, errLen, "corrupted compressed variable");
            else if (!read_element(raw, len, &inner) || inner.type != miMATRIX)
                set_err(err, errLen, "corrupted compressed variable");
            else
                r = add_var(mf, inner.data, inner.nbytes, err, errLen);
            free(raw);
            if (r) {
                mat_free(mf);
                return -1;
            }
#else
            set_err(err, errLen, "compressed MAT-file (-v7): built without zlib, save with -v6 "
                    "(cfg.io.matFormat = '-v6')");
            mat_free(mf);
            return -1;
#endif
        }
        pos += e.total;
    }
    return 0;
}

int mat_load(const char *path, MatFile *mf, char *err, size_t errLen)
{
    FILE *f = fopen(path, "rb");
    unsigned char *data;
    long size;
    int r;

    memset(mf, 0, sizeof *mf);
    if (!f) {
        set_err(err, errLen, "cannot open %s", path);
        return -1;
    }
    if (fseek(f, 0, SEEK_END) != 0 || (size = ftell(f)) < 0 || fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        set_err(err, errLen, "cannot read %s", path);
        return -1;
    }
    data = malloc(size > 0 ? (size_t)size : 1);
    if (!data || fread(data, 1, (size_t)size, f) != (size_t)size) {
        free(data);
        fclose(f);
        set_err(err, errLen, "cannot read %s", path);
        return -1;
    }
    fclose(f);
    r = mat_parse(data, (size_t)size, mf, err, errLen);
    free(data);
    return r;
}

void mat_free(MatFile *mf)
{
    for (int i = 0; i < mf->nVars; i++)
        free_var(&mf->vars[i]);
    free(mf->vars);
    memset(mf, 0, sizeof *mf);
}

const MatVar *mat_get(const MatFile *mf, const char *name)
{
    for (int i = 0; i < mf->nVars; i++)
        if (!strcmp(mf->vars[i].name, name))
            return &mf->vars[i];
    return NULL;
}

const MatVar *mat_field(const MatVar *s, const char *name)
{
    if (!s)
        return NULL;
    for (int i = 0; i < s->nFields; i++)
        if (!strcmp(s->fields[i].name, name))
            return &s->fields[i];
    return NULL;
}

int mat_has_zlib(void)
{
#ifdef RADAR_HAVE_ZLIB
    return 1;
#else
    return 0;
#endif
}

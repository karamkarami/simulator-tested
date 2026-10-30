/*
 * mat_file.h - minimal reader of MATLAB MAT-files (level 5: save -v6 and,
 * when built with zlib, the compressed save -v7).
 *
 * Reads the numeric variables (as double) and 1x1 structs of numeric fields,
 * which is what the DDC files of the simulator (ddc_<n>.mat: I, Q, hdr)
 * hold. Numeric data may be stored in a smaller type than its class (MATLAB
 * does this, e.g. int16 data in uint8 when the values fit); it is converted.
 * Complex parts, cells, chars, sparse and object arrays are skipped.
 * MAT-files -v7.3 (HDF5) and big-endian files are not supported.
 */
#ifndef MAT_FILE_H
#define MAT_FILE_H

#include <stddef.h>

enum {
    MAT_CLASS_STRUCT = 2,
    MAT_CLASS_DOUBLE = 6,
    MAT_CLASS_INT16  = 10
};

typedef struct MatVar {
    char           name[64];
    int            mxClass;      /* MATLAB class (MAT_CLASS_..., 0 when skipped)   */
    size_t         rows, cols;   /* first dimension, product of the others        */
    size_t         numel;
    double        *re;           /* numeric: the values in column-major order     */
    int            nFields;      /* struct: the fields of its first element       */
    struct MatVar *fields;
} MatVar;

typedef struct {
    int     nVars;
    MatVar *vars;
} MatFile;

/* Reads the whole file. Returns 0 on success; on failure err holds the reason. */
int mat_load(const char *path, MatFile *mf, char *err, size_t errLen);
/* Same from a file image in memory. */
int mat_parse(const unsigned char *data, size_t size, MatFile *mf, char *err, size_t errLen);
void mat_free(MatFile *mf);

/* The variable or field of that name, or NULL. */
const MatVar *mat_get(const MatFile *mf, const char *name);
const MatVar *mat_field(const MatVar *s, const char *name);

/* 1 if the library was built with zlib (compressed -v7 files). */
int mat_has_zlib(void);

#endif

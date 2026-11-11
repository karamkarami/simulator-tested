/*
 * mat_write.h - minimal writer of MATLAB MAT-files (level 5, as MATLAB's
 * save -v6, or -v7 with every variable in a compressed element).
 *
 * Enough for the files of the simulator: int16 matrices (DDC output I, Q, raw
 * ADC data) and a 1x1 struct of real double fields (the PRI header hdr).
 * Compressed files (-v7) use zlib when the program is built with it
 * (RADAR_HAVE_ZLIB), else a valid zlib stream of stored (uncompressed)
 * blocks, which MATLAB, Octave and mat_file.c read as well.
 */
#ifndef MAT_WRITE_H
#define MAT_WRITE_H

#include <stddef.h>
#include <stdint.h>

typedef struct MatWriter MatWriter;

/* compressed: 0 = -v6, 1 = -v7. Returns NULL when the file cannot be created. */
MatWriter *mat_write_open(const char *path, int compressed);
/* int16 matrix [rows x cols], column-major data. */
void mat_write_int16(MatWriter *w, const char *name, const int16_t *data, size_t rows, size_t cols);
/* A 1x1 struct of real double fields: field k has n[k] values (a row vector). */
void mat_write_struct(MatWriter *w, const char *name, int nFields, const char *const *fieldNames,
                      const double *const *values, const int *n);
/* Writes the file and frees the writer; returns 0 on success. */
int  mat_write_close(MatWriter *w);
/* 1 if built with zlib. */
int  mat_write_has_zlib(void);

#endif

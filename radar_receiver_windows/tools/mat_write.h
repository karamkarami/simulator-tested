/*
 * mat_write.h - minimal writer of MATLAB MAT-files (level 5, as MATLAB's
 * save -v6, or -v7 with every variable in a compressed element).
 *
 * Enough for the files of the simulator: int16 matrices (DDC output I, Q, raw
 * ADC data), a 1x1 struct of real double fields (the PRI header hdr), and for
 * the logs of the receiver single (real or complex, N-D) and double matrices
 * and strings.
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
/* single matrix of nd dimensions (column-major data), complex when im is not NULL. */
void mat_write_single(MatWriter *w, const char *name, const float *re, const float *im, const size_t *dims, int nd);
/* double matrix [rows x cols], column-major data. */
void mat_write_double(MatWriter *w, const char *name, const double *data, size_t rows, size_t cols);
/* char row vector (a string). */
void mat_write_char(MatWriter *w, const char *name, const char *text);
/* Writes the file and frees the writer; returns 0 on success. */
int  mat_write_close(MatWriter *w);

/* Nested structs (radar_parameters.mat): a 1x1 struct is built field by field (the values are
 * copied), a struct field takes another struct (ownership passes to the parent), then
 * mat_write_tree writes it as variable 'name' and frees it. Field names: at most 63 characters. */
typedef struct MatStruct MatStruct;
MatStruct *mat_struct_new(void);
void mat_struct_double(MatStruct *s, const char *field, const double *v, size_t rows, size_t cols);
void mat_struct_scalar(MatStruct *s, const char *field, double v);
void mat_struct_single(MatStruct *s, const char *field, const float *re, const float *im, const size_t *dims, int nd);
void mat_struct_char(MatStruct *s, const char *field, const char *text);
void mat_struct_struct(MatStruct *s, const char *field, MatStruct *child);
void mat_struct_complex(MatStruct *s, const char *field, const double *re, const double *im, size_t rows, size_t cols);
void mat_struct_logical(MatStruct *s, const char *field, const double *v, size_t rows, size_t cols);
/* char matrix [rows x cols], column-major text (no terminating zero needed) */
void mat_struct_char_matrix(MatStruct *s, const char *field, const char *text, size_t rows, size_t cols);
/* A cell array [rows x cols]: its rows*cols values are added in column-major order with the
 * mat_struct_* functions above (the field name is ignored), then it is put into a struct with
 * mat_struct_struct. A struct array [rows x cols]: the fields of every element are added element
 * after element, always the same fields in the same order. */
MatStruct *mat_cell_new(size_t rows, size_t cols);
/* moves the fields of src to the end of dst and frees src */
void mat_struct_append(MatStruct *dst, MatStruct *src);
MatStruct *mat_struct_array_new(size_t rows, size_t cols);
void mat_struct_free(MatStruct *s);
void mat_write_tree(MatWriter *w, const char *name, MatStruct *s);
/* 1 if built with zlib. */
int  mat_write_has_zlib(void);

#endif

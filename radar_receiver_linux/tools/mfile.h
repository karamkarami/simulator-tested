/*
 * mfile.h - evaluates a MATLAB configuration function such as radar_config.m
 * without MATLAB (used by radar_export_config).
 *
 * The file is run as a script of assignments. Supported: the function line
 * (its output variable is the result, e.g. "function cfg = radar_config()"),
 * assignments to variables and struct fields (cfg.a.b = ...), comments (%,
 * %{ ... %}), continuation lines (...), statements separated by newlines, ';'
 * and ','; real numbers, 'strings' and "strings", [matrices] (rows with ';'
 * or newlines, elements with ',' or spaces, also concatenation of strings),
 * {cell arrays}, the operators + - * / ^ .* ./ .^ ' .' == ~= < > <= >= & | &&
 * || ~ and the colon a:b, a:s:b; indexing x(i), x(i,j), c{i} and field
 * access s.f; the constants Inf NaN pi eps true false; and the functions
 * fullfile fileparts mfilename struct sqrt abs exp log log10 log2 sin cos tan
 * sind cosd tand asind acosd atand deg2rad rad2deg round floor ceil fix mod
 * rem min max sum prod mean numel length zeros ones linspace logical double
 * isempty strcmp strcmpi lower upper.
 *
 * Other functions (e.g. target_scenario, addpath) give an "unevaluated"
 * value: no error unless the exporter reads it. if/for/while/switch blocks
 * and subfunctions are not supported (an error names the line).
 */
#ifndef MFILE_H
#define MFILE_H

#include <stddef.h>

enum { MV_NUM = 1, MV_STR, MV_CELL, MV_STRUCT, MV_UNEVAL };

typedef struct MValue {
    int     type;
    int     rows, cols;          /* MV_NUM, MV_STR (1 x n), MV_CELL                  */
    double *re;                  /* MV_NUM: values, column-major                     */
    char   *str;                 /* MV_STR: text; MV_UNEVAL: why it was not evaluated */
    struct MValue **cells;       /* MV_CELL: rows*cols values                        */
    int     nFields;             /* MV_STRUCT (1x1): fields in the order of creation  */
    char  **names;
    struct MValue **fields;
} MValue;

/* Runs the file, then the extra statements (each a line of MATLAB, e.g.
 * "cfg.prf.mode = 'staggered';"; may be NULL). Returns the output variable
 * of the function line (or 'cfg' without one), NULL on an error (message in
 * err, with the line). */
MValue *mfile_eval(const char *path, const char *const *extra, int nExtra, char *err, size_t errLen);
void    mvalue_free(MValue *v);

/* Field of a struct by a dotted path ("rx.maxRange_m"), NULL if missing. */
const MValue *mvalue_get(const MValue *s, const char *path);

#endif

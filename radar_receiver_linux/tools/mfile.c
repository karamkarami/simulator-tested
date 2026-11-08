/*
 * mfile.c - evaluates a MATLAB configuration function such as radar_config.m
 * (see mfile.h for what is supported).
 *
 * The text is split into tokens (lexer), which a recursive-descent parser
 * evaluates statement by statement: assignments store the value of the right
 * side in the variable or struct field, other statements (addpath(...)) are
 * evaluated and dropped. Operator precedence as in MATLAB: || && | & ,
 * comparisons, colon, + -, * / .* ./ , unary + - ~, ^ .^ , transpose. Inside
 * [ ] and { } a space separates elements ([1 -2 1] has three, [1 - 2] one),
 * a newline or ';' separates rows.
 *
 * Errors end the evaluation (longjmp back to mfile_eval, which frees what it
 * can; the program then stops anyway).
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "mfile.h"

#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846

static int strcasecmp_portable(const char *a, const char *b)
{
    while (*a && tolower((unsigned char)*a) == tolower((unsigned char)*b)) { a++; b++; }
    return tolower((unsigned char)*a) - tolower((unsigned char)*b);
}

/* ------------------------------------------------------------------------ */
/* Values                                                                   */
/* ------------------------------------------------------------------------ */

static void *xmalloc(size_t n)
{
    void *p = calloc(1, n ? n : 1);
    if (!p) { fprintf(stderr, "out of memory\n"); exit(3); }
    return p;
}

static char *xstrdup(const char *s)
{
    char *d = xmalloc(strlen(s) + 1);
    strcpy(d, s);
    return d;
}

static MValue *mv_new(int type)
{
    MValue *v = xmalloc(sizeof *v);
    v->type = type;
    return v;
}

static MValue *mv_num(int rows, int cols)
{
    MValue *v = mv_new(MV_NUM);
    v->rows = rows;
    v->cols = cols;
    v->re = xmalloc(sizeof(double) * (size_t)(rows * cols > 0 ? rows * cols : 1));
    return v;
}

static MValue *mv_scalar(double x)
{
    MValue *v = mv_num(1, 1);
    v->re[0] = x;
    return v;
}

static MValue *mv_str(const char *s)
{
    MValue *v = mv_new(MV_STR);
    v->str = xstrdup(s);
    v->rows = 1;
    v->cols = (int)strlen(s);
    return v;
}

static MValue *mv_uneval(const char *why)
{
    MValue *v = mv_new(MV_UNEVAL);
    v->str = xstrdup(why);
    return v;
}

void mvalue_free(MValue *v)
{
    if (!v) return;
    free(v->re);
    free(v->str);
    if (v->cells) {
        for (int i = 0; i < v->rows * v->cols; i++) mvalue_free(v->cells[i]);
        free(v->cells);
    }
    for (int i = 0; i < v->nFields; i++) {
        free(v->names[i]);
        mvalue_free(v->fields[i]);
    }
    free(v->names);
    free(v->fields);
    free(v);
}

static MValue *mv_copy(const MValue *v)
{
    MValue *c = mv_new(v->type);
    c->rows = v->rows;
    c->cols = v->cols;
    if (v->re) {
        c->re = xmalloc(sizeof(double) * (size_t)(v->rows * v->cols > 0 ? v->rows * v->cols : 1));
        memcpy(c->re, v->re, sizeof(double) * (size_t)(v->rows * v->cols));
    }
    if (v->str) c->str = xstrdup(v->str);
    if (v->cells) {
        c->cells = xmalloc(sizeof(MValue *) * (size_t)(v->rows * v->cols > 0 ? v->rows * v->cols : 1));
        for (int i = 0; i < v->rows * v->cols; i++) c->cells[i] = mv_copy(v->cells[i]);
    }
    if (v->nFields) {
        c->nFields = v->nFields;
        c->names = xmalloc(sizeof(char *) * (size_t)v->nFields);
        c->fields = xmalloc(sizeof(MValue *) * (size_t)v->nFields);
        for (int i = 0; i < v->nFields; i++) {
            c->names[i] = xstrdup(v->names[i]);
            c->fields[i] = mv_copy(v->fields[i]);
        }
    }
    return c;
}

static MValue **struct_slot(MValue *s, const char *name, int create)
{
    for (int i = 0; i < s->nFields; i++)
        if (!strcmp(s->names[i], name)) return &s->fields[i];
    if (!create) return NULL;
    s->names = realloc(s->names, sizeof(char *) * (size_t)(s->nFields + 1));
    s->fields = realloc(s->fields, sizeof(MValue *) * (size_t)(s->nFields + 1));
    if (!s->names || !s->fields) { fprintf(stderr, "out of memory\n"); exit(3); }
    s->names[s->nFields] = xstrdup(name);
    s->fields[s->nFields] = NULL;
    return &s->fields[s->nFields++];
}

const MValue *mvalue_get(const MValue *s, const char *path)
{
    char part[128];
    while (s && *path) {
        const char *dot = strchr(path, '.');
        size_t n = dot ? (size_t)(dot - path) : strlen(path);
        if (n >= sizeof part || s->type != MV_STRUCT) return NULL;
        memcpy(part, path, n);
        part[n] = '\0';
        MValue **slot = struct_slot((MValue *)s, part, 0);
        s = slot ? *slot : NULL;
        path += n + (dot ? 1 : 0);
    }
    return s;
}

/* ------------------------------------------------------------------------ */
/* Lexer                                                                    */
/* ------------------------------------------------------------------------ */

enum { T_EOF, T_NUM, T_STR, T_ID, T_OP, T_NL };

typedef struct {
    int    type;
    char   text[8];      /* operator, or the first characters of an identifier */
    char  *s;            /* identifier or string                               */
    double num;
    int    space;        /* whitespace before the token                        */
    int    line;
} Token;

typedef struct {
    Token  *tok;
    int     n, cap, pos;
    const char *source;  /* file name or "--set" for messages                  */
    int     lineBase;
    /* evaluation */
    int     nVars;
    char  **varNames;
    MValue **vars;
    char   outName[64];
    const char *filePath;
    int     ctx[64];     /* 1 inside [ ] or { } (spaces separate elements)       */
    int     depth;
    int     lexLine;     /* line of the lexer (errors before the tokens exist)    */
    jmp_buf jb;
    char   *err;
    size_t  errLen;
} Parser;

static void fail(Parser *p, const char *fmt, ...)
{
    char msg[400];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(msg, sizeof msg, fmt, ap);
    va_end(ap);
    const int line = p->tok && p->pos < p->n ? p->tok[p->pos].line : p->lexLine;
    if (p->err && p->errLen) snprintf(p->err, p->errLen, "%s:%d: %s", p->source, line, msg);
    longjmp(p->jb, 1);
}

static void push_token(Parser *p, Token t)
{
    if (p->n == p->cap) {
        p->cap = p->cap ? 2 * p->cap : 1024;
        p->tok = realloc(p->tok, sizeof(Token) * (size_t)p->cap);
        if (!p->tok) { fprintf(stderr, "out of memory\n"); exit(3); }
    }
    p->tok[p->n++] = t;
}

static int value_end(const Parser *p)
{
    if (p->n == 0) return 0;
    const Token *t = &p->tok[p->n - 1];
    if (t->type == T_NUM || t->type == T_ID || t->type == T_STR) return 1;
    return t->type == T_OP && (!strcmp(t->text, ")") || !strcmp(t->text, "]") || !strcmp(t->text, "}") ||
                               !strcmp(t->text, "'") || !strcmp(t->text, ".'"));
}

static void lex(Parser *p, const char *src)
{
    char stack[256];
    int sp = 0, line = p->lineBase;
    const char *c = src;
    int space = 0, lineStart = 1;
    while (*c) {
        if (lineStart) {                                  /* block comment %{ ... %} on lines of their own */
            const char *q = c;
            while (*q == ' ' || *q == '\t') q++;
            if (q[0] == '%' && q[1] == '{') {
                const char *e = q + 2;
                while (*e == ' ' || *e == '\t' || *e == '\r') e++;
                if (*e == '\n' || !*e) {
                    int nest = 1;
                    c = e;
                    while (*c && nest > 0) {
                        if (*c == '\n') {
                            line++;
                            const char *r = c + 1;
                            while (*r == ' ' || *r == '\t') r++;
                            if (r[0] == '%' && (r[1] == '{' || r[1] == '}')) {
                                const char *t = r + 2;
                                while (*t == ' ' || *t == '\t' || *t == '\r') t++;
                                if (*t == '\n' || !*t) nest += r[1] == '{' ? 1 : -1;
                            }
                        }
                        c++;
                    }
                    while (*c && *c != '\n') c++;
                    continue;
                }
            }
            lineStart = 0;
        }
        if (*c == ' ' || *c == '\t' || *c == '\r') { space = 1; c++; continue; }
        Token t;
        memset(&t, 0, sizeof t);
        t.space = space;
        t.line = line;
        p->lexLine = line;
        space = 0;
        if (*c == '\n') {
            line++;
            c++;
            lineStart = 1;
            if (sp == 0 || stack[sp - 1] != '(') {        /* a statement end, or a row in [ ] { } */
                t.type = T_NL;
                push_token(p, t);
            }
            space = 1;
            continue;
        }
        if (*c == '%') { while (*c && *c != '\n') c++; continue; }
        if (c[0] == '.' && c[1] == '.' && c[2] == '.') {     /* continuation: the rest of the line is ignored */
            while (*c && *c != '\n') c++;
            if (*c) { c++; line++; }
            space = 1;
            continue;
        }
        if (isdigit((unsigned char)c[0]) || (c[0] == '.' && isdigit((unsigned char)c[1]))) {
            char *end;
            t.type = T_NUM;
            t.num = strtod(c, &end);
            if (end == c) fail(p, "bad number");
            if (*end == 'i' || *end == 'j') { p->pos = p->n; fail(p, "complex numbers are not supported"); }
            c = end;
            push_token(p, t);
            continue;
        }
        if (isalpha((unsigned char)*c) || *c == '_') {
            const char *s = c;
            while (isalnum((unsigned char)*c) || *c == '_') c++;
            t.type = T_ID;
            t.s = xmalloc((size_t)(c - s) + 1);
            memcpy(t.s, s, (size_t)(c - s));
            push_token(p, t);
            continue;
        }
        if (*c == '"' || (*c == '\'' && !(value_end(p) && (!t.space || sp == 0)))) {   /* a string */
            const char q = *c++;
            size_t cap = 64, n = 0;
            char *s = xmalloc(cap);
            for (;;) {
                if (!*c || *c == '\n') { p->pos = p->n; fail(p, "unterminated string"); }
                if (*c == q) {
                    if (c[1] == q) { c++; } else { c++; break; }
                }
                if (n + 2 >= cap) { cap *= 2; s = realloc(s, cap); if (!s) exit(3); }
                s[n++] = *c++;
            }
            s[n] = '\0';
            t.type = T_STR;
            t.s = s;
            push_token(p, t);
            continue;
        }
        static const char *const ops[] = { "==", "~=", "!=", "<=", ">=", "&&", "||", ".*", "./", ".^", ".'", ".\\",
                                           "+", "-", "*", "/", "\\", "^", "<", ">", "&", "|", "~", "!", "=", "(", ")",
                                           "[", "]", "{", "}", ".", ":", ",", ";", "'", "@" };
        size_t k;
        for (k = 0; k < sizeof ops / sizeof ops[0]; k++)
            if (!strncmp(c, ops[k], strlen(ops[k]))) break;
        if (k == sizeof ops / sizeof ops[0]) { p->pos = p->n; fail(p, "unexpected character '%c'", *c); }
        t.type = T_OP;
        strcpy(t.text, ops[k]);
        c += strlen(ops[k]);
        if (!strcmp(t.text, "(") || !strcmp(t.text, "[") || !strcmp(t.text, "{")) {
            if (sp < (int)sizeof stack) stack[sp++] = t.text[0];
        } else if ((!strcmp(t.text, ")") || !strcmp(t.text, "]") || !strcmp(t.text, "}")) && sp > 0) {
            sp--;
        }
        if (!strcmp(t.text, "!=")) strcpy(t.text, "~=");
        if (!strcmp(t.text, "!")) strcpy(t.text, "~");
        push_token(p, t);
    }
    Token e;
    memset(&e, 0, sizeof e);
    e.type = T_NL;
    e.line = line;
    push_token(p, e);
    e.type = T_EOF;
    push_token(p, e);
}

/* ------------------------------------------------------------------------ */
/* Helpers of the evaluation                                                */
/* ------------------------------------------------------------------------ */

static Token *cur(Parser *p) { return &p->tok[p->pos]; }
static int is_op(const Token *t, const char *op) { return t->type == T_OP && !strcmp(t->text, op); }
static int accept(Parser *p, const char *op)
{
    if (is_op(cur(p), op)) { p->pos++; return 1; }
    return 0;
}
static void expect(Parser *p, const char *op)
{
    if (!accept(p, op)) fail(p, "'%s' expected", op);
}

static MValue *get_var(Parser *p, const char *name)
{
    for (int i = 0; i < p->nVars; i++)
        if (!strcmp(p->varNames[i], name)) return p->vars[i];
    return NULL;
}

static MValue **var_slot(Parser *p, const char *name)
{
    for (int i = 0; i < p->nVars; i++)
        if (!strcmp(p->varNames[i], name)) return &p->vars[i];
    p->varNames = realloc(p->varNames, sizeof(char *) * (size_t)(p->nVars + 1));
    p->vars = realloc(p->vars, sizeof(MValue *) * (size_t)(p->nVars + 1));
    if (!p->varNames || !p->vars) exit(3);
    p->varNames[p->nVars] = xstrdup(name);
    p->vars[p->nVars] = NULL;
    return &p->vars[p->nVars++];
}

/* A string as a numeric row of character codes (MATLAB arithmetic on chars). */
static MValue *to_num(Parser *p, MValue *v)
{
    if (v->type == MV_NUM || v->type == MV_UNEVAL) return v;
    if (v->type == MV_STR) {
        MValue *n = mv_num(1, v->cols);
        for (int i = 0; i < v->cols; i++) n->re[i] = (unsigned char)v->str[i];
        mvalue_free(v);
        return n;
    }
    const int type = v->type;
    mvalue_free(v);
    fail(p, "a number or string was expected (not a %s)", type == MV_CELL ? "cell array" : "struct");
    return NULL;
}

static int numel(const MValue *v) { return v->rows * v->cols; }

/* Element-wise operation with scalar and row/column expansion. */
static MValue *elementwise(Parser *p, MValue *a, MValue *b, double (*f)(double, double))
{
    if (a->type == MV_UNEVAL) { mvalue_free(b); return a; }
    if (b->type == MV_UNEVAL) { mvalue_free(a); return b; }
    a = to_num(p, a);
    b = to_num(p, b);
    const int r = a->rows == b->rows ? a->rows : a->rows == 1 ? b->rows : b->rows == 1 ? a->rows : -1;
    const int c = a->cols == b->cols ? a->cols : a->cols == 1 ? b->cols : b->cols == 1 ? a->cols : -1;
    if (r < 0 || c < 0) {
        int ar = a->rows, ac = a->cols, br = b->rows, bc = b->cols;
        mvalue_free(a);
        mvalue_free(b);
        fail(p, "sizes %dx%d and %dx%d do not match", ar, ac, br, bc);
    }
    MValue *y = mv_num(r, c);
    for (int j = 0; j < c; j++)
        for (int i = 0; i < r; i++)
            y->re[i + j * r] = f(a->re[(a->rows == 1 ? 0 : i) + (a->cols == 1 ? 0 : j) * a->rows],
                                 b->re[(b->rows == 1 ? 0 : i) + (b->cols == 1 ? 0 : j) * b->rows]);
    mvalue_free(a);
    mvalue_free(b);
    return y;
}

static double f_add(double a, double b) { return a + b; }
static double f_sub(double a, double b) { return a - b; }
static double f_mul(double a, double b) { return a * b; }
static double f_div(double a, double b) { return a / b; }
static double f_ldiv(double a, double b) { return b / a; }
static double f_pow(double a, double b) { return pow(a, b); }
static double f_eq(double a, double b) { return a == b; }
static double f_ne(double a, double b) { return a != b; }
static double f_lt(double a, double b) { return a < b; }
static double f_gt(double a, double b) { return a > b; }
static double f_le(double a, double b) { return a <= b; }
static double f_ge(double a, double b) { return a >= b; }
static double f_and(double a, double b) { return a != 0 && b != 0; }
static double f_or(double a, double b) { return a != 0 || b != 0; }
static double f_mod(double a, double b) { return b == 0 ? a : a - floor(a / b) * b; }
static double f_rem(double a, double b) { return b == 0 ? a : a - trunc(a / b) * b; }
static double f_min(double a, double b) { return isnan(a) ? b : isnan(b) ? a : a < b ? a : b; }
static double f_max(double a, double b) { return isnan(a) ? b : isnan(b) ? a : a > b ? a : b; }

static MValue *transpose(Parser *p, MValue *a)
{
    if (a->type == MV_UNEVAL) return a;
    if (a->type == MV_CELL) {
        MValue *y = mv_new(MV_CELL);
        y->rows = a->cols;
        y->cols = a->rows;
        y->cells = xmalloc(sizeof(MValue *) * (size_t)(numel(a) > 0 ? numel(a) : 1));
        for (int i = 0; i < a->rows; i++)
            for (int j = 0; j < a->cols; j++) { y->cells[j + i * y->rows] = a->cells[i + j * a->rows]; a->cells[i + j * a->rows] = NULL; }
        mvalue_free(a);
        return y;
    }
    a = to_num(p, a);
    MValue *y = mv_num(a->cols, a->rows);
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < a->cols; j++) y->re[j + i * a->cols] = a->re[i + j * a->rows];
    mvalue_free(a);
    return y;
}

static MValue *matmul(Parser *p, MValue *a, MValue *b)
{
    if (a->type == MV_UNEVAL) { mvalue_free(b); return a; }
    if (b->type == MV_UNEVAL) { mvalue_free(a); return b; }
    a = to_num(p, a);
    b = to_num(p, b);
    if (numel(a) == 1 || numel(b) == 1) return elementwise(p, a, b, f_mul);
    if (a->cols != b->rows) {
        mvalue_free(a);
        mvalue_free(b);
        fail(p, "inner matrix dimensions do not agree");
    }
    MValue *y = mv_num(a->rows, b->cols);
    for (int i = 0; i < a->rows; i++)
        for (int j = 0; j < b->cols; j++) {
            double s = 0;
            for (int k = 0; k < a->cols; k++) s += a->re[i + k * a->rows] * b->re[k + j * b->rows];
            y->re[i + j * a->rows] = s;
        }
    mvalue_free(a);
    mvalue_free(b);
    return y;
}

/* Horizontal (dim 2) or vertical (dim 1) concatenation of a list of values. */
static MValue *concat(Parser *p, MValue **v, int n, int dim, int cell)
{
    int k = 0;
    for (int i = 0; i < n; i++)                           /* [] disappears */
        if (!(v[i]->type == MV_NUM && numel(v[i]) == 0 && !cell)) v[k++] = v[i];
        else mvalue_free(v[i]);
    n = k;
    if (n == 0) return cell ? mv_new(MV_CELL) : mv_num(0, 0);
    if (n == 1) return v[0];
    for (int i = 0; i < n; i++)
        if (v[i]->type == MV_UNEVAL) {
            MValue *u = v[i];
            for (int j = 0; j < n; j++) if (j != i) mvalue_free(v[j]);
            return u;
        }
    int allStr = !cell;
    for (int i = 0; i < n && allStr; i++) allStr = v[i]->type == MV_STR;
    if (allStr && dim == 2) {
        size_t len = 0;
        for (int i = 0; i < n; i++) len += strlen(v[i]->str);
        char *s = xmalloc(len + 1);
        for (int i = 0; i < n; i++) { strcat(s, v[i]->str); mvalue_free(v[i]); }
        MValue *y = mv_str(s);
        free(s);
        return y;
    }
    if (cell) {
        for (int i = 0; i < n; i++)
            if (v[i]->type != MV_CELL) { for (int j = 0; j < n; j++) mvalue_free(v[j]); fail(p, "cannot concatenate"); }
    } else {
        for (int i = 0; i < n; i++) {
            if (v[i]->type == MV_CELL || v[i]->type == MV_STRUCT) {
                for (int j = 0; j < n; j++) mvalue_free(v[j]);
                fail(p, "cannot concatenate cells or structs in [ ]");
            }
            v[i] = to_num(p, v[i]);
        }
    }
    int rows = 0, cols = 0;
    for (int i = 0; i < n; i++) {
        if (dim == 2) {
            if (i > 0 && v[i]->rows != v[0]->rows) { for (int j = 0; j < n; j++) mvalue_free(v[j]); fail(p, "dimensions of the concatenated values do not agree"); }
            cols += v[i]->cols;
            rows = v[0]->rows;
        } else {
            if (i > 0 && v[i]->cols != v[0]->cols) { for (int j = 0; j < n; j++) mvalue_free(v[j]); fail(p, "dimensions of the concatenated rows do not agree"); }
            rows += v[i]->rows;
            cols = v[0]->cols;
        }
    }
    MValue *y = cell ? mv_new(MV_CELL) : mv_num(rows, cols);
    y->rows = rows;
    y->cols = cols;
    if (cell) y->cells = xmalloc(sizeof(MValue *) * (size_t)(rows * cols > 0 ? rows * cols : 1));
    int off = 0;
    for (int i = 0; i < n; i++) {
        for (int r = 0; r < v[i]->rows; r++)
            for (int c = 0; c < v[i]->cols; c++) {
                const int rr = dim == 1 ? off + r : r, cc = dim == 2 ? off + c : c;
                if (cell) { y->cells[rr + cc * rows] = v[i]->cells[r + c * v[i]->rows]; v[i]->cells[r + c * v[i]->rows] = NULL; }
                else y->re[rr + cc * rows] = v[i]->re[r + c * v[i]->rows];
            }
        off += dim == 1 ? v[i]->rows : v[i]->cols;
        mvalue_free(v[i]);
    }
    return y;
}

/* ------------------------------------------------------------------------ */
/* Functions                                                                */
/* ------------------------------------------------------------------------ */

static MValue *arg_num(Parser *p, MValue **a, int n, int i, const char *f)
{
    if (i >= n) fail(p, "%s: not enough arguments", f);
    if (a[i]->type == MV_UNEVAL) return a[i];
    a[i] = to_num(p, a[i]);
    return a[i];
}

static MValue *map1(Parser *p, MValue *x, double (*f)(double))
{
    if (x->type == MV_UNEVAL) return mv_copy(x);
    MValue *y = mv_num(x->rows, x->cols);
    for (int i = 0; i < numel(x); i++) y->re[i] = f(x->re[i]);
    (void)p;
    return y;
}

static double d_sind(double x) { double r = fmod(x, 360.0); return (r == 0 || r == 180 || r == -180) ? 0.0 : sin(x / 180.0 * PI); }
static double d_cosd(double x) { double r = fmod(fabs(x), 360.0); return (r == 90 || r == 270) ? 0.0 : cos(x / 180.0 * PI); }
static double d_tand(double x) { return d_sind(x) / d_cosd(x); }
static double d_asind(double x) { return asin(x) * 180.0 / PI; }
static double d_acosd(double x) { return acos(x) * 180.0 / PI; }
static double d_atand(double x) { return atan(x) * 180.0 / PI; }
static double d_deg2rad(double x) { return PI / 180.0 * x; }
static double d_rad2deg(double x) { return 180.0 / PI * x; }
static double d_round(double x) { return round(x); }
static double d_fix(double x) { return trunc(x); }
static double d_logical(double x) { return x != 0; }
static double d_not(double x) { return x == 0; }
static double d_neg(double x) { return -x; }

/* Reduction over a vector (scalar result) or over the columns of a matrix. */
static MValue *reduce(Parser *p, MValue *x, int op)
{
    if (x->type == MV_UNEVAL) return mv_copy(x);
    const int vec = x->rows == 1 || x->cols == 1;
    const int len = vec ? numel(x) : x->rows, nOut = vec ? 1 : x->cols;
    if (numel(x) == 0) return op == 1 ? mv_scalar(1) : op <= 2 ? mv_scalar(op == 0 ? 0 : 1) : mv_num(0, 0);
    MValue *y = mv_num(1, nOut);
    for (int j = 0; j < nOut; j++) {
        const double *v = x->re + (size_t)j * len;
        double r = op == 1 ? 1.0 : op >= 3 ? NAN : 0.0;
        for (int i = 0; i < len; i++) {
            switch (op) {
            case 0: case 2: r += v[i]; break;          /* sum, mean */
            case 1: r *= v[i]; break;                  /* prod */
            case 3: r = f_min(r, v[i]); break;
            case 4: r = f_max(r, v[i]); break;
            }
        }
        y->re[j] = op == 2 ? r / len : r;
    }
    (void)p;
    return y;
}

static void strip_ext(char *s)
{
    char *dot = strrchr(s, '.');
    char *sep1 = strrchr(s, '/'), *sep2 = strrchr(s, '\\');
    char *sep = sep1 > sep2 ? sep1 : sep2;
    if (dot && (!sep || dot > sep)) *dot = '\0';
}

static MValue *str_arg(Parser *p, MValue **a, int n, int i, const char *f)
{
    if (i >= n) fail(p, "%s: not enough arguments", f);
    if (a[i]->type != MV_STR && a[i]->type != MV_UNEVAL) fail(p, "%s: argument %d must be a string", f, i + 1);
    return a[i];
}

static MValue *call(Parser *p, const char *f, MValue **a, int n)
{
    /* constants and constructors */
    if (!strcmp(f, "pi")) return mv_scalar(PI);
    if (!strcmp(f, "eps")) return mv_scalar(2.220446049250313e-16);
    if (!strcmp(f, "true") || !strcmp(f, "false") || !strcmp(f, "Inf") || !strcmp(f, "inf") || !strcmp(f, "NaN") ||
        !strcmp(f, "nan") || !strcmp(f, "zeros") || !strcmp(f, "ones")) {
        const double x = !strcmp(f, "true") || !strcmp(f, "ones") ? 1.0 : !strcmp(f, "false") || !strcmp(f, "zeros") ? 0.0
                       : (f[0] == 'I' || f[0] == 'i') ? INFINITY : NAN;
        int r = 1, c = 1;
        if (n >= 1) { MValue *d = arg_num(p, a, n, 0, f); if (d->type == MV_UNEVAL) return mv_copy(d); r = c = (int)d->re[0]; if (numel(d) >= 2) c = (int)d->re[1]; }
        if (n >= 2) { MValue *d = arg_num(p, a, n, 1, f); if (d->type == MV_UNEVAL) return mv_copy(d); c = (int)d->re[0]; }
        if (r < 0) r = 0;
        if (c < 0) c = 0;
        MValue *y = mv_num(r, c);
        for (int i = 0; i < r * c; i++) y->re[i] = x;
        return y;
    }
    for (int i = 0; i < n; i++)                           /* an unevaluated argument: so is the result */
        if (a[i]->type == MV_UNEVAL && strcmp(f, "struct") && strcmp(f, "isempty")) return mv_copy(a[i]);
    /* files and strings */
    if (!strcmp(f, "mfilename")) {
        char *s = xstrdup(p->filePath);
        strip_ext(s);
        if (!(n >= 1 && a[0]->type == MV_STR && !strcmp(a[0]->str, "fullpath"))) {
            char *sep1 = strrchr(s, '/'), *sep2 = strrchr(s, '\\'), *sep = sep1 > sep2 ? sep1 : sep2;
            if (sep) memmove(s, sep + 1, strlen(sep + 1) + 1);
        }
        MValue *y = mv_str(s);
        free(s);
        return y;
    }
    if (!strcmp(f, "fileparts")) {
        char *s = xstrdup(str_arg(p, a, n, 0, f)->str);
        char *sep1 = strrchr(s, '/'), *sep2 = strrchr(s, '\\'), *sep = sep1 > sep2 ? sep1 : sep2;
        if (sep) *sep = '\0'; else s[0] = '\0';
        MValue *y = mv_str(s);
        free(s);
        return y;
    }
    if (!strcmp(f, "fullfile")) {
        size_t len = 1;
        for (int i = 0; i < n; i++) len += strlen(str_arg(p, a, n, i, f)->str) + 1;
        char *s = xmalloc(len);
        for (int i = 0; i < n; i++) {
            const char *part = a[i]->str;
            if (!*part) continue;
            const size_t l = strlen(s);
            if (l > 0 && s[l - 1] != '/' && s[l - 1] != '\\') strcat(s, "/");
            strcat(s, part);
        }
        MValue *y = mv_str(s);
        free(s);
        return y;
    }
    if (!strcmp(f, "lower") || !strcmp(f, "upper")) {
        MValue *y = mv_copy(str_arg(p, a, n, 0, f));
        for (char *c = y->str; *c; c++) *c = (char)(f[0] == 'l' ? tolower((unsigned char)*c) : toupper((unsigned char)*c));
        return y;
    }
    if (!strcmp(f, "strcmp") || !strcmp(f, "strcmpi")) {
        if (n < 2) fail(p, "%s: two arguments expected", f);
        if (a[0]->type != MV_STR || a[1]->type != MV_STR) return mv_scalar(0);
        return mv_scalar(f[6] ? !strcasecmp_portable(a[0]->str, a[1]->str) : !strcmp(a[0]->str, a[1]->str));
    }
    if (!strcmp(f, "struct")) {
        if (n % 2) fail(p, "struct: name/value pairs expected");
        MValue *s = mv_new(MV_STRUCT);
        s->rows = s->cols = 1;
        for (int i = 0; i < n; i += 2) {
            if (a[i]->type != MV_STR) { mvalue_free(s); fail(p, "struct: field names must be strings"); }
            MValue *v = a[i + 1];
            MValue **slot = struct_slot(s, a[i]->str, 1);
            if (v->type == MV_CELL && numel(v) == 1) *slot = mv_copy(v->cells[0]);
            else if (v->type == MV_CELL) *slot = mv_uneval("struct arrays (struct with cell values) are not supported");
            else *slot = mv_copy(v);
        }
        return s;
    }
    if (!strcmp(f, "isempty")) {
        if (n < 1) fail(p, "isempty: one argument expected");
        if (a[0]->type == MV_UNEVAL) return mv_copy(a[0]);
        return mv_scalar(a[0]->type == MV_STRUCT ? 0 : numel(a[0]) == 0);
    }
    if (!strcmp(f, "numel") || !strcmp(f, "length")) {
        if (n < 1) fail(p, "%s: one argument expected", f);
        if (a[0]->type == MV_STRUCT) return mv_scalar(1);
        return mv_scalar(f[0] == 'n' ? numel(a[0]) : numel(a[0]) ? (a[0]->rows > a[0]->cols ? a[0]->rows : a[0]->cols) : 0);
    }
    if (!strcmp(f, "double")) return mv_copy(arg_num(p, a, n, 0, f));
    /* element-wise functions of one argument */
    static const struct { const char *name; double (*f)(double); } fn1[] = {
        { "sqrt", sqrt }, { "abs", fabs }, { "exp", exp }, { "log", log }, { "log10", log10 }, { "log2", log2 },
        { "sin", sin }, { "cos", cos }, { "tan", tan }, { "asin", asin }, { "acos", acos }, { "atan", atan },
        { "sind", d_sind }, { "cosd", d_cosd }, { "tand", d_tand }, { "asind", d_asind }, { "acosd", d_acosd },
        { "atand", d_atand }, { "deg2rad", d_deg2rad }, { "rad2deg", d_rad2deg }, { "round", d_round },
        { "floor", floor }, { "ceil", ceil }, { "fix", d_fix }, { "logical", d_logical } };
    for (size_t k = 0; k < sizeof fn1 / sizeof fn1[0]; k++)
        if (!strcmp(f, fn1[k].name)) return map1(p, arg_num(p, a, n, 0, f), fn1[k].f);
    if (!strcmp(f, "mod") || !strcmp(f, "rem") || ((!strcmp(f, "min") || !strcmp(f, "max")) && n >= 2)) {
        MValue *x = mv_copy(arg_num(p, a, n, 0, f)), *y = mv_copy(arg_num(p, a, n, 1, f));
        return elementwise(p, x, y, !strcmp(f, "mod") ? f_mod : !strcmp(f, "rem") ? f_rem : !strcmp(f, "min") ? f_min : f_max);
    }
    if (!strcmp(f, "sum") || !strcmp(f, "prod") || !strcmp(f, "mean") || !strcmp(f, "min") || !strcmp(f, "max"))
        return reduce(p, arg_num(p, a, n, 0, f), !strcmp(f, "sum") ? 0 : !strcmp(f, "prod") ? 1 : !strcmp(f, "mean") ? 2
                                                  : !strcmp(f, "min") ? 3 : 4);
    if (!strcmp(f, "linspace")) {
        const double x0 = arg_num(p, a, n, 0, f)->re[0], x1 = arg_num(p, a, n, 1, f)->re[0];
        const int k = n >= 3 ? (int)floor(arg_num(p, a, n, 2, f)->re[0]) : 100;
        MValue *y = mv_num(1, k > 0 ? k : 0);
        for (int i = 0; i < k; i++) y->re[i] = k == 1 ? x1 : x0 + i * (x1 - x0) / (k - 1);
        if (k > 1) y->re[k - 1] = x1;
        return y;
    }
    if (!strcmp(f, "addpath") || !strcmp(f, "rmpath") || !strcmp(f, "warning") || !strcmp(f, "disp") ||
        !strcmp(f, "fprintf") || !strcmp(f, "rng"))
        return mv_uneval("(no value)");
    char why[160];
    snprintf(why, sizeof why, "%s(...) is not evaluated by radar_export_config (MATLAB function)", f);
    return mv_uneval(why);
}

/* ------------------------------------------------------------------------ */
/* Parser                                                                   */
/* ------------------------------------------------------------------------ */

static MValue *expr(Parser *p);

static int in_matrix(const Parser *p) { return p->depth > 0 && p->ctx[p->depth - 1]; }
static void push_ctx(Parser *p, int matrix)
{
    if (p->depth >= 64) fail(p, "nested too deeply");
    p->ctx[p->depth++] = matrix;
}

/* Arguments of a call or an index, up to the closing ')' or '}'. */
static int parse_args(Parser *p, MValue **args, int max, const char *close)
{
    int n = 0;
    push_ctx(p, 0);
    while (cur(p)->type == T_NL) p->pos++;
    if (!accept(p, close)) {
        for (;;) {
            if (n == max) fail(p, "too many arguments");
            if (is_op(cur(p), ":") && (is_op(&p->tok[p->pos + 1], ",") || is_op(&p->tok[p->pos + 1], close))) {
                p->pos++;
                args[n++] = mv_str(":");                  /* x(:) */
            } else {
                args[n++] = expr(p);
            }
            if (accept(p, close)) break;
            expect(p, ",");
        }
    }
    p->depth--;
    return n;
}

/* x(i), x(i,j), c{i}: 1-based indices (no 'end'). */
static MValue *index_value(Parser *p, MValue *x, MValue **ix, int n, int brace)
{
    if (x->type == MV_UNEVAL) return mv_copy(x);
    if (x->type == MV_STRUCT) fail(p, "struct arrays cannot be indexed");
    if (n < 1 || n > 2) fail(p, "one or two indices expected");
    const int dims[2] = { n == 1 ? numel(x) : x->rows, n == 1 ? 1 : x->cols };
    int cnt[2] = { 1, 1 };
    int *idx[2] = { NULL, NULL };
    for (int d = 0; d < n; d++) {
        if (ix[d]->type == MV_STR && !strcmp(ix[d]->str, ":")) {
            cnt[d] = dims[d];
            idx[d] = xmalloc(sizeof(int) * (size_t)(cnt[d] ? cnt[d] : 1));
            for (int i = 0; i < cnt[d]; i++) idx[d][i] = i;
            continue;
        }
        ix[d] = to_num(p, ix[d]);
        if (ix[d]->type == MV_UNEVAL) { free(idx[0]); return mv_copy(ix[d]); }
        cnt[d] = numel(ix[d]);
        idx[d] = xmalloc(sizeof(int) * (size_t)(cnt[d] ? cnt[d] : 1));
        for (int i = 0; i < cnt[d]; i++) {
            const double v = ix[d]->re[i];
            if (v != floor(v) || v < 1 || v > dims[d]) { free(idx[0]); free(idx[1]); fail(p, "index %g out of bounds (1 ... %d)", v, dims[d]); }
            idx[d][i] = (int)v - 1;
        }
    }
    MValue *y;
    const int linear = n == 1;
    const int rows = linear ? (x->rows == 1 ? 1 : cnt[0]) : cnt[0], cols = linear ? (x->rows == 1 ? cnt[0] : 1) : cnt[1];
    if (x->type == MV_CELL) {
        if (brace) {
            if (cnt[0] * cnt[1] != 1) { free(idx[0]); free(idx[1]); fail(p, "c{...} must select one cell"); }
            y = mv_copy(x->cells[linear ? idx[0][0] : idx[0][0] + idx[1][0] * x->rows]);
        } else {
            y = mv_new(MV_CELL);
            y->rows = rows;
            y->cols = cols;
            y->cells = xmalloc(sizeof(MValue *) * (size_t)(rows * cols > 0 ? rows * cols : 1));
            for (int j = 0; j < cnt[1]; j++)
                for (int i = 0; i < cnt[0]; i++)
                    y->cells[i + j * cnt[0]] = mv_copy(x->cells[linear ? idx[0][i] : idx[0][i] + idx[1][j] * x->rows]);
        }
    } else {
        if (brace) { free(idx[0]); free(idx[1]); fail(p, "{ } indexing of a non-cell"); }
        MValue *src = x->type == MV_STR ? to_num(p, mv_copy(x)) : x;
        y = mv_num(rows, cols);
        for (int j = 0; j < cnt[1]; j++)
            for (int i = 0; i < cnt[0]; i++)
                y->re[i + j * cnt[0]] = src->re[linear ? idx[0][i] : idx[0][i] + idx[1][j] * x->rows];
        if (x->type == MV_STR) {
            char *s = xmalloc((size_t)numel(y) + 1);
            for (int i = 0; i < numel(y); i++) s[i] = (char)y->re[i];
            mvalue_free(y);
            y = mv_str(s);
            free(s);
            mvalue_free(src);
        }
    }
    free(idx[0]);
    free(idx[1]);
    return y;
}

static MValue *matrix_literal(Parser *p, const char *close, int cell)
{
    MValue **rows = NULL, **elems = NULL;
    int nRows = 0, nElems = 0, capR = 0, capE = 0;
    push_ctx(p, 1);
    for (;;) {
        while (cur(p)->type == T_NL || is_op(cur(p), ";")) {
            if (nElems > 0) {                            /* end of a row */
                if (nRows == capR) { capR = capR ? 2 * capR : 8; rows = realloc(rows, sizeof(MValue *) * (size_t)capR); }
                rows[nRows++] = concat(p, elems, nElems, 2, cell);
                nElems = 0;
            }
            p->pos++;
        }
        if (accept(p, close)) break;
        if (cur(p)->type == T_EOF) fail(p, "'%s' expected", close);
        MValue *e = expr(p);
        if (cell) {                                      /* every element is a cell */
            MValue *c = mv_new(MV_CELL);
            c->rows = c->cols = 1;
            c->cells = xmalloc(sizeof(MValue *));
            c->cells[0] = e;
            e = c;
        }
        if (nElems == capE) { capE = capE ? 2 * capE : 8; elems = realloc(elems, sizeof(MValue *) * (size_t)capE); }
        elems[nElems++] = e;
        if (accept(p, ",")) continue;
        const Token *t = cur(p);
        if (t->type == T_NL || is_op(t, ";") || is_op(t, close)) continue;
        if (!t->space) fail(p, "unexpected '%s' in [ ] or { }", t->type == T_OP ? t->text : t->s ? t->s : "token");
    }
    if (nElems > 0) {
        if (nRows == capR) { capR = capR ? 2 * capR : 8; rows = realloc(rows, sizeof(MValue *) * (size_t)capR); }
        rows[nRows++] = concat(p, elems, nElems, 2, cell);
    }
    p->depth--;
    MValue *y = concat(p, rows, nRows, 1, cell);
    free(rows);
    free(elems);
    return y;
}

static MValue *primary(Parser *p)
{
    Token *t = cur(p);
    MValue *v;
    if (t->type == T_NUM) { p->pos++; v = mv_scalar(t->num); }
    else if (t->type == T_STR) { p->pos++; v = mv_str(t->s); }
    else if (accept(p, "(")) {
        push_ctx(p, 0);
        v = expr(p);
        expect(p, ")");
        p->depth--;
    }
    else if (accept(p, "[")) v = matrix_literal(p, "]", 0);
    else if (accept(p, "{")) v = matrix_literal(p, "}", 1);
    else if (t->type == T_ID) {
        p->pos++;
        MValue *var = get_var(p, t->s);
        if (var) {
            v = mv_copy(var);
        } else {
            MValue *args[32];
            int n = 0;
            if (is_op(cur(p), "(") && !(in_matrix(p) && cur(p)->space)) {
                p->pos++;
                n = parse_args(p, args, 32, ")");
            }
            v = call(p, t->s, args, n);
            for (int i = 0; i < n; i++) mvalue_free(args[i]);
        }
    }
    else fail(p, "unexpected %s", t->type == T_NL || t->type == T_EOF ? "end of the statement" : t->type == T_OP ? t->text : "token");
    /* postfix: .field, (index), {index} */
    for (;;) {
        Token *q = cur(p);
        if (in_matrix(p) && q->space) break;
        if (is_op(q, ".") && p->tok[p->pos + 1].type == T_ID) {
            p->pos++;
            const char *name = cur(p)->s;
            p->pos++;
            if (v->type == MV_UNEVAL) continue;
            if (v->type != MV_STRUCT) fail(p, "field '%s' of a value that is not a struct", name);
            MValue **slot = struct_slot(v, name, 0);
            if (!slot) fail(p, "reference to the non-existent field '%s'", name);
            MValue *f = mv_copy(*slot);
            mvalue_free(v);
            v = f;
        } else if (is_op(q, "(") || is_op(q, "{")) {
            const int brace = is_op(q, "{");
            MValue *args[2] = { NULL, NULL };
            p->pos++;
            const int n = parse_args(p, args, 2, brace ? "}" : ")");
            MValue *y = index_value(p, v, args, n, brace);
            for (int i = 0; i < n; i++) mvalue_free(args[i]);
            mvalue_free(v);
            v = y;
        } else {
            break;
        }
    }
    return v;
}

static MValue *postfix(Parser *p)
{
    MValue *v = primary(p);
    while ((is_op(cur(p), "'") || is_op(cur(p), ".'")) && !(in_matrix(p) && cur(p)->space)) {
        p->pos++;
        v = transpose(p, v);
    }
    return v;
}

static MValue *unary(Parser *p);

static MValue *power(Parser *p)
{
    MValue *v = postfix(p);
    while ((is_op(cur(p), "^") || is_op(cur(p), ".^")) && !(in_matrix(p) && cur(p)->space && !p->tok[p->pos + 1].space)) {
        const int dot = is_op(cur(p), ".^");
        p->pos++;
        MValue *e;
        if (is_op(cur(p), "-") || is_op(cur(p), "+") || is_op(cur(p), "~")) e = unary(p);   /* 2^-1 */
        else e = postfix(p);
        if (!dot && v->type == MV_NUM && e->type == MV_NUM && (numel(v) != 1 || numel(e) != 1)) {
            mvalue_free(v);
            mvalue_free(e);
            fail(p, "matrix power: use .^");
        }
        v = elementwise(p, v, e, f_pow);
    }
    return v;
}

static MValue *unary(Parser *p)
{
    if (accept(p, "-")) { MValue *v = unary(p); if (v->type == MV_UNEVAL) return v; v = to_num(p, v); MValue *y = map1(p, v, d_neg); mvalue_free(v); return y; }
    if (accept(p, "+")) return to_num(p, unary(p));
    if (accept(p, "~")) { MValue *v = unary(p); if (v->type == MV_UNEVAL) return v; v = to_num(p, v); MValue *y = map1(p, v, d_not); mvalue_free(v); return y; }
    return power(p);
}

static MValue *mul(Parser *p)
{
    MValue *v = unary(p);
    for (;;) {
        Token *t = cur(p);
        if (in_matrix(p) && t->space && t->type == T_OP && !p->tok[p->pos + 1].space &&
            (is_op(t, "*") || is_op(t, "/"))) break;
        if (accept(p, "*")) v = matmul(p, v, unary(p));
        else if (accept(p, ".*")) v = elementwise(p, v, unary(p), f_mul);
        else if (accept(p, "./")) v = elementwise(p, v, unary(p), f_div);
        else if (accept(p, ".\\")) v = elementwise(p, v, unary(p), f_ldiv);
        else if (accept(p, "/")) {
            MValue *b = unary(p);
            if (b->type == MV_NUM && numel(b) != 1) { mvalue_free(v); mvalue_free(b); fail(p, "matrix division: use ./"); }
            v = elementwise(p, v, b, f_div);
        } else if (accept(p, "\\")) {
            MValue *b = unary(p);
            if (v->type == MV_NUM && numel(v) != 1) { mvalue_free(v); mvalue_free(b); fail(p, "matrix division: use .\\"); }
            v = elementwise(p, v, b, f_ldiv);
        } else break;
    }
    return v;
}

static MValue *add(Parser *p)
{
    MValue *v = mul(p);
    for (;;) {
        Token *t = cur(p);
        if (!is_op(t, "+") && !is_op(t, "-")) break;
        if (in_matrix(p) && t->space && !p->tok[p->pos + 1].space) break;   /* [1 -2]: a new element */
        p->pos++;
        MValue *b = mul(p);
        v = elementwise(p, v, b, t->text[0] == '+' ? f_add : f_sub);
    }
    return v;
}

static MValue *range(Parser *p)
{
    MValue *a = add(p);
    if (!is_op(cur(p), ":")) return a;
    p->pos++;
    MValue *b = add(p), *c = NULL;
    if (accept(p, ":")) c = add(p);
    MValue *lo = a, *step = c ? b : NULL, *hi = c ? c : b;
    MValue *vals[3] = { lo, step, hi };
    for (int i = 0; i < 3; i++)
        if (vals[i] && vals[i]->type == MV_UNEVAL) {
            MValue *u = mv_copy(vals[i]);
            mvalue_free(a); mvalue_free(b); mvalue_free(c);
            return u;
        }
    lo = to_num(p, lo);
    hi = to_num(p, hi);
    if (step) step = to_num(p, step);
    const double x0 = lo->re[0], x1 = hi->re[0], d = step ? step->re[0] : 1.0;
    mvalue_free(lo); mvalue_free(hi); mvalue_free(step);
    if (d == 0 || !isfinite(x0) || !isfinite(x1)) return mv_num(1, 0);
    const double nf = floor((x1 - x0) / d + 1e-10);
    const int n = nf < 0 ? 0 : (int)nf + 1;
    MValue *y = mv_num(1, n);
    for (int i = 0; i < n; i++) y->re[i] = x0 + i * d;
    return y;
}

static MValue *compare(Parser *p)
{
    MValue *v = range(p);
    static const char *const ops[] = { "==", "~=", "<", ">", "<=", ">=" };
    static double (*const fn[])(double, double) = { f_eq, f_ne, f_lt, f_gt, f_le, f_ge };
    for (;;) {
        int k = -1;
        for (int i = 0; i < 6; i++) if (is_op(cur(p), ops[i])) k = i;
        if (k < 0) break;
        if (in_matrix(p) && cur(p)->space && !p->tok[p->pos + 1].space && k >= 2) break;
        p->pos++;
        if (k <= 1 && v->type == MV_STR) {                /* 'abc' == 'abc' as in MATLAB: per character */
            MValue *b = range(p);
            v = elementwise(p, v, b, fn[k]);
            continue;
        }
        v = elementwise(p, v, range(p), fn[k]);
    }
    return v;
}

static MValue *and_(Parser *p)
{
    MValue *v = compare(p);
    while (is_op(cur(p), "&")) { p->pos++; v = elementwise(p, v, compare(p), f_and); }
    return v;
}

static MValue *or_(Parser *p)
{
    MValue *v = and_(p);
    while (is_op(cur(p), "|")) { p->pos++; v = elementwise(p, v, and_(p), f_or); }
    return v;
}

static MValue *andand(Parser *p)
{
    MValue *v = or_(p);
    while (accept(p, "&&")) v = elementwise(p, v, or_(p), f_and);
    return v;
}

static MValue *expr(Parser *p)
{
    MValue *v = andand(p);
    while (accept(p, "||")) v = elementwise(p, v, andand(p), f_or);
    return v;
}

static void skip_statement(Parser *p)
{
    while (cur(p)->type != T_NL && cur(p)->type != T_EOF && !is_op(cur(p), ";")) p->pos++;
}

/* Runs the statements of the token list. stopAtFunction: a second 'function'
 * line (a subfunction) ends the script. */
static void run(Parser *p)
{
    int functions = 0;
    for (;;) {
        Token *t = cur(p);
        if (t->type == T_EOF) return;
        if (t->type == T_NL || is_op(t, ";") || is_op(t, ",")) { p->pos++; continue; }
        if (t->type == T_ID && !strcmp(t->s, "function")) {
            if (functions++ > 0) return;                  /* a subfunction: the main function ended */
            p->pos++;
            /* function out = name(...)  or  function name(...) */
            if (cur(p)->type == T_ID && is_op(&p->tok[p->pos + 1], "=")) {
                snprintf(p->outName, sizeof p->outName, "%s", cur(p)->s);
            } else if (is_op(cur(p), "[")) {
                p->pos++;
                if (cur(p)->type == T_ID) snprintf(p->outName, sizeof p->outName, "%s", cur(p)->s);
            }
            skip_statement(p);
            continue;
        }
        if (t->type == T_ID && !strcmp(t->s, "end")) { p->pos++; continue; }
        if (t->type == T_ID && (!strcmp(t->s, "if") || !strcmp(t->s, "for") || !strcmp(t->s, "while") ||
                                !strcmp(t->s, "switch") || !strcmp(t->s, "try") || !strcmp(t->s, "parfor") ||
                                !strcmp(t->s, "global") || !strcmp(t->s, "persistent")))
            fail(p, "'%s' is not supported by radar_export_config: only assignments (cfg.x = value;)", t->s);
        /* assignment: name(.field)* = expr */
        int k = p->pos, nPath = 0;
        const char *path[16];
        if (t->type == T_ID) {
            path[nPath++] = t->s;
            k++;
            while (is_op(&p->tok[k], ".") && p->tok[k + 1].type == T_ID && nPath < 16) {
                path[nPath++] = p->tok[k + 1].s;
                k += 2;
            }
        }
        if (nPath > 0 && is_op(&p->tok[k], "=")) {
            p->pos = k + 1;
            MValue *v = expr(p);
            MValue **slot = var_slot(p, path[0]);
            for (int i = 1; i < nPath; i++) {
                if (!*slot) { *slot = mv_new(MV_STRUCT); (*slot)->rows = (*slot)->cols = 1; }
                if ((*slot)->type != MV_STRUCT) { mvalue_free(v); fail(p, "'%s' is not a struct", path[i - 1]); }
                slot = struct_slot(*slot, path[i], 1);
            }
            mvalue_free(*slot);
            *slot = v;
        } else if (nPath > 0 && (is_op(&p->tok[k], "(") || is_op(&p->tok[k], "{")) && get_var(p, path[0])) {
            /* x(i) = ... ? */
            int d = 0, j = k;
            do {
                if (is_op(&p->tok[j], "(") || is_op(&p->tok[j], "{") || is_op(&p->tok[j], "[")) d++;
                if (is_op(&p->tok[j], ")") || is_op(&p->tok[j], "}") || is_op(&p->tok[j], "]")) d--;
                j++;
            } while (d > 0 && p->tok[j].type != T_EOF);
            if (is_op(&p->tok[j], "=")) fail(p, "indexed assignment (x(i) = ...) is not supported: assign the whole value");
            mvalue_free(expr(p));
        } else if (is_op(t, "[")) {
            fail(p, "multiple outputs ([a, b] = ...) are not supported");
        } else {
            mvalue_free(expr(p));                         /* e.g. addpath(...) */
        }
        const Token *e = cur(p);
        if (e->type != T_NL && e->type != T_EOF && !is_op(e, ";") && !is_op(e, ","))
            fail(p, "unexpected '%s' after the statement", e->type == T_OP ? e->text : e->s ? e->s : "token");
    }
}

static void free_tokens(Parser *p)
{
    for (int i = 0; i < p->n; i++) free(p->tok[i].s);
    free(p->tok);
    p->tok = NULL;
    p->n = p->cap = p->pos = 0;
}

MValue *mfile_eval(const char *path, const char *const *extra, int nExtra, char *err, size_t errLen)
{
    Parser p;
    MValue *volatile result = NULL;
    memset(&p, 0, sizeof p);
    p.err = err;
    p.errLen = errLen;
    p.filePath = path;
    strcpy(p.outName, "cfg");
    if (err && errLen) err[0] = '\0';

    FILE *f = fopen(path, "rb");
    if (!f) {
        if (err && errLen) snprintf(err, errLen, "cannot open %s", path);
        return NULL;
    }
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *text = xmalloc((size_t)(size > 0 ? size : 0) + 1);
    if (size > 0 && fread(text, 1, (size_t)size, f) != (size_t)size) {
        fclose(f);
        free(text);
        if (err && errLen) snprintf(err, errLen, "cannot read %s", path);
        return NULL;
    }
    fclose(f);
    text[size > 0 ? size : 0] = '\0';

    if (setjmp(p.jb) == 0) {
        p.source = path;
        p.lineBase = 1;
        lex(&p, text);
        run(&p);
        free_tokens(&p);
        for (int i = 0; i < nExtra; i++) {
            char src[32];
            snprintf(src, sizeof src, "--set %d", i + 1);
            p.source = src;
            p.lineBase = 1;
            p.depth = 0;
            lex(&p, extra[i]);
            run(&p);
            free_tokens(&p);
            p.source = path;
        }
        MValue **slot = var_slot(&p, p.outName);
        if (!*slot) fail(&p, "the function does not set its output '%s'", p.outName);
        result = *slot;
        *slot = NULL;
    }
    free_tokens(&p);
    free(text);
    for (int i = 0; i < p.nVars; i++) {
        free(p.varNames[i]);
        mvalue_free(p.vars[i]);
    }
    free(p.varNames);
    free(p.vars);
    return result;
}

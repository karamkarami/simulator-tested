/*
 * test_mfile.c - the MATLAB subset of radar_export_config (tools/mfile.c):
 * syntax of radar_config.m (function line, comments, continuation lines,
 * matrices with spaces and unary minus, multi-line cells, struct(), fullfile,
 * fileparts(mfilename('fullpath'))), operators and precedence, --set
 * statements, unevaluated MATLAB functions and the error messages.
 */
#include "mfile.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
#define CHECK(c) do { if (!(c)) { printf("FAIL line %d: %s\n", __LINE__, #c); failures++; } } while (0)

static double scalar(const MValue *cfg, const char *path)
{
    const MValue *v = mvalue_get(cfg, path);
    return v && v->type == MV_NUM && v->rows * v->cols >= 1 ? v->re[0] : NAN;
}

static const char *text(const MValue *cfg, const char *path)
{
    const MValue *v = mvalue_get(cfg, path);
    return v && v->type == MV_STR ? v->str : "";
}

static void write_file(const char *path, const char *s)
{
    FILE *f = fopen(path, "wb");
    fputs(s, f);
    fclose(f);
}

int main(void)
{
    const char *path = "test_mfile_config.m";
    write_file(path,
        "function cfg = test_mfile_config()\n"
        "%TEST  help text\n"
        "root = fileparts(mfilename('fullpath'));\n"
        "addpath(fullfile(root, 'src'));\n"
        "%{\n"
        "cfg.blocked = 1;\n"
        "%}\n"
        "cfg.files.a = fullfile(root, 'inputs', 'a.mat');   % comment\n"
        "cfg.m = [    0 40;  1000 40; ...\n"
        "         2000 -32];\n"
        "cfg.row = [1 -2 1];  cfg.one = [1 - 2 1]; cfg.empty = [];\n"
        "cfg.p = -2^2 + 3*2^-1;          % -4 + 1.5\n"
        "cfg.c = 299792458;\n"
        "cfg.r = 2*cfg.c/cfg.c;\n"
        "cfg.cmp = (3 > 2) && ~(1 == 2) || false;\n"
        "cfg.s = 'it''s'; cfg.d = \"dq\";\n"
        "cfg.t = struct('a', 1, 'b', [], 'c', 'x');\n"
        "cfg.cells = {\n"
        "    'mountain', [ 30  75], [42e3 62e3], 2500   % first\n"
        "    'hill',     [100 160], [ 3e3 12e3],  200\n"
        "};\n"
        "cfg.col = (1:3)';\n"
        "cfg.idx = cfg.col(2) + cfg.m(2, 1);\n"
        "cfg.fn = round(sqrt(2)*100)/100 + floor(-0.5) + mod(-1, 3) + max([3 7 5]) + numel(1:4);\n"
        "cfg.k = Inf; cfg.n = NaN; cfg.lam = 1e3 * .5e-3;\n"
        "cfg.targets = target_scenario(cfg);\n"
        "end\n"
        "function x = helper()\n"
        "x = 1;\n"
        "end\n");
    const char *sets[] = { "cfg.set1 = 'staggered';", "cfg.c = 3e8;" };
    char err[512];
    MValue *cfg = mfile_eval(path, sets, 2, err, sizeof err);
    CHECK(cfg != NULL);
    if (!cfg) { printf("%s\n", err); return 1; }
    CHECK(!mvalue_get(cfg, "blocked"));
    CHECK(strstr(text(cfg, "files.a"), "inputs/a.mat") != NULL);
    const MValue *m = mvalue_get(cfg, "m");
    CHECK(m && m->rows == 3 && m->cols == 2 && m->re[2] == 2000 && m->re[5] == -32);
    const MValue *row = mvalue_get(cfg, "row"), *one = mvalue_get(cfg, "one");
    CHECK(row && row->cols == 3 && row->re[1] == -2 && one && one->cols == 2 && one->re[0] == -1);
    CHECK(mvalue_get(cfg, "empty")->rows == 0);
    CHECK(scalar(cfg, "p") == -2.5 && scalar(cfg, "r") == 2 && scalar(cfg, "cmp") == 1);
    CHECK(!strcmp(text(cfg, "s"), "it's") && !strcmp(text(cfg, "d"), "dq"));
    CHECK(scalar(cfg, "t.a") == 1 && mvalue_get(cfg, "t.b")->rows == 0 && !strcmp(text(cfg, "t.c"), "x"));
    const MValue *cells = mvalue_get(cfg, "cells");
    CHECK(cells && cells->type == MV_CELL && cells->rows == 2 && cells->cols == 4 &&
          !strcmp(cells->cells[1]->str, "hill") && cells->cells[7]->re[0] == 200);
    const MValue *col = mvalue_get(cfg, "col");
    CHECK(col && col->rows == 3 && col->cols == 1 && scalar(cfg, "idx") == 1002);
    CHECK(fabs(scalar(cfg, "fn") - (1.41 - 1 + 2 + 7 + 4)) < 1e-12);
    CHECK(isinf(scalar(cfg, "k")) && isnan(scalar(cfg, "n")) && fabs(scalar(cfg, "lam") - 0.5) < 1e-15);
    CHECK(mvalue_get(cfg, "targets")->type == MV_UNEVAL);
    CHECK(!strcmp(text(cfg, "set1"), "staggered") && scalar(cfg, "c") == 3e8);
    mvalue_free(cfg);

    /* errors name the line */
    write_file(path, "cfg.a = 1;\ncfg.b = [1 2; 3];\n");
    CHECK(mfile_eval(path, NULL, 0, err, sizeof err) == NULL && strstr(err, ":2:") != NULL);
    write_file(path, "cfg.a = 1;\nif cfg.a\n cfg.b = 2;\nend\n");
    CHECK(mfile_eval(path, NULL, 0, err, sizeof err) == NULL && strstr(err, "'if' is not supported") != NULL);
    write_file(path, "cfg.a = [1 2 3];\ncfg.a(2) = 5;\n");
    CHECK(mfile_eval(path, NULL, 0, err, sizeof err) == NULL && strstr(err, "indexed assignment") != NULL);
    CHECK(mfile_eval("does_not_exist.m", NULL, 0, err, sizeof err) == NULL && strstr(err, "cannot open") != NULL);
    remove(path);
    printf(failures ? "%d FAILURE(S)\n" : "PASS\n", failures);
    return failures ? 1 : 0;
}

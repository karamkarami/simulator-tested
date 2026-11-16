function info = make_phase_code_file(outFile, varargin)
%MAKE_PHASE_CODE_FILE  Phase-code file (cfg.files.phaseCode) from binary codes and mismatch filters.
%
%   MAKE_PHASE_CODE_FILE(outFile, 'shortCodes', codeFiles, 'shortMismatch', filterFiles, ...
%                                 'longCodes',  codeFiles, 'longMismatch',  filterFiles)
%   writes outFile (e.g. inputs/phase_code.mat) in the format of the phase-code file of a radar
%   mode (make_input_files.m, setup_radar.m):
%     shortCode_<k>, longCode_<k>          [nChips x 1] chip values +-1 (the first chip is sent first)
%     shortMismatch_<k>, longMismatch_<k>  [K x 1] mismatch filter of code k, one coefficient per chip
%   Pair k is the k-th code file with the k-th filter file (k = 1 ... n). Give the files of one or
%   both sub-pulses (a cell array of names, or one name for n = 1). The file has no variable
%   'madeFor', so make_input_files / setup_radar use it as it is and never rebuild it.
%
%   Code files (.xlsx / .xls, or .csv / .txt): one row of binary chips 0 / 1 (the first non-empty
%     row of the first sheet; a column works too). 0 -> -1 and 1 -> +1 (phases pi and 0); the
%     opposite mapping changes only the sign of the whole code, which the compression does not see.
%     A code already written as -1 / +1 is taken as it is.
%   Filter files (.mat): the mismatch coefficients as a row (or column) vector, real or complex.
%     The file may hold one numeric variable (any name); with several, name it with 'mismatchVar'.
%     Convention of the simulator: output = conv(code, filter), i.e. a convolution filter (for a
%     correlator h the filter is conj(flip(h))). The filter must not be shorter than its code.
%
%   Options (name, value):
%     'mismatchVar'     ''     name of the variable in the filter files ('' = the only numeric one)
%     'reverseMismatch' false  true: the files hold correlator weights; the filter is conj(flip(h))
%     'chips'           []     [nShort nLong]: expected chips of the codes (width_s/chip_s of the
%                              mode); [] = not checked. A different number stops with a message
%     'verbose'         true   print every pair: chips, filter length, peak sidelobe level (PSL)
%
%   info = MAKE_PHASE_CODE_FILE(...) returns for every sub-pulse the chips, filter lengths and the
%   PSL [dB] of every pair (and of the reversed filter, to check the convention).
%
%   In the mode file:  cfg.files.phaseCode = fullfile(root, 'inputs', 'phase_code_binary.mat');
%                      cfg.wave.<pulse>.type = 'phasecode'; compression = 'mismatch';
%                      cfg.wave.<pulse>.availableCodes = 1:n; codeIndex = one of them
%   (cfg.wave.<pulse>.mismatchLength is not used with a file of your own.)
%
%   Example:
%     make_phase_code_file(fullfile('inputs', 'phase_code_binary.mat'), ...
%         'shortCodes',    {'codes/short_1.xlsx', 'codes/short_2.xlsx'}, ...
%         'shortMismatch', {'codes/short_mis_1.mat', 'codes/short_mis_2.mat'}, ...
%         'longCodes',     {'codes/long_1.xlsx', 'codes/long_2.xlsx'}, ...
%         'longMismatch',  {'codes/long_mis_1.mat', 'codes/long_mis_2.mat'}, ...
%         'chips', [25 125]);

opt = struct('shortCodes', {{}}, 'shortMismatch', {{}}, 'longCodes', {{}}, 'longMismatch', {{}}, ...
             'mismatchVar', '', 'reverseMismatch', false, 'chips', [], 'verbose', true);
if mod(numel(varargin), 2) ~= 0, error('make_phase_code_file: options come in name / value pairs.'); end
for i = 1:2:numel(varargin)
    name = varargin{i};
    f = fieldnames(opt);
    j = find(strcmpi(f, name), 1);
    if isempty(j), error('make_phase_code_file: unknown option ''%s''.', name); end
    opt.(f{j}) = varargin{i+1};
end

names = {'short', 'long'};
out = struct();
info = struct();
for p = 1:2
    nm = names{p};
    cf = as_list(opt.([nm 'Codes']));
    mf = as_list(opt.([nm 'Mismatch']));
    if isempty(cf) && isempty(mf), continue; end
    if numel(cf) ~= numel(mf)
        error('make_phase_code_file: %d %s code files but %d %s mismatch files (one filter per code).', ...
              numel(cf), nm, numel(mf), nm);
    end
    n = numel(cf);
    s = struct('chips', zeros(1, n), 'filterLength', zeros(1, n), 'psl_dB', zeros(1, n), ...
               'pslReversed_dB', zeros(1, n));
    for k = 1:n
        code = read_code(cf{k});
        h = read_filter(mf{k}, opt.mismatchVar);
        if opt.reverseMismatch, h = conj(flipud(h)); end
        N = numel(code);
        if ~isempty(opt.chips) && N ~= opt.chips(p)
            error('make_phase_code_file: %s has %d chips, but the %s pulse needs %d (width_s/chip_s).', ...
                  cf{k}, N, nm, opt.chips(p));
        end
        if k > 1 && N ~= s.chips(1)
            error('make_phase_code_file: %s has %d chips, %s %d: all %s codes need the same length.', ...
                  cf{k}, N, cf{1}, s.chips(1), nm);
        end
        if numel(h) < N
            error('make_phase_code_file: the filter of %s has %d coefficients, fewer than the %d chips of %s.', ...
                  mf{k}, numel(h), N, cf{k});
        end
        s.chips(k) = N;
        s.filterLength(k) = numel(h);
        s.psl_dB(k) = psl(code, h);
        s.pslReversed_dB(k) = psl(code, conj(flipud(h)));
        out.(sprintf('%sCode_%d', nm, k)) = code;
        out.(sprintf('%sMismatch_%d', nm, k)) = h;
        if opt.verbose
            fprintf('%s pair %d: %d chips (%s), filter %d coefficients (%s): PSL %.1f dB\n', ...
                    nm, k, N, short_name(cf{k}), numel(h), short_name(mf{k}), s.psl_dB(k));
        end
        if s.pslReversed_dB(k) < s.psl_dB(k) - 3
            warning('make_phase_code_file:convention', ['%s pair %d: the reversed filter conj(flip(h)) ' ...
                    'gives a much lower PSL (%.1f dB instead of %.1f dB): the file seems to hold correlator ' ...
                    'weights; use ''reverseMismatch'', %s.'], nm, k, s.pslReversed_dB(k), s.psl_dB(k), ...
                    mat2str(~opt.reverseMismatch));
        end
    end
    info.(nm) = s;
end
if isempty(fieldnames(out))
    error('make_phase_code_file: no files given (shortCodes / shortMismatch, longCodes / longMismatch).');
end
folder = fileparts(outFile);
if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
save(outFile, '-struct', 'out', '-v7');
if opt.verbose, fprintf('Phase-code file written: %s\n', outFile); end
end

function c = as_list(x)
if isempty(x), c = {}; elseif ischar(x) || isstring(x), c = cellstr(x); else, c = x(:).'; end
end

function code = read_code(file)
% the first non-empty row (or a single column) of the first sheet: chips 0 / 1 (or -1 / +1)
if ~exist(file, 'file'), error('make_phase_code_file: code file %s not found.', file); end
[~, ~, ext] = fileparts(file);
switch lower(ext)
    case {'.xlsx', '.xls', '.xlsm'}
        if exist('readmatrix', 'file') || exist('readmatrix', 'builtin')
            M = readmatrix(file);
        else
            M = xlsread(file);
        end
    case {'.csv', '.txt'}
        M = dlmread(file);
    otherwise
        error('make_phase_code_file: %s: the code files must be .xlsx, .xls, .csv or .txt.', file);
end
M = double(M);
if isempty(M), error('make_phase_code_file: %s holds no numbers.', file); end
if size(M, 1) > 1 && size(M, 2) > 1                       % a sheet: its first row that holds numbers
    r = find(any(~isnan(M), 2), 1);
    M = M(r, :);
end
v = M(:);
v = v(~isnan(v));                                          % empty cells at the end of the row
if all(v == 0 | v == 1)
    code = 2*v - 1;                                        % 0 -> -1, 1 -> +1
elseif all(v == -1 | v == 1)
    code = v;
else
    bad = v(~(v == 0 | v == 1 | v == -1));
    error('make_phase_code_file: %s holds %g: a binary code has only 0 and 1 (or -1 and +1).', file, bad(1));
end
end

function h = read_filter(file, var)
if ~exist(file, 'file'), error('make_phase_code_file: mismatch file %s not found.', file); end
s = load(file);
f = fieldnames(s);
if isempty(var)
    num = f(cellfun(@(x) isnumeric(s.(x)) && numel(s.(x)) > 1, f));
    if numel(num) ~= 1
        error(['make_phase_code_file: %s holds %d numeric vectors (%s): name the one of the filter ' ...
               'with ''mismatchVar''.'], file, numel(num), strjoin(num.', ', '));
    end
    var = num{1};
elseif ~isfield(s, var)
    error('make_phase_code_file: %s has no variable ''%s'' (it holds: %s).', file, var, strjoin(f.', ', '));
end
h = s.(var);
if ~isvector(h) || ~isnumeric(h) || any(~isfinite(h(:)))
    error('make_phase_code_file: ''%s'' of %s must be a finite numeric vector.', var, file);
end
h = double(h(:));
end

function d = psl(code, h)
% peak sidelobe level of conv(code, h) relative to its main peak [dB] (one sample per chip)
y = abs(conv(code, h));
[pk, i] = max(y);
y(i) = 0;
d = 20*log10(max(max(y), eps*pk) / pk);
end

function s = short_name(f)
[~, a, b] = fileparts(f);
s = [a b];
end

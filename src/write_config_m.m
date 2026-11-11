function write_config_m(cfg, file)
%WRITE_CONFIG_M  Write a configuration struct as a MATLAB file of assignments (cfg.a.b = value;).
%   WRITE_CONFIG_M(cfg, file) writes every field of cfg (numbers, logicals,
%   strings, cell arrays of them, nested structs) so that the C programs of
%   the receiver project read it like radar_config.m: radar_export_config
%   --config-m FILE and radar_generate_ddc --config-m FILE use a cfg changed
%   in MATLAB (e.g. by a test or a script) without retyping the changes as
%   --set statements. cfg.targets (a struct array) is not written: give the
%   targets to radar_generate_ddc with --targets (export_targets_csv.m).
%   Example:
%     cfg = radar_config;  cfg.sim.sectors_deg = [20 45];
%     write_config_m(cfg, 'my_cfg.m');  export_targets_csv(cfg.targets, 'my_targets.csv');
%     % radar_generate_ddc --config-m my_cfg.m --targets my_targets.csv
fid = fopen(file, 'w');
if fid < 0, error('write_config_m: cannot write %s', file); end
c = onCleanup(@() fclose(fid));
fprintf(fid, '%% Configuration written by write_config_m.m\n');
emit(fid, 'cfg', cfg);
end

function emit(fid, path, v)
if isstruct(v) && isscalar(v)
    f = fieldnames(v);
    for i = 1:numel(f)
        if strcmp(path, 'cfg') && strcmp(f{i}, 'targets'), continue; end
        emit(fid, [path '.' f{i}], v.(f{i}));
    end
    if isempty(f), fprintf(fid, '%s = struct();\n', path); end
else
    fprintf(fid, '%s = %s;\n', path, value_text(v));
end
end

function s = value_text(v)
if ischar(v)
    s = ['''' strrep(v, '''', '''''') ''''];
elseif islogical(v) && isscalar(v)
    if v, s = 'true'; else, s = 'false'; end
elseif isnumeric(v) || islogical(v)
    v = double(v);
    if isempty(v)
        s = '[]';
    elseif isscalar(v)
        s = num_text(v);
    else
        rows = cell(1, size(v, 1));
        for r = 1:size(v, 1)
            rows{r} = strjoin(arrayfun(@num_text, v(r, :), 'UniformOutput', false), ' ');
        end
        s = ['[' strjoin(rows, '; ') ']'];
    end
elseif iscell(v)
    rows = cell(1, size(v, 1));
    for r = 1:size(v, 1)
        rows{r} = strjoin(cellfun(@value_text, v(r, :), 'UniformOutput', false), ', ');
    end
    s = ['{' strjoin(rows, '; ') '}'];
elseif isstruct(v)
    error('write_config_m: struct arrays and structs inside cells are not supported');
else
    error('write_config_m: unsupported value of class %s', class(v));
end
end

function s = num_text(x)
if isnan(x), s = 'NaN';
elseif isinf(x), if x > 0, s = 'Inf'; else, s = '-Inf'; end
elseif x == round(x) && abs(x) < 1e15, s = sprintf('%d', x);
else, s = sprintf('%.17g', x);
end
end

function ok = test_new_options()
%TEST_NEW_OPTIONS  Options of radar_config: sub-pulses on/off and on their own carriers,
%   CPI length auto/manual, the range limit of the clutter filter, the recording of the
%   matched filter and Doppler outputs, and the same receiver configuration from
%   MATLAB (export_receiver_config) and from the C program radar_export_config.
%
%   1) cfg.proc.cpiMode: 'auto' gives the beam dwell (114 pulses), 'manual' cfg.proc.cpiLength
%      (and its FFT length), an unknown mode stops with a message
%   2) cfg.mti.maxRange_m: a stationary echo before the limit is cancelled, one beyond it
%      stays at zero Doppler; the noise floor is flat on both sides (own noise gain)
%   3) cfg.wave.long.enable = false: the short-pulse window runs to cfg.rx.maxRange_m, a
%      target at 50 km gives a short-pulse plot, there is no long-pulse plot
%   4) cfg.wave.short.enable = false: a target at 12 km is not seen (blind range of the long
%      pulse), one at 45 km gives a long-pulse plot, no short-pulse plot
%   5) cfg.wave.short.fc_Hz = 165 MHz, cfg.wave.long.fc_Hz = 175 MHz: both targets (15 km,
%      +380 m/s; 50 km, -400 m/s) with the radial velocity of their own carrier (within one
%      Doppler bin; the wavelength of cfg.fc would be 3 - 4 bins off), the beams of the DBF
%      point to sind(el) * fc / fc_Hz
%   6) cfg.io.saveMf / saveDoppler: one mf_<n>.mat per PRI and one doppler_<n>.mat per CPI
%      update, with the sizes of the blocks ([] for a sub-pulse that is off)
%   7) the C program radar_export_config (when built) writes the same records as
%      export_receiver_config for every option above
%   8) the noise floor of the PPI video (colour scale) stays within 0.5 dB while the antenna
%      crosses a 45 deg sector whose rows are 95% clutter (the mountains)

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
base = radar_config();
base.ant.nElements = 1;
base.clutter.enable = false;
base.io.saveRaw = false;
base.display.enable = false;
base.console.enable = false;
base.track.enable = false;
base.io.outFolder = fullfile(root, 'output', 'test_new_options');

%% 1) CPI length
c = base;
s = setup_radar(c);
ok = check(s.cpi.length == 114 && strcmp(s.cpi.mode, 'auto'), sprintf('cpiMode auto: %d pulses', s.cpi.length)) && ok;
c.proc.cpiMode = 'manual';  c.proc.cpiLength = 64;
s = setup_radar(c);
ok = check(s.cpi.length == 64 && s.dop.nfft == 64, sprintf('cpiMode manual, cpiLength 64: %d pulses, FFT %d', ...
           s.cpi.length, s.dop.nfft)) && ok;
c.proc.cpiMode = 'fixed';
ok = check(~isempty(strfind(error_of(@() setup_radar(c)), 'cpiMode')), 'cpiMode ''fixed'' stops with a message') && ok;

%% 2) range limit of the clutter filter
c = base;
c.mti.maxRange_m = 50e3;
s = setup_radar(c);
wv = s.wave.long;
n = s.cpi.length;  nR = wv.nRange;
rng(3);
x = (randn(1, nR, n) + 1j*randn(1, nR, n)) / sqrt(2);
x(1, [100 nR-100], :) = x(1, [100 nR-100], :) + 100;        % stationary echoes (40 dB per pulse)
P = doppler_process(mti_filter(x, s.mtiCoeffs, ones(1, n), wv.mtiCells, s.mtiBypass), s, wv.mtiCells);
z = floor(s.dop.nfft/2) + 1;
out = abs((1:s.dop.nfft) - z) > s.dop.notchBins + 4;
nBefore = mean(mean(P(1, 200:wv.mtiCells-10, out)));
nBeyond = mean(mean(P(1, wv.mtiCells+10:nR-200, out)));
eBefore = 10*log10(P(1, 100, z) / nBefore);  eBeyond = 10*log10(P(1, nR-100, z) / nBeyond);
ok = check(wv.mtiCells < nR && eBefore < 10 && eBeyond > 50 && abs(10*log10(nBeyond / nBefore)) < 0.3, ...
           sprintf(['MTI up to 50 km (%d of %d long-pulse cells): echo before %.1f dB, beyond %.1f dB over the ', ...
                    'noise, noise beyond / before %.2f dB'], wv.mtiCells, nR, eBefore, eBeyond, ...
                   10*log10(nBeyond / nBefore))) && ok;

%% 3) long pulse off
c = base;
c.wave.long.enable = false;
c.targets = make_target('A', 5, 50e3, 40, 3000, -150, 0);
[plots, s] = run_case(c);
onT = plots_near(plots, c.targets, c.ant.height_m, 300);
ok = check(~s.wave.long.enable && s.wave.short.range_m(end) > 79e3 && any(onT & strcmp({plots.pulse}, 'short')) && ...
           ~any(strcmp({plots.pulse}, 'long')), sprintf(['long pulse off: short pulse up to %.1f km, %d plot(s), ', ...
           '%d short-pulse plot(s) on the target at 50 km'], s.wave.short.range_m(end)/1e3, numel(plots), sum(onT))) && ok;

%% 4) short pulse off
c = base;
c.wave.short.enable = false;
c.targets = [make_target('N', 5, 12e3, 40, 1000, -150, 0), make_target('F', 5, 45e3, 60, 3000, 150, 0)];
[plots, s] = run_case(c);
near = plots_near(plots, c.targets(1), c.ant.height_m, 300);
far  = plots_near(plots, c.targets(2), c.ant.height_m, 300);
ok = check(~any(near) && any(far & strcmp({plots.pulse}, 'long')) && ~any(strcmp({plots.pulse}, 'short')), ...
           sprintf('short pulse off: %d plot(s), %d at 12 km (blind range), %d long-pulse plot(s) at 45 km', ...
                   numel(plots), sum(near), sum(far))) && ok;

%% 5) own carriers
c = base;
c.wave.short.fc_Hz = 165e6;
c.wave.long.fc_Hz  = 175e6;
c.targets = [make_target('S', 5, 15e3, 40, 1000, 380, 0), make_target('L', 5, 50e3, 60, 3000, -400, 0)];
c.io.saveMf = true;  c.io.mfFolder = fullfile(c.io.outFolder, 'mf');
c.io.saveDoppler = true;  c.io.dopplerFolder = fullfile(c.io.outFolder, 'doppler');
c.proc.updateEvery = 25;
for f = {c.io.mfFolder, c.io.dopplerFolder}
    if exist(f{1}, 'dir'), rmdir(f{1}, 's'); end
end
[plots, s, res] = run_case(c);
bin = abs(diff(s.dop.f_Hz(1:2))) * s.lambda / 2;               % Doppler bin at cfg.fc [m/s]
good = true;  txt = '';
for i = 1:2
    pulse = {'short', 'long'};
    m = find(plots_near(plots, c.targets(i), c.ant.height_m, 300) & strcmp({plots.pulse}, pulse{i}));
    if isempty(m), good = false;  txt = [txt sprintf(' %s: no plot', c.targets(i).name)]; continue; end %#ok<AGROW>
    [~, ~, ~, vr] = target_geometry(c.targets(i), plots(m(1)).time_s, c.ant.height_m);
    err = abs(plots(m(1)).vr_mps - vr);
    wrong = abs(plots(m(1)).vr_mps * s.lambda / s.wave.(pulse{i}).lambda_m - vr);   % with the wavelength of cfg.fc
    good = good && err < bin;
    txt = [txt sprintf(' %s (%s, %.0f MHz): vr %.1f m/s, truth %.1f (error %.1f, with cfg.fc %.1f m/s);', ...
                       c.targets(i).name, pulse{i}, s.wave.(pulse{i}).fc_Hz/1e6, plots(m(1)).vr_mps, vr, err, wrong)]; %#ok<AGROW>
end
ok = check(good && s.wave.short.fIF_Hz == 30e6 && s.wave.long.fIF_Hz == 30e6 && s.wave.short.fLO_Hz == 135e6 ...
           && s.wave.long.fLO_Hz == 145e6 && isscalar(s.ddc.fNco_Hz) && s.ddc.fNco_Hz == 30e6, ...
           sprintf('own carriers (IF 30 MHz, LO 135 / 145 MHz, NCO fixed, Doppler bin %.1f m/s):%s', bin, txt)) && ok;
% at the ADC both echoes are on the 30 MHz IF; the short-pulse echo of the far target, received
% with the LO of the long pulse, on 30 + (165 - 175) = 20 MHz
sq = s;  sq.cfg.rx.addThermalNoise = false;  sq.cfg.rx.addLeakage = false;  sq.cfg.rx.stc.enable = false;
fPk = zeros(1, 3);  pul = {'short', 'long', 'short'};  tgt = [1 2 2];
for i = 1:3
    [R, azT] = target_geometry(c.targets(tgt(i)), 0, c.ant.height_m);
    x = double(generate_rx_signal(sq, c.targets(tgt(i)), 0, azT, []));
    wv = s.wave.(pul{i});
    n = find(s.tAdc >= wv.tTx_s + 2*R/c.const.c & s.tAdc < wv.tTx_s + 2*R/c.const.c + wv.width_s);
    seg = x(n, 1);
    nf = 2^nextpow2(16 * numel(seg));
    X = abs(fft(seg .* hanning(numel(seg)), nf));
    [~, k] = max(X(1:nf/2));
    fPk(i) = (k - 1) * s.fs / nf;
end
ok = check(all(abs(fPk - [30e6 30e6 20e6]) < 0.3e6), sprintf(['own carriers: echo spectra at the ADC: short %.2f MHz, ' ...
           'long %.2f MHz (30 MHz IF), short-pulse echo in the long window %.2f MHz (20 MHz)'], fPk / 1e6)) && ok;
c48 = c;  c48.ant.nElements = 48;
s48 = setup_radar(c48);
uExp = sind(s48.dbf.beamEl_deg - s48.ant.tilt_deg) * c.fc / 175e6;
ok = check(max(abs(s48.wave.long.beamU - uExp)) < 1e-12 && max(abs(s48.wave.short.beamU - uExp * 175/165)) < 1e-12, ...
           sprintf('48 elements: beam 21 at %.2f deg on fc, %.2f deg (short, 165 MHz), %.2f deg (long, 175 MHz)', ...
                   s48.dbf.beamEl_deg(end), s48.wave.short.beamEl_deg(end), s48.wave.long.beamEl_deg(end))) && ok;

%% 6) recordings of the matched filter and Doppler outputs
mf = dir(fullfile(c.io.mfFolder, 'mf_*.mat'));
dp = dir(fullfile(c.io.dopplerFolder, 'doppler_*.mat'));
good = numel(mf) == res.nProcessed && ~isempty(dp);
if good
    a = load(fullfile(c.io.mfFolder, mf(end).name));
    b = load(fullfile(c.io.dopplerFolder, dp(end).name));
    good = isequal(size(a.short), [s.nBeams s.wave.short.nRange]) && isequal(size(a.long), [s.nBeams s.wave.long.nRange]) && ...
           isequal(size(b.long), [s.nBeams s.wave.long.nRange s.dop.nfft]) && isfield(a, 'hdr') && isfield(b, 'hdr') && ...
           numel(b.longV_mps) == s.dop.nfft;
end
ok = check(good, sprintf('saveMf / saveDoppler: %d matched filter files (%d PRIs), %d Doppler files (update every %d PRIs)', ...
           numel(mf), res.nProcessed, numel(dp), c.proc.updateEvery)) && ok;
c.wave.long.enable = false;  c.io.saveDoppler = false;  c.sim.sectors_deg = [38 39];
if exist(c.io.mfFolder, 'dir'), rmdir(c.io.mfFolder, 's'); end
run_case(c, false);
mf = dir(fullfile(c.io.mfFolder, 'mf_*.mat'));
a = load(fullfile(c.io.mfFolder, mf(1).name));
ok = check(isempty(a.long) && size(a.short, 2) == setup_radar(c).wave.short.nRange, ...
           'saveMf with the long pulse off: its output is [] in the file') && ok;
rmdir(c.io.mfFolder, 's');
if exist(c.io.dopplerFolder, 'dir'), rmdir(c.io.dopplerFolder, 's'); end

%% 8) video noise floor of the PPI over a mountain sector (95% of the bins 30 dB above the noise, 45 deg)
s = setup_radar(c48);
V = ppi_video_init(s, 80e3);
rng(1);  lo = inf;  hi = -inf;
for az = 60:0.036:300                                          % 1 ms PRI at 6 rpm
    row = mean(-log(rand(48, V.nR)) * 500, 1);                 % 48 elements, noise power 500
    if az > 180 && az < 225, k = 1:floor(0.95*V.nR);  row(k) = row(k) * 1000; end
    V = ppi_video_add_row(V, row, az, az / 36);
    if az > 150, lo = min(lo, V.noise_dB);  hi = max(hi, V.noise_dB); end
end
ok = check(abs(V.noise_dB - 10*log10(500)) < 0.5 && hi - lo < 0.5, sprintf(['PPI video noise floor %.1f dB (noise ', ...
           '27.0 dB), %.2f dB of variation while the antenna crosses a mountain sector'], V.noise_dB, hi - lo)) && ok;

%% 7) the same receiver configuration from MATLAB and from C
exe = fullfile(receiver_dir(), 'build', 'radar_export_config');
if ispc, exe = [exe '.exe']; end
if exist(exe, 'file')
    cases = {'', {}
             'cfg.proc.cpiMode = ''manual''; cfg.proc.cpiLength = 64;', {'proc.cpiMode', 'manual'; 'proc.cpiLength', 64}
             'cfg.mti.maxRange_m = 50e3;', {'mti.maxRange_m', 50e3}
             'cfg.wave.long.enable = false;', {'wave.long.enable', false}
             'cfg.wave.short.enable = false;', {'wave.short.enable', false}
             'cfg.wave.short.fc_Hz = 165e6; cfg.wave.long.fc_Hz = 175e6;', {'wave.short.fc_Hz', 165e6; 'wave.long.fc_Hz', 175e6}
             'cfg.prf.mode = ''staggered''; cfg.wave.long.fc_Hz = 180e6; cfg.track.method = ''tomht'';', ...
                 {'prf.mode', 'staggered'; 'wave.long.fc_Hz', 180e6; 'track.method', 'tomht'}};
    for k = 1:size(cases, 1)
        c = radar_config();
        for j = 1:size(cases{k, 2}, 1)
            parts = strsplit(cases{k, 2}{j, 1}, '.');
            c = setfield(c, parts{:}, cases{k, 2}{j, 2});
        end
        fm = [tempname() '.bin'];  fc = [tempname() '.bin'];
        evalc('export_receiver_config(c, fm)');
        cmd = sprintf('"%s" --quiet --out "%s"', exe, fc);
        if ~isempty(cases{k, 1}), cmd = sprintf('%s --set "%s"', cmd, cases{k, 1}); end
        [status, msg] = system(cmd);
        [same, why] = same_records(fm, fc);
        ok = check(status == 0 && same, sprintf('radar_export_config = export_receiver_config: %s%s', ...
                   ternary(isempty(cases{k, 1}), 'default', cases{k, 1}), ternary(status == 0, why, [' ' msg]))) && ok;
        delete(fm);  if exist(fc, 'file'), delete(fc); end
    end
else
    fprintf('  radar_export_config not built (%s): MATLAB / C export comparison skipped\n', exe);
end
fprintf('New options (sub-pulses, carriers, CPI, MTI range, recordings): %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [plots, s, res] = run_case(c, sectors)
% Run the simulator over the sectors of the targets.
s = setup_radar(c);
if nargin < 2 || sectors
    c.sim.sectors_deg = target_sectors(s, c.targets);
end
res = run_radar_sim(c);
plots = res.plots;
if isempty(plots), plots = struct('pulse', {}, 'range_m', {}, 'time_s', {}, 'vr_mps', {}); end
end

function on = plots_near(plots, tgt, h, tol)
on = false(1, numel(plots));
for i = 1:numel(plots)
    R = target_geometry(tgt, plots(i).time_s, h);
    on(i) = abs(plots(i).range_m - R) < tol;
end
end

function msg = error_of(f)
msg = '';
try
    f();
catch err
    msg = err.message;
end
end

function [same, why] = same_records(a, b)
% The records of two receiver configuration files agree (relative 1e-9; NaN = NaN).
A = read_records(a);  B = read_records(b);
same = isequal(sort(fieldnames(A)), sort(fieldnames(B)));
why = '';
if ~same, why = sprintf(' (%d / %d records)', numel(fieldnames(A)), numel(fieldnames(B))); return; end
f = fieldnames(A);
for i = 1:numel(f)
    x = A.(f{i});  y = B.(f{i});
    d = abs(x(:) - y(:)) ./ max(max(abs(x(:)), abs(y(:))), 1e-300);
    d(x(:) == y(:) | (isnan(x(:)) & isnan(y(:)))) = 0;
    if ~isequal(size(x), size(y)) || any(d > 1e-9)
        same = false;  why = sprintf(' (record %s differs)', f{i});  return
    end
end
why = sprintf(' (%d records)', numel(f));
end

function R = read_records(file)
fid = fopen(file, 'r', 'ieee-le');
fread(fid, 8, 'char=>char');
h = fread(fid, 2, 'uint32');
R = struct();
for i = 1:h(2)
    name = deblank(char(fread(fid, 32, 'uint8=>char').'));
    name = name(name ~= 0);
    t = fread(fid, 3, 'uint32');
    v = fread(fid, t(2) * t(3) * (1 + t(1)), 'double');
    if t(1), v = complex(v(1:2:end), v(2:2:end)); end
    R.(name) = reshape(v, t(2), t(3));
end
fclose(fid);
end

function ok = check(cond, msg)
ok = logical(cond);
fprintf('  %-100s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

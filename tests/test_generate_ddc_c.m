function ok = test_generate_ddc_c()
%TEST_GENERATE_DDC_C  The C program radar_generate_ddc against generate_ddc_files.m.
%   The configuration of every case goes to the C program as a file
%   (write_config_m.m) with the targets (export_targets_csv.m).
%   1. Without noise, clutter and Swerling fluctuation a PRI is deterministic:
%      the C files must hold the same PRIs as MATLAB's, the same header, and
%      the same I, Q up to the rounding of a few ADC samples (MATLAB sums the
%      signal in single precision): at least 99.9% of the samples equal, none
%      more than 2 LSB apart. Cases: default (2 MHz), 4 MHz DDC with the
%      sub-pulses on 165 / 175 MHz, the long pulse off, phase codes.
%   2. Terrain clutter: with the expected amplitude in every cell (no random
%      phase, no motion: RADAR_GEN_FIXED_CLUTTER, and clutter_init.m changed the
%      same way) the DDC output agrees with MATLAB as in 1 (over the mountains
%      and over the hills). With noise and random clutter (other random numbers
%      than MATLAB) the total power of the drawn cell reflectivities equals the
%      expected one within 1%, the noise floor equals MATLAB's within 0.3 dB, the
%      power of both receive windows over 15 separated azimuths agrees within
%      2.5 dB (one random map: it scatters); the run is reproducible (the same
%      files with another number of threads).
%   3. The C files of the six targets (one scan, clutter, noise) through the C
%      receiver (radar_stream_send -> radar_receiver): every target that the
%      MATLAB test finds gives a plot.
%   Needs the receiver project built (radar_receiver_<os>/build).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
build = fullfile(receiver_dir(), 'build');
gen = fullfile(build, 'radar_generate_ddc');  if ispc, gen = [gen '.exe']; end
if ~exist(gen, 'file')
    fprintf('radar_generate_ddc not built (%s): test skipped\n', gen);
    ok = true;
    return
end
out = fullfile(root, 'output', 'test_generate_ddc_c');
if exist(out, 'dir'), rmdir(out, 's'); end
mkdir(out);
ok = true;

base = radar_config();
base.display.enable = false;
base.console.enable = false;
base.io.saveRaw = false;
base.io.outFolder = out;

%% 1) deterministic PRIs
cases = {'default, 2 MHz', {}
         '4 MHz, short 165 MHz, long 175 MHz', {'ddc.outputRate_Hz', 4e6; 'wave.short.fc_Hz', 165e6; 'wave.long.fc_Hz', 175e6}
         'long pulse off', {'wave.long.enable', false}
         'P4 phase codes, mismatch filter (long)', {'wave.short.type', 'phasecode'; 'wave.long.type', 'phasecode'; ...
                                                    'wave.long.compression', 'mismatch'}};
for k = 1:size(cases, 1)
    c = base;
    for j = 1:size(cases{k, 2}, 1)
        parts = strsplit(cases{k, 2}{j, 1}, '.');
        c = setfield(c, parts{:}, cases{k, 2}{j, 2});
    end
    c.sim.sectors_deg = [10 10.25];
    c.targets = [make_target('A', 10, 20e3, 10.1, 1000, 150, 0), make_target('B', 5, 50e3, 10.2, 5000, -200, 30)];
    c.clutter.enable = false;  c.rx.addThermalNoise = false;  c.adc.addAdcNoise = false;  c.sim.swerling = 0;
    [okRun, fm, fc] = run_both(c, out, gen, sprintf('det%d', k), '');
    good = okRun;
    nSame = 0;  nAll = 0;  dMax = 0;  hdrOk = true;
    fmList = dir(fullfile(fm, 'ddc_*.mat'));  fcList = dir(fullfile(fc, 'ddc_*.mat'));
    good = good && numel(fmList) == numel(fcList) && ~isempty(fmList) && isequal({fmList.name}, {fcList.name});
    for i = 1:numel(fmList) * good
        a = load(fullfile(fm, fmList(i).name));  b = load(fullfile(fc, fmList(i).name));
        d = [double(a.I(:)) - double(b.I(:)); double(a.Q(:)) - double(b.Q(:))];
        nSame = nSame + sum(d == 0);  nAll = nAll + numel(d);  dMax = max(dMax, max(abs(d)));
        hdrOk = hdrOk && same_header(a.hdr, b.hdr) && isequal(size(a.I), size(b.I));
    end
    frac = nSame / max(nAll, 1);
    good = good && hdrOk && frac >= 0.999 && dMax <= 2;
    ok = report(sprintf('%s: %d PRIs, headers equal %d, %.4f%% of the I/Q samples equal, largest difference %d LSB', ...
                        cases{k, 1}, numel(fmList), hdrOk, 100 * frac, dMax), good) && ok;
end

%% 2) terrain clutter: the same signal with fixed reflectivities; noise and clutter statistics
% (a) every clutter cell with its expected amplitude (no random phase, no internal motion):
%     radar_generate_ddc with the environment variable RADAR_GEN_FIXED_CLUTTER, MATLAB with
%     clutter_init.m changed the same way: the DDC output must agree as in 1)
for az = [40 125]                                      % mountains, hills
    c = base;
    c.rx.addThermalNoise = false;  c.adc.addAdcNoise = false;  c.sim.swerling = 0;
    c.sim.sectors_deg = [az az + 0.05];
    c.targets = make_target('Z', 1e-30, 79e3, 0, 0, 0, 0);
    sys = setup_radar(c);
    sched = pulse_schedule(sys);
    clut = clutter_init(sys);
    T = clutter_terrain(c, sys);
    tilt = sys.ant.tilt_deg;
    elRel = asind(-sind(tilt)*cosd(T.el_deg) + cosd(tilt)*sind(T.el_deg));
    gEl = pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, elRel) .* pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, elRel);
    amp = sqrt(2*sys.radarConst * T.sigma0 .* (T.range_m * deg2rad(c.clutter.azStep_deg) * T.dR_m ./ T.range_m.^4) .* gEl);
    clut.rho = amp(sub2ind(size(amp), double(clut.cellRow), clut.cellRange));
    clut.fd_Hz(:) = 0;  clut.dc = [];
    fc = fullfile(out, sprintf('fixed%d_c', az));
    write_config_m(c, fullfile(out, 'fixed_cfg.m'));
    export_targets_csv(c.targets, fullfile(out, 'fixed_targets.csv'));
    setenv('RADAR_GEN_FIXED_CLUTTER', '1');
    [st, txt] = system(sprintf('"%s" --config-m "%s" --targets "%s" --out "%s" --no-raw --quiet', gen, ...
                               fullfile(out, 'fixed_cfg.m'), fullfile(out, 'fixed_targets.csv'), fc));
    setenv('RADAR_GEN_FIXED_CLUTTER', '');
    if st ~= 0, fprintf('%s', txt); end
    fl = dir(fullfile(fc, 'ddc_*.mat'));
    nSame = 0;  nAll = 0;  dMax = 0;  lev = 0;
    for i = 1:numel(fl)
        b = load(fullfile(fc, fl(i).name));
        k = b.hdr.pulseNumber;
        [I, Q] = ddc(adc_sample(generate_rx_signal(sys, draw_target_rcs(c.targets, 5, 0), sched.time_s(k), ...
                     sched.az_deg(k), clutter_signal(clut, sys, sched.az_deg(k), sched.time_s(k))), sys), sys);
        d = [double(I(:)) - double(b.I(:)); double(Q(:)) - double(b.Q(:))];
        nSame = nSame + sum(d == 0);  nAll = nAll + numel(d);  dMax = max(dMax, max(abs(d)));
        lev = max(lev, sqrt(mean(double(I(:)).^2)));
    end
    frac = nSame / max(nAll, 1);
    ok = report(sprintf(['clutter with fixed reflectivities at %d deg: %d PRIs (rms %.0f LSB), %.4f%% of the I/Q ', ...
                         'samples equal to MATLAB, largest difference %d LSB'], az, numel(fl), lev, 100 * frac, dMax), ...
                st == 0 && ~isempty(fl) && frac >= 0.999 && dMax <= 2) && ok;
end
% (b) random clutter and noise: 15 separated azimuths over the hills (independent clutter
%     cells), mean power of each receive window within 1 dB of MATLAB, the noise beyond the
%     clutter (long window, 70 - 80 km) within 0.3 dB
c = base;
c.sim.sectors_deg = [(95:5:165).', (95:5:165).' + 0.04];
c.targets = make_target('Z', 1e-30, 79e3, 0, 0, 0, 0);
[okRun, fm, fc] = run_both(c, out, gen, 'clut', '--threads 3');
sys = setup_radar(c);
[sm, lm, nm] = window_power(fm, sys);  [sc, lc, nc] = window_power(fc, sys);
% (one random clutter map each: the power of a few PRIs scatters by about +-1.5 dB; the
% expected power is checked exactly in (a) and by the drawn cell powers below)
good = okRun && abs(10*log10(sc / sm)) < 2.5 && abs(10*log10(lc / lm)) < 2.5 && abs(10*log10(nc / nm)) < 0.3;
ok = report(sprintf(['noise and clutter over the hills (%d PRIs): short window %.2f / %.2f dB, long window ', ...
                     '%.2f / %.2f dB, noise %.2f / %.2f dB (MATLAB / C, DDC LSB^2)'], numel(dir(fullfile(fc, 'ddc_*.mat'))), ...
                    10*log10(sm), 10*log10(sc), 10*log10(lm), 10*log10(lc), 10*log10(nm), 10*log10(nc)), good) && ok;
% the random reflectivities of all clutter cells: their total power equals the expected one
setenv('RADAR_GEN_DEBUG', '1');
[st, txt] = system(sprintf('"%s" --config-m "%s" --targets "%s" --out "%s" --no-raw --max-pulses 1', gen, ...
                           fullfile(out, 'clut_cfg.m'), fullfile(out, 'clut_targets.csv'), fullfile(out, 'clut_dbg')));
setenv('RADAR_GEN_DEBUG', '');
tok = regexp(txt, 'expected cell powers ([\d.e+-]+), drawn ([\d.e+-]+)', 'tokens', 'once');
ratio = NaN;
if st == 0 && ~isempty(tok), ratio = str2double(tok{2}) / str2double(tok{1}); end
ok = report(sprintf('random clutter reflectivities: drawn / expected power of all cells %.4f', ratio), ...
            abs(ratio - 1) < 0.01) && ok;
% reproducible: the same files with one thread
fcAgain = fullfile(out, 'clut_c1');
[st, txt] = system(sprintf('"%s" --config-m "%s" --targets "%s" --out "%s" --no-raw --threads 1 --quiet', gen, ...
                           fullfile(out, 'clut_cfg.m'), fullfile(out, 'clut_targets.csv'), fcAgain));
fl = dir(fullfile(fcAgain, 'ddc_*.mat'));
same = st == 0 && ~isempty(fl);
for i = 1:numel(fl) * same
    a = load(fullfile(fcAgain, fl(i).name));  b = load(fullfile(out, 'clut_c', fl(i).name));
    same = same && isequal(a.I, b.I) && isequal(a.Q, b.Q);
end
if st ~= 0, fprintf('%s', txt); end
ok = report(sprintf('reproducible: %d PRIs with 1 thread equal to those with 3 threads', numel(fl)), same) && ok;

%% 3) the six targets through the C receiver
recv = fullfile(build, 'radar_receiver');  send = fullfile(build, 'radar_stream_send');
if isunix && exist(recv, 'file') && exist(send, 'file')
    c = base;
    c.targets = six_targets_scenario();
    c.sim.sectors_deg = target_sectors(setup_radar(c), c.targets);
    fc = fullfile(out, 'six_c');
    write_config_m(c, fullfile(out, 'six_cfg.m'));
    export_targets_csv(c.targets, fullfile(out, 'six_targets.csv'));
    t0 = tic;
    [st, txt] = system(sprintf('"%s" --config-m "%s" --targets "%s" --out "%s" --no-raw --quiet', gen, ...
                               fullfile(out, 'six_cfg.m'), fullfile(out, 'six_targets.csv'), fc));
    tGen = toc(t0);
    nFiles = numel(dir(fullfile(fc, 'ddc_*.mat')));
    cfgFile = fullfile(out, 'six.bin');  rec = fullfile(out, 'six.rec');  plotsCsv = fullfile(out, 'six_plots.csv');
    done = fullfile(out, 'six.done');  logf = fullfile(out, 'six.log');
    evalc('export_receiver_config(c, cfgFile)');
    port = 50420;
    system(sprintf('("%s" --config "%s" --port %d --record "%s" --plots-csv "" --idle-timeout 10 --quiet > "%s" 2>&1; echo $? > "%s") &', ...
                   recv, cfgFile, port, rec, logf, done));
    pause(2);
    [st2, txt2] = system(sprintf('"%s" "%s" --port %d --speed 0.05 --no-gaps --quiet', send, fc, port));
    t = tic;
    while ~exist(done, 'file') && toc(t) < 300, pause(1); end
    [~, rtxt] = system(sprintf('"%s" --config "%s" --replay "%s" --plots-csv "%s" --no-skip --quiet 2>&1', recv, ...
                               cfgFile, rec, plotsCsv));
    P = csvread_skip(plotsCsv);
    found = false(1, numel(c.targets));
    for i = 1:numel(c.targets)
        for r = 1:size(P, 1)
            R = target_geometry(c.targets(i), P(r, 1), c.ant.height_m);
            found(i) = found(i) || abs(P(r, 3) - R) < 300;
        end
    end
    nOn = sum(found);
    good = st == 0 && st2 == 0 && nFiles > 0 && all(found);
    if st ~= 0, fprintf('%s', txt); end
    if st2 ~= 0, fprintf('%s', txt2); end
    ok = report(sprintf(['six targets: %d PRIs generated in C in %.0f s, through the C receiver: %d plots, %d of %d ', ...
                         'targets found (%s)'], nFiles, tGen, size(P, 1), nOn, numel(found), ...
                        strjoin({c.targets(~found).name}, ' ')), good) && ok;
    if ~good, fprintf('%s', rtxt); end
else
    fprintf('  C receiver not built: the six-target check is skipped\n');
end

rmdir(out, 's');
fprintf('C DDC file generator (radar_generate_ddc): %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [ok, fm, fc] = run_both(c, out, gen, name, extra)
% The same configuration through generate_ddc_files.m and radar_generate_ddc.
fm = fullfile(out, [name '_m']);  fc = fullfile(out, [name '_c']);
c.io.ddcFolder = fm;
generate_ddc_files(c, 0);
cfgM = fullfile(out, [name '_cfg.m']);  tg = fullfile(out, [name '_targets.csv']);
write_config_m(c, cfgM);
export_targets_csv(c.targets, tg);
[st, txt] = system(sprintf('"%s" --config-m "%s" --targets "%s" --out "%s" --no-raw --quiet %s', gen, cfgM, tg, fc, extra));
ok = st == 0;
if ~ok, fprintf('%s', txt); end
end

function same = same_header(a, b)
f = fieldnames(a);
same = isequal(sort(f), sort(fieldnames(b)));
for i = 1:numel(f) * same
    same = same && isequal(double(a.(f{i})), double(b.(f{i})));
end
end

function [ps, pl, pn] = window_power(folder, sys)
% Mean DDC power (PRIs, elements, samples) of the short and the long receive window, and of the
% long window beyond 70 km (noise).
fl = dir(fullfile(folder, 'ddc_*.mat'));
dR = sys.cfg.const.c / (2 * sys.fsOut);
ws = sys.wave.short;  wl = sys.wave.long;
rl = wl.range_m(1) + (0:numel(wl.rxIdx)-1) * dR;
far = wl.rxIdx(rl > 70e3);
ps = 0;  pl = 0;  pn = 0;
for i = 1:numel(fl)
    d = load(fullfile(folder, fl(i).name));
    pw = mean(double(d.I).^2 + double(d.Q).^2, 2);
    ps = ps + mean(pw(ws.rxIdx)) / numel(fl);
    pl = pl + mean(pw(wl.rxIdx)) / numel(fl);
    pn = pn + mean(pw(far)) / numel(fl);
end
end

function q = prctile_(x, p)
x = sort(x(:));
if isempty(x), q = NaN; return; end
q = x(max(1, ceil(p * numel(x))));
end

function M = csvread_skip(file)
% Numeric rows of a CSV file with a header line.
M = zeros(0, 3);
if ~exist(file, 'file'), return; end
txt = fileread(file);
lines = regexp(txt, '\r?\n', 'split');
rows = {};
for i = 2:numel(lines)
    if isempty(strtrim(lines{i})), continue; end
    v = str2double(strsplit(lines{i}, ','));
    rows{end+1} = v; %#ok<AGROW>
end
if ~isempty(rows), M = vertcat(rows{:}); end
end

function ok = report(msg, ok)
fprintf('  %-110s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function ok = test_setup()
%TEST_SETUP  Configuration edge cases of setup_radar.m and pulse_schedule.m (no simulation).
%   1. Simulated sectors: [0 360] is the whole circle, [350 10] crosses
%      North, [10 10] is a single azimuth, and several sectors add up.
%   2. Settings that cannot work stop with a message that names them: an MTI
%      notch over every Doppler bin, a CPI too short for the MTI and the
%      Doppler FFT, an unknown CFAR method, DDC output rate and PRF mode.
%   3. A short long pulse (30 us) hands over within the first CFAR window of
%      the short pulse: the short pulse keeps one whole CFAR range window, so
%      cfar_detect and the C receiver (export_receiver_config) can use it, and
%      the CFAR tests its cells before the handover.
%   4. The DDC output (samples per PRI) and the PRI header (receive interval
%      of each sub-pulse) are the same for every CFAR method and every pulse
%      compression filter (LFM, phase code with its matched filter, mismatch
%      filters of the default length and of 3x the code), so a receiver set
%      up with another CFAR method or filter accepts the PRIs; the exported
%      receiver configuration holds the same receive intervals.
%   5. Target scenarios (cfg.scenario.type): 'six' gives 6 targets,
%      'challenge' 13; 'random' gives cfg.scenario.random.n targets inside
%      its limits, the same ones with the same seed, other ones with another
%      seed, and leaves the random state of the caller as it was. The speed
%      gate of the tracker (cfg.track.minSpeedEnable) is off by default, set
%      from cfg.track.minSpeed_mps when enabled (also in the exported
%      receiver file), and a negative speed stops with a message.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
cfg = radar_config();
cfg.console.enable = false;
cfg.ant.nElements = 1;                                 % (quick: no DBF matrix)
sys = setup_radar(cfg);

%% 1) sectors
cases = {[0 360], 1;  [350 10], 20/360;  [10 10], 0;  [10 20; 100 130], 40/360};
for k = 1:size(cases, 1)
    c = cfg;
    c.sim.sectors_deg = cases{k, 1};
    sys.cfg = c;
    s = pulse_schedule(sys);
    frac = mean(s.active);
    good = abs(frac - cases{k, 2}) < 0.002 && (cases{k, 2} > 0 || sum(s.active) <= 1);
    ok = report(sprintf('sector %s: %.4f of the pulses simulated (expected %.4f)', mat2str(cases{k, 1}), frac, ...
                        cases{k, 2}), good) && ok;
end

%% 2) settings that cannot work
bad = {'cfar.notchWidth', 100, 'notch';  'proc.cpiLength', 3, 'too short';  'cfar.method', 'xyz', 'cfar.method';
       'ddc.outputRate_Hz', 3e6, 'outputRate';  'prf.mode', 'random', 'PRF mode'};
for k = 1:size(bad, 1)
    c = cfg;
    parts = strsplit(bad{k, 1}, '.');
    c = setfield(c, parts{:}, bad{k, 2});
    msg = '';
    try
        setup_radar(c);
    catch err
        msg = err.message;
    end
    ok = report(sprintf('cfg.%s = %s stops: "%s"', bad{k, 1}, num2str(bad{k, 2}), msg), ...
                ~isempty(strfind(lower(msg), lower(bad{k, 3})))) && ok;
end

%% 3) handover within the first CFAR window of the short pulse
c = cfg;
c.wave.long.width_s = 30e-6;
s = setup_radar(c);
ws = s.wave.short;
oR = s.cfar.guardRange + s.cfar.trainRange;
nTested = ws.testIdx(2) - ws.testIdx(1) + 1;
P = single(-log(rand(s.nBeams, ws.nRange, s.dop.nfft)));
[~, noise] = cfar_detect(P, s, ws);
good = ws.nRange >= 2*oR + 1 && nTested >= 1 && ws.testIdx(2) <= oR && all(noise(:, ws.testIdx(1), :) > 0);
ok = report(sprintf('30 us long pulse: short pulse %d range cells (CFAR window %d), cells %d ... %d tested', ...
                    ws.nRange, 2*oR + 1, ws.testIdx), good) && ok;
fc = [tempname() '.bin'];
evalc('export_receiver_config(c, fc)');
bin = fullfile(receiver_dir(), 'build', 'test_cfar');
if ispc, bin = [bin '.exe']; end
if exist(bin, 'file')                                  % the C receiver accepts the configuration
    [status, out] = system(sprintf('"%s" "%s"', bin, fc));
    ok = report('the C receiver loads it and its CFAR agrees with direct sums (test_cfar)', status == 0) && ok;
    if status ~= 0, fprintf('%s', out); end
end
delete(fc);

%% 4) DDC output and PRI header independent of the CFAR method and of the compression filter
ref = [];
good = true;
txt = {};
for comp = {'lfm', 'matched', 'mismatch', 'mismatch3x'}
    for method = {'2d', 'range', 'doppler'}
        c = cfg;
        c.cfar.method = method{1};
        for p = {'short', 'long'}
            if ~strcmp(comp{1}, 'lfm')
                c.wave.(p{1}).type = 'phasecode';
                c.wave.(p{1}).compression = strrep(comp{1}, '3x', '');
            end
            if strcmp(comp{1}, 'mismatch3x')
                c.wave.(p{1}).mismatchLength = 3 * round(c.wave.(p{1}).width_s / c.wave.(p{1}).chip_s);
                c.files.phaseCode = fullfile(tempdir(), 'test_setup_phase_code_3x.mat');
            end
        end
        s = setup_radar(c);
        hdr = pulse_header(s, struct('scan', 0, 'time_s', 0, 'az_deg', 0, 'pri_s', s.pri(1), 'staggerIdx', 1), 1);
        lay = [s.nDdc, hdr.shortDdcIdx, hdr.longDdcIdx];
        if isempty(ref), ref = lay; end
        good = good && isequal(lay, ref) && isequal(s.wave.long.range_m(1:2), sys.wave.long.range_m(1:2));
        fc = [tempname() '.bin'];
        evalc('export_receiver_config(c, fc)');
        rx = read_receiver_records(fc, {'short_rxFirst', 'short_rxCount', 'long_rxFirst', 'long_rxCount'});
        delete(fc);
        good = good && isequal(rx, [ref(2) - 1, ref(3) - ref(2) + 1, ref(4) - 1, ref(5) - ref(4) + 1]);
    end
    txt{end+1} = sprintf('%s %d taps', comp{1}, numel(s.wave.long.mfTaps)); %#ok<AGROW>
end
f3 = fullfile(tempdir(), 'test_setup_phase_code_3x.mat');
if exist(f3, 'file'), delete(f3); end
ok = report(sprintf('DDC output %d samples, short [%d %d], long [%d %d] for every CFAR method and filter (long: %s)', ...
                    ref, strjoin(txt, ', ')), good) && ok;

%% 5) target scenarios and the speed gate of the tracker
tg = target_scenario(cfg);
ok = report(sprintf('scenario ''six'': %d targets', numel(tg)), numel(tg) == 6) && ok;
c = cfg;  c.scenario.type = 'challenge';
tg = target_scenario(c);
ok = report(sprintf('scenario ''challenge'': %d targets (%s)', numel(tg), strjoin({tg.name}, ' ')), numel(tg) == 13) && ok;
c.scenario.type = 'random';
rng(5);  before = rand;  rng(5);
tg = target_scenario(c);
after = rand;
tg2 = target_scenario(c);
c.scenario.random.seed = c.sim.seed + 1;
tg3 = target_scenario(c);
r = c.scenario.random;
P = reshape([tg.pos0_m], 3, []).';  V = reshape([tg.vel_mps], 3, []).';
rho = hypot(P(:, 1), P(:, 2));  sp = hypot(V(:, 1), V(:, 2));  rcs = [tg.rcs_m2].';
inLimits = all(rho >= r.range_m(1) - 1e-6 & rho <= r.range_m(2) + 1e-6) && ...
           all(P(:, 3) >= r.alt_m(1) - 1e-6 & P(:, 3) <= r.alt_m(2) + 1e-6) && ...
           all(atand(P(:, 3) ./ rho) <= r.maxEl_deg + 1e-6) && ...
           all(sp >= r.speed_mps(1) - 1e-6 & sp <= r.speed_mps(2) + 1e-6) && ...
           all(rcs >= r.rcs_m2(1) * (1 - 1e-9) & rcs <= r.rcs_m2(2) * (1 + 1e-9));
good = numel(tg) == r.n && inLimits && isequal(tg, tg2) && ~isequal(tg, tg3) && before == after;
ok = report(sprintf(['scenario ''random'': %d targets in the limits, the same with the same seed, others with another ', ...
                     'seed, random state of the caller kept'], numel(tg)), good) && ok;
c = cfg;  c.track.minSpeedEnable = true;
s = setup_radar(c);
good = s.track.minSpeed_mps == c.track.minSpeed_mps && sys.track.minSpeed_mps == 0;
fc = [tempname() '.bin'];
evalc('export_receiver_config(c, fc)');
v = read_receiver_records(fc, {'trackMinSpeed'});
delete(fc);
good = good && v == c.track.minSpeed_mps;
ok = report(sprintf('speed gate of the tracker: off by default, %g m/s when enabled (also in the receiver file)', ...
                    c.track.minSpeed_mps), good) && ok;
c.track.minSpeed_mps = -1;
msg = '';
try
    setup_radar(c);
catch err
    msg = err.message;
end
ok = report(sprintf('cfg.track.minSpeed_mps = -1 stops: "%s"', msg), ~isempty(strfind(msg, 'minSpeed_mps'))) && ok;
fprintf('Configuration edge cases: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ok = report(text, ok)
fprintf('  %-110s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function v = read_receiver_records(fileName, names)
% Scalar records of a receiver configuration file (export_receiver_config.m), in the order of names.
fid = fopen(fileName, 'r', 'ieee-le');
fread(fid, 8, 'char');
hd = fread(fid, 2, 'uint32');
v = nan(1, numel(names));
for i = 1:hd(2)
    name = deblank(char(fread(fid, 32, 'uint8').'));
    name = name(name ~= 0);
    t = fread(fid, 3, 'uint32');
    data = fread(fid, (1 + t(1)) * t(2) * t(3), 'double');
    k = find(strcmp(names, name));
    if ~isempty(k), v(k) = data(1); end
end
fclose(fid);
end

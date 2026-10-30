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

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
cfg = radar_config();
cfg.sim.verbose = false;
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
bin = fullfile(root, 'radar_receiver', 'build', 'test_cfar');
if ispc, bin = [bin '.exe']; end
if exist(bin, 'file')                                  % the C receiver accepts the configuration
    [status, out] = system(sprintf('"%s" "%s"', bin, fc));
    ok = report('the C receiver loads it and its CFAR agrees with direct sums (test_cfar)', status == 0) && ok;
    if status ~= 0, fprintf('%s', out); end
end
delete(fc);
fprintf('Configuration edge cases: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ok = report(text, ok)
fprintf('  %-110s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function [ok, results] = test_six_targets(prfMode, clutterOn, keepRaw)
%TEST_SIX_TARGETS  End-to-end test: six targets around the radar must be detected.
%
%   [ok, results] = TEST_SIX_TARGETS(prfMode, clutterOn, keepRaw) runs one antenna scan
%   with the six targets of six_targets_scenario.m (RCS 0.5 ... 10 m^2,
%   8 ... 46 km, radial speeds up to 290 m/s) and checks that every target
%   gives a plot with the right range, azimuth, elevation and (ambiguous)
%   radial velocity. prfMode is 'fixed' (default) or 'staggered'; clutterOn
%   (default true) switches the ground clutter on; keepRaw = false (default
%   true) deletes the raw files after they have been checked.
%
%   To save time only the azimuth sectors around the targets are simulated
%   (cfg.sim.sectors_deg); time and antenna rotation still run over the whole
%   scan. The raw ADC data of every simulated PRI is written to
%   output/test_<mode>/raw_pulses and the final PPI to output/test_<mode>/ppi.png.

if nargin < 1, prfMode = 'fixed'; end
if nargin < 2, clutterOn = true; end
if nargin < 3, keepRaw = true; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));

%% Configuration
cfg = radar_config();
cfg.prf.mode       = prfMode;
cfg.clutter.enable = clutterOn;
cfg.rx.maxRange_m  = 50e3;                       % targets are within 50 km
cfg.targets        = six_targets_scenario();
cfg.io.outFolder   = fullfile(root, 'output', ['test_' prfMode]);
cfg.io.rawFolder   = fullfile(cfg.io.outFolder, 'raw_pulses');
cfg.display.visible = ~isempty(get(0, 'ScreenSize')) && usejava('awt');
cfg.sim.sectors_deg = target_sectors(setup_radar(cfg), cfg.targets);
if exist(cfg.io.rawFolder, 'dir')
    delete(fullfile(cfg.io.rawFolder, 'pulse_*.mat'));    % only count the files of this run
end

%% Run
results = run_radar_sim(cfg);
sys   = results.sys;
plots = results.plots;

%% Associate plots with targets
gate.range_m = 300;      % includes the range-Doppler coupling of the LFM (~130 m at 290 m/s)
gate.az_deg  = 1.5;
gate.el_deg  = 2;
gate.vrBins  = 1.5;
used = false(1, numel(plots));
ok = true;
fprintf('\n%-4s %6s %6s | %8s %8s | %7s %7s | %6s %6s | %7s %7s | %6s %6s | %s\n', 'Tgt', 'RCS', 'pulse', ...
    'R true', 'R meas', 'az tru', 'az mea', 'el tru', 'el mea', 'v fold', 'v meas', 'SNRexp', 'SNRmea', 'result');
for i = 1:numel(cfg.targets)
    tgt = cfg.targets(i);
    best = 0;
    for j = 1:numel(plots)
        p = plots(j);
        [R, az, el, vr] = target_geometry(tgt, p.time_s, cfg.ant.height_m);
        wv = sys.wave.(p.pulse);
        vSpan = wv.lambda_m / (2*sys.dop.rowPri_s);                 % unambiguous velocity span
        vFold = fold_velocity(vr, vSpan);
        dv = abs(mod(p.vr_mps - vFold + vSpan/2, vSpan) - vSpan/2);
        elOk = isnan(p.el_deg) || abs(p.el_deg - el) < gate.el_deg;   % no elevation without DBF
        match = abs(p.range_m - R) < gate.range_m && abs(wrap180(p.az_deg - az)) < gate.az_deg ...
             && elOk && dv < gate.vrBins*vSpan/sys.dop.nfft;
        if match
            used(j) = true;
            if best == 0 || p.snr_dB > plots(best).snr_dB
                best = j;
            end
        end
    end
    if best == 0
        ok = false;
        fprintf('%-4s  no plot found  <-- FAIL\n', tgt.name);
        continue
    end
    p = plots(best);
    [R, az, el, vr] = target_geometry(tgt, p.time_s, cfg.ant.height_m);
    wv = sys.wave.(p.pulse);
    vFold = fold_velocity(vr, wv.lambda_m / (2*sys.dop.rowPri_s));
    % Expected SNR with the RCS this target had during the scan (Swerling 1)
    tgtScan = results.targets(i);
    tgtScan.rcs_m2 = target_rcs(tgtScan, target_look(sys, az, p.time_s));
    b = snr_budget(sys, tgtScan, p.time_s);
    fprintf('%-4s %6.2f %6s | %8.3f %8.3f | %7.2f %7.2f | %6.2f %6.2f | %7.2f %7.2f | %6.1f %6.1f | ok\n', ...
        tgt.name, tgtScan.rcs_m2, p.pulse, R/1e3, p.range_m/1e3, az, p.az_deg, el, p.el_deg, vFold, p.vr_mps, ...
        b.doppler, p.snr_dB);
end

nFalse = sum(~used);
fprintf('\nPlots: %d, associated with targets: %d, false plots: %d\n', numel(plots), sum(used), nFalse);
for j = find(~used)
    p = plots(j);
    fprintf('   false plot: t %.2f s, az %.1f deg, R %.2f km, el %.1f deg, vr %.1f m/s, SNR %.1f dB, %s\n', ...
        p.time_s, p.az_deg, p.range_m/1e3, p.el_deg, p.vr_mps, p.snr_dB, p.pulse);
end
if nFalse > 2
    ok = false;
end

%% Raw data files: one per simulated PRI, named by the pulse number
files = dir(fullfile(cfg.io.rawFolder, 'pulse_*.mat'));
raw = load(fullfile(cfg.io.rawFolder, files(end).name));
rawOk = numel(files) == results.nProcessed && isa(raw.adcData, 'int16') ...
     && isequal(size(raw.adcData), [sys.nAdc sys.nCh]) ...
     && strcmp(files(end).name, sprintf('pulse_%07d.mat', raw.hdr.pulseNumber)) ...
     && all(isfield(raw.hdr, {'scanNumber', 'antennaAz_deg', 'time_s', 'pri_s'}));
fprintf('Raw data: %d files in %s (%s)\n', numel(files), cfg.io.rawFolder, ternary(rawOk, 'ok', 'FAIL'));
ok = ok && rawOk;
if ~keepRaw
    delete(fullfile(cfg.io.rawFolder, 'pulse_*.mat'));
end
fprintf('Result (%s PRF, clutter %s): %s\n', prfMode, ternary(clutterOn, 'on', 'off'), ternary(ok, 'PASS', 'FAIL'));
end

function v = fold_velocity(vr, vSpan)
% Radial velocity folded into the unambiguous interval [-vSpan/2, vSpan/2).
v = mod(vr + vSpan/2, vSpan) - vSpan/2;
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function [ok, out] = test_scenarios(scenario, varargin)
%TEST_SCENARIOS  Several scans of a target scenario in the terrain clutter: detection and tracking.
%   [ok, out] = TEST_SCENARIOS(scenario, name, value, ...) runs run_radar_sim
%   for the scenario cfg.scenario.type = scenario ('six', 'challenge' or
%   'random') over the azimuth sectors of the targets and checks, per target,
%   the plots and the tracks of the MATLAB processing. Options (default):
%     'nScans'    3            antenna scans
%     'prf'       'fixed'      cfg.prf.mode
%     'terrain'   'regions'    cfg.clutter.terrain ('uniform', 'regions'; 'off' = no clutter)
%     'cnrOffset' 0            dB added to cfg.clutter.cnr_dB
%     'nElements' 48           cfg.ant.nElements (1: no DBF, 2-D tracking)
%     'seed'      1            cfg.sim.seed (and the random scenario)
%   A target is expected to be seen when snr_budget gives it at least
%   3 dB of margin (Pd 0.9 at the CFAR Pfa), it is not at a blind speed and
%   it lies in the instrumented range: it must give a plot (within 1 km and
%   3 deg) in at least nScans - 1 scans (1 of 1) and, with at least
%   cfg.track.confirmHits scans, end with a confirmed track on it (within
%   500 m in range and 2 deg in azimuth, updated in the last scan).
%   The other targets are only reported. Over all: at most 2 false plots per
%   scan (more than 1 km / 3 deg from every target), no false confirmed
%   track, and the ADC overload below cfg.adc.overloadWarn.
%   With the defaults one scan takes about 35 - 60 min in Octave (48 elements).

p = struct('nScans', 3, 'prf', 'fixed', 'terrain', 'regions', 'cnrOffset', 0, 'nElements', 48, 'seed', 1);
for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
cfg = radar_config();
cfg.console.plots = false;  cfg.display.enable = false;  cfg.io.saveRaw = false;
cfg.sim.seed = p.seed;
cfg.sim.nScans = p.nScans;
cfg.prf.mode = p.prf;
cfg.ant.nElements = p.nElements;
if strcmp(p.terrain, 'off')
    cfg.clutter.enable = false;
else
    cfg.clutter.terrain = p.terrain;
end
cfg.clutter.cnr_dB = cfg.clutter.cnr_dB + p.cnrOffset;
cfg.scenario.type = scenario;
cfg.targets = target_scenario(cfg);
sys = setup_radar(cfg);
cfg.sim.sectors_deg = target_sectors(sys, cfg.targets);
label = sprintf('%s scenario, %s PRF, clutter %s (CNR %.1f dB), %d element(s), %d scans', scenario, p.prf, ...
                p.terrain, cfg.clutter.cnr_dB, p.nElements, p.nScans);
fprintf('%s\n', label);
res = run_radar_sim(cfg);
sys = res.sys;
h = cfg.ant.height_m;
T = sys.track.scanTime_s;
P = res.plots;
L = res.tracker.log;                                  % [time id confirmed x y z vx vy vz hits plotIndex]
nT = numel(cfg.targets);
rCover = [sys.wave.short.range_m(1), sys.wave.long.range_m(end)];
ok = true;
out = struct('label', label, 'expected', false(1, nT), 'scansSeen', zeros(1, nT), 'tracked', false(1, nT));

% plots per target and scan, false plots
seen = false(nT, p.nScans);
nFalse = 0;
for i = 1:numel(P)
    near = false;
    s = min(floor(P(i).time_s / T), p.nScans - 1) + 1;
    for k = 1:nT
        [R, az] = target_geometry(cfg.targets(k), P(i).time_s, h);
        if abs(P(i).range_m - R) < 1000 && abs(wrap180(P(i).az_deg - az)) < 3
            near = true;  seen(k, s) = true;
        end
    end
    nFalse = nFalse + ~near;
end

% confirmed tracks: on a target (updated in the last scan) or false
tEnd = max([L(:, 1); 0]);
nFalseTracks = 0;
for id = unique(L(:, 2)).'
    rows = L(L(:, 2) == id, :);
    if ~any(rows(:, 3)), continue; end
    last = rows(end, :);
    pos = last(4:6);
    Rtrk = norm(pos);  azTrk = mod(atan2(pos(1), pos(2)) * 180/pi, 360);
    match = 0;
    for k = 1:nT
        [R, az] = target_geometry(cfg.targets(k), last(1), h);   % slant range (2-D tracks lie at slant range too)
        if abs(Rtrk - R) < 500 && abs(wrap180(azTrk - az)) < 2, match = k; break; end
    end
    if match == 0
        if ~any_target_near(cfg.targets, last(1), pos, h), nFalseTracks = nFalseTracks + 1; end
    elseif last(1) > tEnd - T
        out.tracked(match) = true;
    end
end

% expectations from the SNR budget
for k = 1:nT
    tg = cfg.targets(k);
    t1 = illumination_time(sys, tg);
    b = snr_budget(sys, tg, t1);
    R = target_geometry(tg, t1, h);
    out.expected(k) = ~b.blind && b.margin >= 3 && R >= rCover(1) && R <= rCover(2);
    out.scansSeen(k) = nnz(seen(k, :));
    needTrack = p.nScans >= sys.track.confirmHits;    % a track is confirmed after confirmHits looks
    good = ~out.expected(k) || (out.scansSeen(k) >= max(p.nScans - 1, 1) && (out.tracked(k) || ~needTrack));
    why = '';
    if b.blind, why = ' (blind speed)'; elseif b.margin < 3, why = sprintf(' (margin %.1f dB)', b.margin); end
    ok = report(sprintf('%-4s %5.1f m^2 %5.1f km az %5.1f: margin %5.1f dB%s, plots in %d of %d scans, %s', ...
                        tg.name, tg.rcs_m2, R/1e3, mod(atan2(tg.pos0_m(1), tg.pos0_m(2))*180/pi, 360), b.margin, ...
                        ternary(out.expected(k), '', [' not expected' why]), out.scansSeen(k), p.nScans, ...
                        ternary(out.tracked(k), 'confirmed track', 'no confirmed track')), good) && ok;
end
out.nFalsePlots = nFalse;  out.nFalseTracks = nFalseTracks;  out.adcClippedMax = res.adcClipped.max;
ok = report(sprintf('%d plots, %d false (at most %d); %d false confirmed track(s); ADC clipped at most %.3f%% of a PRI', ...
                    numel(P), nFalse, 2*p.nScans, nFalseTracks, 100*res.adcClipped.max), ...
            nFalse <= 2*p.nScans && nFalseTracks == 0 && res.adcClipped.max <= cfg.adc.overloadWarn) && ok;
fprintf('Scenario test (%s): %s\n', label, ternary(ok, 'PASS', 'FAIL'));
end

function near = any_target_near(targets, t, pos, h)
% A confirmed track within 3 km of a target (e.g. still converging) is not counted as false.
near = false;
for k = 1:numel(targets)
    d = targets(k).pos0_m + targets(k).vel_mps * t - [0 0 h];
    near = near || norm(pos(1:2) - d(1:2)) < 3000;
end
end

function ok = report(text, ok)
fprintf('  %-118s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

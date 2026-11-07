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
%     'gainOffset' 0           dB added to cfg.rx.gain_dB (e.g. -cnrOffset: no ADC overload)
%     'nElements' 48           cfg.ant.nElements (1: no DBF, 2-D tracking)
%     'seed'      1            cfg.sim.seed (and the random scenario)
%     'method'    'ekf'        cfg.track.method: 'ekf' (EKF, GNN) or 'pmbm' (PMBM, IMM-UKF)
%   The sectors around every look at every target are simulated (a crossing
%   target moves several degrees per scan). A target is expected in a scan
%   when snr_budget gives it at least 3 dB of margin at that look (Pd 0.9 at
%   the CFAR Pfa), it is not at a blind speed and it lies in the
%   instrumented range: it must give a plot (within 1 km and 3 deg) in all
%   but one of these scans and, when expected in at least
%   cfg.track.confirmHits scans including the last, end with a confirmed
%   track on it (within 500 m in range and 2 deg in azimuth, updated in the
%   last scan). Not expected in a scan: a blind speed, less than 3 dB of
%   margin, out of range, or another target at least as strong (-3 dB) in
%   its CFAR window (CA-CFAR masking: C6/C7 of the challenge scenario). A
%   target above the highest DBF beam is detected with a wrong elevation (and
%   its elevation sidelobes as extra plots near it): no 3-D track required.
%   A target whose azimuth changes by more than 15 deg between
%   looks (passing close to the radar) is not required to have a confirmed
%   track: a new track knows only the radial velocity, and the next plot
%   lies outside its gate (a limit of the track initiation; reported).
%   With one element (2-D tracking at slant range) a crossing target at high
%   elevation can lose its track as well: the slant-range projection bends
%   its path (reported as a failure).
%   A maneuvering target (C11 - C13 of the challenge scenario: more than
%   2 m/s^2 between its looks) needs its track with the PMBM tracker (IMM
%   with turn and acceleration models); the EKF (one motion model) is only
%   reported for it.
%   The other targets are only reported. Over all: at most 2 false plots per
%   scan (more than 1 km / 3 deg from every target), no false confirmed
%   track, and the ADC overload below cfg.adc.overloadWarn.
%   With the defaults one scan takes about 35 - 60 min in Octave (48 elements).

p = struct('nScans', 3, 'prf', 'fixed', 'terrain', 'regions', 'cnrOffset', 0, 'gainOffset', 0, 'nElements', 48, 'seed', 1, ...
           'method', 'ekf');
for i = 1:2:numel(varargin), p.(varargin{i}) = varargin{i+1}; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
cfg = radar_config();
cfg.console.plots = false;  cfg.display.enable = false;  cfg.io.saveRaw = false;
cfg.sim.seed = p.seed;
cfg.sim.nScans = p.nScans;
cfg.prf.mode = p.prf;
cfg.ant.nElements = p.nElements;
cfg.track.method = p.method;
if strcmp(p.terrain, 'off')
    cfg.clutter.enable = false;
else
    cfg.clutter.terrain = p.terrain;
end
cfg.clutter.cnr_dB = cfg.clutter.cnr_dB + p.cnrOffset;
cfg.rx.gain_dB = cfg.rx.gain_dB + p.gainOffset;
cfg.scenario.type = scenario;
cfg.targets = target_scenario(cfg);
sys = setup_radar(cfg);
T = sys.track.scanTime_s;
h = cfg.ant.height_m;
nT = numel(cfg.targets);
% look times of every target in every scan, and the simulated sectors around them
% (target_sectors.m around every look: a crossing target moves several degrees per scan)
tLook = zeros(nT, p.nScans);
cpiSpan = sys.cpi.length * mean(sys.pri) * sys.ant.rate_dps;
half = sys.ant.azBeamwidth_deg + cpiSpan/2 + 1;
sectors = zeros(0, 2);
for k = 1:nT
    tLook(k, :) = look_times(cfg.targets(k), cfg.ant.az0_deg, sys.ant.rate_dps, T, p.nScans, h);
    for s = find(isfinite(tLook(k, :)))
        [~, az] = target_geometry(cfg.targets(k), tLook(k, s), h);
        sectors(end+1, :) = mod([az - half, az + half], 360); %#ok<AGROW>
    end
end
cfg.sim.sectors_deg = sectors;
label = sprintf('%s scenario, %s PRF, clutter %s (CNR %.1f dB), gain %g dB, %d element(s), %d scans, %s tracker', ...
                scenario, p.prf, p.terrain, cfg.clutter.cnr_dB, cfg.rx.gain_dB, p.nElements, p.nScans, upper(p.method));
fprintf('%s\n', label);
res = run_radar_sim(cfg);
sys = res.sys;
P = res.plots;
L = res.tracker.log;                                  % [time id confirmed x y z vx vy vz hits plotIndex]
rCover = [sys.wave.short.range_m(1), sys.wave.long.range_m(end)];
ok = true;
out = struct('label', label, 'expected', false(1, nT), 'scansExpected', zeros(1, nT), 'scansSeen', zeros(1, nT), ...
             'tracked', false(1, nT));

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

% expectations from the SNR budget at every look
B = cell(nT, p.nScans);
for k = 1:nT
    for s = find(isfinite(tLook(k, :))), B{k, s} = snr_budget(sys, cfg.targets(k), tLook(k, s)); end
end
cw = sys.cfar;                                        % CFAR window (cells, bins on each side)
binV = abs(diff(sys.dop.v_mps(1:2)));
elTop = Inf;
if sys.dbf.enabled, elTop = max(sys.dbf.beamEl_deg) + abs(diff(sys.dbf.beamEl_deg(end-1:end))) / 2; end
for k = 1:nT
    tg = cfg.targets(k);
    expScan = false(1, p.nScans);
    notes = {};
    above = false;
    for s = 1:p.nScans
        if ~isfinite(tLook(k, s))                   % the antenna does not pass the target in this scan
            notes{end+1} = sprintf('scan %d no look (crosses the start azimuth)', s); %#ok<AGROW>
            continue
        end
        b = B{k, s};
        [R, ~, el, vr] = target_geometry(tg, tLook(k, s), h);
        inRange = R >= rCover(1) && R <= rCover(2);
        % CA-CFAR masking: a target at least as strong (-3 dB) inside the CFAR window of this one
        cell = abs(diff(sys.wave.(b.pulse).range_m(1:2)));
        masker = '';
        for j = [1:k-1, k+1:nT]
            if isempty(B{j, s}), continue; end
            [Rj, ~, ~, vrj] = target_geometry(cfg.targets(j), tLook(k, s), h);
            dBin = abs(mod(vrj - vr + sys.track.vrSpan/2, sys.track.vrSpan) - sys.track.vrSpan/2) / binV;
            if strcmp(B{j, s}.pulse, b.pulse) && B{j, s}.doppler >= b.doppler - 3 && ...
                    abs(Rj - R) / cell <= cw.guardRange + cw.trainRange + 1 && dBin <= cw.guardDoppler + cw.trainDoppler + 1
                masker = cfg.targets(j).name;
            end
        end
        expScan(s) = ~b.blind && b.margin >= 3 && inRange && isempty(masker);
        if b.blind, notes{end+1} = sprintf('scan %d blind speed', s); %#ok<AGROW>
        elseif b.margin < 3, notes{end+1} = sprintf('scan %d margin %.1f dB', s, b.margin); %#ok<AGROW>
        elseif ~inRange, notes{end+1} = sprintf('scan %d out of range %.1f km', s, R/1e3); %#ok<AGROW>
        elseif ~isempty(masker), notes{end+1} = sprintf('scan %d in the CFAR window of %s', s, masker); %#ok<AGROW>
        end
        if expScan(s) && el > elTop
            above = true;
            notes{end+1} = sprintf('scan %d el %.1f deg above the DBF beams (%.1f)', s, el, elTop); %#ok<AGROW>
        end
    end
    nExp = nnz(expScan);
    needTrack = nExp >= sys.track.confirmHits && expScan(end);
    % azimuth step between looks: a target that passes close to the radar
    % turns by tens of degrees per scan; a new track (velocity unknown but
    % the radial one) does not gate its next plot then (track initiation limit)
    tk = tLook(k, isfinite(tLook(k, :)));
    azL = zeros(1, numel(tk));
    for s = 1:numel(tk), [~, azL(s)] = target_geometry(tg, tk(s), h); end
    azStep = max([0, abs(wrap180(diff(azL)))]);
    if azStep > 15 && needTrack
        needTrack = false;
        notes{end+1} = sprintf('fast crosser %.0f deg/scan: track not required', azStep); %#ok<AGROW>
    end
    % maneuvering target: acceleration between its looks
    aMax = 0;
    for tt = linspace(min(tk), max(tk), 60)
        [~, ~, acc] = target_state(tg, tt);
        aMax = max(aMax, norm(acc));
    end
    if aMax > 2 && needTrack && ~strcmp(sys.track.method, 'pmbm')
        needTrack = false;
        notes{end+1} = sprintf('maneuver %.0f m/s^2: track required only with the PMBM tracker', aMax); %#ok<AGROW>
    elseif aMax > 2
        notes{end+1} = sprintf('maneuver %.0f m/s^2', aMax); %#ok<AGROW>
    end
    if above && sys.track.dim == 3 && needTrack
        needTrack = false;
        notes{end+1} = 'wrong elevation: 3-D track not required'; %#ok<AGROW>
    end
    out.expected(k) = nExp > 0;
    out.scansExpected(k) = nExp;
    out.scansSeen(k) = nnz(seen(k, expScan));
    good = out.scansSeen(k) >= nExp - (nExp > 1) && (out.tracked(k) || ~needTrack);
    why = '';
    if ~isempty(notes), why = [' (' strjoin(notes, ', ') ')']; end
    ok = report(sprintf('%-4s %5.1f m^2 %5.1f km az %5.1f: margin %5.1f dB, expected in %d of %d scans, plots in %d, %s%s', ...
                        tg.name, tg.rcs_m2, target_geometry(tg, tk(1), h)/1e3, ...
                        mod(atan2(tg.pos0_m(1), tg.pos0_m(2))*180/pi, 360), B{k, find(isfinite(tLook(k, :)), 1)}.margin, nExp, p.nScans, ...
                        nnz(seen(k, :)), ternary(out.tracked(k), 'confirmed track', 'no confirmed track'), why), good) && ok;
end
out.nFalsePlots = nFalse;  out.nFalseTracks = nFalseTracks;  out.adcClippedMax = res.adcClipped.max;
ok = report(sprintf('%d plots, %d false (at most %d); %d false confirmed track(s); ADC clipped at most %.3f%% of a PRI', ...
                    numel(P), nFalse, 2*p.nScans, nFalseTracks, 100*res.adcClipped.max), ...
            nFalse <= 2*p.nScans && nFalseTracks == 0 && res.adcClipped.max <= cfg.adc.overloadWarn) && ok;
fprintf('Scenario test (%s): %s\n', label, ternary(ok, 'PASS', 'FAIL'));
end

function tl = look_times(tg, az0, rate, T, nScans, h)
% Times when the antenna passes the target, one per scan (NaN: none in that
% scan, e.g. a target that crosses the start azimuth against the rotation,
% or along it: its look moves to the next scan, after the end of the run).
tt = 0:0.02:nScans*T;
f = zeros(size(tt));
for k = 1:numel(tt)
    [~, az] = target_geometry(tg, tt(k), h);
    f(k) = wrap180(az0 + rate*tt(k) - az);           % antenna ahead of the target: > 0
end
tl = nan(1, nScans);
if f(1) >= 0 && f(1) < rate * 0.02, tl(1) = 0; end   % the antenna starts at the target
for k = find(f(1:end-1) < 0 & f(2:end) >= 0)
    a = tt(k);  b = tt(k+1);
    for it = 1:30                                    % bisection
        c = (a + b) / 2;
        [~, az] = target_geometry(tg, c, h);
        if wrap180(az0 + rate*c - az) < 0, a = c; else, b = c; end
    end
    s = floor(a / T) + 1;
    if s <= nScans && ~isfinite(tl(s)), tl(s) = a; end
end
end

function near = any_target_near(targets, t, pos, h)
% A confirmed track within 3 km of a target (e.g. still converging) is not counted as false.
near = false;
for k = 1:numel(targets)
    d = target_state(targets(k), t) - [0 0 h];
    near = near || norm(pos(1:2) - d(1:2)) < 3000;
end
end

function ok = report(text, ok)
fprintf('  %-118s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

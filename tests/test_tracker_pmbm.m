function ok = test_tracker_pmbm()
%TEST_TRACKER_PMBM  PMBM tracker with IMM-UKF (cfg.track.method = 'pmbm') on synthetic plots.
%   Every case runs in 3-D (DBF: range, azimuth, elevation) and in 2-D (one
%   element: range and azimuth, slant range), unless noted. The plots come
%   from the true target positions (target_state.m) when the antenna points
%   at the target, with the errors of cfg.track, the radial velocity folded
%   into the Doppler span as the radar measures it, and reach the tracker
%   0.1 s later; the tracker is called every 0.1 s as well (as every PRI in
%   run_radar_sim). The EKF (cfg.track.method = 'ekf', CV model) runs on
%   the same plots for comparison.
%   1. Six targets (six_targets_scenario.m) with 8 random false plots per
%      scan, 5 scans: every target ends with one confirmed track close to its
%      position and velocity, no false plot gives a confirmed track. Also with
%      the staggered PRF (radial velocity folded at +-112 m/s: unfolded).
%   2. Maneuvering target: 250 m/s, straight for 3 scans, then a 4 deg/s
%      coordinated turn (1.8 g) for 30 s (120 deg), then straight again,
%      12 scans: one track follows it all the time (a single confirmed track,
%      updated by every plot after its confirmation), its position error is
%      smaller than the one of the EKF, and the turn model (CT) dominates the
%      IMM during the turn.
%   3. Constant acceleration: 6 m/s^2 along the track for 50 s (100 ->
%      400 m/s), crossing the line of sight at 45 deg, 10 scans: one track
%      follows it, the CA model dominates the IMM, the velocity and
%      acceleration errors at the end are small, the velocity error smaller
%      than the one of the EKF.
%   4. Low detection probability and dense false plots: the six targets with
%      Pd = 0.7 and 30 false plots per scan (cfg.track.pmbm set to these
%      values), 10 scans (T1 receding): every target detected in the last 3
%      scans ends with a confirmed track (the others are reported), at most
%      one false confirmed track. The EKF is reported for comparison.
%   5. Two targets 200 m apart in range on the same azimuth (as the GNN case
%      of test_tracker.m), and two targets crossing each other: two tracks,
%      every one on its own target (no swap).
%   6. Track termination and the speed gate: a target that disappears loses
%      its track after the missed scans (existence below rDelete); with
%      cfg.track.minSpeedEnable a 3 m/s target is deleted, a 100 m/s one kept.
%   7. The thirteen targets of challenge_targets_scenario.m (with the weave
%      C11, the constant acceleration C12, the turning, slowing, climbing
%      C13): a target is detected at a look with probability 0.9 when
%      snr_budget.m gives it at least 3 dB of margin and it is in the
%      instrumented range and not at a blind speed (else not detected), 6
%      scans, 10 false plots per scan, fixed PRF. Every target detected in at
%      least the last 3 scans ends with a confirmed track on it, there is no
%      false confirmed track; C11 - C13 also with the EKF for comparison.
%   When the C receiver is built (radar_receiver_<os>/build/test_tracker), the
%   C PMBM tracker (tracker_pmbm.c) runs on the same calls of every PMBM case:
%   every track update must match MATLAB (plot, track number, confirmation,
%   updates, state).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
c_compare('reset');
for dim = [3 2]
    ok = six_case(dim, 'fixed') && ok;
end
ok = six_case(3, 'staggered') && ok;
for dim = [3 2]
    ok = turn_case(dim) && ok;
    ok = accel_case(dim) && ok;
end
for dim = [3 2]
    ok = clutter_case(dim) && ok;
end
ok = close_case() && ok;
ok = crossing_case() && ok;
ok = termination_case() && ok;
ok = min_speed_case() && ok;
for dim = [3 2]
    ok = challenge_case(dim) && ok;
end
[nC, nBad] = c_compare('get');
if nC > 0
    fprintf('  C tracker (tracker_pmbm.c) on the same calls: %d runs compared, %d different from MATLAB: %s\n', nC, nBad, ...
            ternary(nBad == 0, 'ok', 'FAIL'));
    ok = ok && nBad == 0;
else
    fprintf('  C tracker not built: comparison skipped\n');
end
fprintf('PMBM tracker (IMM-UKF): %s\n', ternary(ok, 'PASS', 'FAIL'));
end

%% ------------------------------------------------------------------ cases
function ok = six_case(dim, prf)
[cfg, sys] = make_cfg(dim, 'pmbm', prf, struct('falsePlotsPerScan', 8));
tgts = six_targets_scenario();
nScans = 5;                                          % (T1 passes over the radar after 50 s)
pl = make_plots(cfg, sys, tgts, nScans, struct('nFalse', 8, 'seed', 7));
tr = run_tracker(sys.track, pl, nScans * sys.track.scanTime_s, cfg);
[good, txt] = check_targets(cfg, tr, tgts, nScans * sys.track.scanTime_s, 500, 25);   % 0.3 deg at 46 km: 240 m
nF = false_tracks(cfg, tr, tgts);
ok = good && nF == 0;
fprintf('  Six targets, %d-D, %s PRF, 8 false plots per scan: %s; %d false confirmed track(s): %s\n', dim, prf, ...
        txt, nF, ternary(ok, 'ok', 'FAIL'));
end

function ok = turn_case(dim)
% 250 m/s, straight 3 scans, 4 deg/s turn for 30 s, straight again.
tg = make_target('M', 2, 25e3, 60, 3000, -150, 200, [30 4 0 0; 60 0 0 0]);
nScans = 12;
[good, e, txt] = single_target_case(dim, tg, nScans, [30 60]);
okE = e.pmbm.pos < e.ekf.pos;
okM = e.modeDuring(3) > max(e.modeDuring(1:2));              % CT dominates during the turn
ok = good && okE && okM;
fprintf(['  Turn 4 deg/s (1.8 g), %d-D: %s; position error (rms over the updates) PMBM %.0f m, EKF %.0f m ' ...
         '(EKF tracks: %d); model probabilities in the turn CV %.2f CA %.2f CT %.2f: %s\n'], dim, txt, ...
        e.pmbm.pos, e.ekf.pos, e.ekf.nTracks, e.modeDuring, ternary(ok, 'ok', 'FAIL'));
end

function ok = accel_case(dim)
% 6 m/s^2 along the track for 50 s, crossing the line of sight at 45 deg.
tg = make_target('A', 2, 20e3, 300, 2000, -100/sqrt(2), 100/sqrt(2), [0 0 6 0; 50 0 0 0]);
nScans = 10;
[good, e, txt] = single_target_case(dim, tg, nScans, [10 50]);
okA = e.modeDuring(2) > max(e.modeDuring([1 3])) && e.pmbm.velEnd < 15 && e.pmbm.accEnd < 2.5;
okE = e.pmbm.vel < e.ekf.vel;
ok = good && okA && okE;
fprintf(['  Constant acceleration 6 m/s^2, %d-D: %s; velocity error (rms) PMBM %.1f m/s, EKF %.1f m/s; at the ' ...
         'end %.1f m/s, acceleration error %.2f m/s^2; model probabilities CV %.2f CA %.2f CT %.2f: %s\n'], dim, ...
        txt, e.pmbm.vel, e.ekf.vel, e.pmbm.velEnd, e.pmbm.accEnd, e.modeDuring, ternary(ok, 'ok', 'FAIL'));
end

function [ok, e, txt] = single_target_case(dim, tg, nScans, tMan)
% One target, PMBM and EKF on the same plots: errors over the confirmed
% updates and the model probabilities while it maneuvers (tMan).
e = struct();
T = [];
for method = {'pmbm', 'ekf'}
    [cfg, sys] = make_cfg(dim, method{1}, 'fixed');
    T = sys.track.scanTime_s;
    pl = make_plots(cfg, sys, tg, nScans, struct('nFalse', 0, 'seed', 3));
    [tr, modes] = run_tracker(sys.track, pl, nScans * T, cfg);
    L = tr.log;
    ids = unique(L(L(:, 3) == 1, 2));
    s = struct('nTracks', numel(unique(L(:, 2))), 'nConfirmed', numel(ids), 'pos', Inf, 'vel', Inf, ...
               'velEnd', Inf, 'accEnd', Inf, 'updates', 0);
    if ~isempty(ids)
        rows = L(ismember(L(:, 2), ids) & L(:, 3) == 1, :);
        err = zeros(size(rows, 1), 2);
        for k = 1:size(rows, 1)
            [pt, vt] = truth_xyz(tg, rows(k, 1), cfg, dim);
            err(k, :) = [norm(rows(k, 4:3+dim) - pt(1:dim)), norm(rows(k, 7:8) - vt(1:2))];
        end
        s.pos = sqrt(mean(err(:, 1).^2));
        s.vel = sqrt(mean(err(:, 2).^2));
        s.updates = size(rows, 1);
        rep = track_report(tr);
        rep = rep([rep.confirmed]);
        if ~isempty(rep)
            [~, ~, at] = target_state(tg, rep(1).t_s);
            [~, vt] = truth_xyz(tg, rep(1).t_s, cfg, dim);
            s.velEnd = norm(rep(1).vel_mps(1:2) - vt(1:2));
            s.accEnd = norm(rep(1).acc_mps2(1:2) - at(1:2));
        end
    end
    e.(method{1}) = s;
    if strcmp(method{1}, 'pmbm')
        nPlots = numel(pl);
        inMan = modes(:, 1) >= tMan(1) + T & modes(:, 1) <= tMan(2) + T;   % updates after the start
        e.modeDuring = mean(modes(inMan, 2:4), 1);
        firstConf = find(L(:, 3) == 1, 1);
        follow = numel(ids) == 1 && all(L(firstConf:end, 2) == ids) && ...
                 size(L, 1) - firstConf + 1 >= nPlots - tr.p.confirmHits;
        ok = follow && e.pmbm.nTracks == 1;
        txt = sprintf('%d plots, %d track(s), %d confirmed, %d confirmed updates', nPlots, e.pmbm.nTracks, ...
                      e.pmbm.nConfirmed, e.pmbm.updates);
    end
end
end

function ok = clutter_case(dim)
% Pd 0.7, 30 false plots per scan, 10 scans.
[cfg, sys] = make_cfg(dim, 'pmbm', 'fixed', struct('pd', 0.7, 'falsePlotsPerScan', 30));
tgts = six_targets_scenario();
tgts(1).vel_mps = -tgts(1).vel_mps;                  % T1 recedes (inbound it passes over the radar after 50 s)
nScans = 10;
[pl, detected] = make_plots(cfg, sys, tgts, nScans, struct('nFalse', 30, 'pd', 0.7, 'seed', 21));
tEnd = nScans * sys.track.scanTime_s;
tr = run_tracker(sys.track, pl, tEnd, cfg);
[cfgE, sysE] = make_cfg(dim, 'ekf', 'fixed');
trE = run_tracker(sysE.track, pl, tEnd);
need = all(detected(:, end-2:end), 2).';            % detected in the last 3 scans: a track is required
has = arrayfun(@(k) on_target(cfg, tr, tgts(k), tEnd, dim), 1:numel(tgts));
hasE = arrayfun(@(k) on_target(cfgE, trE, tgts(k), tEnd, dim), 1:numel(tgts));
nF = false_tracks(cfg, tr, tgts);
nFE = false_tracks(cfgE, trE, tgts);
ok = all(has | ~need) && nF <= 1;
fprintf(['  Pd 0.7, 30 false plots per scan, %d-D, 10 scans: detections per target %s, %d detected in the ' ...
         'last 3 scans; confirmed tracks on targets PMBM %d (%s), EKF %d (%s); false confirmed tracks PMBM %d, ' ...
         'EKF %d; %d global hypotheses, %d Bernoulli components at the end: %s\n'], dim, mat2str(sum(detected, 2).'), ...
        sum(need), sum(has), strjoin({tgts(has).name}, ' '), sum(hasE), strjoin({tgts(hasE).name}, ' '), nF, nFE, ...
        size(tr.G, 1), numel(tr.T), ternary(ok, 'ok', 'FAIL'));
end

function ok = close_case()
% Two targets 200 m apart in range, same azimuth, receding at 100 m/s (2-D).
[cfg, sys] = make_cfg(2, 'pmbm', 'fixed');
tg(1) = make_target('A', 2, 20e3, 30, 0, 100, 0);
tg(2) = make_target('B', 2, 20.2e3, 30, 0, 100, 0);
nScans = 6;
pl = make_plots(cfg, sys, tg, nScans, struct('nFalse', 0, 'seed', 4, 'noise', 0));
tr = run_tracker(sys.track, pl, nScans * sys.track.scanTime_s, cfg);
[good, txt] = check_targets(cfg, tr, tg, nScans * sys.track.scanTime_s, 50, 5);
rep = track_report(tr);
ok = good && sum([rep.confirmed]) == 2;
fprintf('  Two targets 200 m apart (2-D): %s: %s\n', txt, ternary(ok, 'ok', 'FAIL'));
end

function ok = crossing_case()
% Two targets crossing each other at 90 deg at 25 km (3-D, other altitudes; 2-D).
ok = true;
for dim = [3 2]
    [cfg, sys] = make_cfg(dim, 'pmbm', 'fixed');
    T = sys.track.scanTime_s;
    tc = 4.5 * T;                                    % they cross between two looks
    c = [25e3*sind(80), 25e3*cosd(80)];
    v1 = 200 * [sind(45) cosd(45)];  v2 = 200 * [sind(135) cosd(135)];
    tg(1) = target_at(c - v1*tc, [v1 0], 3000, 'X1'); %#ok<AGROW>
    tg(2) = target_at(c - v2*tc, [v2 0], 3600, 'X2'); %#ok<AGROW>
    nScans = 9;
    pl = make_plots(cfg, sys, tg, nScans, struct('nFalse', 4, 'seed', 5));
    tr = run_tracker(sys.track, pl, nScans * T, cfg);
    % the track of each target before the crossing must be the one after it
    L = tr.log(tr.log(:, 3) == 1, :);
    owner = zeros(2, 2);                             % track of target k before / after
    for k = 1:2
        for phase = 1:2
            if phase == 1, sel = L(:, 1) < tc - T/2; else, sel = L(:, 1) > tc + T; end
            rows = L(sel, :);
            best = 0;  bd = Inf;
            for r = 1:size(rows, 1)
                pt = truth_xyz(tg(k), rows(r, 1), cfg, dim);
                dd = norm(rows(r, 4:5) - pt(1:2));
                if dd < bd, bd = dd; best = rows(r, 2); end
            end
            if bd < 1000, owner(k, phase) = best; end
        end
    end
    good = all(owner(:) > 0) && owner(1, 1) == owner(1, 2) && owner(2, 1) == owner(2, 2) && owner(1, 1) ~= owner(2, 1);
    fprintf('  Two targets crossing at 90 deg, %d-D: tracks before %s, after %s: %s\n', dim, mat2str(owner(:, 1).'), ...
            mat2str(owner(:, 2).'), ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
end

function ok = termination_case()
% A target seen for 4 scans, then gone (3-D): its track must be deleted.
[cfg, sys] = make_cfg(3, 'pmbm', 'fixed');
T = sys.track.scanTime_s;
tg = make_target('G', 2, 30e3, 100, 5000, -200, 0);
pl = make_plots(cfg, sys, tg, 4, struct('nFalse', 0, 'seed', 8));
tr = track_init(sys.track);
tr = feed(tr, pl, 4*T);
nConf = sum([tr.tracks.confirmed]);
del = NaN;
for t = 4*T:0.1:12*T                                  % no more plots
    tr = track_update(tr, [], t);
    if isempty(tr.tracks), del = t; break; end
end
scans = (del - pl(end).time_s) / T;
ok = nConf == 1 && isfinite(del) && scans <= sys.track.maxMiss + 1.5;
fprintf('  Track termination: confirmed after 4 scans %d, deleted %.1f scans after the last plot: %s\n', nConf, ...
        scans, ternary(ok, 'ok', 'FAIL'));
end

function ok = min_speed_case()
% Speed gate 5 m/s: a 3 m/s and a 100 m/s target (2-D, plots without errors).
ok = true;
for gate = [false true]
    [cfg, sys] = make_cfg(2, 'pmbm', 'fixed', struct(), gate);
    tg(1) = make_target('S', 2, 20e3, 30, 0, 3, 0);
    tg(2) = make_target('F', 2, 20e3, 120, 0, 100, 0);
    nScans = 4;
    pl = make_plots(cfg, sys, tg, nScans, struct('nFalse', 0, 'seed', 1, 'noise', 0));
    tr = run_tracker(sys.track, pl, nScans * sys.track.scanTime_s, cfg);
    rep = track_report(tr);
    v = arrayfun(@(r) norm(r.vel_mps(1:2)), rep);
    nSlow = sum(v < 50);  nFast = sum(v >= 50 & [rep.confirmed]);
    good = nFast == 1 && nSlow == ~gate;
    fprintf('  Speed gate %s: %d slow and %d fast confirmed track(s) at the end: %s\n', ternary(gate, '5 m/s', 'off'), ...
            nSlow, nFast, ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
end

function ok = challenge_case(dim)
% The challenge scenario, detections from the SNR budget.
[cfg, sys] = make_cfg(dim, 'pmbm', 'fixed', struct('falsePlotsPerScan', 10));
cfg.scenario.type = 'challenge';
tgts = target_scenario(cfg);
nScans = 6;
T = sys.track.scanTime_s;
[pl, detected] = make_plots(cfg, sys, tgts, nScans, struct('nFalse', 10, 'seed', 11, 'budget', true));
tEnd = nScans * T;
tr = run_tracker(sys.track, pl, tEnd, cfg);
[cfgE, sysE] = make_cfg(dim, 'ekf', 'fixed');
trE = run_tracker(sysE.track, pl, tEnd);
ok = true;
names = {};
for k = 1:numel(tgts)
    need = all(detected(k, end-2:end));
    has = on_target(cfg, tr, tgts(k), tEnd, dim);
    hasE = on_target(cfgE, trE, tgts(k), tEnd, dim);
    good = has || ~need;
    if ~good || k >= 11
        names{end+1} = sprintf('%s %s/%s%s', tgts(k).name, ternary(has, 'yes', 'no'), ternary(hasE, 'yes', 'no'), ...
                               ternary(good, '', ' FAIL')); %#ok<AGROW>
    end
    ok = ok && good;
end
nF = false_tracks(cfg, tr, tgts);  nFE = false_tracks(cfgE, trE, tgts);
nNeed = sum(all(detected(:, end-2:end), 2));
nHas = sum(arrayfun(@(k) on_target(cfg, tr, tgts(k), tEnd, dim), 1:numel(tgts)));
nHasE = sum(arrayfun(@(k) on_target(cfgE, trE, tgts(k), tEnd, dim), 1:numel(tgts)));
ok = ok && nF == 0;
fprintf(['  Challenge scenario (13 targets), %d-D, 6 scans: detected in the last 3 scans %d, confirmed tracks on ' ...
         'targets PMBM %d, EKF %d; false confirmed PMBM %d, EKF %d; confirmed track PMBM/EKF: %s: %s\n'], dim, ...
        nNeed, nHas, nHasE, nF, nFE, strjoin(names, ', '), ternary(ok, 'ok', 'FAIL'));
end

%% ------------------------------------------------------------------ helpers
function [cfg, sys] = make_cfg(dim, method, prf, pmbm, minSpeed)
cfg = radar_config();
cfg.console.enable = false;
cfg.track.method = method;
cfg.track.dim = sprintf('%dd', dim);
cfg.prf.mode = prf;
if dim == 2, cfg.ant.nElements = 1; end              % no elevation
if nargin > 3
    f = fieldnames(pmbm);
    for i = 1:numel(f), cfg.track.pmbm.(f{i}) = pmbm.(f{i}); end
end
if nargin > 4 && minSpeed
    cfg.track.minSpeedEnable = true;
    cfg.track.minSpeed_mps = 5;
end
persistent cache
key = sprintf('%d_%s', dim, prf);
if isempty(cache), cache = struct(); end
if isfield(cache, key)                               % the radar setup is the same; only the tracker differs
    sys = cache.(key);
    sys.track = track_setup(cfg, sys);
else
    sys = setup_radar(cfg);
    cache.(key) = sys;
end
end

function tg = target_at(xy, vel, alt, name)
tg = make_target(name, 2, 1, 0, alt, 0, 0);
tg.pos0_m = [xy alt];
tg.vel_mps = vel;
end

function [pl, detected] = make_plots(cfg, sys, tgts, nScans, o)
% Plots of the targets at their looks (and false plots), in time order.
if ~isfield(o, 'pd'), o.pd = 1; end
if ~isfield(o, 'noise'), o.noise = 1; end
if ~isfield(o, 'budget'), o.budget = false; end
p = sys.track;
T = p.scanTime_s;
h = cfg.ant.height_m;
rate = sys.ant.rate_dps;
dim = p.dim;
saved = rng;
rng(o.seed);
pl = struct('time_s', {}, 'az_deg', {}, 'range_m', {}, 'el_deg', {}, 'vr_mps', {});
detected = false(numel(tgts), nScans);
rCover = [sys.wave.short.range_m(1), sys.wave.long.range_m(end)];
looks = cell(1, numel(tgts));
for i = 1:numel(tgts), looks{i} = look_times(tgts(i), cfg, rate, nScans * T); end
for s = 0:nScans-1
    for i = 1:numel(tgts)
        t = looks{i}(s+1);
        if ~isfinite(t), continue; end
        [R, az, el, vr] = target_geometry(tgts(i), t, h);
        pd = o.pd;
        if o.budget
            b = snr_budget(sys, tgts(i), t);
            pd = 0.9 * (~b.blind && b.margin >= 3 && R >= rCover(1) && R <= rCover(2));
        end
        e = randn(1, 4) * o.noise;
        if rand >= pd, continue; end
        detected(i, s+1) = true;
        vm = vr + sqrt(p.Rvr) * e(4);
        vm = mod(vm + p.vrSpan/2, p.vrSpan) - p.vrSpan/2;   % folded as the radar measures it
        elm = el + cfg.track.sigmaEl_deg * e(3);
        if dim == 2, elm = NaN; end
        pl(end+1) = struct('time_s', t, 'az_deg', mod(az + cfg.track.sigmaAz_deg * e(2), 360), ...
                           'range_m', R + cfg.track.sigmaRange_m * e(1), 'el_deg', elm, 'vr_mps', vm); %#ok<AGROW>
    end
    for f = 1:o.nFalse
        az = 360*rand;
        el = 20*rand;  if dim == 2, el = NaN; end
        pl(end+1) = struct('time_s', s*T + mod(az - cfg.ant.az0_deg, 360)/rate, 'az_deg', az, ...
                           'range_m', 5e3 + 60e3*rand, 'el_deg', el, 'vr_mps', p.vrSpan*(rand - 0.5)); %#ok<AGROW>
    end
end
rng(saved);
[~, o2] = sort([pl.time_s]);
pl = pl(o2);
end

function tl = look_times(tg, cfg, rate, tEnd)
% Times when the antenna points at the target, one per scan (NaN: none in that scan).
T = 360 / rate;
h = cfg.ant.height_m;
tt = 0:0.02:tEnd;
f = zeros(size(tt));
for k = 1:numel(tt)
    [~, az] = target_geometry(tg, tt(k), h);
    f(k) = wrap180(cfg.ant.az0_deg + rate * tt(k) - az);      % antenna ahead of the target: > 0
end
tl = nan(1, round(tEnd / T));
for k = find(f(1:end-1) < 0 & f(2:end) >= 0)                   % the antenna passes the target
    a = tt(k);  b = tt(k+1);
    for it = 1:30                                              % bisection
        c = (a + b) / 2;
        [~, az] = target_geometry(tg, c, h);
        if wrap180(cfg.ant.az0_deg + rate * c - az) < 0, a = c; else, b = c; end
    end
    s = floor(a / T) + 1;
    if s <= numel(tl) && ~isfinite(tl(s)), tl(s) = a; end
end
end

function [tr, modes] = run_tracker(p, pl, tEnd, cfg)
% The tracker on the plots; with cfg (PMBM) also the C tracker on the same calls (when built).
tr = track_init(p);
[tr, modes, calls] = feed(tr, pl, tEnd);
if nargin > 3 && strcmp(p.method, 'pmbm')
    c_compare('add', compare_with_c(cfg, calls, tr.log));
end
end

function [tr, modes, calls] = feed(tr, pl, tEnd)
% The plots reach the tracker 0.1 s after their time; the tracker is called every 0.1 s.
% modes: [time, model probabilities] of the confirmed tracks after every update (PMBM).
% calls: one row per plot and call [callTime flush time az range el vr] (time NaN: no plot).
modes = zeros(0, 4);
calls = zeros(0, 7);
k = 1;
nLog = 0;
for t = 0:0.1:tEnd + 0.05
    new = [];
    while k <= numel(pl) && pl(k).time_s + 0.1 <= t + 1e-9
        new = [new, pl(k)]; %#ok<AGROW>
        k = k + 1;
    end
    if isempty(new)
        calls(end+1, :) = [t 0 NaN NaN NaN NaN NaN]; %#ok<AGROW>
    end
    for q = 1:numel(new)
        calls(end+1, :) = [t 0 new(q).time_s new(q).az_deg new(q).range_m new(q).el_deg new(q).vr_mps]; %#ok<AGROW>
    end
    tr = track_update(tr, new, t);
    if isfield(tr, 'T') && size(tr.log, 1) > nLog
        ids = tr.log(nLog+1:end, 2);
        for j = 1:numel(tr.tracks)
            if tr.tracks(j).confirmed && any(ids == tr.tracks(j).id) && numel(tr.tracks(j).mu) == 3
                modes(end+1, :) = [tr.tracks(j).t, tr.tracks(j).mu]; %#ok<AGROW>
            end
        end
        nLog = size(tr.log, 1);
    end
end
tr = track_update(tr, [], tEnd + 0.1, true);
calls(end+1, :) = [tEnd + 0.1, 1, NaN NaN NaN NaN NaN];
end

function ok = compare_with_c(cfg, calls, log)
% The C tracker (radar_receiver_<os>/build/test_tracker, tracker_pmbm.c) on the
% same calls: every track update must match the MATLAB log (when built).
ok = true;
bin = fullfile(receiver_dir(), 'build', 'test_tracker');
if ispc, bin = [bin '.exe']; end
if ~exist(bin, 'file'), ok = []; return; end
folder = tempname();
mkdir(folder);
fp = fullfile(folder, 'calls.csv');  fl = fullfile(folder, 'log.csv');  fc = fullfile(folder, 'config.bin');
write_csv(fp, calls, [repmat('%.17g,', 1, 6) '%.17g\n']);
write_csv(fl, log, [repmat('%.17g,', 1, 10) '%.17g\n']);
evalc('export_receiver_config(cfg, fc)');
[status, out] = system(sprintf('"%s" "%s" "%s" "%s" %.17g', bin, fc, fp, fl, calls(end, 1)));
ok = status == 0;
if ~ok, fprintf('    C tracker: %s', regexprep(out, '\n(?=.)', '\n    ')); end
c_compare('text', out);
delete(fp);  delete(fl);  delete(fc);  rmdir(folder);
end

function varargout = c_compare(action, v)
% Results of the comparisons with the C tracker: 'reset', 'add' (ok, [] = not built),
% 'text' (its last output), 'get' -> [nCompared, nFailed, lastText].
persistent n nBad last
if isempty(n), n = 0;  nBad = 0;  last = ''; end
switch action
    case 'reset', n = 0;  nBad = 0;  last = '';
    case 'add'
        if ~isempty(v), n = n + 1;  nBad = nBad + ~v; end
    case 'text', last = v;
    case 'get', varargout = {n, nBad, last};
end
end

function write_csv(file, M, fmt)
fid = fopen(file, 'w');
fprintf(fid, fmt, M.');
fclose(fid);
end

function [pt, vt] = truth_xyz(tg, t, cfg, dim)
% True position relative to the antenna (2-D: at slant range in the horizontal plane, as the plots).
[pos, vt] = target_state(tg, t);
pt = pos - [0 0 cfg.ant.height_m];
if dim == 2
    [R, az] = target_geometry(tg, t, cfg.ant.height_m);
    pt = [R*sind(az), R*cosd(az), 0];
end
end

function [ok, txt] = check_targets(cfg, tr, tgts, tEnd, limPos, limVel)
% Every target: a confirmed track within limPos [m] and limVel [m/s] (2-D: slant range).
rep = track_report(tr);
conf = rep([rep.confirmed]);
dim = tr.p.dim;
ok = true;
eP = zeros(1, numel(tgts));  eV = eP;
for i = 1:numel(tgts)
    best = Inf;  bv = Inf;
    for j = 1:numel(conf)
        [pt, vt] = truth_xyz(tgts(i), conf(j).t_s, cfg, dim);
        dd = norm(conf(j).pos_m(1:2) - pt(1:2));
        if dd < best, best = dd; bv = norm(conf(j).vel_mps(1:2) - vt(1:2)); end
    end
    eP(i) = best;  eV(i) = bv;
    ok = ok && best < limPos && bv < limVel;
end
txt = sprintf('%d confirmed tracks for %d targets, largest position error %.0f m, velocity error %.1f m/s', ...
              numel(conf), numel(tgts), max(eP), max(eV));
ok = ok && numel(conf) >= numel(tgts);
if tEnd < 0, ok = false; end
end

function n = false_tracks(cfg, tr, tgts)
% Confirmed tracks (in the log) whose last confirmed update is more than 2 km from every target.
L = tr.log;
n = 0;
for id = unique(L(L(:, 3) == 1, 2)).'
    rows = L(L(:, 2) == id & L(:, 3) == 1, :);
    last = rows(end, :);
    d = Inf;
    for k = 1:numel(tgts)
        pt = truth_xyz(tgts(k), last(1), cfg, tr.p.dim);
        d = min(d, norm(last(4:5) - pt(1:2)));
    end
    n = n + (d > 2000);
end
end

function has = on_target(cfg, tr, tg, tEnd, dim)
% A confirmed track within 1 km of the target, updated in the last scan.
has = false;
rep = track_report(tr);
for j = find([rep.confirmed])
    if rep(j).t_s < tEnd - tr.p.scanTime_s, continue; end
    pt = truth_xyz(tg, rep(j).t_s, cfg, dim);
    has = has || norm(rep(j).pos_m(1:2) - pt(1:2)) < 1000;
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

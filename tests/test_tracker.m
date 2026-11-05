function ok = test_tracker()
%TEST_TRACKER  EKF tracker (track_update.m) on synthetic plots, 3-D and 2-D, CV and CA models.
%   1. The six targets of six_targets_scenario.m are seen once per scan (when
%      the antenna passes them) for 5 scans, with the plot errors of
%      cfg.track (range, azimuth, elevation), plus 8 random false plots per
%      scan. Every target must end with one confirmed track close to its true
%      position and velocity, and no false plot may give a confirmed track.
%      3-D: range, azimuth, elevation; 2-D: range and azimuth (the elevation
%      is unknown, as with a single element, so the track lies at slant
%      range). Both motion models (constant velocity and constant
%      acceleration; from 5 plots the CA model estimates the velocity of
%      these unaccelerated targets less precisely: wider limits). The plots
%      carry the radial velocity with its error (cfg.track.sigmaVr_mps; false
%      plots a random one; drawn from a stream of their own, so the position
%      errors are those of the case without it). Every case runs with and
%      without the radial velocity measurement (cfg.track.useRadialVelocity):
%      with it the radial velocity of the tracks after 2 and after 5 scans and
%      their velocity after 2 scans must be closer to the truth, and it must
%      not confirm more false tracks. With the staggered PRF the plot
%      velocity is folded (+-112 m/s at 170 MHz): the tracks must still have
%      the true radial velocity (unfolded with the predicted range rate).
%   2. GNN: two targets 200 m apart in range at the same azimuth. At the
%      second scan the first plot is closer to the other (still tentative)
%      track: a nearest-neighbour association would swap them; the global
%      assignment must give every track its own target.
%   3. A target accelerating at 4 m/s^2 for 12 scans: the constant-
%      acceleration model must follow it (and estimate the acceleration),
%      closer than the constant-velocity model. With the radial velocity the
%      CV track must be closer than without it (this target accelerates
%      nearly along the line of sight) and CA still better in velocity.
%   4. Speed gate (cfg.track.minSpeedEnable, minSpeed_mps = 5 m/s): a target
%      at 3 m/s and one at 100 m/s, plots without errors, 4 scans. Without
%      the gate both have a track; with it the slow track is deleted at its
%      second update (a track with one update is kept: its velocity is not
%      measured yet), the fast one is kept (as the C test_tracker).
%   When the C receiver is built, radar_receiver_<os>/build/test_tracker runs the
%   C tracker on the plots of case 1: every track update must match MATLAB.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
for model = {'cv', 'ca'}
    for dim = [3 2]
        [okVr, eVr] = run_case(dim, model{1}, true);
        [okNo, eNo] = run_case(dim, model{1}, false);
        better = eVr.vel(1) < eNo.vel(1) && all(eVr.vr < eNo.vr) && eVr.nFalse <= eNo.nFalse;
        fprintf(['  %d-D %s, radial velocity used / not used: velocity error after 2 scans %.1f / %.1f m/s, ', ...
                 'after 5 scans %.1f / %.1f m/s; radial velocity error after 2 scans %.1f / %.1f m/s, after 5 ', ...
                 'scans %.1f / %.1f m/s; false confirmed tracks %d / %d: %s\n'], dim, upper(model{1}), ...
                eVr.vel(1), eNo.vel(1), eVr.vel(2), eNo.vel(2), eVr.vr(1), eNo.vr(1), eVr.vr(2), eNo.vr(2), ...
                eVr.nFalse, eNo.nFalse, ternary(better, 'better', 'FAIL'));
        ok = okVr && okNo && better && ok;
    end
end
% staggered PRF: the plot velocity is folded into +-vrSpan/2 (225 m/s at 170 MHz); the tracker unfolds it
[okS, eS] = run_case(3, 'cv', true, 'staggered');
good = okS && eS.vr(2) < 5;
fprintf(['  3-D CV, staggered PRF (plot velocity folded): radial velocity error of the tracks after 5 scans ', ...
         '%.1f m/s (unambiguous): %s\n'], eS.vr(2), ternary(good, 'ok', 'FAIL'));
ok = good && ok;
ok = gnn_case() && ok;
ok = maneuver_case() && ok;
ok = min_speed_case() && ok;
fprintf('EKF tracker: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [ok, e] = run_case(dim, model, useVr, prfMode)
% e.vel, e.vr: mean velocity / radial velocity error of the tracks on the
% targets after 2 and after 5 scans [m/s]; e.nFalse: confirmed tracks not on a target
cfg = radar_config();
cfg.console.enable = false;
cfg.track.dim = sprintf('%dd', dim);
cfg.track.model = model;
cfg.track.useRadialVelocity = useVr;
if nargin > 3, cfg.prf.mode = prfMode; end
if dim == 2, cfg.ant.nElements = 1; end              % no elevation
sys = setup_radar(cfg);
tgts = six_targets_scenario();
p = sys.track;
rng(99);                                              % radial velocity errors: a stream of their own
vrNoise = randn(5, numel(tgts));
vrFalse = rand(5, 8);
rng(7);
tr = track_init(p);
T = p.scanTime_s;
allPlots = zeros(0, 5);                               % plots in processing order (for the C tracker)
rate = sys.ant.rate_dps;
nScans = 5;
sigmaVr = sqrt(p.Rvr);
e = struct('vel', [0 0], 'vr', [0 0], 'nFalse', 0);
for scan = 0:nScans-1
    pl = [];
    for i = 1:numel(tgts)                             % time when the antenna points at the target
        t = scan*T;
        for it = 1:3
            [~, az] = target_geometry(tgts(i), t, cfg.ant.height_m);
            t = scan*T + mod(az - cfg.ant.az0_deg, 360) / rate;
        end
        [R, az, el, vr] = target_geometry(tgts(i), t, cfg.ant.height_m);
        q = make_plot(t, R + cfg.track.sigmaRange_m*randn, az + cfg.track.sigmaAz_deg*randn, ...
                      el + cfg.track.sigmaEl_deg*randn, dim, ...
                      mod(vr + sigmaVr*vrNoise(scan+1, i) + p.vrSpan/2, p.vrSpan) - p.vrSpan/2);   % as measured
        pl = [pl, q]; %#ok<AGROW>
    end
    for f = 1:8                                       % false plots
        az = 360*rand;
        pl = [pl, make_plot(scan*T + az/rate, 5e3 + 45e3*rand, az, 20*rand, dim, p.vrSpan*(vrFalse(scan+1, f) - 0.5))]; %#ok<AGROW>
    end
    [~, order] = sort([pl.time_s]);
    for k = order                                     % plots in time order, maintenance every plot
        tr = track_update(tr, pl(k), pl(k).time_s);
        allPlots(end+1, :) = [pl(k).time_s, pl(k).az_deg, pl(k).range_m, pl(k).el_deg, pl(k).vr_mps]; %#ok<AGROW>
    end
    if scan == 1, [e.vel(1), e.vr(1)] = velocity_error(track_report(tr), tgts, cfg); end
end
tEnd = nScans*T;
tr = track_update(tr, [], tEnd, true);              % associate the last batch
rep = track_report(tr);
conf = rep([rep.confirmed]);
ok = numel(conf) == numel(tgts);
[e.vel(2), e.vr(2), e.nFalse] = velocity_error(rep, tgts, cfg);
fprintf('  %d-D %s, radial velocity %s: %d tracks, %d confirmed (expected %d), %d updates\n', dim, upper(model), ...
        ternary(useVr, 'used', 'not used'), numel(rep), numel(conf), numel(tgts), size(tr.log, 1));
for i = 1:numel(tgts)
    best = inf;  bi = 0;
    for j = 1:numel(conf)
        pos0 = tgts(i).pos0_m + tgts(i).vel_mps*conf(j).t_s - [0 0 cfg.ant.height_m];
        dist = norm(conf(j).pos_m(1:2) - pos0(1:2));
        if dist < best, best = dist; bi = j; end
    end
    if bi == 0, ok = false; continue; end
    tk = conf(bi);
    pos0 = tgts(i).pos0_m + tgts(i).vel_mps*tk.t_s;
    if dim == 2                       % without elevation the track is at slant range (as the plots)
        [R, az] = target_geometry(tgts(i), tk.t_s, cfg.ant.height_m);
        best = norm(tk.pos_m(1:2) - R*[sind(az) cosd(az)]);
    end
    dv = norm(tk.vel_mps(1:2) - tgts(i).vel_mps(1:2));
    dh = NaN;
    if dim == 3, dh = tk.alt_m - pos0(3); end
    if strcmp(model, 'cv'), lim = [250 25]; else, lim = [300 45]; end   % CA: 3 unknowns per axis
    good = best < lim(1) && dv < lim(2) && (dim == 2 || abs(dh) < 400);
    fprintf('    %s: track %d, position error %5.0f m, velocity error %5.1f m/s, height error %5.0f m %s\n', ...
            tgts(i).name, tk.id, best, dv, dh, ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
ok = compare_with_c(cfg, allPlots, tr.log, tEnd) && ok;
end

function [eVel, eVr, nFalse] = velocity_error(rep, tgts, cfg)
% Mean over the targets of the horizontal velocity error and of the radial
% velocity error of the track nearest each target; confirmed tracks more than
% 2 km from every target.
eVel = 0;  eVr = 0;
for i = 1:numel(tgts)
    best = inf;  dv = NaN;  dr = NaN;
    for j = 1:numel(rep)
        pos0 = tgts(i).pos0_m + tgts(i).vel_mps*rep(j).t_s - [0 0 cfg.ant.height_m];
        dist = norm(rep(j).pos_m(1:2) - pos0(1:2));
        if dist < best
            best = dist;
            dv = norm(rep(j).vel_mps(1:2) - tgts(i).vel_mps(1:2));
            [~, ~, ~, vr] = target_geometry(tgts(i), rep(j).t_s, cfg.ant.height_m);
            dr = abs(rep(j).vr_mps - vr);
        end
    end
    eVel = eVel + dv / numel(tgts);
    eVr = eVr + dr / numel(tgts);
end
nFalse = 0;
for j = find([rep.confirmed])
    d = arrayfun(@(i) norm(rep(j).pos_m(1:2) - (tgts(i).pos0_m(1:2) + tgts(i).vel_mps(1:2)*rep(j).t_s)), 1:numel(tgts));
    nFalse = nFalse + (min(d) > 2000);
end
end

function ok = compare_with_c(cfg, plots, log, tEnd)
% The C tracker on the same plots (when built): it must give the same track updates.
ok = true;
root = fileparts(fileparts(mfilename('fullpath')));
bin = fullfile(receiver_dir(), 'build', 'test_tracker');
if ispc, bin = [bin '.exe']; end
if ~exist(bin, 'file')
    fprintf('    C tracker not built: comparison skipped\n');
    return
end
folder = tempname();
mkdir(folder);
fp = fullfile(folder, 'plots.csv');  fl = fullfile(folder, 'log.csv');  fc = fullfile(folder, 'config.bin');
write_csv(fp, plots, '%.17g,%.17g,%.17g,%.17g,%.17g\n');
write_csv(fl, log, [repmat('%.17g,', 1, 10) '%.17g\n']);
evalc('export_receiver_config(cfg, fc)');            % dimension, model, gate, batch of this case
[status, out] = system(sprintf('"%s" "%s" "%s" "%s" %.17g', bin, fc, fp, fl, tEnd));
fprintf('    C tracker: %s', regexprep(out, '\n(?=.)', '\n    '));
ok = status == 0;
delete(fp);  delete(fl);  delete(fc);  rmdir(folder);
end

function ok = gnn_case()
% Two targets 200 m apart in range, receding at 100 m/s, plots without errors.
% Scan 2: plot 1 is 800 m from track 2 and 1000 m from track 1 (both still at
% the first plot, velocity unknown): nearest neighbour would give it to track 2.
cfg = radar_config();
cfg.console.enable = false;
cfg.ant.nElements = 1;
cfg.track.dim = '2d';
sys = setup_radar(cfg);
p = sys.track;
T = p.scanTime_s;
az = 30;
R0 = [20e3, 20.2e3];
v = 100;
tr = track_init(p);
for scan = 0:5
    t = scan*T;
    pl = [make_plot(t, R0(1) + v*t, az, NaN, 2, v), make_plot(t, R0(2) + v*t, az, NaN, 2, v)];
    tr = track_update(tr, pl, t, true);
end
rep = track_report(tr);
ok = numel(rep) == 2 && all([rep.confirmed]) && all([rep.hits] == 6);
for k = 1:numel(rep)
    R = norm(rep(k).pos_m(1:2));
    i = 1 + (abs(R - R0(2) - v*rep(k).t_s) < abs(R - R0(1) - v*rep(k).t_s));
    dv = norm(rep(k).vel_mps(1:2) - v*[sind(az) cosd(az)]);
    ok = ok && abs(R - R0(i) - v*rep(k).t_s) < 5 && dv < 1 && rep(k).id == i;
end
fprintf('  GNN: two targets 200 m apart, %d tracks, updates %s, every track on its own target: %s\n', ...
        numel(rep), mat2str([rep.hits]), ternary(ok, 'ok', 'FAIL'));
end

function ok = maneuver_case()
% A target accelerating along its track: constant-velocity and constant-acceleration
% models, without and with the radial velocity measurement (same plot errors).
ok = true;
errNo = [];
for useVr = [false true]
    rng(99);                                          % radial velocity errors: a stream of their own
    vrNoise = randn(2, 12);
    rng(11);
    err = zeros(2, 3);
    names = {'cv', 'ca'};
    for mi = 1:2
        cfg = radar_config();
        cfg.console.enable = false;
        cfg.ant.nElements = 1;
        cfg.track.dim = '2d';
        cfg.track.model = names{mi};
        cfg.track.useRadialVelocity = useVr;
        sys = setup_radar(cfg);
        p = sys.track;
        T = p.scanTime_s;
        p0 = [5e3 20e3];  v0 = [150 0];  a0 = [0 4];    % m, m/s, m/s^2
        tr = track_init(p);
        nScans = 12;
        for scan = 0:nScans-1
            t = scan*T;
            pos = p0 + v0*t + a0*t^2/2;
            R = norm(pos) + cfg.track.sigmaRange_m*randn;
            az = atan2(pos(1), pos(2))*180/pi + cfg.track.sigmaAz_deg*randn;
            vr = pos*(v0 + a0*t).'/norm(pos) + sqrt(p.Rvr)*vrNoise(mi, scan+1);
            tr = track_update(tr, make_plot(t, R, az, NaN, 2, vr), t, true);
        end
        rep = track_report(tr);
        [~, k] = max([rep.hits]);
        t = rep(k).t_s;
        err(mi, :) = [norm(rep(k).pos_m(1:2) - (p0 + v0*t + a0*t^2/2)), norm(rep(k).vel_mps(1:2) - (v0 + a0*t)), ...
                      norm(rep(k).acc_mps2(1:2) - a0)];
    end
    if ~useVr
        % without the radial velocity the CV model lags behind: CA is closer in position and velocity
        good = err(2, 1) < err(1, 1) && err(2, 2) < 10 && err(2, 2) < err(1, 2) && err(2, 3) < 1.5;
        errNo = err;
    else
        % the measured radial velocity keeps also the CV track close (this target accelerates nearly
        % along the line of sight); CA still estimates velocity and acceleration better
        good = err(2, 2) < 10 && err(2, 2) < err(1, 2) && err(2, 3) < 1.5 && err(2, 1) < 200 && ...
               err(1, 1) < errNo(1, 1) && err(1, 2) < errNo(1, 2);
    end
    fprintf(['  Maneuver (4 m/s^2), radial velocity %s: CV position error %.0f m, velocity error %.1f m/s; ' ...
             'CA %.0f m, %.1f m/s, acceleration error %.2f m/s^2: %s\n'], ternary(useVr, 'used', 'not used'), ...
            err(1, 1:2), err(2, :), ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
end

function ok = min_speed_case()
% A slow (3 m/s) and a fast (100 m/s) target receding along the line of sight, without
% and with the speed gate of 5 m/s.
ok = true;
speed = [3 100];  az = [30 120];  R0 = 20e3;
for gate = [false true]
    cfg = radar_config();
    cfg.console.enable = false;
    cfg.ant.nElements = 1;
    cfg.track.dim = '2d';
    cfg.track.minSpeedEnable = gate;
    cfg.track.minSpeed_mps = 5;
    sys = setup_radar(cfg);
    T = sys.track.scanTime_s;
    tr = track_init(sys.track);
    for scan = 0:3
        t = scan*T;
        pl = [make_plot(t, R0 + speed(1)*t, az(1), NaN, 2, speed(1)), make_plot(t, R0 + speed(2)*t, az(2), NaN, 2, speed(2))];
        tr = track_update(tr, pl, t, true);
        if scan == 0, afterFirst = numel(tr.tracks); end
    end
    rep = track_report(tr);
    v = arrayfun(@(r) norm(r.vel_mps(1:2)), rep);
    nSlow = sum(v < 50);  nFast = sum(v >= 50);
    good = afterFirst == 2 && nFast == 1 && nSlow == ~gate;
    fprintf('  Speed gate %s: %d tracks after the first scan, at the end %d slow and %d fast track(s): %s\n', ...
            ternary(gate, '5 m/s', 'off'), afterFirst, nSlow, nFast, ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
end

function write_csv(file, M, fmt)
fid = fopen(file, 'w');
fprintf(fid, fmt, M.');                               % NaN is written as NaN (read by strtod)
fclose(fid);
end

function q = make_plot(t, R, az, el, dim, vr)
if dim == 2, el = NaN; end
q = struct('time_s', t, 'az_deg', mod(az, 360), 'range_m', R, 'el_deg', el, 'vr_mps', vr, 'snr_dB', 20, ...
           'power', 1, 'beam', 1, 'rIdx', 1, 'dIdx', 1, 'pulse', 'long', 'hits', 60);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

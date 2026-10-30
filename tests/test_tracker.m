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
%      these unaccelerated targets less precisely: wider limits).
%   2. GNN: two targets 200 m apart in range at the same azimuth. At the
%      second scan the first plot is closer to the other (still tentative)
%      track: a nearest-neighbour association would swap them; the global
%      assignment must give every track its own target.
%   3. A target accelerating at 4 m/s^2 for 12 scans: the constant-
%      acceleration model must follow it (and estimate the acceleration),
%      closer than the constant-velocity model.
%   When the C receiver is built, radar_receiver/build/test_tracker runs the
%   C tracker on the plots of case 1: every track update must match MATLAB.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
for model = {'cv', 'ca'}
    for dim = [3 2]
        ok = run_case(dim, model{1}) && ok;
    end
end
ok = gnn_case() && ok;
ok = maneuver_case() && ok;
fprintf('EKF tracker: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ok = run_case(dim, model)
cfg = radar_config();
cfg.sim.verbose = false;
cfg.track.dim = sprintf('%dd', dim);
cfg.track.model = model;
if dim == 2, cfg.ant.nElements = 1; end              % no elevation
sys = setup_radar(cfg);
tgts = six_targets_scenario();
p = sys.track;
rng(7);
tr = track_init(p);
T = p.scanTime_s;
allPlots = zeros(0, 4);                               % plots in processing order (for the C tracker)
rate = sys.ant.rate_dps;
nScans = 5;
for scan = 0:nScans-1
    pl = [];
    for i = 1:numel(tgts)                             % time when the antenna points at the target
        t = scan*T;
        for it = 1:3
            [~, az] = target_geometry(tgts(i), t, cfg.ant.height_m);
            t = scan*T + mod(az - cfg.ant.az0_deg, 360) / rate;
        end
        [R, az, el] = target_geometry(tgts(i), t, cfg.ant.height_m);
        q = make_plot(t, R + cfg.track.sigmaRange_m*randn, az + cfg.track.sigmaAz_deg*randn, ...
                      el + cfg.track.sigmaEl_deg*randn, dim);
        pl = [pl, q]; %#ok<AGROW>
    end
    for f = 1:8                                       % false plots
        az = 360*rand;
        pl = [pl, make_plot(scan*T + az/rate, 5e3 + 45e3*rand, az, 20*rand, dim)]; %#ok<AGROW>
    end
    [~, order] = sort([pl.time_s]);
    for k = order                                     % plots in time order, maintenance every plot
        tr = track_update(tr, pl(k), pl(k).time_s);
        allPlots(end+1, :) = [pl(k).time_s, pl(k).az_deg, pl(k).range_m, pl(k).el_deg]; %#ok<AGROW>
    end
end
tEnd = nScans*T;
tr = track_update(tr, [], tEnd, true);              % associate the last batch
rep = track_report(tr);
conf = rep([rep.confirmed]);
ok = numel(conf) == numel(tgts);
fprintf('  %d-D %s: %d tracks, %d confirmed (expected %d), %d updates\n', dim, upper(model), numel(rep), ...
        numel(conf), numel(tgts), size(tr.log, 1));
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

function ok = compare_with_c(cfg, plots, log, tEnd)
% The C tracker on the same plots (when built): it must give the same track updates.
ok = true;
root = fileparts(fileparts(mfilename('fullpath')));
bin = fullfile(root, 'radar_receiver', 'build', 'test_tracker');
if ispc, bin = [bin '.exe']; end
if ~exist(bin, 'file')
    fprintf('    C tracker not built: comparison skipped\n');
    return
end
folder = tempname();
mkdir(folder);
fp = fullfile(folder, 'plots.csv');  fl = fullfile(folder, 'log.csv');  fc = fullfile(folder, 'config.bin');
write_csv(fp, plots, '%.17g,%.17g,%.17g,%.17g\n');
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
cfg.sim.verbose = false;
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
    pl = [make_plot(t, R0(1) + v*t, az, NaN, 2), make_plot(t, R0(2) + v*t, az, NaN, 2)];
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
% A target accelerating along its track: constant-velocity and constant-acceleration models.
rng(11);
err = zeros(2, 3);
names = {'cv', 'ca'};
for mi = 1:2
    cfg = radar_config();
    cfg.sim.verbose = false;
    cfg.ant.nElements = 1;
    cfg.track.dim = '2d';
    cfg.track.model = names{mi};
    sys = setup_radar(cfg);
    p = sys.track;
    T = p.scanTime_s;
    p0 = [5e3 20e3];  v0 = [150 0];  a0 = [0 4];        % m, m/s, m/s^2
    tr = track_init(p);
    nScans = 12;
    for scan = 0:nScans-1
        t = scan*T;
        pos = p0 + v0*t + a0*t^2/2;
        R = norm(pos) + cfg.track.sigmaRange_m*randn;
        az = atan2(pos(1), pos(2))*180/pi + cfg.track.sigmaAz_deg*randn;
        tr = track_update(tr, make_plot(t, R, az, NaN, 2), t, true);
    end
    rep = track_report(tr);
    [~, k] = max([rep.hits]);
    t = rep(k).t_s;
    err(mi, :) = [norm(rep(k).pos_m(1:2) - (p0 + v0*t + a0*t^2/2)), norm(rep(k).vel_mps(1:2) - (v0 + a0*t)), ...
                  norm(rep(k).acc_mps2(1:2) - a0)];
end
ok = err(2, 1) < err(1, 1) && err(2, 2) < 10 && err(2, 2) < err(1, 2) && err(2, 3) < 1.5;
fprintf(['  Maneuver (4 m/s^2): CV position error %.0f m, velocity error %.1f m/s; ' ...
         'CA %.0f m, %.1f m/s, acceleration error %.2f m/s^2: %s\n'], err(1, 1:2), err(2, :), ternary(ok, 'ok', 'FAIL'));
end

function write_csv(file, M, fmt)
fid = fopen(file, 'w');
fprintf(fid, fmt, M.');                               % NaN is written as NaN (read by strtod)
fclose(fid);
end

function q = make_plot(t, R, az, el, dim)
if dim == 2, el = NaN; end
q = struct('time_s', t, 'az_deg', mod(az, 360), 'range_m', R, 'el_deg', el, 'vr_mps', 0, 'snr_dB', 20, ...
           'power', 1, 'beam', 1, 'rIdx', 1, 'dIdx', 1, 'pulse', 'long', 'hits', 60);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

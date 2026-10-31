function ok = test_vhf_challenges()
%TEST_VHF_CHALLENGES  Demanding cases of the 170 MHz radar (one CPI each, clutter on).
%   At 170 MHz the echoes are 25 dB stronger than at 3 GHz (wavelength
%   squared): strong targets, strong clutter and a coarse Doppler axis
%   (6.9 m/s per bin, +-441 m/s unambiguous with the fixed PRF).
%   1. ADC overload: with the 3 GHz receiver gain (60 dB) adc_overload.m finds
%      far more clipped receive samples than cfg.adc.overloadWarn (and
%      run_radar_sim warns); the default gain stays below it.
%   2. Strong targets near the radar (5, 12, 20 km, up to 10 m^2): exactly
%      one detection, on the target cell, and no long-pulse detection in the
%      zone that their eclipsed long-pulse echo reaches (eclipse_blanking.m;
%      the number of ghosts it removed is printed).
%   3. A target behind a strong near target, on the same azimuth (other
%      Doppler): still detected in the long pulse.
%   4. Fixed PRF, velocities 5 (inside the MTI notch: blind), 30, 200 and
%      430 m/s on one azimuth: all but the blind one are detected, with the
%      radial velocity within one Doppler bin (no ambiguity at 170 MHz).
%   5. Staggered PRF: the de-interleaved rows are sampled every sum(PRI), so
%      a target at that row blind speed (225 m/s) is in the notch: snr_budget
%      flags it and it is not detected; a target at 150 m/s is (at its
%      velocity folded into the row interval).
%   6. Two targets 300 m apart on the same azimuth (other Doppler): two detections.
%   7. Pruning under load: a map with about 33000 CFAR detections is pruned
%      in well under a second, with exactly the result of the direct
%      pairwise rules.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
cfg = radar_config();
cfg.sim.verbose = false;  cfg.display.enable = false;  cfg.io.saveRaw = false;
cfg.rx.maxRange_m = 50e3;
fprintf('fc = %.0f MHz, receiver gain %g dB\n', cfg.fc/1e6, cfg.rx.gain_dB);

%% 1) ADC overload
for g = [60, cfg.rx.gain_dB]
    c = cfg;  c.rx.gain_dB = g;
    sys = setup_radar(c);
    rng(1);
    clut = clutter_init(sys);
    frac = 0;
    for az = 0:45:315
        frac = max(frac, adc_overload(adc_sample(generate_rx_signal(sys, [], 0, az, clutter_signal(clut, sys, az, 0)), sys), sys));
    end
    good = (g == 60) == (frac > c.adc.overloadWarn);
    ok = report(sprintf('receiver gain %2g dB: up to %.2f%% of the receive samples clipped (warning above %.0f%%)', ...
                        g, 100*frac, 100*c.adc.overloadWarn), good) && ok;
end

sys = setup_radar(cfg);
h = cfg.ant.height_m;

%% 2) strong targets near the radar: one detection, no eclipse ghosts
cases = {10, 5e3, -150;  1, 12e3, -236;  1, 20e3, 129};          % RCS, range, radial velocity
for k = 1:size(cases, 1)
    tg = make_target('N', cases{k, 1}, cases{k, 2}, 40, 1000, cases{k, 3}, 0);
    [r, nGhost] = cpi_reports(cfg, sys, tg, 40, true);
    [Rt, ~, ~, vr] = target_geometry(tg, r.t0, h);
    s = r.short;  l = r.long;
    zone = [l.range_m] <= Rt + sys.eclipseExtent_m;
    good = numel(s) == 1 && abs(s(1).range_m - Rt) < 75 && abs(s(1).vr_mps - vr) < 7 && ~any(zone);
    ok = report(sprintf(['%4.1f m^2 at %4.1f km, %4.0f m/s: %d short-pulse detection(s), %d long-pulse detection(s) in its ', ...
                         'eclipse zone (%d ghosts removed by the eclipse blanking)'], cases{k, 1}, Rt/1e3, vr, ...
                        numel(s), nnz(zone), nGhost), good) && ok;
end

%% 3) a target behind a strong near target, same azimuth
tg = [make_target('N', 1, 12e3, 40, 1000, -236, 0), make_target('F', 1, 30e3, 40, 3000, 150, 0)];
r = cpi_reports(cfg, sys, tg, 40, true);
[Rf, ~, ~, vf] = target_geometry(tg(2), r.t0, h);
hit = abs([r.long.range_m] - Rf) < 75 & abs([r.long.vr_mps] - vf) < 7;
ok = report(sprintf('1 m^2 at 30 km behind 1 m^2 at 12 km (same azimuth): detected in the long pulse: %d', any(hit)), ...
            any(hit)) && ok;

%% 4) fixed PRF: blind speed and unambiguous velocities
vel = [5 30 200 430];
tg = [];
for k = 1:numel(vel)
    tg = [tg, make_target(sprintf('V%d', k), 1, (12 + 7*k)*1e3, 40, 1500, vel(k), 0)]; %#ok<AGROW>
end
r = cpi_reports(cfg, sys, tg, 40, true);
txt = '';  good = true;
for k = 1:numel(tg)
    [Rt, ~, ~, vr] = target_geometry(tg(k), r.t0, h);
    det = [r.short, r.long];
    found = any(abs([det.range_m] - Rt) < 75 & abs([det.vr_mps] - vr) < abs(diff(sys.dop.v_mps(1:2))));
    expected = abs(vr) > (sys.dop.notchBins + 0.5) * abs(diff(sys.dop.v_mps(1:2)));
    good = good && found == expected;
    txt = [txt, sprintf(' %.0f m/s: %s', vr, ternary(found, 'detected', 'not detected'))]; %#ok<AGROW>
end
ok = report(sprintf('fixed PRF (bin %.1f m/s, notch +-%d bins):%s', abs(diff(sys.dop.v_mps(1:2))), sys.dop.notchBins, txt), ...
            good) && ok;

%% 5) staggered PRF: blind speed of the de-interleaved rows
cs = cfg;  cs.prf.mode = 'staggered';
ss = setup_radar(cs);
vBlind = ss.lambda / (2 * ss.dop.rowPri_s);
tg = [make_target('B', 1, 20e3, 40, 1500, vBlind, 0), make_target('S', 1, 30e3, 40, 1500, 150, 0)];
r = cpi_reports(cs, ss, tg, 40, true);
det = [r.short, r.long];
bud = [snr_budget(ss, tg(1), r.t0), snr_budget(ss, tg(2), r.t0)];
[Rb, ~, ~, vb] = target_geometry(tg(1), r.t0, h);
[Rs, ~, ~, vs] = target_geometry(tg(2), r.t0, h);
vSpan = 2 * ss.lambda / (4 * ss.dop.rowPri_s);
fold = @(v) mod(v + vSpan/2, vSpan) - vSpan/2;
foundB = any(abs([det.range_m] - Rb) < 75);
foundS = any(abs([det.range_m] - Rs) < 75 & abs(fold([det.vr_mps] - vs)) < abs(diff(ss.dop.v_mps(1:2))));
good = bud(1).blind && ~bud(2).blind && ~foundB && foundS;
ok = report(sprintf(['staggered PRF: %.0f m/s (row blind speed %.0f m/s) flagged blind %d, detected %d; ', ...
                     '%.0f m/s detected %d'], vb, vBlind, bud(1).blind, foundB, vs, foundS), good) && ok;

%% 6) two targets 300 m apart
tg = [make_target('A', 1, 30e3, 40, 3000, 100, 0), make_target('B', 1, 30.3e3, 40, 3000, -60, 0)];
r = cpi_reports(cfg, sys, tg, 40, true);
n = 0;
for k = 1:2
    [Rt, ~, ~, vr] = target_geometry(tg(k), r.t0, h);
    n = n + any(abs([r.long.range_m] - Rt) < 75 & abs([r.long.vr_mps] - vr) < 7);
end
ok = report(sprintf('two targets 300 m apart on one azimuth: %d of 2 detected', n), n == 2) && ok;

%% 7) pruning under load
rng(5);
wv = sys.wave.long;
P = single(-log(rand(sys.nBeams, wv.nRange, sys.dop.nfft)));
P(3, 100, 40) = 1e7;  P(3, 60:140, 40) = max(P(3, 60:140, 40), 1e4);   % strong target with range sidelobes
isDet = P > 3.4;                                                        % about 20000 cells
noise = ones(size(P), 'single');
t = tic;  d = prune_peaks(isDet, P, noise, sys, wv);  dt = toc(t);
ref = prune_direct(isDet, P, sys, wv);
good = isequal(d.beam, ref(:, 1)) && isequal(d.rIdx, ref(:, 2)) && isequal(d.dIdx, ref(:, 3)) && dt < 2;
ok = report(sprintf('pruning of %d CFAR detections: %d peaks kept in %.2f s, same as the direct rules: %d', ...
                    nnz(isDet), numel(d.beam), dt, isequal(d.beam, ref(:, 1))), good) && ok;

fprintf('170 MHz challenges: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [r, nGhost] = cpi_reports(cfg, sys, targets, az0, clutterOn)
% One CPI centred on azimuth az0: detection reports of both sub-pulses (as
% run_radar_sim) and the number of long-pulse detections the eclipse
% blanking removed.
rng(11);
targets = draw_target_rcs(targets, 2, 0);
clut = [];
if clutterOn, clut = clutter_init(sys); end
n = sys.cpi.length;
r.t0 = (az0 - cfg.ant.az0_deg) / sys.ant.rate_dps;
names = {'short', 'long'};
for i = 1:2, buf.(names{i}) = cpi_buffer_init(sys.nBeams, sys.wave.(names{i}).nRange, n); end
st = mod(0:n-1, sys.nStagger) + 1;
dt = reshape(sys.pri(st(1:end-1)), 1, []);                     % PRI before every pulse
tp = r.t0 + [0, cumsum(dt)] - sum(dt) / 2;                      % pulse times, centred on t0
for k = 1:n
    az = mod(cfg.ant.az0_deg + sys.ant.rate_dps * tp(k), 360);
    ci = [];
    if clutterOn, ci = clutter_signal(clut, sys, az, tp(k)); end
    [I, Q] = ddc(adc_sample(generate_rx_signal(sys, targets, tp(k), az, ci), sys), sys);
    bb = complex(single(I), single(Q));
    for i = 1:2
        wv = sys.wave.(names{i});
        seg = ddc_segment(bb, wv);
        if sys.dbf.enabled, bm = digital_beamforming(seg, sys.dbf.W); else, bm = seg.'; end
        buf.(names{i}) = cpi_buffer_push(buf.(names{i}), matched_filter(bm, wv), st(k), az, tp(k));
    end
end
nGhost = 0;
for i = 1:2
    wv = sys.wave.(names{i});
    [cpi, info] = cpi_buffer_read(buf.(names{i}));
    rd = doppler_process(mti_filter(cpi, sys.mtiCoeffs, info.staggerIdx), sys);
    [isDet, noise] = cfar_detect(rd, sys, wv);
    dets = prune_peaks(isDet, rd, noise, sys, wv);
    if i == 1
        shortDets = dets;
    else
        n0 = numel(dets.beam);
        dets = eclipse_blanking(dets, shortDets, sys);
        nGhost = n0 - numel(dets.beam);
    end
    r.(names{i}) = detection_reports(dets, rd, sys, wv, info);
end
end

function kept = prune_direct(isDet, P, sys, wv)
% The pruning rules of prune_peaks.m applied pair by pair: [beam rIdx dIdx] of the kept peaks.
[nB, nR, nD] = size(P);
idx = find(isDet);
[b, r, d] = ind2sub([nB nR nD], idx);
peak = false(size(idx));
for i = 1:numel(idx)
    blk = P(max(b(i)-1, 1):min(b(i)+1, nB), max(r(i)-1, 1):min(r(i)+1, nR), mod(d(i)-2:d(i), nD) + 1);
    peak(i) = P(idx(i)) >= max(blk(:));
end
b = b(peak);  r = r(peak);  d = d(peak);  pw = double(P(idx(peak)));
ratio = 10^(sys.cfg.prune.sidelobeMargin_dB/10);
keep = true(size(pw));
for i = 1:numel(pw)
    for j = find(pw > ratio*pw(i)).'
        dDop = abs(d(i) - d(j));  sameD = min(dDop, nD - dDop) <= 1;
        sameB = abs(b(i) - b(j)) <= 1;  sameR = abs(r(i) - r(j)) <= 1;
        if (sameB && sameD && abs(r(i) - r(j)) <= numel(wv.mfTaps)) || (sameR && sameD) || (sameR && sameB)
            keep(i) = false;
            break
        end
    end
end
kept = [b(keep), r(keep), d(keep)];
end

function ok = report(text, ok)
fprintf('  %-118s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

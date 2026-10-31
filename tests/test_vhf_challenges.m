function ok = test_vhf_challenges()
%TEST_VHF_CHALLENGES  Demanding cases of the 170 MHz radar (one CPI each, clutter on).
%   At 170 MHz the echoes are 25 dB stronger than at 3 GHz (wavelength
%   squared): strong targets, strong clutter and a coarse Doppler axis
%   (6.9 m/s per bin, +-441 m/s unambiguous with the fixed PRF).
%   1. ADC overload: at the receiver gain of 60 dB without STC adc_overload.m
%      finds far more clipped receive samples than cfg.adc.overloadWarn (and
%      run_radar_sim warns); with the STC of radar_config (stc_gain.m) no
%      sample is clipped, nor at 30 dB without STC.
%   2. Strong targets near the radar (5, 8, 12, 20 km, up to 10 m^2):
%      exactly one detection, on the target cell, and no long-pulse
%      detection at the target velocity in the zone that their eclipsed
%      long-pulse echo reaches (up to the target range + the long-pulse
%      blind range): the long-window STC (cfg.rx.stc.curveLong) attenuates
%      the eclipsed echo. The number of such ghosts without the STC (at 30 dB)
%      is printed for comparison.
%   3. Targets behind a strong near target, on the same azimuth (other
%      Doppler), at 30 km and at 23.5 km (just behind the long-pulse blind
%      range, where the long-window STC is deepest): detected in the long
%      pulse, and no other long-pulse detection at the velocity of one of
%      the three targets (eclipse ghost or range sidelobe).
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
cfg.console.enable = false;  cfg.display.enable = false;  cfg.io.saveRaw = false;
cfg.rx.maxRange_m = 50e3;
fprintf('fc = %.0f MHz, receiver gain %g dB\n', cfg.fc/1e6, cfg.rx.gain_dB);

%% 1) ADC overload
gains = {60, false;  cfg.rx.gain_dB, cfg.rx.stc.enable;  30, false};    % receiver gain, STC
for k = 1:size(gains, 1)
    c = cfg;  c.rx.gain_dB = gains{k, 1};  c.rx.stc.enable = gains{k, 2};
    sys = setup_radar(c);
    rng(1);
    clut = clutter_init(sys);
    frac = 0;
    for az = 0:45:315
        frac = max(frac, adc_overload(adc_sample(generate_rx_signal(sys, [], 0, az, clutter_signal(clut, sys, az, 0)), sys), sys));
    end
    good = (k == 1) == (frac > c.adc.overloadWarn);
    ok = report(sprintf('receiver gain %2g dB, STC %-3s: up to %.3f%% of the receive samples clipped (warning above %.0f%%)', ...
                        gains{k, 1}, ternary(gains{k, 2}, 'on', 'off'), 100*frac, 100*c.adc.overloadWarn), good) && ok;
end

sys = setup_radar(cfg);
h = cfg.ant.height_m;
c0 = cfg;  c0.rx.gain_dB = 30;  c0.rx.stc.enable = false;           % for comparison: no STC
sys0 = setup_radar(c0);
blindLong = cfg.const.c * sys.wave.long.width_s / 2;               % long-pulse blind range

%% 2) strong targets near the radar: one detection, no eclipse ghosts
cases = {10, 5e3, -150;  2, 7.7e3, -40;  1, 12e3, -236;  0.5, 20.6e3, 129};   % RCS, range, radial velocity
for k = 1:size(cases, 1)
    tg = make_target('N', cases{k, 1}, cases{k, 2}, 40, 1000, cases{k, 3}, 0);
    r = cpi_reports(cfg, sys, tg, 40, true);
    [Rt, ~, ~, vr] = target_geometry(tg, r.t0, h);
    s = r.short;
    nGhost = eclipse_ghosts(r.long, Rt, vr, blindLong);
    r0 = cpi_reports(c0, sys0, tg, 40, true);
    nGhost0 = eclipse_ghosts(r0.long, Rt, vr, blindLong);
    good = numel(s) == 1 && abs(s(1).range_m - Rt) < 75 && abs(s(1).vr_mps - vr) < 7 && nGhost == 0;
    ok = report(sprintf(['%4.1f m^2 at %4.1f km, %4.0f m/s: %d short-pulse detection(s), %d long-pulse ghost(s) of its ', ...
                         'eclipsed echo (%d without STC)'], cases{k, 1}, Rt/1e3, vr, numel(s), nGhost, nGhost0), good) && ok;
end

%% 3) targets behind a strong near target, same azimuth
tg = [make_target('N', 1, 12e3, 40, 1000, -236, 0), make_target('F', 1, 30e3, 40, 3000, 150, 0), ...
      make_target('G', 1, 23.5e3, 40, 2000, 80, 0)];
r = cpi_reports(cfg, sys, tg, 40, true);
other = false(size(r.long));                 % at the velocity of a target, not on a target cell
found = false(1, 2);
for k = 1:3
    [Rk, ~, ~, vk] = target_geometry(tg(k), r.t0, h);
    atV = abs([r.long.vr_mps] - vk) < 7;
    hit = abs([r.long.range_m] - Rk) < 75 & atV;
    if k > 1, found(k - 1) = any(hit); end
    other = (other | atV) & ~hit;
end
ok = report(sprintf(['1 m^2 at 30 km and at 23.5 km (deepest long-window STC) behind 1 m^2 at 12 km (same azimuth): ', ...
                     'detected in the long pulse %d %d, other long-pulse detections at their velocities %d'], ...
                    found, nnz(other)), all(found) && ~any(other)) && ok;

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

function n = eclipse_ghosts(long, Rt, vr, blindLong)
% Long-pulse detections at the velocity of a target at range Rt < blindLong
% within the reach of its eclipsed long-pulse echo.
n = 0;
if isempty(long), return; end
n = nnz([long.range_m] <= Rt + blindLong + 1e3 & abs([long.vr_mps] - vr) < 7);
end

function r = cpi_reports(cfg, sys, targets, az0, clutterOn)
% One CPI centred on azimuth az0: detection reports of both sub-pulses (as run_radar_sim).
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
for i = 1:2
    wv = sys.wave.(names{i});
    [cpi, info] = cpi_buffer_read(buf.(names{i}));
    rd = doppler_process(mti_filter(cpi, sys.mtiCoeffs, info.staggerIdx), sys);
    [isDet, noise] = cfar_detect(rd, sys, wv);
    r.(names{i}) = detection_reports(prune_peaks(isDet, rd, noise, sys, wv), rd, sys, wv, info);
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

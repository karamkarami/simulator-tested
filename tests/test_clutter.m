function ok = test_clutter()
%TEST_CLUTTER  Clutter power (CNR) and terrain types (no full scan).
%   1. CNR: uniform land, one element, STC off: the measured clutter-to-noise
%      ratio at the matched filter output (7 km, 60 azimuths) is the CNR that
%      clutter_cnr.m expects for cfg.clutter.cnr_dB = 60, 72.7 and 85 dB
%      (within 1 dB), and cfg.clutter.cnr_dB = [] keeps cfg.clutter.sigma0_dB.
%   2. Terrain map (cfg.clutter.terrain = 'regions'): desert and sea are
%      15 dB below land (measured as in 1, at 7 km, within 2 dB: Rayleigh
%      clutter over 20 - 40 azimuths), sea and desert end at the radar
%      horizon, the mountains are seen beyond it, their back slopes are
%      mostly shadowed, and every elevated region returns from its own
%      elevation (one steering vector per group).
%   3. Internal motion of every type: land and hills Gaussian 0.1 m/s,
%      mountains 0.05 m/s, desert 0.02 m/s, sea 1 m/s around the radial
%      component of its 2 m/s drift toward the East.
%   4. ADC: no receive sample clipped over 360 deg with the default CNR and
%      STC; 15 dB more CNR clips the mountains and hills above
%      cfg.adc.overloadWarn (run_radar_sim warns).
%   5. One CPI (48 elements, fixed PRF) over every terrain type with a target
%      in it: the target is detected (range within 75 m, radial velocity
%      within 7 m/s) and nothing else.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;

%% 1) CNR of uniform land
cfg = radar_config();
cfg.console.enable = false;  cfg.display.enable = false;  cfg.io.saveRaw = false;
c1 = cfg;
c1.ant.nElements = 1;  c1.rx.stc.enable = false;  c1.clutter.terrain = 'uniform';
c1.adc.quantize = false;  c1.adc.addAdcNoise = false;  c1.ddc.fixedPoint = false;  c1.rx.addLeakage = false;
for cnr = [60 72.7 85]
    c = c1;  c.clutter.cnr_dB = cnr;
    [meas, expd] = measure_cnr(c, 7000, 0:6:354);
    ok = report(sprintf('uniform land, cfg.clutter.cnr_dB = %5.1f dB (at 5 km): CNR at 7 km %.2f dB, expected %.2f dB', ...
                        cnr, meas, expd), abs(meas - expd) < 1) && ok;
end
c = c1;  c.clutter.cnr_dB = [];  c.clutter.sigma0_dB = -25;
s = setup_radar(c);
ok = report(sprintf('cfg.clutter.cnr_dB = []: sigma0 %.2f dB (cfg.clutter.sigma0_dB -25), CNR %.1f dB', ...
                    s.clutter.sigma0_dB, s.clutter.cnr_dB), abs(s.clutter.sigma0_dB + 25) < 1e-9) && ok;

%% 2) terrain map
c = c1;  c.clutter.terrain = 'regions';
[landM, landE] = measure_cnr(c, 7000, [345:3:357, 0:3:96]);
for tc = {'desert', 185:3:245;  'sea', 284:3:338}.'
    [m, e] = measure_cnr(c, 7000, tc{2});
    rel = c.clutter.types.(tc{1}).sigma0Rel_dB;
    ok = report(sprintf('%-6s CNR at 7 km %.1f dB (expected %.1f), land %.1f dB (expected %.1f): %+.1f dB (type %+g dB)', ...
                        tc{1}, m, e, landM, landE, m - landM, rel), abs(m - e) < 1.5 && abs(landM - landE) < 1.5 && ...
                abs(m - landM - rel) < 2) && ok;
end
sys = setup_radar(cfg);
T = clutter_terrain(cfg, sys);
horizon = sqrt(2 * 4/3 * 6371e3 * cfg.ant.height_m);
isT = @(name) T.type == find(strcmp(T.typeNames, name));
R = repmat(T.range_m, numel(T.az_deg), 1);
flatFar = (isT('sea') | isT('desert') | isT('land')) & R > horizon + 100;
mtn = isT('mountain');
ok = report(sprintf('horizon %.1f km: %d cells of flat terrain beyond it seen; mountains (%.0f - %.0f km) %d of %d cells seen', ...
                    horizon/1e3, nnz(flatFar & T.visible), min(R(mtn))/1e3, max(R(mtn))/1e3, nnz(mtn & T.visible), ...
                    nnz(mtn)), ~any(flatFar(:) & T.visible(:)) && nnz(mtn & T.visible) > 0.2*nnz(mtn)) && ok;
rr = cfg.clutter.regions{strcmp(cfg.clutter.regions(:, 1), 'mountain'), 3};
front = mtn & R < mean(rr);  back = mtn & R > mean(rr);
fF = nnz(front & T.visible) / nnz(front);  fB = nnz(back & T.visible) / nnz(back);
ok = report(sprintf('mountain range: %.0f%% of the front slopes seen, %.0f%% of the back slopes (shadow behind the ridge)', ...
                    100*fF, 100*fB), fF > 0.6 && fB < 0.3) && ok;
el = sys.clutter.groupEl_deg;
ok = report(sprintf('elevation groups %s deg, steering vectors %d x %d', mat2str(el.', 3), size(sys.clutterSteer)), ...
            numel(el) == 3 && el(1) == cfg.clutter.elevation_deg && all(el(2:end) > 0) && ...
            isequal(size(sys.clutterSteer), [sys.nCh, numel(el)])) && ok;

%% 3) internal motion
rng(3);
cl = clutter_init(sys);
v = cl.fd_Hz * sys.lambda / 2;
az = cl.az_deg(cl.cellRow);
expect = {'land', 0, 0.1;  'hill', 0, 0.1;  'mountain', 0, 0.05;  'desert', 0, 0.02};
for k = 1:size(expect, 1)
    m = cl.cellType == find(strcmp(cl.typeNames, expect{k, 1}));
    ok = report(sprintf('%-8s %6d cells: radial velocity mean %+.3f m/s, std %.3f m/s (expected %+g, %g)', expect{k, 1}, ...
                        nnz(m), mean(v(m)), std(v(m)), expect{k, 2}, expect{k, 3}), ...
                abs(mean(v(m)) - expect{k, 2}) < 0.01 && abs(std(v(m)) / expect{k, 3} - 1) < 0.05) && ok;
end
sea = cl.cellType == find(strcmp(cl.typeNames, 'sea'));
st = cfg.clutter.types.sea;
d = v(sea) - st.meanVel_mps * cosd(az(sea) - st.meanDir_deg);
w = sea & abs(az - 310) < 5;
vw = mean(st.meanVel_mps * cosd(az(w) - st.meanDir_deg));
ok = report(sprintf(['sea %d cells: drift removed mean %+.3f m/s, std %.3f m/s (expected 0, %g); at az 305 - 315 ', ...
                     'mean %+.2f m/s (expected %+.2f: toward the radar)'], nnz(sea), mean(d), std(d), st.icmStd_mps, ...
                    mean(v(w)), vw), abs(mean(d)) < 0.02 && abs(std(d) / st.icmStd_mps - 1) < 0.05 && ...
            abs(mean(v(w)) - vw) < 0.05 && vw < -1) && ok;

%% 4) ADC overload
for extra = [0 15]
    c = cfg;  c.ant.nElements = 1;  c.clutter.cnr_dB = cfg.clutter.cnr_dB + extra;
    s = setup_radar(c);
    rng(1);
    cl = clutter_init(s);
    azs = 0:4:356;
    frac = zeros(size(azs));
    for k = 1:numel(azs)
        frac(k) = adc_overload(adc_sample(generate_rx_signal(s, [], 0, azs(k), clutter_signal(cl, s, azs(k), 0)), s), s);
    end
    [fMax, kMax] = max(frac);
    if extra == 0
        good = fMax == 0;
    else
        good = fMax > c.adc.overloadWarn;
    end
    ok = report(sprintf('CNR %.1f dB: up to %.3f%% of the receive samples clipped (az %d; warning above %.0f%%)', ...
                        c.clutter.cnr_dB, 100*fMax, azs(kMax), 100*c.adc.overloadWarn), good) && ok;
end

%% 5) one CPI over every terrain type
cases = {'land',     make_target('L', 1, 6e3,   0, 1000,  120, 0);
         'hill',     make_target('H', 2, 8e3, 125,  500,  -60, 0);
         'mountain', make_target('M', 2, 48e3, 50, 4500, -150, 0);
         'desert',   make_target('D', 1, 7e3, 215,  800,  100, 0);
         'sea',      make_target('S', 1, 7e3, 310,  300,  -80, 0)};
h = cfg.ant.height_m;
for k = 1:size(cases, 1)
    tg = cases{k, 2};
    [~, azT] = target_geometry(tg, 0, h);
    r = cpi_reports(cfg, sys, tg, azT);
    [Rt, ~, ~, vr] = target_geometry(tg, r.t0, h);
    all_ = [r.short(:); r.long(:)];
    hit = false(size(all_));
    for i = 1:numel(all_)
        hit(i) = abs(all_(i).range_m - Rt) < 75 && abs(all_(i).vr_mps - vr) < 7;
    end
    ok = report(sprintf('over %-8s: %g m^2 at %.1f km, %+.0f m/s: %d detection(s) on the target, %d other', ...
                        cases{k, 1}, tg.rcs_m2, Rt/1e3, vr, nnz(hit), nnz(~hit)), any(hit) && all(hit)) && ok;
end
fprintf('Clutter and terrain: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function [meas, expd] = measure_cnr(cfg, R, azs)
% Clutter and noise power at the matched filter output of the short pulse at range R (one element).
sys = setup_radar(cfg);
c0 = cfg;  c0.rx.addThermalNoise = false;
s0 = setup_radar(c0);
rng(1);
clut = clutter_init(s0);
wv = sys.wave.short;
[~, k] = min(abs(wv.range_m - R));
ks = k + (-5:5);
pc = 0;  pn = 0;
for az = azs
    rp = range_profile(s0, generate_rx_signal(s0, [], 0, az, clutter_signal(clut, s0, az, 0)), wv);
    pc = pc + mean(abs(rp(ks)).^2);
    rp = range_profile(sys, generate_rx_signal(sys, [], 0, az, []), wv);
    pn = pn + mean(abs(rp).^2);
end
meas = 10*log10(pc / pn);
T = clutter_terrain(cfg, sys);
[~, a] = min(abs(T.az_deg - azs(1)));
[~, r] = min(abs(T.range_m - R));
rel = 10*log10(T.sigma0(a, r) / 10^(sys.clutter.sigma0_dB/10));     % type of the measured cells
expd = 10*log10(clutter_cnr(sys, 10^(sys.clutter.sigma0_dB/10), R, cfg.clutter.elevation_deg)) + rel;
end

function rp = range_profile(sys, x, wv)
[I, Q] = ddc(adc_sample(x, sys), sys);
b = double(I) + 1j*double(Q);
rp = matched_filter(b(wv.rxIdx).', wv);
end

function r = cpi_reports(cfg, sys, targets, az0)
% One CPI centred on azimuth az0 with clutter: detection reports of both sub-pulses (as run_radar_sim).
rng(11);
targets = draw_target_rcs(targets, 2, 0);
clut = clutter_init(sys);
n = sys.cpi.length;
r.t0 = (az0 - cfg.ant.az0_deg) / sys.ant.rate_dps;
names = {'short', 'long'};
for i = 1:2, buf.(names{i}) = cpi_buffer_init(sys.nBeams, sys.wave.(names{i}).nRange, n); end
st = mod(0:n-1, sys.nStagger) + 1;
dt = reshape(sys.pri(st(1:end-1)), 1, []);
tp = r.t0 + [0, cumsum(dt)] - sum(dt) / 2;
for k = 1:n
    az = mod(cfg.ant.az0_deg + sys.ant.rate_dps * tp(k), 360);
    ci = clutter_signal(clut, sys, az, tp(k));
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

function ok = report(text, ok)
fprintf('  %-128s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

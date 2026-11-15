function ok = test_modes()
%TEST_MODES  Radar modes, single/dual sub-modes, the scan-to-scan PRF and the CFAR limit cfg.cfar.maxPfa.
%   1. radar_config: modes 0, 1, 2 available, mode 0 'medium_mode' active
%      with both sub-pulses phase codes and mismatch filters, 'dual'.
%   2. apply_radar_mode: mode 1 'lfm_mode' (LFM, matched filters); its
%      'single' sub-mode keeps only cfg.wave.singlePulse (the long pulse); back
%      to mode 0 'dual' restores both sub-pulses and the phase codes.
%   3. setup_radar refuses a mode that is not available, more than 6 modes,
%      repeated or out-of-range mode numbers, an unknown sub-mode and an
%      unknown cfg.wave.singlePulse, with a message that names the setting.
%   4. Scan-to-scan PRF: one set of slow-time settings per PRF of
%      cfg.files.scan2scan (CPI = the pulses in the 3 dB beam, FFT, Doppler
%      axis); pulse_schedule gives every scan its PRF (the pulses of a scan
%      t0 + n*PRI), the PRI header carries the number of the PRF in its
%      stagger index; select_prf_profile puts the settings of a PRF into sys.
%   5. A simulated target near North in two scans (one element, no
%      clutter): plots in both scans with the right radial velocity, each on
%      the Doppler axis of the PRF of its scan.
%   6. cfg.cfar.maxPfa = 1e-4: the lowest threshold has Pfa 1e-4.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;

%% 1) defaults
cfg = radar_config();
ok = report('modes 0, 1, 2 available, mode 0 ''medium_mode'' active, dual', isequal(cfg.modes.available, [0 1 2]) && ...
            cfg.modes.active == 0 && strcmp(cfg.modes.name, 'medium_mode') && strcmp(cfg.modes.subMode, 'dual')) && ok;
ok = report('mode 0: both sub-pulses phase codes with mismatch filters', ...
            strcmp(cfg.wave.short.type, 'phasecode') && strcmp(cfg.wave.long.type, 'phasecode') && ...
            strcmp(cfg.wave.short.compression, 'mismatch') && strcmp(cfg.wave.long.compression, 'mismatch') && ...
            cfg.wave.short.enable && cfg.wave.long.enable) && ok;

%% 2) another mode and the single-pulse sub-mode
c = cfg;  c.modes.active = 1;  c = apply_radar_mode(c);
ok = report('mode 1 ''lfm_mode'': both sub-pulses LFM, matched filters', strcmp(c.modes.name, 'lfm_mode') && ...
            strcmp(c.wave.short.type, 'lfm') && strcmp(c.wave.long.type, 'lfm') && ...
            strcmp(c.wave.long.compression, 'matched') && c.wave.short.enable && c.wave.long.enable) && ok;
c.modes.subMode = 'single';  c = apply_radar_mode(c);
ok = report('mode 1, single: only cfg.wave.singlePulse (long) on', ~c.wave.short.enable && c.wave.long.enable && ...
            strcmp(c.wave.singlePulse, 'long')) && ok;
c.console.enable = false;  c.ant.nElements = 1;
s = setup_radar(c);
ok = report('setup of the single sub-mode: no short-pulse cells, the long pulse from its blind range', ...
            ~s.wave.short.enable && s.wave.short.nRange == 0 && s.wave.long.nRange > 0 && ...
            s.mode.active == 1 && strcmp(s.mode.subMode, 'single')) && ok;
c.modes.active = 0;  c.modes.subMode = 'dual';  c = apply_radar_mode(c);
ok = report('back to mode 0, dual: both sub-pulses and the phase codes again', c.wave.short.enable && ...
            c.wave.long.enable && strcmp(c.wave.short.type, 'phasecode')) && ok;

%% 3) invalid modes
base = cfg;  base.console.enable = false;  base.ant.nElements = 1;
bad = {'modes.active', 4, 'cfg.modes.active';  'modes.available', 0:6, 'cfg.modes.available';
       'modes.available', [0 0 1], 'cfg.modes.available';  'modes.available', [0 7], 'cfg.modes.available';
       'modes.subMode', 'triple', 'cfg.modes.subMode';  'wave.singlePulse', 'middle', 'cfg.wave.singlePulse'};
for i = 1:size(bad, 1)
    c = base;
    f = strsplit(bad{i, 1}, '.');
    c.(f{1}).(f{2}) = bad{i, 2};
    msg = '';
    try
        setup_radar(c);
    catch err
        msg = err.message;
    end
    v = bad{i, 2};
    if ischar(v), v = ['''' v '''']; else, v = mat2str(v); end
    ok = report(sprintf('refused: %s = %s', bad{i, 1}, v), ~isempty(strfind(msg, bad{i, 3}))) && ok;
end

%% 4) scan-to-scan PRF: settings per PRF, schedule, header
c = base;
c.prf.mode = 'scan2scan';
c.sim.nScans = 3;
s = setup_radar(c);
st = load(c.files.scan2scan);
n = numel(st.prf_Hz);
cpiOk = true;
for j = 1:n
    q = s.prfProfiles(j);
    cpiOk = cpiOk && abs(q.pri - 1/st.prf_Hz(j)) < 1e-15 && q.cpi.length == round(s.cpi.time_s * st.prf_Hz(j)) && ...
            q.dop.nfft == 2^nextpow2(q.cpi.length) && abs(q.dop.f_Hz(2) - q.dop.f_Hz(1) - st.prf_Hz(j)/q.dop.nfft) < 1e-9;
end
ok = report(sprintf('one set of slow-time settings per PRF (%d PRFs): PRI, CPI in the beam, FFT, Doppler axis', n), ...
            s.scan2scan && numel(s.prfProfiles) == n && cpiOk && s.prfProfile == 1) && ok;
sc = pulse_schedule(s);
tScan = 360 / s.ant.rate_dps;
schedOk = true;
for k = 1:c.sim.nScans
    j = mod(k - 1, n) + 1;
    in = sc.time_s >= (k-1)*tScan - 1e-9 & sc.time_s < k*tScan - 1e-9;
    schedOk = schedOk && all(sc.prfIdx(in) == j) && all(abs(sc.pri_s(in) - 1/st.prf_Hz(j)) < 1e-15) && ...
              abs(sum(in) - ceil(tScan * st.prf_Hz(j))) <= 1;
end
hdrA = pulse_header(s, sc, 1);  hdrB = pulse_header(s, sc, find(sc.prfIdx == 2, 1));
ok = report('schedule: every scan with its PRF; header stagger index = number of the PRF', schedOk && ...
            all(sc.staggerIdx == 1) && hdrA.staggerIndex == 1 && hdrB.staggerIndex == 2 && ...
            abs(hdrB.pri_s - 1/st.prf_Hz(2)) < 1e-15) && ok;
s2 = select_prf_profile(s, 2);
ok = report('select_prf_profile: PRI, CPI, FFT, axes and Doppler spans of PRF 2', s2.prfProfile == 2 && ...
            s2.pri == s.prfProfiles(2).pri && s2.cpi.length == s.prfProfiles(2).cpi.length && ...
            isequal(s2.dop.f_Hz, s.prfProfiles(2).dop.f_Hz) && ...
            abs(s2.track.vrSpanPulse.long - s2.wave.long.lambda_m * st.prf_Hz(2) / 2) < 1e-9 && ...
            s2.track.vrUnambiguous == (c.track.maxSpeed_mps < s.track.vrSpanMin / 2)) && ok;

%% 5) a target crossing North in two scans
c = base;
c.prf.mode = 'scan2scan';
c.sim.nScans = 2;
c.sim.sectors_deg = [355 10];
c.clutter.enable = false;
c.display.enable = false;
c.io.saveRaw = false;
c.track.enable = false;
c.console.enable = false;
vr = -120;
c.targets = make_target('A', 10, 30e3, 4, 2000, vr, 0);   % (4 deg: a whole CPI after the start of scan 1)
r = run_radar_sim(c);
p = r.plots;
t = [p.time_s];
inScan = {t < tScan, t >= tScan};
found = false(1, 2);
for k = 1:2
    q = p(inScan{k});
    bin = abs(diff(s.prfProfiles(k).dop.f_Hz(1:2))) * r.sys.wave.long.lambda_m / 2;
    found(k) = any(abs([q.range_m] - (30e3 + vr * [q.time_s])) < 500 & abs([q.vr_mps] - vr) < 1.5 * bin);   % (it closes in)
end
ok = report(sprintf('target near North: a plot in scan 1 (PRF %.0f Hz) and in scan 2 (PRF %.0f Hz), right velocity', ...
            st.prf_Hz(1), st.prf_Hz(2)), all(found) && r.sys.prfProfile == 2) && ok;

%% 6) lowest CFAR threshold
s = setup_radar(base);
ok = report(sprintf('cfg.cfar.maxPfa = %.0e: the lowest threshold %.2f dB has Pfa 1e-4', base.cfar.maxPfa, ...
            10*log10(s.cfar.alphaMin)), base.cfar.maxPfa == 1e-4 && ...
            abs(cfar_pfa(s.cfar.alphaMin, s.cfar, s.dop.nLooks) / 1e-4 - 1) < 1e-3) && ok;
end

function ok = report(text, ok)
fprintf('  %-110s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

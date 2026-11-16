function ok = test_modes()
%TEST_MODES  Radar modes (8), one sub-pulse per PRI, the scan-to-scan PRF and the CFAR limit cfg.cfar.maxPfa.
%   1. radar_config: modes 0, 1, 2 available, mode 0 'medium_mode' active
%      with both sub-pulses phase codes and mismatch filters.
%   2. apply_radar_mode: mode 1 'lfm_mode' (LFM, matched filters); all 8
%      modes 0 ... 7 load and set up; mode 7 with the short pulse switched
%      off (a mode with one sub-pulse per PRI); back to mode 0 restores both
%      sub-pulses and the phase codes.
%   3. setup_radar refuses a mode that is not available, more than 8 modes,
%      repeated or out-of-range mode numbers, the removed cfg.modes.subMode /
%      cfg.wave.singlePulse and phase-code chips shorter than a DDC sample,
%      with a message that names the setting.
%   4. Scan-to-scan PRF: one set of slow-time settings per PRF of
%      cfg.files.scan2scan (CPI = the pulses in the 3 dB beam, FFT, Doppler
%      axis); pulse_schedule gives every scan its PRF (the pulses of a scan
%      t0 + n*PRI), the PRI header carries the number of the PRF in its
%      stagger index; select_prf_profile puts the settings of a PRF into sys.
%   5. A simulated target near North in two scans (one element, no
%      clutter): plots in both scans with the right radial velocity, each on
%      the Doppler axis of the PRF of its scan.
%   6. cfg.cfar.maxPfa = 1e-4: the lowest threshold has Pfa 1e-4.
%   7. The whole processing chain belongs to the mode: every mode file sets
%      the same fields (waveforms, timing, PRF, files, STC, pulse compression,
%      CPI, Doppler window, MTI, CFAR, pruning, reports, extraction, tracking),
%      radar_config.m none of them; apply_radar_mode replaces changed
%      processing settings by those of the mode; a mode of its own (a mode
%      file earlier on the path) with another CFAR, tracker, Doppler window,
%      staggered PRF from its own file, CPI length and STC gives a system with
%      exactly these.
%   8. Transmit times and receive windows (mode 2): short pulse at 11 us,
%      received 40 ... 248 us, long pulse at 321 us, received 533 ... 1780 us:
%      the DDC samples of each window lie in it, the range cells run from
%      c*(rxStart - tTx)/2 to c*(rxEnd - tTx - width)/2 (4.35 ... 31.8 km,
%      31.8 ... 200 km, no gap), the record ends after the long window, the
%      PRI header carries the windows; invalid windows stop with a message.
%   9. A simulation in mode 2 (one element, no clutter): a target at 20 km
%      (short pulse sent at 11 us) and one at 120 km give plots at their
%      true range and radial velocity.
%  10. cfg.proc.dopplerWindow: Hann, Blackman and rectangular windows with
%      their MTI notch widths; an unknown window stops with a message.
%  11. Manual DDC record (cfg.rx.manualDdc): nDdc = 1800 samples, the first
%      at ddcT0_s = 0 (the record starts before the PRI start): the sizes, the
%      time of the first sample (rounded to the ADC grid), the PRI header, the
%      range cells of both windows as with the automatic record (within one
%      DDC sample), a simulation finds targets in both windows at their range;
%      a missing nDdc, a record beyond the PRI and a set window end after the
%      record stop with a message.
%  12. Code agility: cfg.wave.<pulse>.codeIndex selects the code / mismatch
%      pair of the file among cfg.wave.<pulse>.availableCodes (default: the
%      first); more available pairs rebuild the default file (the first pairs
%      unchanged); a file of its own with any pair numbers; a file with only
%      <pulse>Code / <pulse>Mismatch is pair 1; a codeIndex not available, an
%      available pair missing from the file, an invalid list and the removed
%      cfg.wave.maxCodes stop with a message.
%  13. Stagger pairs: cfg.prf.staggerSet selects the pair (T, w) of the file
%      among cfg.prf.availableStaggerSets (its PRIs and length(T) weight rows,
%      CPI of whole periods); invalid choices, files (an available pair
%      missing, w without length(T) rows) and the removed cfg.prf.maxStaggerSets
%      stop with a message; a bad pair the mode may not use does not matter.
%  14. Time-varying MTI weights of the staggered PRF: every output uses the
%      row of the stagger position of its oldest pulse (direct sum); on slowly
%      moving clutter (0.3 m/s) both null DC, but the time-varying weights
%      leave a much smaller residue than [1 -2 1] on the staggered pulses.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;

%% 1) defaults
cfg = radar_config();
ok = report('modes 0, 1, 2 available, mode 0 ''medium_mode'' active', isequal(cfg.modes.available, [0 1 2]) && ...
            cfg.modes.active == 0 && strcmp(cfg.modes.name, 'medium_mode') && ~isfield(cfg.modes, 'subMode') && ...
            ~isfield(cfg.wave, 'singlePulse')) && ok;
ok = report('mode 0: both sub-pulses phase codes with mismatch filters', ...
            strcmp(cfg.wave.short.type, 'phasecode') && strcmp(cfg.wave.long.type, 'phasecode') && ...
            strcmp(cfg.wave.short.compression, 'mismatch') && strcmp(cfg.wave.long.compression, 'mismatch') && ...
            cfg.wave.short.enable && cfg.wave.long.enable) && ok;

%% 2) the other modes, all 8 of them, and a mode with one sub-pulse per PRI
c = cfg;  c.modes.active = 1;  c = apply_radar_mode(c);
ok = report('mode 1 ''lfm_mode'': both sub-pulses LFM, matched filters', strcmp(c.modes.name, 'lfm_mode') && ...
            strcmp(c.wave.short.type, 'lfm') && strcmp(c.wave.long.type, 'lfm') && ...
            strcmp(c.wave.long.compression, 'matched') && c.wave.short.enable && c.wave.long.enable) && ok;
allOk = true;
for k = 0:7
    c = cfg;  c.modes.available = 0:7;  c.modes.active = k;  c = apply_radar_mode(c);
    c.console.enable = false;  c.ant.nElements = 1;
    s = setup_radar(c);
    allOk = allOk && s.mode.active == k && isequal(s.mode.available, 0:7) && s.wave.short.nRange > 0 && ...
            s.wave.long.nRange > 0 && (k < 3 || strcmp(c.modes.name, sprintf('mode_%d', k)));
end
ok = report('modes 0 ... 7 (cfg.modes.available = 0:7): every mode file loads and sets up', allOk) && ok;
c = cfg;  c.modes.available = 0:7;  c.modes.active = 7;  c = apply_radar_mode(c);
c.wave.short.enable = false;  c.console.enable = false;  c.ant.nElements = 1;
s = setup_radar(c);
ok = report('mode 7 with the short pulse off: no short-pulse cells, the long pulse from its blind range', ...
            ~s.wave.short.enable && s.wave.short.nRange == 0 && s.wave.long.nRange > 0 && s.mode.active == 7) && ok;
c.modes.active = 0;  c = apply_radar_mode(c);
ok = report('back to mode 0: both sub-pulses and the phase codes again', c.wave.short.enable && ...
            c.wave.long.enable && strcmp(c.wave.short.type, 'phasecode')) && ok;

%% 3) invalid modes
base = cfg;  base.console.enable = false;  base.ant.nElements = 1;
bad = {'modes.active', 4, 'cfg.modes.active';  'modes.available', 0:8, 'cfg.modes.available';
       'modes.available', [0 0 1], 'cfg.modes.available';  'modes.available', [0 8], 'cfg.modes.available';
       'modes.subMode', 'single', 'cfg.modes.subMode';  'wave.singlePulse', 'long', 'cfg.wave.singlePulse'};
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
% phase-code chips shorter than a DDC sample (0.25 us at 2 MHz; 25 chips as in the default file)
c = base;  c.wave.short.chip_s = 0.25e-6;  c.wave.short.width_s = 6.25e-6;
msg = '';
try
    setup_radar(c);
catch err
    msg = err.message;
end
ok = report('refused: short-pulse chips of 0.25 us (shorter than the 0.5 us DDC sample)', ...
            ~isempty(strfind(msg, 'cfg.wave.short.chip_s'))) && ok;

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

%% 7) the whole processing chain per mode
procFields = {'wave', 'prf', 'proc', 'mti', 'cfar', 'prune', 'report', 'extract', 'track'};
common = struct('modes', cfg.modes);
fieldsOk = true;
setsOf = cell(1, 6);
for k = 0:5
    m = feval(sprintf('radar_config_mode_%d', k), common);
    setsOf{k+1} = flat_fields(m);
    fieldsOk = fieldsOk && all(isfield(m, procFields)) && isfield(m.rx, 'stc') && isfield(m.rx, 'maxRange_m') && ...
               all(isfield(m.files, {'stagger', 'scan2scan', 'phaseCode', 'mtiFir'})) && ...
               isequal(sort(setsOf{k+1}), sort(setsOf{1}));
end
txt = fileread(fullfile(root, 'radar_config.m'));
noneInConfig = isempty(regexp(txt, '(?m)^cfg\.(proc|mti|cfar|prune|report|extract|track|prf|wave)\.', 'once')) && ...
               isempty(regexp(txt, '(?m)^cfg\.rx\.(stc|maxRange_m)', 'once'));
ok = report(sprintf('every mode file sets the same %d settings (processing chain, timing, files); radar_config.m none', ...
            numel(setsOf{1})), fieldsOk && noneInConfig) && ok;
c = cfg;
c.cfar.type = 'os';  c.track.method = 'pmbm';  c.proc.cpiMode = 'manual';  c.rx.stc.enable = false;
c = apply_radar_mode(c);
ok = report('apply_radar_mode replaces changed processing settings (CFAR, tracker, CPI, STC) by those of the mode', ...
            strcmp(c.cfar.type, 'ca') && strcmp(c.track.method, 'ekf') && strcmp(c.proc.cpiMode, 'auto') && ...
            c.rx.stc.enable) && ok;
% a mode of its own: mode 5 from a folder earlier on the path
tmp = tempname();  mkdir(tmp);
src = fileread(fullfile(root, 'modes', 'radar_config_mode_5.m'));
src = strrep(src, 'root = fileparts(fileparts(mfilename(''fullpath'')));', sprintf('root = ''%s'';', root));
ownStagger = fullfile(tmp, 'stagger_mode_5.mat');
priVector_s = [900 1000 1100 950]' * 1e-6;  mtiCoeffs = repmat([1 -2 1], 4, 1);   %#ok<NASGU>
save(ownStagger, 'priVector_s', 'mtiCoeffs');
edits = {'cfg.cfar.type   = ''ca'';', 'cfg.cfar.type   = ''go'';';
         'cfg.cfar.method = ''2d'';', 'cfg.cfar.method = ''range'';';
         'cfg.track.method          = ''ekf'';', 'cfg.track.method          = ''tomht'';';
         'cfg.track.model           = ''cv'';', 'cfg.track.model           = ''ca'';';
         'cfg.proc.dopplerWindow = ''hamming'';', 'cfg.proc.dopplerWindow = ''hann'';';
         'cfg.prf.mode  = ''fixed'';', 'cfg.prf.mode  = ''staggered'';';
         sprintf('cfg.files.stagger   = fullfile(root, ''inputs'', ''stagger_prf.mat'');'), sprintf('cfg.files.stagger = ''%s'';', ownStagger);
         'cfg.prf.availableStaggerSets = [1 2 3 4];', 'cfg.prf.availableStaggerSets = 1;';
         'cfg.proc.cpiMode      = ''auto'';', 'cfg.proc.cpiMode      = ''manual'';';
         'cfg.proc.cpiLength    = 114;', 'cfg.proc.cpiLength    = 88;';
         'cfg.rx.stc.enable = true;', 'cfg.rx.stc.enable = false;'};
allFound = true;
for i = 1:size(edits, 1)
    allFound = allFound && ~isempty(strfind(src, edits{i, 1}));
    src = strrep(src, edits{i, 1}, edits{i, 2});
end
fid = fopen(fullfile(tmp, 'radar_config_mode_5.m'), 'w');  fwrite(fid, src);  fclose(fid);
addpath(tmp);
c = base;  c.modes.available = [0 5];  c.modes.active = 5;
c = apply_radar_mode(c);
c.console.enable = false;  c.ant.nElements = 1;
s5 = setup_radar(c);
rmpath(tmp);
delete(fullfile(tmp, '*'));  rmdir(tmp);
w = s5.dop.window * norm(hann_ref(s5.cpi.rowLen)) / norm(s5.dop.window);
ok = report(sprintf(['a mode of its own (mode 5): GO-CFAR on the range window, TOMHT with the CA model, its own ', ...
                     'stagger file (%d PRIs), CPI %d pulses, Hann window, STC off'], s5.nStagger, s5.cpi.length), ...
            allFound && strcmp(s5.cfar.type, 'go') && strcmp(s5.cfar.method, 'range') && ...
            strcmp(c.track.method, 'tomht') && strcmp(c.track.model, 'ca') && s5.staggered && s5.nStagger == 4 && ...
            isequal(s5.pri(:).', priVector_s(:).') && s5.cpi.length == 88 && max(abs(w - hann_ref(s5.cpi.rowLen))) < 1e-12 && ...
            all(s5.rx.stcGain == 1)) && ok;

%% 8) transmit times and receive windows (mode 2)
c = base;  c.modes.active = 2;  c = apply_radar_mode(c);
c.console.enable = false;  c.ant.nElements = 1;
s2 = setup_radar(c);
ws = s2.wave.short;  wl = s2.wave.long;
C = c.const.c;  dR = C / (2 * s2.fsOut);
inWin = @(wv) all(s2.tDdc(wv.rxIdx) >= wv.rxStart_s) && all(s2.tDdc(wv.rxIdx) < wv.rxEnd_s);
rangeOk = abs(ws.range_m(1) - C*(ws.rxStart_s - ws.tTx_s)/2) <= dR && ...
          abs(ws.range_m(end) - C*(ws.rxEnd_s - ws.tTx_s - ws.width_s)/2) <= dR && ...
          abs(wl.range_m(1) - C*(wl.rxStart_s - wl.tTx_s)/2) <= dR && ...
          abs(wl.range_m(end) - C*(wl.rxEnd_s - wl.tTx_s - wl.width_s)/2) <= dR && ...
          wl.range_m(1) - ws.range_m(end) <= dR;
hdr = pulse_header(s2, struct('scan', 1, 'time_s', 0, 'az_deg', 0, 'pri_s', s2.pri, 'staggerIdx', 1), 1);
ok = report(sprintf(['mode 2: short pulse at %.0f us, DDC samples of the windows %.1f ... %.1f us and %.1f ... %.1f us, ', ...
                     'range %.2f ... %.2f km and %.2f ... %.2f km'], ws.tTx_s*1e6, s2.tDdc(ws.rxIdx([1 end]))*1e6, ...
                     s2.tDdc(wl.rxIdx([1 end]))*1e6, ws.range_m([1 end])/1e3, wl.range_m([1 end])/1e3), ...
            inWin(ws) && inWin(wl) && rangeOk && s2.tAdc(end) < s2.pri && s2.tAdc(end) >= wl.rxEnd_s && ...
            isequal(hdr.shortDdcIdx, ws.rxIdx([1 end])) && isequal(hdr.longDdcIdx, wl.rxIdx([1 end])) && ...
            isequal(s2.rxWindow_s.short, [ws.rxStart_s ws.rxEnd_s])) && ok;
bad = {'wave.short.rxStart_s', 30e-6, 'cfg.wave.short.rxStart_s';   % before the end of the short pulse (36 us)
       'wave.short.rxEnd_s', 330e-6, 'cfg.wave.short.rxEnd_s';      % after the long pulse is sent
       'wave.long.rxEnd_s', 1.9e-3, 'shortest PRI';                 % after the PRI
       'wave.long.rxEnd_s', 500e-6, 'cfg.wave.long.rxEnd_s';        % before its start
       'wave.short.tTx_s', -1e-6, 'cfg.wave.short.tTx_s'};
for i = 1:size(bad, 1)
    cb = c;
    f = strsplit(bad{i, 1}, '.');
    cb.(f{1}).(f{2}).(f{3}) = bad{i, 2};
    msg = '';
    try
        setup_radar(cb);
    catch err
        msg = err.message;
    end
    ok = report(sprintf('refused: cfg.%s = %g us', bad{i, 1}, bad{i, 2}*1e6), ~isempty(strfind(msg, bad{i, 3}))) && ok;
end

%% 9) a simulation in mode 2
c.sim.nScans = 1;
c.sim.sectors_deg = [5 25];
c.clutter.enable = false;
c.display.enable = false;
c.io.saveRaw = false;
c.track.enable = false;
c.targets = [make_target('S', 10, 20e3, 13, 2000, -150, 0), make_target('L', 200, 120e3, 17, 6000, 200, 0)];
r = run_radar_sim(c);
p = r.plots;
found = '';
for i = 1:2
    tg = c.targets(i);
    pulse = {'short', 'long'};
    q = p(strcmp({p.pulse}, pulse{i}));
    for j = 1:numel(q)
        [pos, vel] = target_state(tg, q(j).time_s);
        R = norm(pos - [0 0 c.ant.height_m]);
        vr = dot(vel, (pos - [0 0 c.ant.height_m]) / R);
        bin = abs(diff(r.sys.dop.f_Hz(1:2))) * r.sys.wave.(pulse{i}).lambda_m / 2;
        if abs(q(j).range_m - R) < 100 && abs(q(j).vr_mps - vr) < 1.5 * bin
            found = [found tg.name]; %#ok<AGROW>
            break
        end
    end
end
ok = report(sprintf('mode 2 simulation: target S at 20 km (short pulse) and L at 120 km (long pulse) found: %s (%d plots)', ...
            found, numel(p)), strcmp(found, 'SL')) && ok;

%% 10) Doppler window
c = base;
ref = {'hann', 2, @hann_ref;  'blackman', 3, @(N) 0.42 - 0.5*cos(2*pi*(0:N-1)'/(N-1)) + 0.08*cos(4*pi*(0:N-1)'/(N-1));
       'rect', 1, @(N) ones(N, 1)};
for i = 1:size(ref, 1)
    c.proc.dopplerWindow = ref{i, 1};
    s = setup_radar(c);
    w0 = ref{i, 3}(s.cpi.rowLen);
    ok = report(sprintf('cfg.proc.dopplerWindow = ''%s'': window of %d samples, MTI notch +-%d bins', ref{i, 1}, ...
                s.cpi.rowLen, s.dop.notchBins), max(abs(s.dop.window - w0 / norm(w0))) < 1e-12 && ...
                s.dop.notchBins == round(c.cfar.notchWidth * ref{i, 2} * s.dop.nfft / s.cpi.rowLen)) && ok;
end
c.proc.dopplerWindow = 'kaiser';
msg = '';
try, setup_radar(c); catch err, msg = err.message; end
ok = report('refused: cfg.proc.dopplerWindow = ''kaiser''', ~isempty(strfind(msg, 'cfg.proc.dopplerWindow'))) && ok;

%% 11) manual DDC record
sA = setup_radar(base);
c = base;  c.rx.manualDdc = true;  c.rx.nDdc = 1800;  c.rx.ddcT0_s = 0;
sM = setup_radar(c);
hdr = pulse_header(sM, struct('scan', 1, 'time_s', 0, 'az_deg', 0, 'pri_s', sM.pri, 'staggerIdx', 1), 1);
dt = 1 / sM.fsOut;
cellsOk = true;
for p = {'short', 'long'}
    a = sA.wave.(p{1});  m = sM.wave.(p{1});
    cellsOk = cellsOk && abs(m.range_m(1) - a.range_m(1)) <= c.const.c*dt/2 + 1 && ...
              all(sM.tDdc(m.rxIdx) >= sA.tDdc(a.rxIdx(1)) - dt) && m.nRange > 0;
end
ok = report(sprintf(['manual DDC record: %d samples (automatic %d), the first at %.4f us (set 0, ADC grid), record ', ...
                     '%.2f ... %.1f us, range %.2f / %.2f km as automatic'], sM.nDdc, sA.nDdc, sM.ddcT0_s*1e6, ...
                     sM.rxStart_s*1e6, sM.tAdc(end)*1e6, sM.wave.short.range_m(1)/1e3, sM.wave.long.range_m(1)/1e3), ...
            sM.nDdc == 1800 && abs(sM.ddcT0_s) <= 0.5/sM.fs + 1e-15 && sM.rxStart_s < 0 && sM.tAdc(end) < sM.pri && ...
            sM.ddcRecord.manual && hdr.ddcT0_s == sM.ddcT0_s && hdr.rxStart_s == sM.rxStart_s && cellsOk) && ok;
c.sim.nScans = 1;  c.sim.sectors_deg = [5 25];  c.clutter.enable = false;  c.display.enable = false;
c.io.saveRaw = false;  c.track.enable = false;
c.targets = [make_target('S', 10, 15e3, 13, 2000, -150, 0), make_target('L', 100, 60e3, 17, 4000, 200, 0)];
r = run_radar_sim(c);
found = find_targets(r, c);
ok = report(sprintf('manual DDC record, simulation: S at 15 km (short) and L at 60 km (long) found: %s', found), ...
            strcmp(found, 'SL')) && ok;
bad = {'rx.nDdc', [], 'cfg.rx.nDdc';  'rx.nDdc', 2100, 'shortest PRI';  'rx.nDdc', 12.5, 'cfg.rx.nDdc'};
for i = 1:size(bad, 1)
    cb = base;  cb.rx.manualDdc = true;  cb.rx.nDdc = 1800;
    f = strsplit(bad{i, 1}, '.');  cb.(f{1}).(f{2}) = bad{i, 2};
    msg = '';
    try, setup_radar(cb); catch err, msg = err.message; end
    ok = report(sprintf('manual record refused: cfg.%s = %s', bad{i, 1}, mat2str(bad{i, 2})), ~isempty(strfind(msg, bad{i, 3}))) && ok;
end
cb = base;  cb.modes.active = 2;  cb = apply_radar_mode(cb);  cb.console.enable = false;  cb.ant.nElements = 1;
cb.rx.manualDdc = true;  cb.rx.nDdc = 3000;
msg = '';
try, setup_radar(cb); catch err, msg = err.message; end
ok = report('manual record refused: mode 2 window end (1780 us) after the record of 3000 samples', ...
            ~isempty(strfind(msg, 'cfg.wave.long.rxEnd_s'))) && ok;

%% 12) code agility
tmpd = tempname();  mkdir(tmpd);
cc = base;
cc.files.phaseCode = fullfile(tmpd, 'codes.mat');            % (a default file of its own: the inputs stay)
f = [];
codesOk = true;
list = cc.wave.short.availableCodes;  nL = numel(list);
for k = list
    cc.wave.short.codeIndex = k;  cc.wave.long.codeIndex = list(nL + 1 - k);
    sk = setup_radar(cc);
    if isempty(f), f = load(cc.files.phaseCode); end
    codesOk = codesOk && isequal(sk.wave.short.code, f.(sprintf('shortCode_%d', k))) && ...
              isequal(sk.wave.long.code, f.(sprintf('longCode_%d', list(nL + 1 - k)))) && ...
              isequal(sk.wave.short.mismatch, f.(sprintf('shortMismatch_%d', k))) && sk.wave.short.codeIndex == k && ...
              isequal(sk.wave.short.availableCodes, list);
end
ok = report(sprintf('code agility: codeIndex %s selects the code / mismatch pair of the file (both sub-pulses)', ...
            mat2str(list)), codesOk) && ok;
c2 = cc;  c2.wave.short.availableCodes = [4 2];  c2.wave.short.codeIndex = [];
s2 = setup_radar(c2);
ok = report('code agility: availableCodes [4 2] without codeIndex uses pair 4 (the first)', ...
            s2.wave.short.codeIndex == 4 && isequal(s2.wave.short.code, f.shortCode_4)) && ok;
c6 = cc;  c6.wave.short.availableCodes = 1:6;  c6.wave.short.codeIndex = 6;
s6 = setup_radar(c6);  f6 = load(c6.files.phaseCode);
same = true;
for k = 1:4, same = same && isequal(f6.(sprintf('shortCode_%d', k)), f.(sprintf('shortCode_%d', k))); end
ok = report('code agility: availableCodes 1 ... 6 rebuilds the default file with 6 pairs (pairs 1 ... 4 unchanged)', ...
            same && isfield(f6, 'shortCode_6') && ~isfield(f6, 'shortCode_7') && f6.madeFor.nCodes == 6 && ...
            isequal(s6.wave.short.code, f6.shortCode_6)) && ok;
u = rmfield(f, 'madeFor');
two = struct('shortCode_1', u.shortCode_1, 'shortMismatch_1', u.shortMismatch_1, 'shortCode_2', u.shortCode_2, ...
             'shortMismatch_2', u.shortMismatch_2, 'longCode_1', u.longCode_1, 'longMismatch_1', u.longMismatch_1, ...
             'longCode_2', u.longCode_2, 'longMismatch_2', u.longMismatch_2);
one = struct('shortCode', u.shortCode_3, 'shortMismatch', u.shortMismatch_3, 'longCode', u.longCode_1, 'longMismatch', u.longMismatch_1);
gap = struct('shortCode_2', u.shortCode_2, 'shortMismatch_2', u.shortMismatch_2, 'shortCode_7', u.shortCode_3, ...
             'shortMismatch_7', u.shortMismatch_3, 'longCode_1', u.longCode_1, 'longMismatch_1', u.longMismatch_1);
save(fullfile(tmpd, 'all.mat'), '-struct', 'u');  save(fullfile(tmpd, 'two.mat'), '-struct', 'two');
save(fullfile(tmpd, 'one.mat'), '-struct', 'one');  save(fullfile(tmpd, 'gap.mat'), '-struct', 'gap');
c1 = base;  c1.files.phaseCode = fullfile(tmpd, 'one.mat');
c1.wave.short.availableCodes = 1;  c1.wave.long.availableCodes = 1;  c1.wave.short.codeIndex = 1;  c1.wave.long.codeIndex = 1;
s1 = setup_radar(c1);
ok = report('code agility: a file with only shortCode / shortMismatch holds pair 1', ...
            isequal(s1.wave.short.code, u.shortCode_3) && s1.wave.short.codeIndex == 1) && ok;
cg = c1;  cg.files.phaseCode = fullfile(tmpd, 'gap.mat');  cg.wave.short.availableCodes = [2 7];  cg.wave.short.codeIndex = 7;
sg = setup_radar(cg);
ok = report('code agility: a file of its own with the pairs 2 and 7 (availableCodes [2 7], codeIndex 7)', ...
            isequal(sg.wave.short.code, u.shortCode_3) && sg.wave.short.codeIndex == 7) && ok;
% {file, short.availableCodes, short.codeIndex, long.availableCodes, part of the message}
bad = {'all.mat', 1:4, 5, 1:4, 'is not one of cfg.wave.short.availableCodes';
       'all.mat', [1 3], 2, 1:4, 'is not one of cfg.wave.short.availableCodes';
       'two.mat', 1:3, 1, [1 2], 'holds pair 3, but';
       'one.mat', [1 2], 1, 1, 'holds pair 2, but';
       'two.mat', [1 2], 1, 1:4, 'cfg.wave.long.availableCodes holds pair 3';
       'all.mat', [0 1], 1, 1:4, 'cfg.wave.short.availableCodes must hold';
       'all.mat', [1 1], 1, 1:4, 'cfg.wave.short.availableCodes must hold'};
for i = 1:size(bad, 1)
    cb = base;  cb.files.phaseCode = fullfile(tmpd, bad{i, 1});
    cb.wave.short.availableCodes = bad{i, 2};  cb.wave.short.codeIndex = bad{i, 3};
    cb.wave.long.availableCodes = bad{i, 4};  cb.wave.long.codeIndex = 1;
    msg = '';
    try, setup_radar(cb); catch err, msg = err.message; end
    ok = report(sprintf('code agility refused: %s, availableCodes %s, codeIndex %d, long %s', bad{i, 1}, ...
                mat2str(bad{i, 2}), bad{i, 3}, mat2str(bad{i, 4})), ~isempty(strfind(msg, bad{i, 5}))) && ok;
end
cb = base;  cb.wave.maxCodes = 4;
msg = '';
try, setup_radar(cb); catch err, msg = err.message; end
ok = report('code agility: the removed cfg.wave.maxCodes is refused', ~isempty(strfind(msg, 'cfg.wave.maxCodes no longer exists'))) && ok;

%% 13) stagger pairs (T, w)
cs = base;  cs.prf.mode = 'staggered';
cs.files.stagger = fullfile(tmpd, 'stagger.mat');            % (a default file of its own)
st = [];
setsOk = true;
for k = cs.prf.availableStaggerSets
    cs.prf.staggerSet = k;
    sk = setup_radar(cs);
    if isempty(st), st = load(cs.files.stagger); end
    T = st.(sprintf('priVector_s_%d', k));  w = st.(sprintf('mtiCoeffs_%d', k));
    setsOk = setsOk && isequal(sk.pri(:), T(:)) && isequal(sk.mtiCoeffs, w) && sk.staggerSet == k && ...
             mod(sk.cpi.length, numel(T)) == 0 && sk.cpi.nRows == numel(T) && size(w, 1) == numel(T) && ...
             isequal(sk.availableStaggerSets, cs.prf.availableStaggerSets);
end
ok = report(sprintf('stagger pairs: staggerSet %s selects T and w (length(T) rows) of the file (PRIs, weights, rows)', ...
            mat2str(cs.prf.availableStaggerSets)), setsOk) && ok;
c3 = cs;  c3.prf.availableStaggerSets = [3 1];  c3.prf.staggerSet = [];
s3 = setup_radar(c3);
ok = report('stagger pairs: availableStaggerSets [3 1] without staggerSet uses pair 3 (the first)', ...
            s3.staggerSet == 3 && isequal(s3.pri(:), st.priVector_s_3(:))) && ok;
v = rmfield(st, 'madeFor');
sv = struct('priVector_s_1', v.priVector_s_1, 'mtiCoeffs_1', v.mtiCoeffs_1, 'priVector_s_2', v.priVector_s_2, ...
            'mtiCoeffs_2', v.mtiCoeffs_2);
save(fullfile(tmpd, 'st_all.mat'), '-struct', 'v');  save(fullfile(tmpd, 'st_two.mat'), '-struct', 'sv');
bw = v;  bw.mtiCoeffs_2 = bw.mtiCoeffs_2(1:2, :);
save(fullfile(tmpd, 'st_badw.mat'), '-struct', 'bw');
c3 = cs;  c3.files.stagger = fullfile(tmpd, 'st_badw.mat');  c3.prf.availableStaggerSets = [1 3];  c3.prf.staggerSet = 3;
s3 = setup_radar(c3);
ok = report('stagger pairs: a pair with bad weights that the mode may not use (2) does not matter', s3.staggerSet == 3) && ok;
% {file, availableStaggerSets, staggerSet, part of the message}
bad = {'st_all.mat', 1:4, 5, 'is not one of cfg.prf.availableStaggerSets';
       'st_two.mat', 1:3, 1, 'holds pair 3, but';
       'st_badw.mat', [1 2], 1, 'needs 4 rows, one weight sequence per stagger position';
       'st_all.mat', [2 -1], 2, 'cfg.prf.availableStaggerSets must hold'};
for i = 1:size(bad, 1)
    cb = cs;  cb.files.stagger = fullfile(tmpd, bad{i, 1});
    cb.prf.availableStaggerSets = bad{i, 2};  cb.prf.staggerSet = bad{i, 3};
    msg = '';
    try, setup_radar(cb); catch err, msg = err.message; end
    ok = report(sprintf('stagger pairs refused: %s, availableStaggerSets %s, staggerSet %d', bad{i, 1}, ...
                mat2str(bad{i, 2}), bad{i, 3}), ~isempty(strfind(msg, bad{i, 4}))) && ok;
end
cb = cs;  cb.prf.maxStaggerSets = 4;
msg = '';
try, setup_radar(cb); catch err, msg = err.message; end
ok = report('stagger pairs: the removed cfg.prf.maxStaggerSets is refused', ...
            ~isempty(strfind(msg, 'cfg.prf.maxStaggerSets no longer exists'))) && ok;

%% 14) time-varying MTI weights
cs.prf.staggerSet = 1;
s1 = setup_radar(cs);
K = s1.nStagger;  n = 4*K + 2;
pos = mod(0:n-1, K) + 1;                                      % stagger position of every pulse
tP = [0, cumsum(s1.pri(pos(1:end-1)).')];                     % pulse times
x = reshape(exp(1j*2*pi*(2*0.3/s1.lambda)*tP), 1, 1, n);      % clutter moving at 0.3 m/s
y = squeeze(mti_filter(x, s1.mtiCoeffs, pos)).';
L = size(s1.mtiCoeffs, 2);
yRef = zeros(1, n - L + 1);
for m = 1:n - L + 1
    yRef(m) = s1.mtiCoeffs(pos(m), :) * squeeze(x(1, 1, m:m+L-1));   % the row of the oldest pulse
end
yFix = squeeze(mti_filter(x, repmat([1 -2 1], K, 1), pos)).';
dc = squeeze(mti_filter(ones(1, 1, n), s1.mtiCoeffs, pos));
ratio = 10*log10(mean(abs(y).^2) / mean(abs(yFix).^2));
ok = report(sprintf(['time-varying MTI: every output with the weights of its oldest pulse''s stagger position; ', ...
                     'clutter at 0.3 m/s: residue %.1f dB below [1 -2 1] on the staggered pulses, DC nulled'], -ratio), ...
            max(abs(y - yRef)) < 1e-12 && max(abs(dc)) < 1e-12 && ratio < -20) && ok;
rmdir(tmpd, 's');
end

function found = find_targets(r, c)
% Names of the targets with a plot of their sub-pulse at their range (100 m) and radial velocity.
found = '';
p = r.plots;
if isempty(p), return; end
pulse = {'short', 'long'};
for i = 1:numel(c.targets)
    tg = c.targets(i);
    q = p(strcmp({p.pulse}, pulse{i}));
    for j = 1:numel(q)
        [pos, vel] = target_state(tg, q(j).time_s);
        R = norm(pos - [0 0 c.ant.height_m]);
        vr = dot(vel, (pos - [0 0 c.ant.height_m]) / R);
        bin = abs(diff(r.sys.dop.f_Hz(1:2))) * r.sys.wave.(pulse{i}).lambda_m / 2;
        if abs(q(j).range_m - R) < 100 && abs(q(j).vr_mps - vr) < 1.5 * bin
            found = [found tg.name]; %#ok<AGROW>
            break
        end
    end
end
end

function w = hann_ref(N)
w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
end

function names = flat_fields(s, prefix)
% Names of all leaf fields of a struct ('a.b.c'), except cfg.modes.
if nargin < 2, prefix = ''; end
names = {};
f = fieldnames(s);
for i = 1:numel(f)
    if isempty(prefix) && strcmp(f{i}, 'modes'), continue; end
    v = s.(f{i});
    if isstruct(v) && isscalar(v)
        names = [names, flat_fields(v, [prefix f{i} '.'])]; %#ok<AGROW>
    else
        names{end+1} = [prefix f{i}]; %#ok<AGROW>
    end
end
end

function ok = report(text, ok)
fprintf('  %-110s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

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

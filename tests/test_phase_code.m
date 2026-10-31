function ok = test_phase_code()
%TEST_PHASE_CODE  Phase codes from cfg.files.phaseCode, matched and mismatch filters.
%   1. The default file holds the P4 codes of both sub-pulses and their
%      mismatch filters (one coefficient per chip).
%   2. A target without noise through signal generation, ADC, DDC and pulse
%      compression (1 element, 2 and 4 MHz DDC output, both sub-pulses): with
%      the matched and with the mismatch filter the peak is on the range cell
%      of the target (the range axis takes the earlier peak of the longer
%      mismatch filter into account); the mismatch filter has lower range
%      sidelobes and loses less than 1.5 dB of SNR.
%   3. A user file (no 'madeFor') with a long-pulse mismatch filter designed
%      at the DDC output rate on the pulse after the DDC, for 4 sampling
%      phases (<pulse>MismatchRate_Hz, mismatch_filter.m): used as it is,
%      peak on the target cell, lower sidelobes than the chip-rate filter.
%   4. Mismatch filters of 3x the code length (cfg.wave.<pulse>.mismatchLength):
%      the same range cells and DDC output as with the matched filter, and a
%      target near the edges of both receive intervals (5 and 23 km, 24 and
%      79 km) is on its cell with lower sidelobes (within a pulse length of the
%      peak: beyond 52.5 km the short-pulse echo of the same target is in the
%      long-pulse interval) than with the matched filter and less than 1.5 dB
%      loss.
%   5. The whole chain (run_radar_sim) with phase codes and mismatch filters:
%      one plot of a target at 40 km at the right range (default and 3x
%      filters).
%   6. When the C receiver is built, its core test (radar_receiver/build/
%      test_core) passes with the exported configuration of the mismatch
%      filters (default and 3x).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
cfg = radar_config();
cfg.console.enable = false;
cfg.display.enable = false;
cfg.io.saveRaw = false;
cfg.ant.nElements = 1;
cfg.clutter.enable = false;
for p = {'short', 'long'}
    cfg.wave.(p{1}).type = 'phasecode';
end

%% 1) default file
f = load(cfg.files.phaseCode);
good = true;
for p = {'short', 'long'}
    N = round(cfg.wave.(p{1}).width_s / cfg.wave.(p{1}).chip_s);
    m = (0:N-1)';
    good = good && max(abs(f.([p{1} 'Code']) - exp(1j*(pi*m.^2/N - pi*m)))) < 1e-12 ...
                && numel(f.([p{1} 'Mismatch'])) > N;
end
ok = report('default file: P4 codes and mismatch filters of both sub-pulses', good) && ok;

%% 2) a target through the chain, matched and mismatch filter
R = struct('short', 15e3, 'long', 45e3);
for rate = [2e6 4e6]
    res = struct();
    for comp = {'matched', 'mismatch'}
        c = cfg;
        c.ddc.outputRate_Hz = rate;
        for p = {'short', 'long'}, c.wave.(p{1}).compression = comp{1}; end
        res.(comp{1}) = compress_target(c, R);
    end
    for p = {'short', 'long'}
        a = res.matched.(p{1});  b = res.mismatch.(p{1});
        loss = 10*log10(a.peak / b.peak);              % unit noise gain filters
        good = a.cellErr <= 0.5 && b.cellErr <= 0.5 && b.psl < a.psl - 1 && loss < 1.5;
        ok = report(sprintf(['%g MHz %-5s: target cell error matched %.2f, mismatch %.2f; PSL %.1f -> %.1f dB, ' ...
                             'loss %.2f dB'], rate/1e6, p{1}, a.cellErr, b.cellErr, a.psl, b.psl, loss), good) && ok;
    end
end

%% 3) a user file: mismatch filter designed at the DDC rate on the received long pulse
cDef = cfg;                                           % mismatch filters of the default file
for p = {'short', 'long'}, cDef.wave.(p{1}).compression = 'mismatch'; end
sys = setup_radar(cDef);
u = rmfield(f, 'madeFor');
u.longMismatch = mismatch_filter(ddc_pulse(sys, sys.wave.long, 0:0.25:0.75), round(5*sys.wave.long.chip_s*sys.fsOut), 3);
u.longMismatchRate_Hz = sys.fsOut;
cUser = cDef;
cUser.files.phaseCode = [tempname() '.mat'];
save(cUser.files.phaseCode, '-struct', 'u');
userRes = compress_target(cUser, R);
chipRes = compress_target(cDef, R);
f2 = load(cUser.files.phaseCode);
delete(cUser.files.phaseCode);
matchedPsl = compress_target(cfg, R);
good = ~isfield(f2, 'madeFor') && userRes.long.cellErr <= 0.5 && userRes.long.psl < chipRes.long.psl - 3;
ok = report(sprintf(['user file, long-pulse filter at the DDC rate: PSL %.1f dB (matched %.1f, chip rate %.1f), ' ...
                     'cell error %.2f'], userRes.long.psl, matchedPsl.long.psl, chipRes.long.psl, userRes.long.cellErr), ...
            good) && ok;

%% 4) mismatch filters of 3x the code length
c3 = cDef;
c3.files.phaseCode = [tempname() '.mat'];             % default file built for this length
for p = {'short', 'long'}
    c3.wave.(p{1}).mismatchLength = 3 * round(c3.wave.(p{1}).width_s / c3.wave.(p{1}).chip_s);
end
sysM = setup_radar(cfg);                              % matched filters
sys3 = setup_radar(c3);
good = sysM.nDdc == sys3.nDdc;
for p = {'short', 'long'}
    good = good && isequal(sys3.wave.(p{1}).range_m, sysM.wave.(p{1}).range_m) ...
                && isequal(sys3.wave.(p{1}).rxIdx, sysM.wave.(p{1}).rxIdx) ...
                && numel(sys3.wave.(p{1}).mfTaps) == 3 * numel(sysM.wave.(p{1}).mfTaps);
end
ok = report(sprintf(['3x mismatch filters (%d and %d taps): same DDC output (%d samples), receive intervals and ' ...
                     'range cells (%.2f ... %.2f km, %.2f ... %.2f km) as the matched filters'], ...
                    numel(sys3.wave.short.mfTaps), numel(sys3.wave.long.mfTaps), sys3.nDdc, ...
                    sys3.wave.short.range_m([1 end])/1e3, sys3.wave.long.range_m([1 end])/1e3), good) && ok;
for edge = {struct('short', 5e3, 'long', 24e3), struct('short', 23e3, 'long', 79e3)}
    a = compress_target(cfg, edge{1});
    b = compress_target(c3, edge{1});
    for p = {'short', 'long'}
        loss = 10*log10(a.(p{1}).peak / b.(p{1}).peak);
        good = b.(p{1}).cellErr <= 0.5 && b.(p{1}).pslNear < a.(p{1}).pslNear - 1 && loss < 1.5;
        ok = report(sprintf('3x mismatch, %-5s target at %4.1f km: cell error %.2f, PSL %.1f dB (matched %.1f dB), loss %.2f dB', ...
                            p{1}, edge{1}.(p{1})/1e3, b.(p{1}).cellErr, b.(p{1}).pslNear, a.(p{1}).pslNear, loss), good) && ok;
    end
end

%% 5) the whole chain with mismatch filters (default length and 3x)
c = cDef;
c.targets = make_target('A', 10, 40e3, 20, 1500, 120, 0);   % 113 m/s radial (outside the MTI notch)
c.sim.sectors_deg = [14 26];
c.track.enable = false;
c.io.outFolder = fullfile(root, 'output', 'test_phase_code');
c3.targets = c.targets;  c3.sim.sectors_deg = c.sim.sectors_deg;
c3.track.enable = false;  c3.io.outFolder = c.io.outFolder;
for run = {{c, 'mismatch filters'}, {c3, '3x mismatch filters'}}
    res = run_radar_sim(run{1}{1});
    [Rt, azt] = target_geometry(c.targets, 20/res.sys.ant.rate_dps, c.ant.height_m);
    good = numel(res.plots) == 1 && abs(res.plots(1).range_m - Rt) < 75 && abs(res.plots(1).az_deg - azt) < 0.5;
    if isempty(res.plots), txt = 'no plot'; else, txt = sprintf('%.1f m from the target', res.plots(1).range_m - Rt); end
    ok = report(sprintf('run_radar_sim, phase codes with %s: %d plot(s), %s', run{1}{2}, numel(res.plots), txt), ...
                good) && ok;
end

%% 6) the C receiver core test with the exported mismatch configurations
bin = fullfile(root, 'radar_receiver', 'build', 'test_core');
if ispc, bin = [bin '.exe']; end
if exist(bin, 'file')
    for run = {{c, 'mismatch filters'}, {c3, '3x mismatch filters'}}
        fc = [tempname() '.bin'];
        evalc('export_receiver_config(run{1}{1}, fc)');
        [status, out] = system(sprintf('"%s" "%s"', bin, fc));
        delete(fc);
        tok = regexp(out, 'plot: long, range ([\d.]+) m', 'tokens', 'once');
        ok = report(sprintf('C receiver (test_core) with the %s: %s', run{1}{2}, ...
                            strtrim(regexprep(out(max(1, end-60):end), '\s+', ' '))), status == 0 && ~isempty(tok)) && ok;
    end
else
    fprintf('  C receiver not built: test_core skipped\n');
end
delete(c3.files.phaseCode);
fprintf('Phase codes, matched and mismatch filters: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function out = compress_target(cfg, R)
% A target at R.short and R.long without noise: range profile of each sub-pulse.
cfg.rx.addThermalNoise = false;  cfg.rx.addLeakage = false;  cfg.adc.addAdcNoise = false;
cfg.adc.quantize = false;        cfg.ddc.fixedPoint = false;
sys = setup_radar(cfg);
for p = {'short', 'long'}
    wv = sys.wave.(p{1});
    tgt = make_target('T', 10, R.(p{1}), 0, cfg.ant.height_m, 0, 0);
    tgt = draw_target_rcs(tgt, 2, 0);
    [Rt, ~, ~] = target_geometry(tgt, 0, cfg.ant.height_m);
    [I, Q] = ddc(adc_sample(generate_rx_signal(sys, tgt, 0, 0, []), sys), sys);
    bb = complex(I, Q);
    rp = abs(matched_filter(ddc_segment(bb, wv).', wv)).^2;
    [pk, k] = max(rp);
    cellSize = abs(diff(wv.range_m(1:2)));
    main = round(1.5 * wv.chip_s * sys.fsOut);       % mainlobe and its first sample around the peak
    sl = rp;  sl(max(1, k-main):min(end, k+main)) = 0;
    Lp = round(wv.width_s * sys.fsOut);              % sidelobes within a pulse length of the peak
    near = sl(max(1, k-Lp):min(end, k+Lp));
    out.(p{1}) = struct('peak', pk, 'cellErr', abs(wv.range_m(k) - Rt) / cellSize, ...
                        'psl', 10*log10(max(sl) / pk), 'pslNear', 10*log10(max(near) / pk));
end
end

function ok = report(text, ok)
fprintf('  %-100s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

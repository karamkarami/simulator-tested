function ok = test_short_pulse_ghost()
%TEST_SHORT_PULSE_GHOST  A strong far target must not create a ghost plot.
%
%   Both sub-pulses share the 30 MHz IF. The short-pulse echo of a target
%   beyond c*(tTxLong + widthLong)/2 = 52.5 km therefore also falls inside the
%   long-pulse receive window. The long-pulse matched filter does not match
%   it: the result is a response spread over about 10 km of range, 30 ... 52 km
%   in front of the target, up to about 24 dB below the target. The 2-D CA-CFAR
%   adapts to it like to a local rise of the noise floor, so no plot should
%   appear there. The test runs one 100 m^2 target at 65 km (about 47 dB SNR,
%   fixed PRF, no clutter, raw data not saved) and checks that the target is
%   detected and that no plot lies in the ghost region.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));

cfg = radar_config();
cfg.clutter.enable = false;
cfg.io.saveRaw = false;
cfg.display.enable = false;
cfg.sim.verbose = false;
cfg.io.outFolder = fullfile(root, 'output', 'test_ghost');
cfg.targets = make_target('G', 100, 65e3, 40, 3500, -170, 0);
sys = setup_radar(cfg);
cfg.sim.sectors_deg = target_sectors(sys, cfg.targets);
ghostShift_m = cfg.const.c * (sys.wave.long.tTx_s - sys.wave.short.tTx_s) / 2;

res = run_radar_sim(cfg);
nTarget = 0;
nGhost  = 0;
for p = res.plots
    R = target_geometry(cfg.targets, p.time_s, cfg.ant.height_m);
    nTarget = nTarget + (abs(p.range_m - R) < 300);
    nGhost  = nGhost + (strcmp(p.pulse, 'long') && ...
              abs(R - ghostShift_m - p.range_m) <= cfg.const.c*sys.wave.long.width_s/2);
end
ok = nTarget >= 1 && nGhost == 0;
fprintf('Short-pulse ghost test: %d plot(s), %d on the target, %d in the ghost region: %s\n', ...
    numel(res.plots), nTarget, nGhost, ternary(ok, 'PASS', 'FAIL'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function ok = test_peak_interpolation()
%TEST_PEAK_INTERPOLATION  Interpolation of the detected peaks (peak_interpolation.m).
%   A peak that is Gaussian in Doppler and across the beams (a parabola in
%   dB) gives its exact fractional Doppler bin (radial velocity) and beam
%   (elevation) with cfg.report.dopplerInterp / elevationInterp; the Doppler
%   bin wraps at the edge of the axis, the lowest beam keeps its centre,
%   and without the options the centres of the bin and of the beam stay. The
%   C receiver's test_cfar checks the same in processor.c.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
cfg = radar_config();
cfg.clutter.enable = false;
cfg.report.dopplerInterp = true;
cfg.report.elevationInterp = true;
sys = setup_radar(cfg);
nB = sys.nBeams;  nfft = sys.dop.nfft;  nR = 20;  r = 10;
df = 1 / (nfft * sys.dop.rowPri_s);
ok = true;

for on = [true false]
    sys.report.dopplerInterp = on;
    sys.report.elevationInterp = on;
    d0 = nfft/4 + 0.3;  b0 = round(nB/2) - 0.35;
    rep = interp_one(sys, d0, b0, nR, r);
    k = round(d0);
    vrTrue = -(sys.dop.f_Hz(k) + (d0 - k) * df) * sys.lambda / 2;
    vrBin  = -sys.dop.f_Hz(k) * sys.lambda / 2;
    elTrue = sys.ant.tilt_deg + asind(interp1(1:nB, sys.dbf.beamU, b0));
    elBeam = sys.ant.tilt_deg + asind(sys.dbf.beamU(round(b0)));
    if on
        ok = report(abs(rep.vr_mps - vrTrue) < 1e-3 && abs(rep.el_deg - elTrue) < 1e-3, ...
                    sprintf('fractional Doppler bin and beam (vr %.3f m/s, el %.3f deg)', rep.vr_mps, rep.el_deg)) && ok;
        rep = interp_one(sys, 1 - 0.2, b0, nR, r);          % the first bin: its neighbour is the last one
        ok = report(abs(rep.vr_mps + (sys.dop.f_Hz(1) - 0.2*df) * sys.lambda / 2) < 1e-3, ...
                    'the Doppler axis wraps at its edge') && ok;
        rep = interp_one(sys, d0, 1.2, nR, r);
        ok = report(abs(rep.el_deg - (sys.ant.tilt_deg + asind(sys.dbf.beamU(1)))) < 1e-9, ...
                    'the lowest beam keeps its centre') && ok;
    else
        ok = report(abs(rep.vr_mps - vrBin) < 1e-9 && abs(rep.el_deg - elBeam) < 1e-9, ...
                    'options off: centres of the Doppler bin and of the beam') && ok;
    end
end
% several reports, one of them at the lowest beam (only the others are interpolated in elevation)
sys.report.dopplerInterp = true;  sys.report.elevationInterp = true;
pk = [nfft/4 + 0.3, round(nB/2) - 0.35, 5;  nfft/2 + 7.2, 1.2, 9;  nfft/4 - 3.4, round(nB/2) + 2.3, 14];
P = zeros(nB, nR, nfft, 'single');
for i = 1:size(pk, 1)
    for dd = 1:nfft
        P(:, pk(i, 3), dd) = single(1e4 * exp(-0.3*(dd - pk(i, 1))^2 - 0.4*((1:nB) - pk(i, 2)).^2));
    end
end
k = round(pk(:, 1));  bb = round(pk(:, 2));
dets = struct('beam', bb.', 'rIdx', pk(:, 3).', 'dIdx', k.', 'power', 1, 'noise', 1);
wv = sys.wave.long;  wv.range_m = (1:nR) * 75;
rep = detection_reports(dets, P, sys, wv, struct('tCenter_s', 0, 'azCenter_deg', 0));
vr = -(sys.dop.f_Hz(k).' + (pk(:, 1) - k) * df) * sys.lambda / 2;
el = sys.ant.tilt_deg + asind(interp1(1:nB, sys.dbf.beamU, pk(:, 2)));
el(2) = sys.ant.tilt_deg + asind(sys.dbf.beamU(1));
ok = report(numel(rep) == 3 && max(abs([rep.vr_mps].' - vr)) < 1e-3 && max(abs([rep.el_deg].' - el)) < 1e-3, ...
            'three reports in different cells (one at the lowest beam)') && ok;
fprintf('Peak interpolation: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function rep = interp_one(sys, d0, b0, nR, r)
% Reports of one Gaussian peak at the fractional Doppler bin d0 and beam b0 (1-based).
nB = sys.nBeams;  nfft = sys.dop.nfft;
P = zeros(nB, nR, nfft, 'single');
for d = 1:nfft
    dd = abs(d - d0);  dd = min(dd, nfft - dd);
    P(:, r, d) = single(1e4 * exp(-0.3*dd^2 - 0.4*((1:nB) - b0).^2));
end
k = mod(round(d0) - 1, nfft) + 1;  b = round(b0);
dets = struct('beam', b, 'rIdx', r, 'dIdx', k, 'power', P(b, r, k), 'noise', 1);
info = struct('tCenter_s', 0, 'azCenter_deg', 0);
wv = sys.wave.long;
wv.range_m = (1:nR) * 75;
rep = detection_reports(dets, P, sys, wv, info);
end

function ok = report(cond, msg)
ok = logical(cond);
fprintf('  %-62s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

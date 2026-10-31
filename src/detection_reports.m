function rep = detection_reports(dets, P, sys, wv, info)
%DETECTION_REPORTS  Convert pruned CFAR detections to physical quantities.
%   rep = DETECTION_REPORTS(dets, P, sys, wv, info) returns a struct array with
%   range, radial velocity (ambiguous, + = away), elevation (interpolated
%   between adjacent beams; NaN without beamforming), azimuth and time of the
%   CPI centre and the SNR.
%   The SNR uses the noise floor of the beam (median of its range-Doppler map
%   scaled to the mean), which, unlike the CFAR noise estimate, is not raised
%   by the range sidelobes of strong targets.

rep = struct('time_s', {}, 'az_deg', {}, 'range_m', {}, 'el_deg', {}, 'vr_mps', {}, ...
             'snr_dB', {}, 'power', {}, 'beam', {}, 'rIdx', {}, 'dIdx', {}, 'pulse', {});
nB = sys.nBeams;
if isempty(dets.beam), return; end
noiseFloor = median(reshape(P, nB, []), 2) / sys.dop.medianToMean;   % per beam
b = dets.beam(:);  rI = dets.rIdx(:);  dI = dets.dIdx(:);
n = numel(b);
% Parabolic interpolation of the beam powers (dB) around the strongest beam
u = sys.dbf.beamU(b);  u = u(:);
mid = find(b > 1 & b < nB);
if ~isempty(mid)
    [nB_, nR_] = size(P(:, :, 1));
    lin = @(bb) bb + nB_*(rI(mid) - 1) + nB_*nR_*(dI(mid) - 1);
    p1 = double(10*log10(P(lin(b(mid) - 1))));
    p2 = double(10*log10(P(lin(b(mid)))));
    p3 = double(10*log10(P(lin(b(mid) + 1))));
    den = p1(:) - 2*p2(:) + p3(:);
    delta = 0.5*(p1(:) - p3(:)) ./ den;
    ok = den < 0 & isfinite(delta);                  % (a zero power gives -Inf dB)
    if any(ok)
        u(mid(ok)) = interp1(1:nB, sys.dbf.beamU, b(mid(ok)) + delta(ok));
    end
end
snr = 10*log10(max(double(dets.power(:) ./ noiseFloor(b)) - 1, eps));
rep = struct('time_s', info.tCenter_s, 'az_deg', info.azCenter_deg, ...
             'range_m', num2cell(reshape(wv.range_m(rI), 1, [])), ...
             'el_deg', num2cell(sys.ant.tilt_deg + asind(u.')), ...
             'vr_mps', num2cell(reshape(-sys.dop.f_Hz(dI) * wv.lambda_m / 2, 1, [])), ...
             'snr_dB', num2cell(snr.'), 'power', num2cell(double(dets.power(:)).'), ...
             'beam', num2cell(b.'), 'rIdx', num2cell(rI.'), 'dIdx', num2cell(dI.'), 'pulse', wv.name);
rep = reshape(rep, 1, n);
end

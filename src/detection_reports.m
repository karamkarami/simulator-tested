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
for i = 1:numel(dets.beam)
    b = dets.beam(i);
    % Parabolic interpolation of the beam powers (dB) around the strongest beam
    u = sys.dbf.beamU(b);
    if b > 1 && b < nB
        pdB = double(10*log10(P(b-1:b+1, dets.rIdx(i), dets.dIdx(i))));
        den = pdB(1) - 2*pdB(2) + pdB(3);
        delta = 0.5*(pdB(1) - pdB(3))/den;
        if den < 0 && isfinite(delta)                % (a zero power gives -Inf dB)
            u = interp1(1:nB, sys.dbf.beamU, b + delta);
        end
    end
    r.time_s  = info.tCenter_s;
    r.az_deg  = info.azCenter_deg;
    r.range_m = wv.range_m(dets.rIdx(i));
    r.el_deg  = sys.ant.tilt_deg + asind(u);
    r.vr_mps  = -sys.dop.f_Hz(dets.dIdx(i)) * wv.lambda_m / 2;
    r.snr_dB  = 10*log10(max(double(dets.power(i)/noiseFloor(b)) - 1, eps));
    r.power   = double(dets.power(i));
    r.beam    = b;
    r.rIdx    = dets.rIdx(i);
    r.dIdx    = dets.dIdx(i);
    r.pulse   = wv.name;
    rep(end+1) = r; %#ok<AGROW>
end
end

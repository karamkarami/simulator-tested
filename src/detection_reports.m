function rep = detection_reports(dets, P, sys, wv, info, interp)
%DETECTION_REPORTS  Convert pruned CFAR detections to physical quantities.
%   rep = DETECTION_REPORTS(dets, P, sys, wv, info) returns a struct array with
%   range, radial velocity (ambiguous, + = away), elevation (NaN without
%   beamforming), azimuth and time of the CPI centre and the SNR. The velocity
%   and the elevation are interpolated around the peak as set in
%   cfg.report.dopplerInterp / elevationInterp (peak_interpolation.m).
%   DETECTION_REPORTS(..., false) returns the centres of the Doppler bin and of
%   the beam instead (run_radar_sim.m then calls peak_interpolation.m itself,
%   to time it as a block of its own).
%   The SNR uses the noise floor of the beam (median of its range-Doppler map
%   scaled to the mean), which, unlike the CFAR noise estimate, is not raised
%   by the range sidelobes of strong targets; with cfg.report.snrNoise =
%   'cfar' it uses the CFAR noise estimate of the cell (dets.noise) instead,
%   and the median is not computed.

rep = struct('time_s', {}, 'az_deg', {}, 'range_m', {}, 'el_deg', {}, 'vr_mps', {}, ...
             'snr_dB', {}, 'power', {}, 'beam', {}, 'rIdx', {}, 'dIdx', {}, 'pulse', {});
nB = sys.nBeams;
if isempty(dets.beam), return; end
b = dets.beam(:);  rI = dets.rIdx(:);  dI = dets.dIdx(:);
if isfield(sys, 'report') && strcmp(sys.report.snrNoise, 'cfar')
    noise = dets.noise(:);                                             % CFAR window of the cell
else
    noiseFloor = median(reshape(P, nB, []), 2) / sys.dop.medianToMean;   % per beam
    noise = noiseFloor(b);
end
n = numel(b);
u = sys.dbf.beamU(b);  u = u(:);                           % beam centres (NaN without DBF)
snr = 10*log10(max(double(dets.power(:) ./ noise) - 1, eps));
rep = struct('time_s', info.tCenter_s, 'az_deg', info.azCenter_deg, ...
             'range_m', num2cell(reshape(wv.range_m(rI), 1, [])), ...
             'el_deg', num2cell(sys.ant.tilt_deg + asind(u.')), ...
             'vr_mps', num2cell(reshape(-sys.dop.f_Hz(dI) * wv.lambda_m / 2, 1, [])), ...
             'snr_dB', num2cell(snr.'), 'power', num2cell(double(dets.power(:)).'), ...
             'beam', num2cell(b.'), 'rIdx', num2cell(rI.'), 'dIdx', num2cell(dI.'), 'pulse', wv.name);
rep = reshape(rep, 1, n);
if nargin < 6 || interp
    rep = peak_interpolation(rep, P, sys);
end
end

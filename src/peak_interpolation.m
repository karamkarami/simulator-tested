function rep = peak_interpolation(rep, P, sys)
%PEAK_INTERPOLATION  Refine the radial velocity and the elevation of the detected peaks.
%   rep = PEAK_INTERPOLATION(rep, P, sys) takes the reports of one CPI update
%   (detection_reports.m: velocity of the centre of the Doppler bin, elevation
%   of the centre of the beam) and the range-Doppler power maps P [nBeams x
%   nRange x nfft] they come from, and interpolates every peak with a parabola
%   through the powers [dB] of the peak cell and its two neighbours:
%
%     sys.report.dopplerInterp   (cfg.report.dopplerInterp): along Doppler,
%         bins d-1, d, d+1 (circular: the Doppler axis wraps). The vertex
%         d + delta, delta = 0.5*(p(d-1) - p(d+1)) / (p(d-1) - 2 p(d) + p(d+1)),
%         gives the Doppler frequency f(d) + delta * bin width and the radial
%         velocity -f * lambda / 2 (still ambiguous: modulo the Doppler span).
%     sys.report.elevationInterp (cfg.report.elevationInterp): across the
%         beams b-1, b, b+1, when the peak beam has both neighbours (DBF, not
%         the lowest or highest beam). The fractional beam b + delta is turned
%         into u = sin(angle from broadside) by linear interpolation of the beam
%         directions sys.dbf.beamU (as interp1), then into the elevation.
%
%   A peak is interpolated only when the parabola has a maximum (negative
%   curvature) and all three powers are positive; otherwise the bin or beam
%   centre is kept. The pruned detections are local peaks, so |delta| <= 0.5.
%   The C receiver does the same (processor.c, peak_interpolation()).

if isempty(rep), return; end
doDop = isfield(sys, 'report') && isfield(sys.report, 'dopplerInterp') && sys.report.dopplerInterp;
doEl  = (~isfield(sys, 'report') || ~isfield(sys.report, 'elevationInterp') || sys.report.elevationInterp) ...
        && sys.dbf.enabled && sys.nBeams > 2;
if ~doDop && ~doEl, return; end
[nB, nR, nfft] = size(P);
b = [rep.beam];  r = [rep.rIdx];  d = [rep.dIdx];
lin = @(bb, dd) bb + nB*(r - 1) + nB*nR*(dd - 1);
if doDop
    [delta, ok] = vertex(P, lin(b, mod(d - 2, nfft) + 1), lin(b, d), lin(b, mod(d, nfft) + 1));
    df = 1 / (nfft * sys.dop.rowPri_s);                        % Doppler bin width [Hz]
    for i = find(ok)
        rep(i).vr_mps = -(sys.dop.f_Hz(d(i)) + delta(i) * df) * sys.lambda / 2;
    end
end
if doEl
    mid = find(b > 1 & b < sys.nBeams);
    [delta, ok] = vertex(P, lin(b(mid) - 1, d(mid)), lin(b(mid), d(mid)), lin(b(mid) + 1, d(mid)));
    for i = mid(ok)
        k = find(mid == i);
        u = interp1(1:sys.nBeams, sys.dbf.beamU, b(i) + delta(k));
        rep(i).el_deg = sys.ant.tilt_deg + asind(u);
    end
end
end

function [delta, ok] = vertex(P, i1, i2, i3)
% Vertex offset of the parabola through the powers [dB] at the linear indices i1, i2 (peak), i3.
p1 = 10*log10(double(P(i1)));
p2 = 10*log10(double(P(i2)));
p3 = 10*log10(double(P(i3)));
den = p1(:).' - 2*p2(:).' + p3(:).';
delta = 0.5*(p1(:).' - p3(:).') ./ den;
ok = den < 0 & isfinite(delta);                             % (a zero power gives -Inf dB)
end

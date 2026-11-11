function [isDet, noise] = cfar_detect(P, sys, wv)
%CFAR_DETECT  CFAR detection on every range-Doppler map: the detector of cfg.cfar.type and window of cfg.cfar.method.
%   [isDet, noise] = CFAR_DETECT(P, sys, wv) tests the cells of P
%   [nBeams x nRange x nDoppler] of sub-pulse wv (only the range cells
%   wv.testIdx(1):wv.testIdx(2)) with the detector sys.cfar.type:
%     'ca'  cfar_ca.m    cell averaging: the mean of the training cells
%     'go'  cfar_goso.m  greatest of: the larger of the means of the two
%                        halves of the window (leading / lagging)
%     'so'  cfar_goso.m  smallest of: the smaller of them
%     'os'  cfar_os.m    order statistic: the osK-th smallest training cell
%   and the window sys.cfar.method ('2d', 'range', 'doppler'; guard and
%   training cells on each side, Doppler circular, the range window shifted
%   inwards at the edges). A cell is detected when P > sys.cfar.alpha times
%   the statistic. noise is the noise estimate of the tested cells: the
%   statistic / sys.cfar.noiseScale (its mean for unit noise: 1 for CA), 0
%   elsewhere (OS: only at its detections). The MTI notch is applied only
%   through cfg.cfar.notchMaxRange_km (sys.cfar.notchMaxRange_m); every Doppler bin is tested beyond it.
nD = size(P, 3);
zeroBin = floor(nD/2) + 1;                         % zero Doppler after fftshift
notch = max(1, zeroBin-sys.dop.notchBins):min(nD, zeroBin+sys.dop.notchBins);
switch sys.cfar.type
    case 'ca'
        [isDet, noise] = cfar_ca(P, sys, wv);
    case {'go', 'so'}
        [isDet, noise] = cfar_goso(P, sys, wv);
    case 'os'
        [isDet, noise] = cfar_os(P, sys, wv);
    otherwise
        error('Unknown CFAR type %s.', sys.cfar.type);
end
if isfield(wv, 'range_m')
    inNotchRange = reshape(wv.range_m, 1, []) <= sys.cfar.notchMaxRange_m;
else
    inNotchRange = true(1, size(P, 2));           % legacy callers without a range axis
end
if any(inNotchRange)
    isDet(:, inNotchRange, notch) = false;
    noise(:, inNotchRange, notch) = 0;
end
end

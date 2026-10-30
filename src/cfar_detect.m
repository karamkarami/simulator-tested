function [isDet, noise] = cfar_detect(P, sys, wv)
%CFAR_DETECT  Cell-averaging CFAR on every range-Doppler map (2-D, 1-D range or 1-D Doppler).
%   [isDet, noise] = CFAR_DETECT(P, sys, wv) tests the cells of P
%   [nBeams x nRange x nDoppler] of sub-pulse wv. The noise of a cell is the
%   mean of the training cells of its window (sys.cfar: guardRange/trainRange
%   and guardDoppler/trainDoppler on each side; the 1-D methods have no window
%   in the other dimension), excluding the guard cells. A cell is detected
%   when P > sys.cfar.alpha * noise.
%   Only the range cells wv.testIdx(1):wv.testIdx(2) are tested (and only
%   the cells around them are used). The Doppler axis is circular. At the
%   first and last range cells the range window is shifted inwards to stay
%   inside the data: it keeps its length, the cell under test is no longer in
%   its middle and its guard cells are cut at the edge (so it has a few more
%   training cells there). The sys.dop.notchBins Doppler bins on each side of
%   zero (MTI notch, clutter residue) are not tested. noise holds the noise
%   estimate of the tested cells (0 elsewhere).

c  = sys.cfar;
gR = c.guardRange;    oR = gR + c.trainRange;      % half sizes of the guard and outer windows
gD = c.guardDoppler;  oD = gD + c.trainDoppler;
[nB, nR, nD] = size(P);
isDet = false(nB, nR, nD);
noise = zeros(nB, nR, nD, class(P));
t = wv.testIdx;
if t(2) < t(1), return; end
% range cells the windows of the tested cells need (also the shifted windows
% of tested cells within oR cells of an edge)
lo = max(1, min(t(1), nR - oR) - oR);
hi = min(nR, max(t(2), oR + 1) + oR);
Q = P(:, lo:hi, :);
sumOuter = box_sum(Q, oR, oD);                     % centred sums (zero padded in range)
sumGuard = box_sum(Q, gR, gD);
r  = t(1):t(2);
rc = min(max(r, oR + 1), nR - oR);                 % centre of the (shifted) outer window
nGuard = (min(r + gR, nR) - max(r - gR, 1) + 1) * (2*gD + 1);   % guard cells inside the data
nTrain = (2*oR + 1) * (2*oD + 1) - nGuard;          % = sys.cfar.nTrain away from the edges
est = (sumOuter(:, rc - lo + 1, :) - sumGuard(:, r - lo + 1, :)) ./ reshape(nTrain, 1, [], 1);
isDet(:, r, :) = P(:, r, :) > c.alpha * est;
noise(:, r, :) = est;
zeroBin = floor(nD/2) + 1;                         % zero Doppler after fftshift
isDet(:, :, max(1, zeroBin-sys.dop.notchBins):min(nD, zeroBin+sys.dop.notchBins)) = false;
end

function S = box_sum(P, hr, hd)
% Sum of P over a (2*hr+1) x (2*hd+1) window centred on every cell (range x Doppler).
S = convn(P, ones(1, 2*hr+1), 'same');                         % range (zero padded)
S = cat(3, S(:, :, end-hd+1:end), S, S(:, :, 1:hd));           % circular padding in Doppler
S = convn(S, ones(1, 1, 2*hd+1), 'valid');
end

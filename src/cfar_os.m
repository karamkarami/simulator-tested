function [isDet, noise] = cfar_os(P, sys, wv, notch)
%CFAR_OS  Order-statistic (OS) CFAR (cfg.cfar.type = 'os'), called by cfar_detect.m.
%   [isDet, noise] = CFAR_OS(P, sys, wv, notch): the statistic of a tested
%   cell is the k-th smallest of the N training cells of its window
%   (sys.cfar, as cfar_ca.m: guard cells removed, the range window shifted
%   inwards at the edges, Doppler circular), k = sys.cfar.osK for the full
%   window of sys.cfar.nTrain cells and round(osK*N/nTrain) at the range
%   edges (N larger: the guard cells cut). A cell is detected when
%   P > alpha * X(k), i.e. when at least k training cells t have alpha*t < P
%   (in double).
%   Only the cells that can survive the pruning are tested: the local maxima
%   of their 3 x 3 x 3 neighbourhood in beam x range x Doppler (as
%   prune_peaks.m, Doppler circular), in the tested range cells and outside
%   the Doppler bins notch (the MTI notch). Every other cell is not detected;
%   the result after the pruning is the same as with every cell tested, for
%   a fraction of the work (the C receiver does the same). sys.cfar.osAllCells
%   = true (tests only) tests every cell.
%   noise = X(k) / sys.cfar.noiseScale at the detected cells (0 elsewhere).
c  = sys.cfar;
gR = c.guardRange;    oR = gR + c.trainRange;
gD = c.guardDoppler;  oD = gD + c.trainDoppler;
[nB, nR, nD] = size(P);
isDet = false(nB, nR, nD);
noise = zeros(nB, nR, nD, class(P));
t = wv.testIdx;
if t(2) < t(1), return; end
% candidates: local maxima of the tested cells outside the notch
M = local_max(P);
cand = false(nB, nR, nD);
if isfield(c, 'osAllCells') && c.osAllCells          % tests (false alarm rate): every tested cell
    cand(:, t(1):t(2), :) = true;
else
    cand(:, t(1):t(2), :) = P(:, t(1):t(2), :) >= M(:, t(1):t(2), :);
end
cand(:, :, notch) = false;
idx = find(cand);
if isempty(idx), return; end
[b, r, d] = ind2sub([nB nR nD], idx);
% window offsets around the window centre rc (range) and the cell (Doppler)
[offR, offD] = ndgrid(-oR:oR, -oD:oD);
offR = offR(:).';  offD = offD(:).';
N = c.nTrain;
alpha = c.alpha;
block = 20000;                                       % candidates per block (memory)
for s = 1:block:numel(idx)
    k = s:min(s + block - 1, numel(idx));
    rr = r(k);  bb = b(k);  dd = d(k);
    rc = min(max(rr, oR + 1), nR - oR);
    ri = rc + offR;                                  % rows of the window [nk x nWin]
    di = mod(dd - 1 + offD, nD) + 1;
    guard = abs(ri - rr) <= gR & abs(offD) <= gD;
    V = double(P(bb + nB*(ri - 1) + nB*nR*(di - 1)));
    pc = double(P(idx(k)));
    nGuard = (min(rr + gR, nR) - max(rr - gR, 1) + 1) * (2*gD + 1);
    Nr = numel(offR) - nGuard;                       % training cells of the cell
    kr = max(1, min(Nr, round(c.osK * Nr / N)));
    cnt = sum(alpha * V < pc & ~guard, 2);
    det = cnt >= kr;
    isDet(idx(k(det))) = true;
    for j = find(det).'                              % X(k) of the detections: the noise estimate
        v = sort(V(j, ~guard(j, :)));
        noise(idx(k(j))) = v(kr(j)) / c.noiseScale;
    end
end
end

function M = local_max(P)
% Maximum over the 3 x 3 x 3 neighbourhood (beam, range: inside the map; Doppler circular).
[nB, nR, nD] = size(P);
M = P;
if nB > 1
    M = max(M, cat(1, M(2:end, :, :), -inf(1, nR, nD, class(P))));
    M = max(M, cat(1, -inf(1, nR, nD, class(P)), M(1:end-1, :, :)));
end
M2 = max(M, cat(2, M(:, 2:end, :), -inf(nB, 1, nD, class(P))));
M = max(M2, cat(2, -inf(nB, 1, nD, class(P)), M(:, 1:end-1, :)));
M = max(M, max(circshift(M, 1, 3), circshift(M, -1, 3)));
end

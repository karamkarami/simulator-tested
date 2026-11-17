function [isDet, noise] = cfar_ce(P, sys, wv)
%CFAR_CE  Censored (CE) CFAR (cfg.cfar.type = 'ce'), called by cfar_detect.m.
%   [isDet, noise] = CFAR_CE(P, sys, wv): the statistic of a tested cell is
%   the mean of its training cells without the largest ones (censored mean
%   level, CMLD): of the N training cells of its window (sys.cfar, as
%   cfar_ca.m: guard cells removed, the range window shifted inwards at the
%   edges, Doppler circular) the N - round(nCens * N / nTrain) smallest are
%   averaged, nCens = nTrain - sys.cfar.ceKeep the cells censored for the full
%   window (cfg.cfar.ceCensor of its training cells; N is larger at the range
%   edges, where the guard cells are cut). Other targets and clutter edges in
%   the window are the largest cells, so they do not raise the threshold. A
%   cell is detected when P > alpha * mean (in double); alpha from the Pfa of
%   the statistic's distribution (cfar_ce_table.m).
%   Only the cells that can survive the pruning are tested (as cfar_os.m: the
%   local maxima of their 3 x 3 x 3 neighbourhood in beam x range x Doppler,
%   in the tested range cells; cfar_detect masks the MTI notch afterwards).
%   sys.cfar.osAllCells = true (tests only) tests every cell.
%   noise = mean / sys.cfar.noiseScale at the detected cells (0 elsewhere).
c  = sys.cfar;
gR = c.guardRange;    oR = gR + c.trainRange;
gD = c.guardDoppler;  oD = gD + c.trainDoppler;
[nB, nR, nD] = size(P);
isDet = false(nB, nR, nD);
noise = zeros(nB, nR, nD, class(P));
t = wv.testIdx;
if t(2) < t(1), return; end
M = local_max(P);
cand = false(nB, nR, nD);
if isfield(c, 'osAllCells') && c.osAllCells
    cand(:, t(1):t(2), :) = true;
else
    cand(:, t(1):t(2), :) = P(:, t(1):t(2), :) >= M(:, t(1):t(2), :);
end
idx = find(cand);
if isempty(idx), return; end
[b, r, d] = ind2sub([nB nR nD], idx);
[offR, offD] = ndgrid(-oR:oR, -oD:oD);
offR = offR(:).';  offD = offD(:).';
N = c.nTrain;
nCens = N - c.ceKeep;
alpha = c.alpha;
block = 20000;
for s = 1:block:numel(idx)
    k = s:min(s + block - 1, numel(idx));
    rr = r(k);  bb = b(k);  dd = d(k);
    rc = min(max(rr, oR + 1), nR - oR);
    ri = rc + offR;
    di = mod(dd - 1 + offD, nD) + 1;
    guard = abs(ri - rr) <= gR & abs(offD) <= gD;
    V = double(P(bb + nB*(ri - 1) + nB*nR*(di - 1)));
    V(guard) = Inf;                                  % the guard cells go to the end of every sorted row
    pc = double(P(idx(k)));
    nGuard = (min(rr + gR, nR) - max(rr - gR, 1) + 1) * (2*gD + 1);
    Nr = numel(offR) - nGuard;                       % training cells of the cell
    keep = max(1, Nr - round(nCens * Nr / N));       % the smallest ones that are averaged
    S = cumsum(sort(V, 2), 2);
    mu = S(sub2ind(size(S), (1:numel(k)).', keep)) ./ keep;
    det = pc > alpha * mu;
    isDet(idx(k(det))) = true;
    noise(idx(k(det))) = mu(det) / c.noiseScale;
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

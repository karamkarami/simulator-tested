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
tested = t(1):t(2);
if isfield(c, 'osAllCells') && c.osAllCells
    cand = true(nB, numel(tested), nD);
else
    % Neighbours one cell beyond the tested interval still participate in peak picking.
    lo = max(1, t(1)-1);  hi = min(nR, t(2)+1);
    M = local_max(P(:, lo:hi, :));
    cand = P(:, tested, :) >= M(:, tested-lo+1, :);
end
idx = find(cand(:));
if isempty(idx), return; end
[b, r, d] = ind2sub([nB numel(tested) nD], idx);
r = r + t(1) - 1;
idx = b + nB*(r-1) + nB*nR*(d-1);
rc = min(max(r, oR+1), nR-oR);
% Cells with the same range shift have identical training offsets and censor count.
[shift, order] = sort(rc - r);
idx = idx(order);
groups = [1; find(diff(shift)) + 1; numel(idx) + 1];
[offR, offD] = ndgrid(-oR:oR, -oD:oD);
offR = offR(:);  offD = offD(:);
if oD > 0
    % Pad the needed range rows in Doppler so even wrapped windows use linear offsets.
    lo = min(rc)-oR;  hi = max(rc)+oR;
    padded = P(:, lo:hi, mod(-oD:nD+oD-1, nD) + 1);
else
    lo = 1;  hi = nR;
    padded = P;
end
strideD = nB*(hi-lo+1);
centre = b + nB*(rc-lo) + strideD*(d-1+oD);
centre = centre(order);
nCens = c.nTrain - c.ceKeep;
for group = 1:numel(groups)-1
    first = groups(group);  last = groups(group+1)-1;
    dr = shift(first);
    train = abs(dr + offR) > gR | abs(offD) > gD;
    offsets = nB*offR(train) + strideD*offD(train);
    nr = numel(offsets);
    keep = max(1, nr - round(nCens * nr / c.nTrain));
    % One window per column, with bounded temporary storage even for wide windows.
    block = max(1, floor(131072 / nr));
    for s = first:block:last
        k = s:min(s+block-1, last);
        cut = idx(k);
        values = reshape(padded(offsets + centre(k).'), nr, []);
        values = sort(values, 1);
        totals = cumsum(double(values(1:keep, :)), 1);
        mu = totals(end, :).' / keep;
        det = reshape(double(P(cut)), [], 1) > c.alpha * mu;
        isDet(cut(det)) = true;
        noise(cut(det)) = mu(det) / c.noiseScale;
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
if nD > 1
    M = max(M, max(circshift(M, 1, 3), circshift(M, -1, 3)));
end
end

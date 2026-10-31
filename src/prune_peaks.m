function dets = prune_peaks(isDet, P, noise, sys, wv)
%PRUNE_PEAKS  Keep only the CFAR detections that are real target peaks.
%   dets = PRUNE_PEAKS(isDet, P, noise, sys, wv) prunes the detections of a
%   range-Doppler-beam map P [nBeams x nRange x nDoppler] in two steps:
%   1) Peaks: a detected cell is kept only if it is the maximum of its
%      3 x 3 x 3 neighbourhood (beam x range x Doppler, Doppler circular), so
%      the cells on the skirts of a target response are dropped.
%   2) Sidelobes: a peak that is more than sys.cfg.prune.sidelobeMargin_dB
%      weaker than another peak on the same line is dropped:
%        range sidelobe     - same beam and Doppler, within one pulse length
%        elevation sidelobe - same range and Doppler, any beam
%        Doppler sidelobe   - same range and beam, any Doppler bin
%   (same = within +-1 cell). dets holds column vectors beam, rIdx, dIdx,
%   power and noise.

[nB, nR, nD] = size(P);

%% 1) Local peaks (compared with the 26 neighbours, Doppler circular)
idx = find(isDet);
[b, r, d] = ind2sub([nB nR nD], idx);
v = P(idx);
keep = true(size(idx));
for db = -1:1
    for dr = -1:1
        for dd = -1:1
            bb = b + db;  rr = r + dr;
            in = bb >= 1 & bb <= nB & rr >= 1 & rr <= nR;
            nb = bb(in) + nB*(rr(in) - 1) + nB*nR*mod(d(in) - 1 + dd, nD);   % linear index of the neighbour
            keep(in) = keep(in) & v(in) >= P(nb);
        end
    end
end
idx = idx(keep);  b = b(keep);  r = r(keep);  d = d(keep);
pw  = double(P(idx));

%% 2) Sidelobe blanking
% Only the peaks more than the margin stronger can blank a peak: with the
% powers sorted, they are the first k(i) (binary search for all peaks at once).
ratio = 10^(sys.cfg.prune.sidelobeMargin_dB/10);
pulseCells = numel(wv.mfTaps);                       % extent of the range sidelobes
n = numel(pw);
[pws, ord] = sort(pw, 'descend');
x  = ratio * pw;
lo = zeros(n, 1);  hi = n * ones(n, 1);              % pws(1:lo) > x, pws(hi+1:end) <= x
act = lo < hi;
while any(act)
    mid = floor((lo(act) + hi(act) + 1) / 2);
    gt  = pws(mid) > x(act);
    a   = find(act);
    lo(a(gt))  = mid(gt);
    hi(a(~gt)) = mid(~gt) - 1;
    act = lo < hi;
end
keep = true(n, 1);
for i = find(lo > 0).'
    j = ord(1:lo(i));
    dDop  = abs(d(i) - d(j));
    sameD = min(dDop, nD - dDop) <= 1;
    sameB = abs(b(i) - b(j)) <= 1;
    dR    = abs(r(i) - r(j));
    sameR = dR <= 1;
    keep(i) = ~any((sameB & sameD & dR <= pulseCells) | (sameR & sameD) | (sameR & sameB));
end
idx = idx(keep);
dets.beam  = b(keep);
dets.rIdx  = r(keep);
dets.dIdx  = d(keep);
dets.power = P(idx);
dets.noise = noise(idx);
end

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

%% 1) Local peaks
idx = find(isDet);
[b, r, d] = ind2sub([nB nR nD], idx);
keep = false(size(idx));
for i = 1:numel(idx)
    bb = max(b(i)-1, 1):min(b(i)+1, nB);
    rr = max(r(i)-1, 1):min(r(i)+1, nR);
    dd = mod(d(i)-2:d(i), nD) + 1;                    % circular Doppler neighbours
    block = P(bb, rr, dd);
    keep(i) = P(idx(i)) >= max(block(:));
end
idx = idx(keep);  b = b(keep);  r = r(keep);  d = d(keep);
pw  = double(P(idx));

%% 2) Sidelobe blanking
ratio = 10^(sys.cfg.prune.sidelobeMargin_dB/10);
pulseCells = numel(wv.mfTaps);                       % extent of the range sidelobes
keep = true(size(idx));
for i = 1:numel(idx)
    for j = find(pw > ratio*pw(i)).'
        dDop  = abs(d(i) - d(j));
        sameD = min(dDop, nD - dDop) <= 1;
        sameB = abs(b(i) - b(j)) <= 1;
        sameR = abs(r(i) - r(j)) <= 1;
        if (sameB && sameD && abs(r(i) - r(j)) <= pulseCells) || (sameR && sameD) || (sameR && sameB)
            keep(i) = false;
            break
        end
    end
end
idx = idx(keep);
dets.beam  = b(keep);
dets.rIdx  = r(keep);
dets.dIdx  = d(keep);
dets.power = P(idx);
dets.noise = noise(idx);
end

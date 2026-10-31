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
%   With many detections (a low CFAR threshold) both steps use maxima over
%   the whole map (max filters) instead of comparing detection by detection:
%   the same result, in a time that depends on the map size and hardly on
%   the number of detections.

[nB, nR, nD] = size(P);

%% 1) Local peaks (compared with the 26 neighbours, Doppler circular; P >= 0)
nDet = nnz(isDet);
if nDet > 3000                                       % many detections: 3x3x3 maximum of the whole map
    idx = find(isDet & P >= range_max(max_doppler3(max_beam3(P)), 1));
    [b, r, d] = ind2sub([nB nR nD], idx);
else
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
end
pw  = double(P(idx));

%% 2) Sidelobe blanking
% A peak is blanked when one of its three neighbourhoods holds a peak more
% than the margin stronger, i.e. when the largest peak power there is:
%   range sidelobe     - beam +-1, Doppler +-1, range within one pulse length
%   elevation sidelobe - range +-1, Doppler +-1, any beam
%   Doppler sidelobe   - range +-1, beam +-1, any Doppler bin
% With many peaks (a low CFAR threshold) the maxima are taken over the whole
% map (max filters), so the time depends on the map size only, not on the
% number of detections; with few peaks each one is compared with the peaks
% more than the margin stronger (same result).
ratio = 10^(sys.cfg.prune.sidelobeMargin_dB/10);
pulseCells = numel(wv.mfTaps);                       % extent of the range sidelobes
n = numel(pw);
if n > 300
    M = zeros(nB, nR, nD, 'single');
    M(idx) = P(idx);                                 % the peaks, 0 elsewhere
    MbD = max_doppler3(max_beam3(M));                % beam +-1, Doppler +-1
    A = range_max(MbD, pulseCells);
    blank = double(A(idx)) > ratio * pw;                                  % range sidelobe
    E = range_max(max_doppler3(max(M, [], 1)), 1);   % any beam, range +-1, Doppler +-1
    blank = blank | double(E(r + nR*(d - 1))) > ratio * pw;               % elevation sidelobe
    D = range_max(max_beam3(max(M, [], 3)), 1);      % any Doppler, range +-1, beam +-1
    blank = blank | double(D(b + nB*(r - 1))) > ratio * pw;               % Doppler sidelobe
    keep = ~blank;
else
    [pws, ord] = sort(pw, 'descend');
    x  = ratio * pw;
    lo = zeros(n, 1);  hi = n * ones(n, 1);          % pws(1:lo) > x, pws(hi+1:end) <= x
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
end
idx = idx(keep);
dets.beam  = b(keep);
dets.rIdx  = r(keep);
dets.dIdx  = d(keep);
dets.power = P(idx);
dets.noise = noise(idx);
end

function A = max_beam3(A)
% Maximum over beam -1 ... +1 (dimension 1, no wrap).
if size(A, 1) > 1
    A = max(A, max(A([2:end end], :, :), A([1 1:end-1], :, :)));
end
end

function A = max_doppler3(A)
% Maximum over Doppler -1 ... +1 (dimension 3, circular).
if size(A, 3) > 1
    A = max(A, max(circshift(A, 1, 3), circshift(A, -1, 3)));
end
end

function R = range_max(A, L)
% Maximum over range -L ... +L (dimension 2, zero outside; A >= 0). Long
% windows: van Herk / Gil-Werman (blocks of the window length, prefix and
% suffix maxima of every block), so the time does not depend on L.
[n1, nR, n3] = size(A);
if L <= 2
    R = A;
    for s = 1:L
        R = max(R, max(cat(2, A(:, 1+s:end, :), zeros(n1, s, n3, class(A))), ...
                       cat(2, zeros(n1, s, n3, class(A)), A(:, 1:end-s, :))));
    end
    return
end
w  = 2*L + 1;
nb = ceil((nR + 2*L) / w);                           % blocks of length w
B  = zeros(n1, nb*w, n3, class(A));
B(:, L+1:L+nR, :) = A;
B  = reshape(B, n1, w, nb, n3);
G  = reshape(cummax(B, 2), n1, nb*w, n3);                     % max from the block start
H  = reshape(flip(cummax(flip(B, 2), 2), 2), n1, nb*w, n3);   % max to the block end
j  = 1:nR;                                           % window j ... j+2L (padded)
R  = max(H(:, j, :), G(:, j + w - 1, :));
end

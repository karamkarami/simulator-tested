function [isDet, noise] = cfar_goso(P, sys, wv)
%CFAR_GOSO  Greatest-of (GO) and smallest-of (SO) CFAR (cfg.cfar.type = 'go' / 'so'), called by cfar_detect.m.
%   [isDet, noise] = CFAR_GOSO(P, sys, wv) splits the training cells of the
%   window of every tested cell (sys.cfar, as cfar_ca.m: guard cells removed,
%   the range window shifted inwards at the edges, Doppler circular) into two
%   halves and takes the mean of each:
%     2-D and 1-D range window: the cells at a smaller range than the cell
%       under test (leading) and at a larger one (lagging); the training cells
%       on its own range (2-D: the Doppler bins outside the guard) are in
%       neither half
%     1-D Doppler window: the bins below and above the cell under test
%   GO takes the larger mean, SO the smaller (at a range edge one half may be
%   empty: the other one). A cell is detected when P > sys.cfar.alpha * that
%   mean; noise = mean / sys.cfar.noiseScale. Sums in double.
c  = sys.cfar;
gR = c.guardRange;    oR = gR + c.trainRange;
gD = c.guardDoppler;  oD = gD + c.trainDoppler;
[nB, nR, nD] = size(P);
isDet = false(nB, nR, nD);
noise = zeros(nB, nR, nD, class(P));
t = wv.testIdx;
if t(2) < t(1), return; end
r = t(1):t(2);
if strcmp(c.method, 'doppler')
    Q = double(P(:, r, :));
    Q = cat(3, Q(:, :, end-oD+1:end), Q, Q(:, :, 1:oD));        % circular padding in Doppler
    C = cat(3, zeros(nB, numel(r)), cumsum(Q, 3));               % C(:,:,i+1) = sum of the padded bins 1 ... i
    d = 1:nD;                                                    % padded index of bin d: d + oD
    sA = C(:, :, d + oD - gD) - C(:, :, d);                      % bins d-oD ... d-gD-1
    sB = C(:, :, d + 2*oD + 1) - C(:, :, d + oD + gD + 1);       % bins d+gD+1 ... d+oD
    nA = oD - gD;  nBh = nA;
    mA = sA / nA;  mB = sB / nBh;
else
    lo = max(1, min(t(1), nR - oR) - oR);                        % range cells the windows need
    hi = min(nR, max(t(2), oR + 1) + oR);
    Q = double(P(:, lo:hi, :));
    DO = box_doppler(Q, oD);                                     % Doppler sums over +-oD, +-gD
    DG = box_doppler(Q, gD);
    CO = cat(2, zeros(nB, 1, nD), cumsum(DO, 2));                % CO(:,i+1,:) = sum of the cells lo ... lo+i-1
    CG = cat(2, zeros(nB, 1, nD), cumsum(DG, 2));
    rc = min(max(r, oR + 1), nR - oR);                           % centre of the (shifted) outer window
    a1 = rc - oR;  b2 = rc + oR;                                 % outer window rows
    g1 = max(r - gR, 1);  g2 = min(r + gR, nR);                  % guard rows inside the data
    I = @(x) x - lo + 1;                                         % row index into CO, CG (+1: the zero row)
    sA = (CO(:, I(r), :) - CO(:, I(a1), :)) - (CG(:, I(r), :) - CG(:, I(g1), :));          % rows a1 ... r-1
    sB = (CO(:, I(b2) + 1, :) - CO(:, I(r) + 1, :)) - (CG(:, I(g2) + 1, :) - CG(:, I(r) + 1, :));   % rows r+1 ... b2
    nA = (r - a1) * (2*oD + 1) - (r - g1) * (2*gD + 1);         % training cells of each half (per range cell)
    nBh = (b2 - r) * (2*oD + 1) - (g2 - r) * (2*gD + 1);
    mA = sA ./ reshape(max(nA, 1), 1, [], 1);
    mB = sB ./ reshape(max(nBh, 1), 1, [], 1);
    mA(:, nA == 0, :) = mB(:, nA == 0, :);                       % an empty half: the other one
    mB(:, nBh == 0, :) = mA(:, nBh == 0, :);
end
if strcmp(c.type, 'go'), z = max(mA, mB); else, z = min(mA, mB); end
isDet(:, r, :) = double(P(:, r, :)) > c.alpha * z;
noise(:, r, :) = z / c.noiseScale;
end

function S = box_doppler(Q, h)
% Sum over the Doppler bins d-h ... d+h (circular) of every cell.
S = Q;
for k = 1:h
    S = S + circshift(Q, k, 3) + circshift(Q, -k, 3);
end
end

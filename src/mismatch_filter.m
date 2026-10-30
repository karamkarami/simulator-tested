function h = mismatch_filter(Y, pad, mainHalf, nIter)
%MISMATCH_FILTER  Mismatch filter with the smallest peak range sidelobe.
%   h = MISMATCH_FILTER(Y, pad) designs a filter of size(Y,1) + 2*pad
%   coefficients for the sampled pulse Y (a phase code with one sample per
%   chip, or the pulse after the DDC from ddc_pulse.m, one column per
%   sampling phase). The output conv(Y(:,k), h) is 1 at the main peak of
%   the first column (the peak of its matched filter, in the middle of the
%   full convolution), and the largest sidelobe over all columns is
%   minimised: iteratively reweighted least squares (Lawson), the weight of
%   every sidelobe grows with its level.
%   h = MISMATCH_FILTER(Y, pad, mainHalf, nIter): samples within mainHalf of
%   a main peak count as main lobe (default 0: one sample per chip; 3 for a
%   pulse at 2 samples per chip after the DDC), nIter iterations (300).
%   Example (a filter at the DDC rate for the long pulse, see README.md):
%     h = mismatch_filter(ddc_pulse(sys, sys.wave.long, 0:0.25:0.75), 10, 3);
if nargin < 3, mainHalf = 0; end
if nargin < 4, nIter = 300; end
[N, nCol] = size(Y);
K = N + 2*pad;
As = [];
for col = 1:nCol
    A = zeros(N+K-1, K);                             % conv(Y(:,col), h) = A*h
    for j = 1:K
        A(j:j+N-1, j) = Y(:, col);
    end
    [~, c] = max(abs(A * conj(flipud([zeros(pad, 1); Y(:, col); zeros(pad, 1)]))));
    side = true(N+K-1, 1);
    side(max(1, c-mainHalf):min(end, c+mainHalf)) = false;
    As = [As; A(side, :)]; %#ok<AGROW>
    if col == 1, a = A(c, :).'; end
end
w = ones(size(As, 1), 1) / size(As, 1);
for it = 1:nIter
    h = (As' * (w .* As) + 1e-9*eye(K)) \ conj(a);
    h = h / (a.' * h);                               % main peak exactly 1
    w = w .* abs(As * h);
    w = w / sum(w);
end
end

function y = mti_filter(x, coeffs, staggerIdx, nFiltered, bypass)
%MTI_FILTER  Clutter cancellation along slow time.
%   y = MTI_FILTER(x, coeffs, staggerIdx) filters x [nBeams x nRange x nPulses]
%   (oldest pulse first) along the pulse dimension. coeffs [nStagger x L] holds
%   one set of L weights per stagger position; output m uses pulses m ... m+L-1
%   and the weights of the stagger position of pulse m. With a fixed PRF there
%   is a single row [1 -2 1] (three-pulse canceller).
%   y = MTI_FILTER(x, coeffs, staggerIdx, nFiltered) filters only the range
%   cells 1 ... nFiltered (cfg.mti.maxRange_m, wv.mtiCells); the other cells
%   pass unfiltered: output m is pulse m + L - 1 (the weights sys.mtiBypass,
%   or bypass when given), so every cell has the same nPulses - L + 1 outputs.
L = size(coeffs, 2);
nOut = size(x, 3) - L + 1;
nR = size(x, 2);
if nargin < 4 || isempty(nFiltered), nFiltered = nR; end
nFiltered = min(max(nFiltered, 0), nR);
if nFiltered < nR
    if nargin < 5 || isempty(bypass)
        bypass = zeros(size(coeffs));  bypass(:, end) = 1;
    end
    y = zeros(size(x, 1), nR, nOut, 'like', x);
    if nFiltered > 0
        y(:, 1:nFiltered, :) = mti_filter(x(:, 1:nFiltered, :), coeffs, staggerIdx);
    end
    y(:, nFiltered+1:end, :) = mti_filter(x(:, nFiltered+1:end, :), bypass, staggerIdx);
    return
end
y = 0;
for l = 1:L
    w = coeffs(staggerIdx(1:nOut), l);                      % weight l of every output
    if all(w == w(1))
        if w(1) == 0, continue; end                         % (bypass: only the last tap)
        w = w(1);                                           % fixed PRF: one weight (faster)
    else
        w = reshape(w, 1, 1, nOut);
    end
    y = y + x(:, :, l:l+nOut-1) .* w;
end
if isequal(y, 0)
    y = zeros(size(x, 1), size(x, 2), nOut, 'like', x);
end
end

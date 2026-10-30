function y = mti_filter(x, coeffs, staggerIdx)
%MTI_FILTER  Clutter cancellation along slow time.
%   y = MTI_FILTER(x, coeffs, staggerIdx) filters x [nBeams x nRange x nPulses]
%   (oldest pulse first) along the pulse dimension. coeffs [nStagger x L] holds
%   one set of L weights per stagger position; output m uses pulses m ... m+L-1
%   and the weights of the stagger position of pulse m. With a fixed PRF there
%   is a single row [1 -2 1] (three-pulse canceller).
L = size(coeffs, 2);
nOut = size(x, 3) - L + 1;
y = 0;
for l = 1:L
    w = coeffs(staggerIdx(1:nOut), l);                      % weight l of every output
    if all(w == w(1))
        w = w(1);                                           % fixed PRF: one weight (faster)
    else
        w = reshape(w, 1, 1, nOut);
    end
    y = y + x(:, :, l:l+nOut-1) .* w;
end
end

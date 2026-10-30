function y = fir_decimate_valid(x, h, D)
%FIR_DECIMATE_VALID  FIR filter (valid time-domain convolution) followed by decimation.
%   y = FIR_DECIMATE_VALID(x, h, D) slides the filter h over every column of
%   x and keeps every D-th output of the valid part of the convolution:
%       y(k) = sum_j h(j) * x((k-1)*D + L + 1 - j),   j = 1 ... L = numel(h)
%   which is the same as conv(x, h, 'valid') followed by y(1:D:end), but only
%   the kept outputs are computed (no FFT).
%   Polyphase form: the input is cut into blocks of D samples and the
%   time-reversed filter into Q = ceil(L/D) blocks of D taps. One matrix
%   product applies every filter block to every input block; output k is the
%   sum of filter block q applied to input block k+q-1 (q = 1 ... Q). With
%   integer samples and taps (the fixed-point DDC) every product and sum is
%   exact in double, so the result does not depend on the order of the sums.

[N, nCol] = size(x);
L = numel(h);
K = floor((N - L)/D) + 1;                           % number of output samples
if K < 1
    y = zeros(0, nCol);
    return
end
Q = ceil(L/D);                                      % filter blocks of D taps
g = zeros(Q*D, 1);
g(1:L) = flipud(h(:));                              % y(k) = sum_i g(i) * x((k-1)*D + i)
M = K + Q - 1;                                      % input blocks that are used
if M*D <= N
    xb = x(1:M*D, :);
else
    xb = [x; zeros(M*D - N, nCol)];                 % (these samples meet the zero taps of g)
end
% A(m, col, q): filter block q applied to input block m of column col
A = reshape(reshape(xb, D, M*nCol).' * reshape(g, D, Q), M, nCol, Q);
y = A(1:K, :, 1);
for q = 2:Q
    y = y + A(q:q+K-1, :, q);
end
end

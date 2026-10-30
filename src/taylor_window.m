function w = taylor_window(N, nbar, sll_dB)
%TAYLOR_WINDOW  Taylor window (same definition as MATLAB's taylorwin).
%   w = TAYLOR_WINDOW(N, nbar, sll_dB) returns an N x 1 window with nbar
%   nearly constant-level sidelobes adjacent to the main lobe and a peak
%   sidelobe level of sll_dB (negative, e.g. -35).

A   = acosh(10^(-sll_dB/20)) / pi;
sp2 = nbar^2 / (A^2 + (nbar - 0.5)^2);
xi  = ((0:N-1)' - 0.5*N + 0.5) / N;
w   = ones(N, 1);
n   = (1:nbar-1)';
for m = 1:nbar-1
    p = n(n ~= m);
    num = prod(1 - (m^2/sp2) ./ (A^2 + (n - 0.5).^2));
    den = prod(1 - m^2 ./ p.^2);
    Fm  = (-1)^(m+1) * num / (2*den);
    w   = w + 2*Fm*cos(2*pi*m*xi);
end
end

function w = hamming_window(N)
%HAMMING_WINDOW  N x 1 Hamming window (no toolbox needed).
if N == 1
    w = 1;
else
    w = 0.54 - 0.46*cos(2*pi*(0:N-1)'/(N-1));
end
end

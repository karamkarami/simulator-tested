function [P, Y] = doppler_process(y, sys)
%DOPPLER_PROCESS  Doppler filter bank with non-coherent integration of rows.
%   P = DOPPLER_PROCESS(y, sys) turns the MTI output y [nBeams x nRange x nOut]
%   (oldest first) into range-Doppler power maps P [nBeams x nRange x nfft]
%   (zero Doppler in the middle, axis sys.dop.f_Hz).
%
%   The slow-time samples of each beam and range cell are arranged as rows:
%     fixed PRF     : consecutive blocks of nfft pulses
%     staggered PRF : de-interleaved by stagger position, so that every row is
%                     uniformly sampled with a period equal to the sum of the PRIs
%   Every row is weighted with a Hamming window and transformed with an
%   nfft-point FFT; |.|^2 is averaged over the rows (non-coherent integration).
%   Finally every Doppler bin is divided by its noise gain (sys.dop.noiseGain)
%   so that white noise at the MTI input gives a flat noise floor.
%   Y [nBeams x nRange x nfft x nRows] holds the complex spectra of the rows.

nRows  = sys.cpi.nRows;
rowLen = sys.cpi.rowLen;
[nB, nR, ~] = size(y);
y = y(:, :, end-nRows*rowLen+1:end);                       % most recent samples
if sys.staggered
    y = reshape(y, nB, nR, nRows, rowLen);                 % (:,:,k,j) = sample (j-1)*nRows + k
    y = permute(y, [1 2 4 3]);                             % rows along dimension 4
else
    y = reshape(y, nB, nR, rowLen, nRows);                 % (:,:,j,k) = sample (k-1)*rowLen + j
end
w = reshape(sys.dop.window, 1, 1, rowLen);
Y = fftshift(fft(y .* w, sys.dop.nfft, 3), 3);
P = mean(abs(Y).^2, 4) ./ sys.dop.noiseGain;
end

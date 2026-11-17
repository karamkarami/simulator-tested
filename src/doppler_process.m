function [P, Y] = doppler_process(y, sys, nFiltered, staggerIdx)
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
%   P = DOPPLER_PROCESS(y, sys, nFiltered): only the range cells 1 ... nFiltered
%   passed the clutter filter (cfg.mti.maxRange_m, wv.mtiCells); the others are
%   divided by the noise gain of the unfiltered pulses (sys.dop.noiseGainBypass).
%   P = DOPPLER_PROCESS(y, sys, nFiltered, staggerIdx): staggerIdx holds the stagger position of
%   every MTI output (that of its oldest pulse; cpi_buffer_read: info.staggerIdx). Staggered PRF with
%   cfg.proc.staggerDoppler = 'nufft' (sys.dop.method): the whole CPI is one non-uniformly sampled
%   sequence; its non-uniform DFT on the velocity grid (setup_radar.m: nufft_setup) uses the matrix of
%   the stagger position of the oldest output (position 1 without staggerIdx), P = |Y|^2 (the matrix
%   already equalises the noise; the cells beyond nFiltered are divided by sys.dop.nufft.bypassGain).

if isfield(sys.dop, 'method') && strcmp(sys.dop.method, 'nufft')
    rowLen = sys.cpi.rowLen;
    [nB, nR, nY] = size(y);
    p = 1;
    if nargin >= 4 && ~isempty(staggerIdx), p = staggerIdx(nY - rowLen + 1); end
    A = sys.dop.nufft.A(:, :, p);
    Y = reshape(y(:, :, end-rowLen+1:end), nB*nR, rowLen) * A.';   % [cells x nPoints]
    P = reshape(abs(Y).^2, nB, nR, []);
    Y = reshape(Y, nB, nR, []);
    if nargin >= 3 && ~isempty(nFiltered) && nFiltered < nR
        nF = max(nFiltered, 0);
        P(:, nF+1:end, :) = P(:, nF+1:end, :) ./ reshape(sys.dop.nufft.bypassGain(:, p), 1, 1, []);
    end
    return
end
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
P = mean(abs(Y).^2, 4);
if nargin >= 3 && ~isempty(nFiltered) && nFiltered < nR
    nF = max(nFiltered, 0);
    P(:, 1:nF, :) = P(:, 1:nF, :) ./ reshape(sys.dop.noiseGain, 1, 1, []);
    P(:, nF+1:end, :) = P(:, nF+1:end, :) ./ reshape(sys.dop.noiseGainBypass, 1, 1, []);
else
    P = P ./ sys.dop.noiseGain;
end
end

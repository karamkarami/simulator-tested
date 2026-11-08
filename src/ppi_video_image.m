function [rgb, alpha] = ppi_video_image(V, tMin_s)
%PPI_VIDEO_IMAGE  The PPI video map as an RGB image (ppi_video_init.m).
%   [rgb, alpha] = PPI_VIDEO_IMAGE(V, tMin_s) returns the nPix x nPix x 3
%   image (row 1 = South) and its transparency (0 outside the video range and
%   where there is no data newer than tMin_s). Colour: the power in dB from
%   3 dB below the noise floor (black) to V.span_dB above it (white).
lo = V.noise_dB - 3;
if isnan(lo), lo = 0; end
M = V.map;                                                    % largest bin under every pixel
K = V.binsPerPix;
if K > 1
    M(:, end+1:V.nRPix*K) = 0;
    M = reshape(max(reshape(M, V.nAz, K, V.nRPix), [], 2), V.nAz, V.nRPix);
end
dB = 10*log10(double(M(V.pixLin)) + eps);
k = round((dB - lo) / (V.span_dB + 3) * 255) + 1;
k = min(max(k, 1), 256);
valid = M(V.pixLin) > 0 & V.time(V.pixAz) > tMin_s;
n = V.nPix^2;
rgb = zeros(n, 3);
rgb(V.pixInside, :) = V.cmap(k, :);
alpha = zeros(V.nPix);
alpha(V.pixInside) = valid;
rgb = reshape(rgb, V.nPix, V.nPix, 3);
end

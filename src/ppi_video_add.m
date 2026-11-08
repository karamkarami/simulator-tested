function V = ppi_video_add(V, I, Q, az_deg, t_s)
%PPI_VIDEO_ADD  Add the video of one PRI to the PPI video map (ppi_video_init.m).
%   V = PPI_VIDEO_ADD(V, I, Q, az_deg, t_s) averages the DDC output I, Q
%   [nDdc x nElements] non-coherently over the elements, cuts out and joins
%   the receive intervals of both sub-pulses and writes the result into the
%   azimuth bin of az_deg: the mean of the PRIs while the antenna is in the
%   bin (bins skipped by a jump of less than 10 deg get the same data).
I = double(I(V.sampleIdx, :));
Q = double(Q(V.sampleIdx, :));
p = mean(I.^2 + Q.^2, 2).';                                   % joined video [LSB^2]
row = accumarray(V.sampleBin(:), p(:), [V.nR 1]).' ./ max(V.binCount, 1);
b = min(floor(mod(az_deg, 360) / V.azRes_deg) + 1, V.nAz);
if b == V.azBin
    V.sum = V.sum + row;  V.n = V.n + 1;
else
    if V.azBin > 0                                            % fill a small gap of the antenna
        gap = mod(b - V.azBin, V.nAz) - 1;
        if gap > 0 && gap * V.azRes_deg < 10
            k = mod(V.azBin + (1:gap) - 1, V.nAz) + 1;
            V.map(k, :) = repmat(single(row), gap, 1);
            V.time(k) = t_s;
        end
    end
    V.azBin = b;  V.sum = row;  V.n = 1;
end
V.map(b, :) = single(V.sum / V.n);
V.time(b) = t_s;
m = median(10*log10(row(V.binCount > 0) + eps));              % noise floor: running median
if isnan(V.noise_dB), V.noise_dB = m; else, V.noise_dB = 0.98*V.noise_dB + 0.02*m; end
end

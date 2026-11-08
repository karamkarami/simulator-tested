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
V = ppi_video_add_row(V, row, az_deg, t_s);
end

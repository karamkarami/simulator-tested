function row = ppi_profile_row(V, kind, x)
%PPI_PROFILE_ROW  Power of the selected beam after DBF or pulse compression, joined in range bins.
%   row = PPI_PROFILE_ROW(V, kind, x) with x = {short, long}: the complex
%   samples of the beam for both sub-pulses, kind 'dbf' (DBF output: the
%   pulse compression input of ddc_segment.m, [1 x numel(wv.ddcIdx)]) or 'mf'
%   (pulse compression output, [1 x nRange]). Returns |x|^2 in the range bins
%   of the PPI map (V, ppi_video_init.m): the short pulse up to the first
%   long-pulse range cell, then the long pulse (the mean of the samples of a
%   bin; 0 where there is none). The C receiver does the same (processor.c,
%   profile_add).
L = V.(kind);
num = zeros(1, V.nR);  cnt = zeros(1, V.nR);
for p = 1:2
    v = double(x{p}(L{p}.idx));
    num = num + accumarray(L{p}.bin(:), abs(v(:)).^2, [V.nR 1]).';
    cnt = cnt + accumarray(L{p}.bin(:), 1, [V.nR 1]).';
end
row = num ./ max(cnt, 1);
end

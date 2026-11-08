function V = ppi_video_reset(V, kind, beam)
%PPI_VIDEO_RESET  Empty the raw data map of the PPI for other data (ppi_video_init.m).
%   V = PPI_VIDEO_RESET(V, kind, beam): kind 'video', 'dbf', 'mf' (or '' =
%   none: the map stays empty) and the beam of 'dbf' / 'mf'. The map, the
%   azimuth bins and the noise floor start again.
if isempty(kind), kind = 'none'; end
V.kind = kind;
V.beam = beam;
V.map(:) = 0;
V.time(:) = -Inf;
V.sum = zeros(1, V.nR);  V.n = 0;  V.azBin = 0;
V.noise_dB = NaN;
end

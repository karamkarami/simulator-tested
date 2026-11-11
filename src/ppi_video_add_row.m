function V = ppi_video_add_row(V, row, az_deg, t_s)
%PPI_VIDEO_ADD_ROW  Add the raw data of one PRI (range bins of V.dR_m) to the PPI map.
%   V = PPI_VIDEO_ADD_ROW(V, row, az_deg, t_s) writes row [1 x V.nR] (power,
%   0 = no sample in the bin) into the azimuth bin of az_deg: the mean of the
%   PRIs while the antenna is in the bin (bins skipped by a jump of less than
%   10 deg get the same data). Used by ppi_video_add.m (video) and for the DBF
%   and MF outputs of the selected beam (ppi_profile_row.m).
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
v = sort(row(row > 0));                                        % noise floor: running 10th percentile (the
if isempty(v), v = 0; end
m = 10*log10(v(max(1, ceil(0.1*numel(v)))) + eps);             % raw video is mostly clutter: long-pulse smear)
if isfield(V, 'kind') && any(strcmp(V.kind, {'dbf', 'mf'}))
    m = m + V.p10Sample_dB;                                    % one complex sample per bin: exponential power,
else                                                           % its 10th percentile is 9.77 dB below the mean;
    m = m + V.p10Video_dB;                                     % video: the mean over the elements (48: 0.86 dB)
end
% The 10th percentile of a row is the noise only while 10% of its bins hold no echo: over the
% mountains the clutter covers more than 90% of a row, and a running mean of the rows raised the
% colour scale while the antenna crossed them (the noise floor changed colour for a moment).
% The noise does not depend on the azimuth: the 20th percentile of the newest estimate of every
% azimuth bin (a clutter sector smaller than 80% of the circle cannot move it); the running
% mean until 8 azimuth bins have data.
if ~isfield(V, 'azNoise') || numel(V.azNoise) ~= V.nAz, V.azNoise = nan(V.nAz, 1); end
V.azNoise(b) = m;
e = sort(V.azNoise(isfinite(V.azNoise)));
if numel(e) >= 8
    V.noise_dB = e(floor((numel(e) - 1) / 5) + 1);
elseif isnan(V.noise_dB)
    V.noise_dB = m;
else
    V.noise_dB = 0.98*V.noise_dB + 0.02*min(m, V.noise_dB + 0.05);
end
end

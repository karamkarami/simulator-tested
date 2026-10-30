function [data, info] = cpi_buffer_read(buf)
%CPI_BUFFER_READ  Buffer content in chronological order (oldest pulse first).
%   [data, info] = CPI_BUFFER_READ(buf) returns data [nBeams x nRange x nPulses]
%   and info with the stagger index, azimuth and time of every pulse, plus the
%   azimuth and time of the CPI centre.
n = size(buf.data, 3);
order = [buf.last+1:n, 1:buf.last];
data = buf.data(:, :, order);
info.staggerIdx = buf.staggerIdx(order);
info.az_deg     = buf.az_deg(order);
info.time_s     = buf.time_s(order);
info.azCenter_deg = mod(atan2d(mean(sind(info.az_deg)), mean(cosd(info.az_deg))), 360);
info.tCenter_s    = mean(info.time_s);
end

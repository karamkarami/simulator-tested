function buf = cpi_buffer_push(buf, rp, staggerIdx, az_deg, time_s)
%CPI_BUFFER_PUSH  Add a new range profile; the oldest pulse is overwritten when full.
n = size(buf.data, 3);
buf.last = mod(buf.last, n) + 1;
buf.data(:, :, buf.last) = rp;
buf.staggerIdx(buf.last) = staggerIdx;
buf.az_deg(buf.last)     = az_deg;
buf.time_s(buf.last)     = time_s;
buf.count   = min(buf.count + 1, n);
buf.nPushed = buf.nPushed + 1;
end

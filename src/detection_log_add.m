function log = detection_log_add(log, reports)
%DETECTION_LOG_ADD  Add the reports of a CPI update (detection_reports.m) to the log.
%   When the buffer is full the oldest detections are overwritten.
n = numel(reports);
if n == 0, return; end
cap = size(log.M, 1);
rows = [[reports.time_s]; [reports.az_deg]; [reports.range_m]; [reports.el_deg]; [reports.vr_mps]; ...
        [reports.snr_dB]; 1 + strcmp({reports.pulse}, 'long'); [reports.beam]; [reports.rIdx]; [reports.dIdx]].';
if n > cap, rows = rows(end-cap+1:end, :); n = cap; end
idx = mod(log.next - 1 + (0:n-1), cap) + 1;
log.M(idx, :) = rows;
log.next = mod(idx(end), cap) + 1;
log.n = min(cap, log.n + n);
end

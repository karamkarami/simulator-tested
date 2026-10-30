function log = detection_log_init(capacity)
%DETECTION_LOG_INIT  Ring buffer of the pruned CFAR detections shown on the PPI.
%   log = DETECTION_LOG_INIT(capacity) keeps the newest 'capacity' detections
%   (the reports that go into the plot extractor), one row each:
%     [time_s az_deg range_m el_deg vr_mps snr_dB pulse beam rIdx dIdx]
%   with pulse = 1 (short) or 2 (long). detection_log_add.m adds the reports
%   of a CPI update; log.M(1:log.n, :) are the stored rows (in no order).
log.M = zeros(capacity, 10);
log.n = 0;                 % rows in use
log.next = 1;              % row written next (oldest row once the buffer is full)
end

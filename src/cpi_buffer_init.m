function buf = cpi_buffer_init(nBeams, nRange, nPulses)
%CPI_BUFFER_INIT  Empty circular buffer holding the range profiles of one CPI.
buf.data       = complex(zeros(nBeams, nRange, nPulses, 'single'));  % [beam x range x pulse]
buf.staggerIdx = zeros(1, nPulses);    % stagger position of every stored pulse
buf.az_deg     = zeros(1, nPulses);    % antenna azimuth of every stored pulse
buf.time_s     = zeros(1, nPulses);    % transmit time of every stored pulse
buf.count      = 0;                    % number of valid pulses in the buffer
buf.last       = 0;                    % slot of the newest pulse
buf.nPushed    = 0;                    % pulses pushed since the last reset
end

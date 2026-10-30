function s = pulse_schedule(sys)
%PULSE_SCHEDULE  Transmit time, antenna azimuth, scan number and stagger position of every PRI.
%   s = PULSE_SCHEDULE(sys) covers cfg.sim.nScans antenna revolutions. s.active
%   marks the pulses inside cfg.sim.sectors_deg (all pulses when it is empty).
cfg   = sys.cfg;
tEnd  = cfg.sim.nScans * 360 / sys.ant.rate_dps;
n     = ceil(tEnd / min(sys.pri)) + 1;
s.staggerIdx = mod(0:n-1, sys.nStagger) + 1;
pri          = sys.pri(s.staggerIdx);
s.pri_s      = pri(:).';
s.time_s     = [0, cumsum(s.pri_s(1:end-1))];
keep = s.time_s < tEnd;
s.staggerIdx = s.staggerIdx(keep);
s.pri_s      = s.pri_s(keep);
s.time_s     = s.time_s(keep);
azUnwrapped  = cfg.ant.az0_deg + sys.ant.rate_dps * s.time_s;
s.az_deg     = mod(azUnwrapped, 360);
s.scan       = floor((azUnwrapped - cfg.ant.az0_deg) / 360) + 1;
s.nPulses    = numel(s.time_s);
s.active     = true(1, s.nPulses);
sectors = cfg.sim.sectors_deg;
if ~isempty(sectors)
    s.active(:) = false;
    for i = 1:size(sectors, 1)
        width = mod(sectors(i, 2) - sectors(i, 1), 360);
        if width == 0 && sectors(i, 2) ~= sectors(i, 1)
            width = 360;                                         % e.g. [0 360]: the whole circle
        end
        s.active = s.active | mod(s.az_deg - sectors(i, 1), 360) <= width;
    end
end
end

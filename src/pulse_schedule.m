function s = pulse_schedule(sys)
%PULSE_SCHEDULE  Transmit time, antenna azimuth, scan number and stagger position of every PRI.
%   s = PULSE_SCHEDULE(sys) covers cfg.sim.nScans antenna revolutions. s.active
%   marks the pulses inside cfg.sim.sectors_deg (all pulses when it is empty).
%   s.prfIdx is the PRF of every pulse with the scan-to-scan PRF (cfg.prf.mode =
%   'scan2scan'): scan s uses sys.priScan_s(mod(s - 1, numel) + 1); a pulse
%   belongs to the scan in which it starts (the PRF changes with the first pulse
%   after the antenna passes cfg.ant.az0_deg). Otherwise s.prfIdx is 1.
cfg   = sys.cfg;
tScan = 360 / sys.ant.rate_dps;
tEnd  = cfg.sim.nScans * tScan;
if isfield(sys, 'scan2scan') && sys.scan2scan
    nP = numel(sys.priScan_s);
    time = cell(1, cfg.sim.nScans);  pri = time;  idx = time;
    t0 = 0;
    for sc = 1:cfg.sim.nScans
        j = mod(sc - 1, nP) + 1;
        tk = t0:sys.priScan_s(j):(sc*tScan - 1e-12);        % pulses starting in this scan
        if isempty(tk), continue; end
        time{sc} = tk;  pri{sc} = repmat(sys.priScan_s(j), size(tk));  idx{sc} = repmat(j, size(tk));
        t0 = tk(end) + sys.priScan_s(j);                     % the last PRI of the scan keeps its PRF
    end
    s.time_s = [time{:}];
    s.pri_s = [pri{:}];
    s.prfIdx = [idx{:}];
    s.staggerIdx = ones(size(s.time_s));
else
    n     = ceil(tEnd / min(sys.pri)) + 1;
    s.staggerIdx = mod(0:n-1, sys.nStagger) + 1;
    pri          = sys.pri(s.staggerIdx);
    s.pri_s      = pri(:).';
    s.time_s     = [0, cumsum(s.pri_s(1:end-1))];
    keep = s.time_s < tEnd;
    s.staggerIdx = s.staggerIdx(keep);
    s.pri_s      = s.pri_s(keep);
    s.time_s     = s.time_s(keep);
    s.prfIdx     = ones(size(s.time_s));
end
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

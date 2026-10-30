function n = target_look(sys, azTarget_deg, t_s)
%TARGET_LOOK  Number of the antenna pass over a target (1, 2, ...).
%   n = TARGET_LOOK(sys, azTarget_deg, t_s) is the pass of the antenna over
%   a target at azimuth azTarget_deg at time t_s. The number changes when
%   the antenna points away from the target (180 deg from it), never while
%   the beam is on it. A Swerling 1 RCS drawn per pass therefore stays
%   constant during the whole dwell, also for a target near the start
%   azimuth of the scans.
rot = sys.ant.rate_dps * t_s;                            % antenna rotation since t = 0 [deg]
rel = mod(azTarget_deg - sys.cfg.ant.az0_deg, 360);      % rotation of the first pass [deg]
n = max(1, floor((rot - rel + 180) / 360) + 1);
end

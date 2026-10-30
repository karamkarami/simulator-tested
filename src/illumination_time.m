function t = illumination_time(sys, tgt)
%ILLUMINATION_TIME  Time [s] at which the antenna boresight crosses a target (first scan).
%   The target moves while the antenna turns, so the crossing time is found by
%   a few fixed-point iterations: t = (target azimuth at t - az0) / rotation rate.
t = 0;
for it = 1:4
    [~, az] = target_geometry(tgt, t, sys.ant.height_m);
    t = mod(az - sys.cfg.ant.az0_deg, 360) / sys.ant.rate_dps;
end
end

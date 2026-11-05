function targets = random_targets_scenario(r, seed)
%RANDOM_TARGETS_SCENARIO  r.n targets with random position, RCS and velocity.
%   targets = RANDOM_TARGETS_SCENARIO(r, seed) draws r.n targets (names R1,
%   R2, ...) with the random seed "seed" (the same seed gives the same
%   targets; the random state of the caller is restored afterwards). The
%   limits come from r (cfg.scenario.random of radar_config.m):
%     rcs_m2    [min max] mean RCS, log-uniform [m^2]
%     range_m   [min max] ground range at t = 0, uniform [m]
%     alt_m     [min max] altitude, uniform [m], at most range*tand(maxEl_deg)
%     maxEl_deg highest elevation of a target at t = 0 [deg]
%     speed_mps [min max] horizontal speed, uniform [m/s]; heading uniform 0 - 360 deg
%   The azimuth at t = 0 is uniform in 0 - 360 deg.
saved = rng;
rng(seed);
targets = make_target('R1', 1, 1e3, 0, 0, 0, 0);
for k = 1:r.n
    rcs   = exp(log(r.rcs_m2(1)) + rand * (log(r.rcs_m2(2)) - log(r.rcs_m2(1))));
    range = r.range_m(1) + rand * diff(r.range_m);
    az    = 360 * rand;
    altMax = min(r.alt_m(2), range * tand(r.maxEl_deg));
    altMin = min(r.alt_m(1), altMax);
    alt   = altMin + rand * (altMax - altMin);
    speed = r.speed_mps(1) + rand * diff(r.speed_mps);
    hdg   = 360 * rand;                              % heading, clockwise from North
    targets(k) = make_target(sprintf('R%d', k), rcs, range, az, alt, ...
                             speed * cosd(hdg - az), speed * sind(hdg - az));   % radial, cross (clockwise)
end
targets = targets(1:r.n);
rng(saved);
end

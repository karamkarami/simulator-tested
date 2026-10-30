function sectors = target_sectors(sys, targets)
%TARGET_SECTORS  Azimuth sectors around targets (to simulate only where targets are).
%   sectors = TARGET_SECTORS(sys, targets) returns one row [azStart azEnd] per
%   target for cfg.sim.sectors_deg. Each sector is centred on the azimuth where
%   the antenna crosses the target and extends +-(beamwidth + half a CPI +
%   1 deg), so the CPI buffer is full before the target enters the beam and
%   its plot is closed after it leaves.
cpiSpan = sys.cpi.length * mean(sys.pri) * sys.ant.rate_dps;   % azimuth covered by one CPI
half = sys.ant.azBeamwidth_deg + cpiSpan/2 + 1;
sectors = zeros(numel(targets), 2);
for i = 1:numel(targets)
    t = illumination_time(sys, targets(i));
    [~, az] = target_geometry(targets(i), t, sys.ant.height_m);
    sectors(i, :) = mod([az - half, az + half], 360);
end
end

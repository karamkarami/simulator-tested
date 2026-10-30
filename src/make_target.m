function tgt = make_target(name, rcs_m2, range_m, az_deg, alt_m, vRadial_mps, vCross_mps)
%MAKE_TARGET  Point target flying with a constant velocity.
%   tgt = MAKE_TARGET(name, rcs, range, az, alt, vRadial, vCross) places the
%   target at ground range "range" [m], azimuth az [deg, clockwise from North]
%   and altitude alt [m] at t = 0. The horizontal velocity has a component
%   vRadial along the line of sight (positive = moving away) and vCross
%   perpendicular to it (positive = clockwise).
%   Coordinates: x = East, y = North, z = up, origin on the ground below the radar.

tgt.name    = name;
tgt.rcs_m2  = rcs_m2;
tgt.pos0_m  = [range_m*sind(az_deg), range_m*cosd(az_deg), alt_m];
tgt.vel_mps = vRadial_mps*[sind(az_deg), cosd(az_deg), 0] ...
            + vCross_mps*[cosd(az_deg), -sind(az_deg), 0];
end

function tgt = make_target(name, rcs_m2, range_m, az_deg, alt_m, vRadial_mps, vCross_mps, maneuvers)
%MAKE_TARGET  Point target: constant velocity, or with maneuvers (turns, accelerations).
%   tgt = MAKE_TARGET(name, rcs, range, az, alt, vRadial, vCross) places the
%   target at ground range "range" [m], azimuth az [deg, clockwise from North]
%   and altitude alt [m] at t = 0. The horizontal velocity has a component
%   vRadial along the line of sight (positive = moving away) and vCross
%   perpendicular to it (positive = clockwise).
%   tgt = MAKE_TARGET(..., maneuvers) adds maneuvers (target_state.m): rows
%   [tStart_s turnRate_dps accelAlong_mps2 accelUp_mps2], each valid from
%   tStart_s until the tStart_s of the next row (the last one until the end):
%     turnRate_dps    horizontal coordinated turn (positive = to the right,
%                     clockwise seen from above): the heading changes at this
%                     rate, the horizontal speed only by accelAlong
%     accelAlong_mps2 acceleration along the horizontal velocity (positive =
%                     faster); without a turn the target flies with a constant
%                     acceleration on a straight line
%     accelUp_mps2    vertical acceleration (positive = up)
%   Before the first row and in a row [t 0 0 0] the velocity is constant.
%   Example: a 20 s turn at 4 deg/s to the right after 5 s, then straight on
%       make_target('M', 2, 30e3, 45, 3000, -200, 0, [5 4 0 0; 25 0 0 0])
%   Coordinates: x = East, y = North, z = up, origin on the ground below the radar.
%   pos0_m and vel_mps are the position and velocity at t = 0; target_state.m
%   gives them at any time.

tgt.name    = name;
tgt.rcs_m2  = rcs_m2;
tgt.pos0_m  = [range_m*sind(az_deg), range_m*cosd(az_deg), alt_m];
tgt.vel_mps = vRadial_mps*[sind(az_deg), cosd(az_deg), 0] ...
            + vCross_mps*[cosd(az_deg), -sind(az_deg), 0];
if nargin < 8 || isempty(maneuvers), maneuvers = zeros(0, 4); end
if size(maneuvers, 2) ~= 4 || any(diff(maneuvers(:, 1)) <= 0) || any(maneuvers(:, 1) < 0) || ~all(isfinite(maneuvers(:)))
    error('make_target: maneuvers must be rows [tStart_s turnRate_dps accelAlong_mps2 accelUp_mps2], tStart >= 0 and increasing.');
end
tgt.maneuvers = maneuvers;
end

function [R, az, el, vr] = target_geometry(tgt, t, antHeight_m)
%TARGET_GEOMETRY  Range, azimuth, elevation and radial velocity of a target at time t.
%   R [m] is measured from the antenna phase center, az [deg] clockwise from
%   North, el [deg] above the horizon (flat earth) and vr [m/s] is the radial
%   velocity (positive = moving away from the radar). The position and
%   velocity come from target_state.m (constant velocity, or the maneuvers of
%   the target: turns, accelerations).
[pos, vel] = target_state(tgt, t);
d  = pos - [0 0 antHeight_m];
R  = norm(d);
az = mod(atan2d(d(1), d(2)), 360);
el = asind(d(3)/R);
vr = dot(vel, d)/R;
end

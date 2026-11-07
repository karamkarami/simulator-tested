function [pos, vel, acc] = target_state(tgt, t)
%TARGET_STATE  Position, velocity and acceleration of a target at time t.
%   [pos, vel, acc] = TARGET_STATE(tgt, t) returns row vectors [x y z] (x =
%   East, y = North, z = up, origin on the ground below the radar) in m, m/s
%   and m/s^2. Without maneuvers (tgt.maneuvers empty or missing) the target
%   flies with the constant velocity tgt.vel_mps from tgt.pos0_m (t = 0).
%   Maneuvers (make_target.m): rows [tStart_s turnRate_dps accelAlong_mps2
%   accelUp_mps2], each valid until the next tStart_s. In a row the heading
%   turns at turnRate (positive = clockwise), the horizontal speed changes by
%   accelAlong and the vertical speed by accelUp; the motion is integrated in
%   closed form (exact for turns and constant accelerations). A horizontal
%   speed that would become negative stops at 0 (the target hovers).
pos = tgt.pos0_m;
vel = tgt.vel_mps;
M = zeros(0, 4);
if isfield(tgt, 'maneuvers'), M = tgt.maneuvers; end
tc = 0;                                       % time of pos, vel
cur = [0 0 0];                                % [turn rate deg/s, along, up] of the current segment
for k = 1:size(M, 1)
    if t <= M(k, 1), break; end
    [pos, vel] = segment(pos, vel, M(k, 1) - tc, cur);
    tc = M(k, 1);
    cur = M(k, 2:4);
end
[pos, vel, acc] = segment(pos, vel, t - tc, cur);
end

function [p, v, a] = segment(p, v, tau, m)
% Motion over tau with the turn rate m(1) [deg/s, clockwise], the along-track
% acceleration m(2) and the vertical acceleration m(3).
w = m(1) * pi/180;  aT = m(2);  aZ = m(3);
s0 = hypot(v(1), v(2));
th0 = atan2(v(1), v(2));                      % heading, clockwise from North
th = th0;
s = s0;
if s0 > 0 || aT > 0
    tauH = tau;                               % horizontal motion lasts tauH
    if aT < 0 && s0 + aT*tau < 0, tauH = -s0 / aT; end   % the speed reaches 0: it stays there
    th = th0 + w*tauH;
    if abs(w*tauH) < 1e-4                     % (nearly) straight: series in w
        Is  = tauH*sin(th0) + w*tauH^2/2*cos(th0);
        Ic  = tauH*cos(th0) - w*tauH^2/2*sin(th0);
        Ius = tauH^2/2*sin(th0) + w*tauH^3/3*cos(th0);
        Iuc = tauH^2/2*cos(th0) - w*tauH^3/3*sin(th0);
    else
        Is  = (cos(th0) - cos(th)) / w;
        Ic  = (sin(th) - sin(th0)) / w;
        Ius = -tauH*cos(th)/w + (sin(th) - sin(th0))/w^2;
        Iuc =  tauH*sin(th)/w + (cos(th) - cos(th0))/w^2;
    end
    s = s0 + aT*tauH;
    p(1:2) = p(1:2) + [s0*Is + aT*Ius, s0*Ic + aT*Iuc];
    if tauH < tau, aT = 0;  w = 0; end        % stopped: no more horizontal acceleration
end
p(3) = p(3) + v(3)*tau + aZ*tau^2/2;
v = [s*sin(th), s*cos(th), v(3) + aZ*tau];
% acceleration: along-track aT plus the centripetal s*w across the velocity
a = [aT*sin(th) + s*w*cos(th), aT*cos(th) - s*w*sin(th), aZ];
end

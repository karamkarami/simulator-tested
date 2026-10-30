function rep = track_report(tracker)
%TRACK_REPORT  The current tracks for the display (one struct per track).
%   Fields: id, confirmed, t_s (last update), pos_m [x y z], vel_mps
%   [vx vy vz] and acc_mps2 [ax ay az] (0 with the constant-velocity model)
%   relative to the antenna (x East, y North, z up; z = 0 in 2-D), alt_m
%   (height above ground, NaN in 2-D), hits, misses, dim.
tk = tracker.tracks;
d = tracker.p.dim;
rep = struct('id', {}, 'confirmed', {}, 't_s', {}, 'pos_m', {}, 'vel_mps', {}, 'acc_mps2', {}, 'alt_m', {}, ...
             'hits', {}, 'misses', {}, 'dim', {});
for i = 1:numel(tk)
    pos = zeros(1, 3);  vel = zeros(1, 3);  acc = zeros(1, 3);
    pos(1:d) = tk(i).x(1:d);
    vel(1:d) = tk(i).x(d+1:2*d);
    if numel(tk(i).x) > 2*d, acc(1:d) = tk(i).x(2*d+1:3*d); end
    alt = NaN;
    if d == 3, alt = pos(3) + tracker.p.height_m; end
    rep(i) = struct('id', tk(i).id, 'confirmed', tk(i).confirmed, 't_s', tk(i).t, 'pos_m', pos, ...
                    'vel_mps', vel, 'acc_mps2', acc, 'alt_m', alt, 'hits', tk(i).hits, 'misses', tk(i).misses, 'dim', d);
end
end

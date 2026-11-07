function rep = track_report(tracker)
%TRACK_REPORT  The current tracks for the display (one struct per track).
%   Fields: id, confirmed, t_s (last update), pos_m [x y z], vel_mps
%   [vx vy vz] and acc_mps2 [ax ay az] (0 with the constant-velocity model)
%   relative to the antenna (x East, y North, z up; z = 0 in 2-D), alt_m
%   (height above ground, NaN in 2-D), vr_mps (radial velocity of the track,
%   unambiguous), hits, misses, dim, method ('EKF' or 'PMBM'), existence
%   (PMBM: existence probability over the global hypotheses, EKF: NaN) and
%   modeProb (PMBM: the IMM model probabilities as text, e.g.
%   'CV 0.90 CA 0.04 CT 0.06'; EKF: '').
tk = tracker.tracks;
d = tracker.p.dim;
isPmbm = isfield(tracker, 'T');
method = 'EKF';
if isPmbm, method = 'PMBM'; end
rep = struct('id', {}, 'confirmed', {}, 't_s', {}, 'pos_m', {}, 'vel_mps', {}, 'acc_mps2', {}, 'alt_m', {}, ...
             'vr_mps', {}, 'hits', {}, 'misses', {}, 'dim', {}, 'method', {}, 'existence', {}, 'modeProb', {});
for i = 1:numel(tk)
    pos = zeros(1, 3);  vel = zeros(1, 3);  acc = zeros(1, 3);
    pos(1:d) = tk(i).x(1:d);
    vel(1:d) = tk(i).x(d+1:2*d);
    if numel(tk(i).x) >= 3*d, acc(1:d) = tk(i).x(2*d+1:3*d); end
    alt = NaN;
    if d == 3, alt = pos(3) + tracker.p.height_m; end
    ex = NaN;  modes = '';
    if isPmbm
        ex = tk(i).r;
        names = upper(tracker.p.pmbm.models);
        parts = cell(1, numel(names));
        for k = 1:numel(names), parts{k} = sprintf('%s %.2f', names{k}, tk(i).mu(k)); end
        modes = strjoin(parts, ' ');
    end
    rep(i) = struct('id', tk(i).id, 'confirmed', tk(i).confirmed, 't_s', tk(i).t, 'pos_m', pos, ...
                    'vel_mps', vel, 'acc_mps2', acc, 'alt_m', alt, 'vr_mps', pos*vel.'/norm(pos), ...
                    'hits', tk(i).hits, 'misses', tk(i).misses, 'dim', d, 'method', method, ...
                    'existence', ex, 'modeProb', modes);
end
end

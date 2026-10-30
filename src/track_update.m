function tr = track_update(tr, plots, t_s, flush)
%TRACK_UPDATE  Extended Kalman filter tracking of the plots, GNN association.
%   tr = TRACK_UPDATE(tr, plots, t_s) takes the new plots of time t_s (call
%   it every PRI, also without plots) and deletes the tracks that were not
%   updated for too long. tr comes from track_init(sys.track).
%   tr = TRACK_UPDATE(tr, plots, t_s, true) also associates the plots still
%   waiting in the batch (end of the data, skipped sector).
%
%   State relative to the antenna (x East, y North, z up), per axis the
%   position and its derivatives, ordered [positions; velocities
%   (; accelerations)]: 2-D [x y vx vy], 3-D [x y z vx vy vz], and with the
%   constant-acceleration model ('ca') [... ax ay (az)].
%   Prediction over dt, per axis with the spectral density q of p.q:
%     'cv' nearly constant velocity: F = [1 dt; 0 1],
%          Q = q*[dt^3/3 dt^2/2; dt^2/2 dt]                (white acceleration)
%     'ca' nearly constant acceleration: F = [1 dt dt^2/2; 0 1 dt; 0 0 1],
%          Q = q*[dt^5/20 dt^4/8 dt^3/6; dt^4/8 dt^3/3 dt^2/2; dt^3/6 dt^2/2 dt]
%                                                         (white jerk)
%   Measurement of a plot: 2-D z = [ground range; azimuth] (ground range =
%   range*cos(elevation) when the elevation is known), 3-D z = [range;
%   azimuth; elevation], with the nonlinear h(x) linearised at the predicted
%   state (EKF; the azimuth innovation is wrapped to +-pi).
%
%   Association: global nearest neighbour (GNN) on batches of plots. A plot
%   waits in the batch until the antenna has turned p.batchGap_s without a
%   new plot, or p.batchMax_s after the first plot of the batch. Then the
%   cost of every plot-track pair is d^2 + ln|S| when the plot lies in the
%   gate of the track (d^2 = nu'*inv(S)*nu <= p.gate, the track predicted to
%   the plot time), and the assignment (every track at most one plot, every
%   plot at most one track) has the largest number of pairs and, among those,
%   the smallest total cost (Hungarian algorithm). The plots of the batch are
%   then taken in their order: an assigned plot updates its track (EKF,
%   Joseph form), a plot without a track starts a tentative track at the
%   measured position with zero velocity (and acceleration) and the
%   uncertainty p.sigmaV0 (p.sigmaA0). A tentative track with p.confirmHits
%   updates is confirmed.
%   Missed scans: the antenna passed the predicted position half a scan ago
%   without an update; a track with more than p.maxMiss (tentative:
%   p.maxMissTentative) missed scans is deleted.
%   The updates are appended to tr.log ([time id confirmed x y z vx vy vz hits
%   plotIndex], z = 0 in 2-D; plotIndex counts the plots given to the tracker).

if nargin < 4, flush = false; end
p = tr.p;
% 1) the batch is complete: the antenna turned on without a new plot, or it is old enough
if ~isempty(tr.pending) && (t_s - tr.pendingLast >= p.batchGap_s || t_s - tr.pendingFirst >= p.batchMax_s)
    tr = associate(tr);
end
% 2) the new plots join the batch
for k = 1:numel(plots)
    tr.nPlots = tr.nPlots + 1;
    if isempty(tr.pending), tr.pendingFirst = t_s; end
    tr.pending(end+1, :) = [plots(k).time_s, plots(k).az_deg, plots(k).range_m, plots(k).el_deg, tr.nPlots];
    tr.pendingLast = t_s;
end
if flush && ~isempty(tr.pending)
    tr = associate(tr);
end

% 3) missed scans and deletion
if isempty(tr.tracks), return; end
miss = max(0, floor((t_s - [tr.tracks.t]) / p.scanTime_s - 0.5));
for j = 1:numel(tr.tracks)
    tr.tracks(j).misses = miss(j);
end
limit = p.maxMissTentative * ones(size(miss));
limit([tr.tracks.confirmed]) = p.maxMiss;
tr.tracks(miss > limit) = [];
end

function tr = associate(tr)
% GNN association of the plots of the batch, then the updates and new tracks.
p = tr.p;
B = tr.pending;
tr.pending = zeros(0, 5);
m = size(B, 1);
nT = numel(tr.tracks);
Z = zeros(p.dim, m);
usable = false(1, m);
for i = 1:m
    [Z(:, i), usable(i)] = measurement(B(i, :), p.dim);
end
cost = inf(nT, m);                                  % d^2 + ln|S| of the gated pairs
for j = 1:nT
    for i = find(usable)
        [~, ~, nu, ~, S] = predict_measure(tr.tracks(j), Z(:, i), B(i, 1), p);
        d2 = nu' * (S \ nu);
        if d2 <= p.gate
            cost(j, i) = d2 + log(det(S));
        end
    end
end
trackOf = gnn_assign(cost);                         % track of every plot (0 = none)

n = p.nd * p.dim;
for i = 1:m
    if ~usable(i), continue; end
    tp = B(i, 1);
    j = trackOf(i);
    if j > 0                                        % EKF update (Joseph form)
        tk = tr.tracks(j);
        [x, P, nu, H, S] = predict_measure(tk, Z(:, i), tp, p);
        K = P * H' / S;
        I_KH = eye(n) - K*H;
        tk.x = x + K*nu;
        tk.P = I_KH*P*I_KH' + K*p.R*K';
        tk.P = (tk.P + tk.P') / 2;
        tk.t = tp;
        tk.hits = tk.hits + 1;
        tk.misses = 0;
        tk.confirmed = tk.confirmed || tk.hits >= p.confirmHits;
        tr.tracks(j) = tk;
    else                                            % new tentative track
        [pos, J] = position(Z(:, i), p.dim);
        P = zeros(n);
        d = p.dim;
        P(1:d, 1:d) = J*p.R*J';
        P(d+1:2*d, d+1:2*d) = p.sigmaV0^2 * eye(d);
        if p.nd == 3
            P(2*d+1:end, 2*d+1:end) = p.sigmaA0^2 * eye(d);
        end
        tk = struct('id', tr.nextId, 'confirmed', p.confirmHits <= 1, 'x', [pos; zeros(n - d, 1)], ...
                    'P', P, 't', tp, 'hits', 1, 'misses', 0);
        tr.nextId = tr.nextId + 1;
        tr.tracks(end+1) = tk;
        j = numel(tr.tracks);
    end
    tk = tr.tracks(j);
    xyz = zeros(1, 6);
    xyz([1:p.dim, 4:3+p.dim]) = tk.x(1:2*p.dim).';
    tr.log(end+1, :) = [tk.t, tk.id, tk.confirmed, xyz, tk.hits, B(i, 5)];
end
end

function trackOf = gnn_assign(cost)
% Assignment of plots (columns) to tracks (rows) with the most finite-cost
% pairs and, among those, the smallest total cost. trackOf(i) = row of
% column i (0 = not assigned).
m = size(cost, 2);
trackOf = zeros(1, m);
rows = find(any(isfinite(cost), 2));
cols = find(any(isfinite(cost), 1));
if isempty(rows), return; end
C = cost(rows, cols);
ok = isfinite(C);
C(ok) = C(ok) - min(C(ok));                         % >= 0
% a forbidden pair costs more than any assignment with one more allowed pair
C(~ok) = (max(C(ok)) + 1) * (min(size(C)) + 1);
if numel(rows) <= numel(cols)
    rowOfCol = hungarian(C);                        % row of every column
    for c = 1:numel(cols)
        r = rowOfCol(c);
        if r > 0 && ok(r, c), trackOf(cols(c)) = rows(r); end
    end
else
    colOfRow = hungarian(C.');                      % column of every row
    for r = 1:numel(rows)
        c = colOfRow(r);
        if c > 0 && ok(r, c), trackOf(cols(c)) = rows(r); end
    end
end
end

function rowOf = hungarian(A)
% Minimum-cost assignment of every row of A (n x m, n <= m) to its own column
% (shortest augmenting paths with row and column potentials, O(n^2 m)).
% rowOf(j) = row assigned to column j, 0 = none. Index 1 of u, v, p, way
% stands for the virtual row/column 0.
[n, m] = size(A);
u = zeros(n+1, 1);  v = zeros(m+1, 1);
p = zeros(m+1, 1);  way = zeros(m+1, 1);
for i = 1:n
    p(1) = i;
    j0 = 0;
    minv = inf(m+1, 1);
    used = false(m+1, 1);
    while true
        used(j0+1) = true;
        i0 = p(j0+1);
        delta = inf;
        j1 = 0;
        for j = 1:m
            if ~used(j+1)
                cur = A(i0, j) - u(i0+1) - v(j+1);
                if cur < minv(j+1)
                    minv(j+1) = cur;  way(j+1) = j0;
                end
                if minv(j+1) < delta
                    delta = minv(j+1);  j1 = j;
                end
            end
        end
        for j = 0:m
            if used(j+1)
                u(p(j+1)+1) = u(p(j+1)+1) + delta;
                v(j+1) = v(j+1) - delta;
            else
                minv(j+1) = minv(j+1) - delta;
            end
        end
        j0 = j1;
        if p(j0+1) == 0, break; end
    end
    while j0 ~= 0                                   % augment along the path
        j1 = way(j0+1);
        p(j0+1) = p(j1+1);
        j0 = j1;
    end
end
rowOf = p(2:end).';
end

function [z, ok] = measurement(b, dim)
% Measurement vector of a plot of the batch ([time az range el plotIndex]).
az = deg2rad(b(2));
el = deg2rad(b(4));
if dim == 3
    z = [b(3); az; el];
    ok = isfinite(el);
else
    rho = b(3);
    if isfinite(el), rho = rho * cos(el); end
    z = [rho; az];
    ok = true;
end
end

function [x, P, nu, H, S] = predict_measure(tk, z, tp, p)
% Track predicted to the plot time, innovation, Jacobian and its covariance.
[x, P] = predict(tk.x, tk.P, tp - tk.t, p);
[nu, H] = innovation(z, x, p);
S = H*P*H' + p.R;
end

function [x, P] = predict(x, P, dt, p)
% Constant-velocity or constant-acceleration prediction over dt (may be slightly negative).
a = abs(dt);
if p.nd == 2
    F1 = [1 dt; 0 1];
    Q1 = [a^3/3 a^2/2; a^2/2 a];
else
    F1 = [1 dt dt^2/2; 0 1 dt; 0 0 1];
    Q1 = [a^5/20 a^4/8 a^3/6; a^4/8 a^3/3 a^2/2; a^3/6 a^2/2 a];
end
F = kron(F1, eye(p.dim));                           % state order [positions; velocities; ...]
x = F*x;
P = F*P*F' + kron(Q1, diag(p.q));
end

function [nu, H] = innovation(z, x, p)
% Innovation z - h(x) and Jacobian H of h at x.
px = x(1);  py = x(2);
rh2 = px^2 + py^2;  rh = sqrt(rh2);
if p.dim == 3
    pz = x(3);
    R2 = rh2 + pz^2;  R = sqrt(R2);
    h = [R; atan2(px, py); atan2(pz, rh)];
    Hp = [px/R,           py/R,           pz/R;
          py/rh2,         -px/rh2,        0;
          -px*pz/(R2*rh), -py*pz/(R2*rh), rh/R2];
else
    h = [rh; atan2(px, py)];
    Hp = [px/rh,  py/rh;
          py/rh2, -px/rh2];
end
H = [Hp, zeros(p.dim, numel(x) - p.dim)];
nu = z - h;
nu(2) = mod(nu(2) + pi, 2*pi) - pi;                 % azimuth wrap
end

function [pos, J] = position(z, dim)
% Position of a measurement and the Jacobian d(position)/d(z).
az = z(2);
if dim == 3
    R = z(1);  el = z(3);
    pos = R * [cos(el)*sin(az); cos(el)*cos(az); sin(el)];
    J = [cos(el)*sin(az),  R*cos(el)*cos(az), -R*sin(el)*sin(az);
         cos(el)*cos(az), -R*cos(el)*sin(az), -R*sin(el)*cos(az);
         sin(el),          0,                  R*cos(el)];
else
    rho = z(1);
    pos = rho * [sin(az); cos(az)];
    J = [sin(az),  rho*cos(az);
         cos(az), -rho*sin(az)];
end
end

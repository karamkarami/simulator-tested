function tr = track_pmbm(tr, plots, t_s, flush)
%TRACK_PMBM  Multi-target Bayesian tracker: PMBM filter with IMM-UKF state estimation.
%   tr = TRACK_PMBM(tr, plots, t_s, flush) is the tracker of track_update.m for
%   cfg.track.method = 'pmbm' (tr from track_init.m, parameters p.pmbm from
%   track_setup.m). Call it every PRI with the new plots of time t_s (also
%   without plots); flush = true also associates the plots still waiting in
%   the batch (end of the data, skipped sector).
%
%   Multi-target density: Poisson multi-Bernoulli mixture (PMBM).
%   - Undetected targets: a Poisson point process, here with a uniform
%     intensity p.pmbm.logLambdaU in the measurement space (range x azimuth
%     (x elevation) (x radial velocity)): its update by a plot is the birth of
%     a new Bernoulli component (measurement-driven birth).
%   - Detected targets: Bernoulli components (tracks); every track holds
%     single-target hypotheses (STH): existence probability r and an IMM
%     density, the history of plots it explains.
%   - Global hypotheses: tr.G (rows) picks one STH per track; tr.logw are
%     their log weights (normalised). They are the data association
%     hypotheses of the multi-Bernoulli mixture.
%   Update with the plots of a batch (the GNN batch of the EKF, p.batchGap_s,
%   p.batchMax_s): for every global hypothesis the cost of plot j to track i
%   is -log(r*ps*pd*l_ij / l0_i) (l_ij the IMM-UKF likelihood of the plot,
%   l0_i the weight of no plot: 1 - r*pd), of plot j as a new target or a
%   false plot -log(lambdaC + pd*lambdaU); Murty's algorithm gives the
%   ceil(maxHyp * w) best assignments of each global hypothesis. A track
%   gated by a plot of the batch, whose antenna look lies inside the batch,
%   is missed when it gets no plot; a track whose look comes later keeps its
%   hypothesis unchanged (its plot may still come).
%   Missed looks: when the antenna has passed the predicted azimuth of an
%   STH (p.pmbm.latency_s ago, so that its plot would be in the tracker)
%   without a plot, r becomes r*(1 - pd)/(1 - r*pd) (after the survival
%   ps per scan) and the global hypotheses holding it get the factor
%   1 - r*pd.
%   Pruning: STH below p.pmbm.rPrune become "no target", global hypotheses
%   below p.pmbm.logHypPrune or beyond p.pmbm.maxHyp are dropped, tracks
%   whose existence (sum over the global hypotheses) is below p.pmbm.rDelete
%   are deleted. Speed gate (p.minSpeed_mps) and maxMiss as with the EKF.
%   Confirmation: the probability that a target made the newest plot of the
%   track (sum over all global hypotheses of the existence of the track that
%   holds this plot as its newest one, whatever its earlier plots) is at
%   least p.pmbm.rConfirm, and the track has at least p.confirmHits plots in
%   the most likely hypothesis (it stays confirmed). The tracks shown are those of the most likely global
%   hypothesis; when it changes, a confirmed track that is no longer in it
%   hands its number (and its confirmation) to the track of the new one that
%   holds its newest plot (the same target, another history: e.g. its first
%   plot now explained as a false plot), so the track numbers stay with the
%   targets.
%
%   State estimation per STH: IMM over the models p.pmbm.models ('cv', 'ca',
%   'ct') on the common state [position; velocity; acceleration; turn rate]
%   relative to the antenna (x East, y North (, z up)); 2-D: 7 values, 3-D:
%   10. Mixing with the Markov chain per scan (switchProb, raised to the
%   power dt/scan time), then per model the UKF: CV and CA are linear (the
%   unscented transform of a linear model is exact: plain Kalman
%   prediction), CT (coordinated turn in the horizontal plane with the turn
%   rate of the state, positive = counter-clockwise in x-y) by sigma
%   points; the measurement (range, azimuth (, elevation), radial velocity
%   as with the EKF, also its unfolding with the staggered PRF) is predicted
%   by sigma points for every model. 2-D: ground range and azimuth (slant
%   range without the elevation, as the 2-D EKF).
%
%   Output as the EKF: tr.tracks (the tracks of the most likely global
%   hypothesis: id, confirmed, x = the IMM estimate [position; velocity;
%   acceleration; turn rate], P, t, hits, misses, r = probability that a
%   target made its newest plot (its existence), mu = model
%   probabilities) and tr.log (one row per plot that updates a track of the
%   most likely hypothesis: [time id confirmed x y z vx vy vz hits plotIndex]).

if nargin < 4, flush = false; end
p = tr.p;
m = p.pmbm;
changed = false;
% 1) the batch is complete: the antenna turned on without a new plot, or it is old enough
if ~isempty(tr.pending) && (t_s - tr.pendingLast >= p.batchGap_s || t_s - tr.pendingFirst >= p.batchMax_s)
    tr = associate(tr);
    tr = maintain(tr, t_s - m.latency_s);
    tr.tMaint = t_s;
    changed = true;
end
% 2) the new plots join the batch
for k = 1:numel(plots)
    tr.nPlots = tr.nPlots + 1;
    if isempty(tr.pending), tr.pendingFirst = t_s; end
    tr.pending(end+1, :) = [plots(k).time_s, plots(k).az_deg, plots(k).range_m, plots(k).el_deg, tr.nPlots, ...
                            plots(k).vr_mps];
    tr.pendingLast = t_s;
end
if flush && ~isempty(tr.pending)
    tr = associate(tr);
    changed = true;
end
% 3) missed looks (only when no plot is waiting: it could be the one of the look)
if isempty(tr.pending) && (flush || t_s - tr.tMaint >= m.maintStep_s)
    [tr, ch] = maintain(tr, t_s - m.latency_s);
    tr.tMaint = t_s;
    changed = changed || ch;
end
if changed, tr = publish(tr); end
end

%% ---------------------------------------------------------------- association
function tr = associate(tr)
% PMBM update with the plots of the batch.
p = tr.p;  m = p.pmbm;
B = tr.pending;
tr.pending = zeros(0, 6);
tr.batch = tr.batch + 1;
mAll = size(B, 1);
Z = zeros(p.dim, mAll);  Zvr = nan(1, mAll);  usable = false(1, mAll);
for i = 1:mAll
    [Z(:, i), usable(i), Zvr(i)] = measurement(B(i, :), p.dim);
end
B = B(usable, :);  Z = Z(:, usable);  Zvr = Zvr(usable);
nPl = size(B, 1);
if nPl == 0, return; end
[~, order] = sort(B(:, 1));                         % plots in time order
B = B(order, :);  Z = Z(:, order);  Zvr = Zvr(order);
T = p.scanTime_s;
tMissLimit = max(B(:, 1)) + p.batchGap_s;           % a look up to here is inside the batch
hasVr = p.useVr & isfinite(Zvr);                    % the radial velocity is part of the plot
% new target (Poisson part) or false plot
logC = m.logLambdaC + hasVr * m.logVrDensity;
logE = log(m.pd) + m.logLambdaU + hasVr * m.logVrDensity;
logNew = logsumexp2(logC, logE);
rNew = exp(logE - logNew);

% single-target hypotheses of the tracks for the plots: weights and updated densities
nTr = numel(tr.T);
L = cell(1, nTr);  L0 = cell(1, nTr);  Upd = cell(1, nTr);  Miss = cell(1, nTr);
for i = 1:nTr
    H = tr.T(i).h;
    nh = numel(H);
    L{i} = -inf(nh, nPl);  L0{i} = zeros(nh, 1);
    Upd{i} = cell(nh, nPl);  Miss{i} = cell(nh, 1);
    used = unique(tr.G(:, i)).';
    for h = used
        sth = H(h);
        if sth.r <= 0, continue; end
        [xm, Pm] = imm_estimate(sth);
        for j = 1:nPl
            tp = B(j, 1);
            if tp < sth.tOpp + T/2, continue; end     % this pass already had its plot (or its miss)
            if ~coarse_gate(xm, Pm, sth.t, Z(:, j), tp, p), continue; end
            [ok, logl, hn] = imm_ukf_update(sth, Z(:, j), Zvr(j), tp, p, hasVr(j));
            if ~ok, continue; end
            rp = sth.r * m.ps^((tp - sth.tOpp) / T);
            L{i}(h, j) = log(rp) + log(m.pd) + logl;
            hn.plot = B(j, 5);  hn.batch = tr.batch;
            hn.hist = [hn.hist(max(1, end-6):end), B(j, 5)];   % its newest plots
            hn.tNext = next_look(hn, p);
            Upd{i}{h, j} = hn;
        end
        if any(isfinite(L{i}(h, :))) && sth.tNext <= tMissLimit   % its look is in this batch: missed
            [Miss{i}{h}, L0{i}(h)] = miss_update(sth, tr.T(i).confirmed, p);
        end
    end
end

% global hypotheses: k best assignments of each
nG = size(tr.G, 1);
w = exp(tr.logw - max(tr.logw));  w = w / sum(w);
newRows = zeros(0, nTr + nPl);
newLogW = zeros(0, 1);
keys = cell(1, nTr);                                 % per track the children [parent STH, plot (0 = none)]
for i = 1:nTr, keys{i} = zeros(0, 2); end
for g = 1:nG
    % columns: the tracks whose STH of this hypothesis gates a plot, then one new target per plot
    part = zeros(1, 0);
    for i = 1:nTr
        if any(isfinite(L{i}(tr.G(g, i), :))), part(end+1) = i; end %#ok<AGROW>
    end
    nP = numel(part);
    C = inf(nPl, nP + nPl);
    base = 0;
    for c = 1:nP
        i = part(c);  h = tr.G(g, i);
        C(:, c) = -(L{i}(h, :) - L0{i}(h)).';
        base = base + L0{i}(h);
    end
    C(sub2ind(size(C), 1:nPl, nP + (1:nPl))) = -logNew;
    k = max(1, ceil(m.maxHyp * w(g)));
    [A, cost] = murty(C, k);
    for a = 1:size(A, 1)
        row = zeros(1, nTr + nPl);
        plotOf = zeros(1, nTr);
        for j = 1:nPl
            if A(a, j) <= nP, plotOf(part(A(a, j))) = j; end
        end
        for i = 1:nTr
            key = [tr.G(g, i), plotOf(i)];
            idx = find(keys{i}(:, 1) == key(1) & keys{i}(:, 2) == key(2), 1);
            if isempty(idx)
                keys{i}(end+1, :) = key;
                idx = size(keys{i}, 1);
            end
            row(i) = idx;
        end
        for j = 1:nPl                               % new track of plot j: 2 = born, 1 = no target
            row(nTr + j) = 1 + (A(a, j) == nP + j);
        end
        newRows(end+1, :) = row; %#ok<AGROW>
        newLogW(end+1, 1) = tr.logw(g) + base - cost(a); %#ok<AGROW>
    end
end

% the children become the STH of the tracks
for i = 1:nTr
    H = tr.T(i).h;
    K = keys{i};
    Hn = repmat(empty_sth(), 1, size(K, 1));
    for c = 1:size(K, 1)
        h = K(c, 1);  j = K(c, 2);
        if j > 0
            Hn(c) = Upd{i}{h, j};
        elseif ~isempty(Miss{i}{h})
            Hn(c) = Miss{i}{h};
        else
            Hn(c) = H(h);                           % not looked at in this batch: unchanged
        end
    end
    tr.T(i).h = Hn;
end
for j = 1:nPl                                       % new tracks, one per plot
    born = new_sth(Z(:, j), Zvr(j), B(j, :), rNew(j), p);
    born.batch = tr.batch;
    born.tNext = next_look(born, p);
    tr.T(end+1) = struct('id', 0, 'confirmed', false, 'h', [empty_sth(), born], 'rMarg', rNew(j), 'rTarget', 0, ...
                         'plot', B(j, 5));
end
tr.G = newRows;
tr.logw = newLogW;
tr = cleanup(tr);
end

function [tr, changed] = maintain(tr, tSafe)
% Missed looks: the antenna passed the predicted azimuth of an STH before
% tSafe (its plot would be in the tracker by now) without a plot.
p = tr.p;
changed = false;
for i = 1:numel(tr.T)
    for h = 1:numel(tr.T(i).h)
        sth = tr.T(i).h(h);
        while sth.r > 0 && sth.tNext <= tSafe
            [sth, l0] = miss_update(sth, tr.T(i).confirmed, p);
            rows = tr.G(:, i) == h;
            tr.logw(rows) = tr.logw(rows) + l0;
            changed = true;
        end
        tr.T(i).h(h) = sth;
    end
end
if changed, tr = cleanup(tr); end
end

function [sth, logl0] = miss_update(sth, confirmed, p)
% No plot at the look of the STH: survival, then the missed detection.
m = p.pmbm;
rp = sth.r * m.ps^((sth.tNext - sth.tOpp) / p.scanTime_s);
l0 = 1 - rp * m.pd;
logl0 = log(l0);
sth.r = rp * (1 - m.pd) / l0;
sth.tOpp = sth.tNext;
sth.misses = sth.misses + 1;
limit = p.maxMissTentative;
if confirmed, limit = p.maxMiss; end
if sth.misses > limit, sth.r = 0; end             % as the EKF: too many missed scans
if sth.r > 0
    sth.tNext = next_look(sth, p);
else
    sth.tNext = Inf;
end
end

function tr = cleanup(tr)
% Pruning of STH, tracks and global hypotheses; track numbers, confirmation.
p = tr.p;  m = p.pmbm;
% STH: too unlikely or too slow -> no target; one "no target" STH per track
for i = 1:numel(tr.T)
    H = tr.T(i).h;
    for h = 1:numel(H)
        if H(h).r > 0 && H(h).r < m.rPrune, H(h).r = 0; end
        if H(h).r > 0 && p.minSpeed_mps > 0 && H(h).hits >= p.minSpeedHits
            xm = imm_estimate(H(h));
            if norm(xm(p.dim+1:2*p.dim)) < p.minSpeed_mps, H(h).r = 0; end
        end
    end
    zero = find([H.r] <= 0);
    if numel(zero) > 1
        col = tr.G(:, i);
        col(ismember(col, zero)) = zero(1);
        tr.G(:, i) = col;
    end
    tr.T(i).h = H;
end
for pass = 1:2
    tr = merge_rows(tr);
    % global hypotheses: weight and number
    keep = tr.logw >= max(tr.logw) + m.logHypPrune;
    tr.G = tr.G(keep, :);  tr.logw = tr.logw(keep);
    [~, o] = sort(tr.logw, 'descend');
    o = o(1:min(numel(o), m.maxHyp));
    tr.G = tr.G(o, :);  tr.logw = tr.logw(o);
    tr.logw = tr.logw - logsumexp(tr.logw);
    % tracks: existence over the global hypotheses
    w = exp(tr.logw);
    nTr = numel(tr.T);
    del = false(1, nTr);
    for i = 1:nTr
        r = [tr.T(i).h.r].';
        tr.T(i).rMarg = sum(w .* reshape(r(tr.G(:, i)), [], 1));
        del(i) = tr.T(i).rMarg < m.rDelete;
    end
    if ~any(del) && pass > 1, break; end
    tr.T(del) = [];
    tr.G(:, del) = [];
end
% unused STH are removed
for i = 1:numel(tr.T)
    [used, ~, idx] = unique(tr.G(:, i));
    tr.T(i).h = tr.T(i).h(used);
    tr.G(:, i) = idx;
end
% probability that a target made a plot: sum over the global hypotheses of
% the existence of the track that holds it as its newest plot (whatever its
% earlier plots: the hypotheses may differ only in the history of a target)
w = exp(tr.logw);
pl = zeros(0, 1);  pw = zeros(0, 1);
for i = 1:numel(tr.T)
    H = tr.T(i).h;
    for g = 1:numel(w)
        sth = H(tr.G(g, i));
        if sth.r > 0
            pl(end+1, 1) = sth.plot; %#ok<AGROW>
            pw(end+1, 1) = w(g) * sth.r; %#ok<AGROW>
        end
    end
end
[up, ~, k] = unique(pl);
pSum = accumarray(k, pw, [numel(up) 1]);
% track numbers (in the order of their first plot) and confirmation
[~, best] = max(tr.logw);
[~, o] = sort([tr.T.plot]);
for i = o
    if tr.T(i).id == 0
        tr.T(i).id = tr.nextId;
        tr.nextId = tr.nextId + 1;
    end
    tr.T(i).rTarget = 0;
    if isempty(best), continue; end
    sth = tr.T(i).h(tr.G(best, i));
    if sth.r <= 0, continue; end
    tr.T(i).rTarget = min(1, pSum(up == sth.plot));
    if tr.T(i).rTarget >= m.rConfirm && sth.hits >= p.confirmHits
        tr.T(i).confirmed = true;
    end
end
end

function tr = merge_rows(tr)
% Equal global hypotheses (the same STH of every track) are merged.
if isempty(tr.T)
    tr.G = zeros(1, 0);  tr.logw = 0;               % one hypothesis: no track
    return
end
[U, ~, idx] = unique(tr.G, 'rows');
if size(U, 1) == size(tr.G, 1), return; end
lw = -inf(size(U, 1), 1);
for k = 1:numel(idx)
    lw(idx(k)) = logsumexp2(lw(idx(k)), tr.logw(k));
end
tr.G = U;  tr.logw = lw;
end

function tr = publish(tr)
% The tracks of the most likely global hypothesis (tr.tracks) and the log of
% the updates by the plots of the newest batch.
p = tr.p;
d = p.dim;
tk = struct('id', {}, 'confirmed', {}, 'x', {}, 'P', {}, 't', {}, 'hits', {}, 'misses', {}, 'r', {}, 'mu', {});
rows = zeros(0, 11);
if ~isempty(tr.logw) && ~isempty(tr.T)
    [~, best] = max(tr.logw);
    tr = keep_numbers(tr, best);
    for i = 1:numel(tr.T)
        sth = tr.T(i).h(tr.G(best, i));
        if sth.r <= 0, continue; end
        [x, P] = imm_estimate(sth);
        tk(end+1) = struct('id', tr.T(i).id, 'confirmed', tr.T(i).confirmed, 'x', x, 'P', P, 't', sth.t, ...
                           'hits', sth.hits, 'misses', sth.misses, 'r', tr.T(i).rTarget, 'mu', sth.mu); %#ok<AGROW>
        if sth.batch == tr.batch && sth.plot > tr.logged   % updated by a plot of this batch
            xyz = zeros(1, 6);
            xyz([1:d, 4:3+d]) = x(1:2*d).';
            rows(end+1, :) = [sth.t, tr.T(i).id, tr.T(i).confirmed, xyz, sth.hits, sth.plot]; %#ok<AGROW>
        end
    end
end
if ~isempty(rows)
    [~, o] = sort(rows(:, 11));
    tr.log = [tr.log; rows(o, :)];
    tr.logged = max(rows(:, 11));
end
[~, o] = sort([tk.id]);
tr.tracks = tk(o);
% the confirmed tracks shown, with their newest plot (keep_numbers)
tr.shown = struct('id', {}, 'plot', {});
if ~isempty(tr.logw) && ~isempty(tr.T)
    for i = 1:numel(tr.T)
        sth = tr.T(i).h(tr.G(best, i));
        if sth.r > 0 && tr.T(i).confirmed
            tr.shown(end+1) = struct('id', tr.T(i).id, 'plot', sth.plot);
        end
    end
end
end

function tr = keep_numbers(tr, best)
% A confirmed track shown before that is not in the most likely global
% hypothesis any more: the track of this hypothesis that holds its newest
% plot takes its number and its confirmation (the two tracks swap numbers).
if isempty(tr.shown), return; end
ids = [tr.T.id];
exists = false(1, numel(tr.T));
for i = 1:numel(tr.T), exists(i) = tr.T(i).h(tr.G(best, i)).r > 0; end
for s = 1:numel(tr.shown)
    P = tr.shown(s);
    iP = find(ids == P.id, 1);
    if ~isempty(iP) && exists(iP), continue; end     % still shown
    for iQ = find(exists)
        if any([tr.shown.id] == tr.T(iQ).id), continue; end   % shown as a confirmed track itself
        if any(tr.T(iQ).h(tr.G(best, iQ)).hist == P.plot)
            if ~isempty(iP)
                tr.T(iP).id = tr.T(iQ).id;
                tr.T(iP).confirmed = tr.T(iQ).confirmed;
            end
            tr.T(iQ).id = P.id;
            tr.T(iQ).confirmed = true;
            ids = [tr.T.id];
            break
        end
    end
end
end

%% ---------------------------------------------------------------- IMM-UKF
function [ok, logl, sth] = imm_ukf_update(sth, z, zvr, tp, p, hasVr)
% IMM-UKF update of an STH with the plot z (radial velocity zvr) of time tp.
% ok: inside the gate of at least one model; logl: log of the likelihood of
% the plot (mixture over the models); sth: the updated STH (r = 1).
m = p.pmbm;
M = numel(m.models);
d = p.dim;
dt = tp - sth.t;
% interaction (mixing) with the transition matrix of |dt| (uniform switching)
lk = m.switchLambda ^ (abs(dt) / p.scanTime_s);
Pi = (1 - lk) / M * ones(M) + lk * eye(M);
c = sth.mu * Pi;                                     % predicted model probabilities
W = (Pi .* sth.mu.') ./ max(c, realmin);             % W(i, j) = P(model i before | model j now)
n = m.n;
xp = zeros(n, M);  Pp = zeros(n, n, M);
for j = 1:M
    x0 = sth.x * W(:, j);
    P0 = zeros(n);
    for i = 1:M
        dx = sth.x(:, i) - x0;
        P0 = P0 + W(i, j) * (sth.P(:, :, i) + dx*dx');
    end
    [xp(:, j), Pp(:, :, j)] = predict_model(x0, (P0 + P0')/2, dt, m.models{j}, p);
end
% measurement prediction of every model (with the radial velocity when it is in the plot)
nz = d + hasVr;
Zh = zeros(nz, M);  S = zeros(nz, nz, M);  Cxz = zeros(n, nz, M);
Rm = p.R;
if hasVr, Rm = blkdiag(p.R, p.Rvr); end
for j = 1:M
    [Zh(:, j), S(:, :, j), Cxz(:, :, j)] = ukf_measure(xp(:, j), Pp(:, :, j), Rm, hasVr, p);
end
useVr = hasVr;
zz = z;
if hasVr                                            % the plot velocity is known modulo vrSpan
    k = 0;
    if ~p.vrUnambiguous
        sd = sqrt(max(squeeze(S(end, end, :)) - p.Rvr, 0));
        hv = Zh(end, :) * c.';
        useVr = max(sd) < p.vrSpan / 6;              % unfolded only with a well predicted range rate
        k = round((hv - zvr) / p.vrSpan);
    end
    if useVr
        zz = [z; zvr + k*p.vrSpan];
    else
        Zh = Zh(1:d, :);  S = S(1:d, 1:d, :);  Cxz = Cxz(:, 1:d, :);
    end
end
gate = m.gate;
if useVr, gate = m.gateVr; end
logLam = zeros(1, M);
d2 = zeros(1, M);
nu = zeros(numel(zz), M);
for j = 1:M
    v = zz - Zh(:, j);
    v(2) = mod(v(2) + pi, 2*pi) - pi;                % azimuth wrap
    nu(:, j) = v;
    Sj = S(:, :, j);
    d2(j) = v' * (Sj \ v);
    logLam(j) = -0.5 * (d2(j) + log(det(2*pi*Sj)));
end
ok = min(d2) <= gate;
lw = log(max(c, realmin)) + logLam;
logl = logsumexp(lw);
if hasVr && ~useVr, logl = logl + m.logVrDensity; end   % the plot velocity is not explained: uniform
if ~ok, return; end
for j = 1:M                                         % UKF update of every model
    Sj = S(:, :, j);
    K = Cxz(:, :, j) / Sj;
    xp(:, j) = xp(:, j) + K * nu(:, j);
    Pj = Pp(:, :, j) - K * Sj * K';
    Pp(:, :, j) = (Pj + Pj') / 2;
end
sth.x = xp;  sth.P = Pp;
sth.mu = exp(lw - max(lw));
sth.mu = sth.mu / sum(sth.mu);
sth.t = tp;  sth.tOpp = tp;
sth.r = 1;
sth.hits = sth.hits + 1;
sth.misses = 0;
end

function [x, P] = predict_model(x, P, dt, model, p)
% Prediction of one model of the IMM over dt (may be slightly negative).
m = p.pmbm;
d = p.dim;
n = m.n;
a = abs(dt);
I = eye(d);
ip = 1:d;  iv = d+1:2*d;  ia = 2*d+1:3*d;  iw = n;
switch model
    case 'cv'
        F = zeros(n);
        F(ip, ip) = I;  F(ip, iv) = dt*I;  F(iv, iv) = I;       % acceleration, turn rate -> 0
        Q = zeros(n);
        Q([ip iv], [ip iv]) = kron([a^3/3 a^2/2; a^2/2 a], diag(m.qCv));
        Q(ia, ia) = m.varUnusedA * I;  Q(iw, iw) = m.varUnusedW;
        x = F*x;  P = F*P*F' + Q;
    case 'ca'
        F = zeros(n);
        F(ip, ip) = I;  F(ip, iv) = dt*I;  F(ip, ia) = dt^2/2*I;
        F(iv, iv) = I;  F(iv, ia) = dt*I;  F(ia, ia) = I;          % turn rate -> 0
        Q = zeros(n);
        Q([ip iv ia], [ip iv ia]) = kron([a^5/20 a^4/8 a^3/6; a^4/8 a^3/3 a^2/2; a^3/6 a^2/2 a], diag(m.qCa));
        Q(iw, iw) = m.varUnusedW;
        x = F*x;  P = F*P*F' + Q;
    case 'ct'                                       % coordinated turn: unscented transform
        X = sigma_points(x, P, m.ukfScale);
        for s = 1:size(X, 2)
            X(:, s) = f_ct(X(:, s), dt, d);
        end
        x = X * m.Wm.';
        D = X - x;
        P = (D .* m.Wc) * D.';
        Q = zeros(n);
        Q([ip iv], [ip iv]) = kron([a^3/3 a^2/2; a^2/2 a], diag(m.qCt));
        Q(ia, ia) = m.varUnusedA * I;  Q(iw, iw) = m.qTurn * a;
        P = P + Q;
end
P = (P + P') / 2;
end

function x = f_ct(x, dt, d)
% Coordinated turn in the horizontal plane (turn rate x(end), counter-clockwise
% in x-y), constant vertical velocity; the acceleration becomes 0.
w = x(end);
vx = x(d+1);  vy = x(d+2);
wt = w * dt;
if abs(wt) < 1e-6
    sw = dt - w^2*dt^3/6;  cw = w*dt^2/2;           % sin(wt)/w, (1 - cos(wt))/w
else
    sw = sin(wt) / w;  cw = (1 - cos(wt)) / w;
end
x(1) = x(1) + sw*vx - cw*vy;
x(2) = x(2) + cw*vx + sw*vy;
x(d+1) = cos(wt)*vx - sin(wt)*vy;
x(d+2) = sin(wt)*vx + cos(wt)*vy;
if d == 3, x(3) = x(3) + dt*x(6); end
x(2*d+1:3*d) = 0;
end

function [zh, S, C] = ukf_measure(x, P, Rm, withVr, p)
% Unscented transform of the measurement: mean zh, covariance S (with the
% measurement noise Rm) and the state-measurement cross covariance C.
m = p.pmbm;
X = sigma_points(x, P, m.ukfScale);
Zs = h_meas(X, p.dim, withVr);
Zs(2, :) = Zs(2, 1) + mod(Zs(2, :) - Zs(2, 1) + pi, 2*pi) - pi;   % azimuth around the centre point
zh = Zs * m.Wm.';
Dz = Zs - zh;
Dx = X - x;
S = (Dz .* m.Wc) * Dz.' + Rm;
S = (S + S') / 2;
C = (Dx .* m.Wc) * Dz.';
end

function Z = h_meas(X, d, withVr)
% Measurement of the states (columns of X): 3-D [range; azimuth; elevation],
% 2-D [ground range; azimuth], and the range rate (2-D: of the ground range).
px = X(1, :);  py = X(2, :);
rh = sqrt(px.^2 + py.^2);
az = atan2(px, py);
if d == 3
    pz = X(3, :);
    R = sqrt(rh.^2 + pz.^2);
    Z = [R; az; atan2(pz, rh)];
else
    R = rh;
    Z = [rh; az];
end
if withVr
    Z(end+1, :) = sum(X(1:d, :) .* X(d+1:2*d, :), 1) ./ R;
end
end

function X = sigma_points(x, P, scale)
% 2n+1 sigma points of N(x, P) (Cholesky factor; P only positive semi-definite:
% a small diagonal loading, raised until the factorization succeeds).
[Lc, flag] = chol(P, 'lower');
load = 1e-12 * max(abs(diag(P)));
while flag ~= 0
    [Lc, flag] = chol(P + load * eye(numel(x)), 'lower');
    load = 10 * load;
end
S = scale * Lc;
X = [x, x + S, x - S];
end

function [x, P] = imm_estimate(sth)
% Combined estimate of the IMM (moment matching over the models).
x = sth.x * sth.mu.';
if nargout > 1
    P = zeros(numel(x));
    for j = 1:numel(sth.mu)
        dx = sth.x(:, j) - x;
        P = P + sth.mu(j) * (sth.P(:, :, j) + dx*dx');
    end
end
end

%% ---------------------------------------------------------------- track birth, looks
function sth = new_sth(z, zvr, b, r, p)
% Bernoulli born from a plot: at the measured position, velocity unknown
% (the velocity prior conditioned on the radial velocity when it is unambiguous).
m = p.pmbm;
d = p.dim;
n = m.n;
[pos, J] = position(z, d);
x = zeros(n, 1);
P = zeros(n);
x(1:d) = pos;
P(1:d, 1:d) = J*p.R*J';
D = p.sigmaV0^2 * eye(d);                          % velocity prior (3-D: vertical up to maxClimb)
if d == 3, D(3, 3) = m.sigmaClimb0^2; end
P(d+1:2*d, d+1:2*d) = D;
P(2*d+1:3*d, 2*d+1:3*d) = p.sigmaA0^2 * eye(d);
P(n, n) = m.sigmaTurn0^2;
if p.useVr && p.vrUnambiguous && isfinite(zvr)      % the prior conditioned on the radial velocity
    u = pos / norm(pos);
    s2 = u'*D*u + p.Rvr;
    x(d+1:2*d) = D*u * zvr / s2;
    Pv = D - (D*u)*(u'*D) / s2;
    P(d+1:2*d, d+1:2*d) = (Pv + Pv') / 2;
end
M = numel(m.models);
sth = empty_sth();
sth.r = r;
sth.x = repmat(x, 1, M);
sth.P = repmat(P, [1 1 M]);
sth.mu = m.mu0;
sth.t = b(1);  sth.tOpp = b(1);
sth.hits = 1;
sth.plot = b(5);
sth.hist = b(5);
end

function s = empty_sth()
% A single-target hypothesis without a target (r = 0).
s = struct('r', 0, 'x', [], 'P', [], 'mu', [], 't', NaN, 'tOpp', NaN, 'tNext', Inf, 'hits', 0, ...
           'misses', 0, 'plot', 0, 'batch', 0, 'hist', zeros(1, 0));
end

function tn = next_look(sth, p)
% Time of the next look of the antenna at the STH after its last one: the
% antenna points at the predicted azimuth between half a scan and 1.5 scans later.
T = p.scanTime_s;
xm = imm_estimate(sth);
d = p.dim;
t0 = sth.tOpp + T/2;
antAz = p.az0_deg + p.rate_dps * t0;
tn = t0 + T/2;
for it = 1:4
    q = xm(1:2) + xm(d+1:d+2) * (tn - sth.t);
    az = atan2(q(1), q(2)) * 180/pi;
    tn = t0 + mod(az - antAz, 360) / p.rate_dps;
end
end

function ok = coarse_gate(xm, Pm, t, z, tp, p)
% Quick test before the UKF: the plot lies within a wide circle around the
% position predicted with constant velocity (and up to 50 m/s^2).
d = p.dim;
dt = abs(tp - t);
q = xm(1:2) + xm(d+1:d+2) * (tp - t);
if d == 3
    rh = z(1) * cos(z(3));
else
    rh = z(1);
end
zp = rh * [sin(z(2)); cos(z(2))];
s2 = trace(Pm(1:2, 1:2)) + dt^2 * trace(Pm(d+1:d+2, d+1:d+2)) + dt^4/4 * trace(Pm(2*d+1:2*d+2, 2*d+1:2*d+2));
sm = sqrt(p.R(1, 1)) + rh * sqrt(p.R(2, 2));
ok = norm(zp - q) <= 6 * sqrt(s2) + 25 * dt^2 + 6 * sm + 1000;
end

%% ---------------------------------------------------------------- assignment
function [A, cost] = murty(C, k)
% The k best assignments of the rows of C (n x m, n <= m, Inf = forbidden)
% to distinct columns (Murty's algorithm). A(a, i) = column of row i in the
% a-th best assignment, cost(a) its total cost.
n = size(C, 1);
A = zeros(0, n);  cost = zeros(0, 1);
[a0, c0] = lap(C);
if isempty(a0), return; end
Q = struct('C', {C}, 'a', {a0}, 'c', {c0});
while ~isempty(Q) && size(A, 1) < k
    [~, b] = min([Q.c]);
    node = Q(b);  Q(b) = [];
    A(end+1, :) = node.a; %#ok<AGROW>
    cost(end+1, 1) = node.c; %#ok<AGROW>
    Cw = node.C;
    for r = 1:n
        Cs = Cw;
        Cs(r, node.a(r)) = Inf;                     % without the pair (r, a(r))
        [as, cs] = lap(Cs);
        if ~isempty(as)
            Q(end+1) = struct('C', Cs, 'a', as, 'c', cs); %#ok<AGROW>
        end
        keep = Cw(r, node.a(r));                    % following subproblems keep (r, a(r))
        Cw(r, :) = Inf;  Cw(:, node.a(r)) = Inf;
        Cw(r, node.a(r)) = keep;
    end
end
end

function [a, c] = lap(C)
% Minimum-cost assignment of every row of C (n <= m) to its own column;
% empty when no assignment avoids the forbidden (Inf) pairs.
[n, m] = size(C);
ok = isfinite(C);
a = [];  c = Inf;
if n == 0, a = zeros(1, 0); c = 0; return; end
if ~all(any(ok, 2)), return; end
F = C(ok);
A = C - min(F);
A(~ok) = (max(F) - min(F) + 1) * (n + 1);
rowOf = hungarian(A);
a = zeros(1, n);
for j = 1:m
    if rowOf(j) > 0, a(rowOf(j)) = j; end
end
idx = sub2ind([n m], 1:n, a);
if any(~ok(idx)), a = []; return; end
c = sum(C(idx));
end

function rowOf = hungarian(A)
% Minimum-cost assignment of every row of A (n x m, n <= m) to its own column
% (shortest augmenting paths with potentials, as in track_ekf.m).
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
        cur = A(i0, :).' - u(i0+1) - v(2:end);
        upd = ~used(2:end) & cur < minv(2:end);
        minv([false; upd]) = cur(upd);
        way([false; upd]) = j0;
        cand = minv(2:end);
        cand(used(2:end)) = Inf;
        [delta, j1] = min(cand);
        uj = find(used);
        u(p(uj)+1) = u(p(uj)+1) + delta;
        v(uj) = v(uj) - delta;
        minv(~used) = minv(~used) - delta;
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

%% ---------------------------------------------------------------- helpers
function [z, ok, zvr] = measurement(b, dim)
% Measurement vector of a plot of the batch ([time az range el plotIndex vr])
% and its radial velocity (2-D: of the ground range when the elevation is known).
az = deg2rad(b(2));
el = deg2rad(b(4));
zvr = b(6);
if dim == 3
    z = [b(3); az; el];
    ok = isfinite(el);
else
    rho = b(3);
    if isfinite(el), rho = rho * cos(el);  zvr = zvr / cos(el); end
    z = [rho; az];
    ok = true;
end
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

function s = logsumexp(x)
mx = max(x(:));
if ~isfinite(mx), s = mx; return; end
s = mx + log(sum(exp(x(:) - mx)));
end

function s = logsumexp2(a, b)
% Element-wise log(exp(a) + exp(b)).
mx = max(a, b);
s = mx + log(exp(a - mx) + exp(b - mx));
s(~isfinite(mx)) = mx(~isfinite(mx));
end

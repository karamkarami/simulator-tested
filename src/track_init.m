function tr = track_init(params)
%TRACK_INIT  Empty EKF tracker (see track_update.m).
tr.p = params;
tr.tracks = struct('id', {}, 'confirmed', {}, 'x', {}, 'P', {}, 't', {}, 'hits', {}, 'misses', {});
tr.nextId = 1;
tr.log = zeros(0, 11);   % one row per update: [time id confirmed x y z vx vy vz hits plotIndex]
tr.nPlots = 0;
tr.pending = zeros(0, 5);      % GNN batch: [time az range el plotIndex] of the plots not yet associated
tr.pendingFirst = -Inf;        % arrival time of the first and of the newest plot of the batch
tr.pendingLast = -Inf;
end

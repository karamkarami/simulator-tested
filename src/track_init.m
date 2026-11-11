function tr = track_init(params)
%TRACK_INIT  Empty tracker (see track_update.m): EKF, or PMBM (params.method = 'pmbm', track_pmbm.m).
tr.p = params;
tr.tracks = struct('id', {}, 'confirmed', {}, 'x', {}, 'P', {}, 't', {}, 'hits', {}, 'misses', {});
tr.nextId = 1;
tr.log = zeros(0, 11);   % one row per update: [time id confirmed x y z vx vy vz hits plotIndex]
tr.nPlots = 0;
tr.pending = zeros(0, 7);      % GNN batch: [time az range el plotIndex vr vrSpan] of the plots not yet associated
tr.pendingFirst = -Inf;        % arrival time of the first and of the newest plot of the batch
tr.pendingLast = -Inf;
if isfield(params, 'method') && strcmp(params.method, 'pmbm')
    tr.tracks = struct('id', {}, 'confirmed', {}, 'x', {}, 'P', {}, 't', {}, 'hits', {}, 'misses', {}, ...
                       'r', {}, 'mu', {});
    tr.T = struct('id', {}, 'confirmed', {}, 'h', {}, 'rMarg', {}, 'rTarget', {}, 'plot', {});   % Bernoulli components
    tr.G = zeros(1, 0);        % global hypotheses: STH of every track (one row each)
    tr.logw = 0;               % their log weights
    tr.batch = 0;              % batches associated
    tr.tMaint = -Inf;          % last check of the missed looks
    tr.logged = 0;             % newest plot in tr.log
    tr.shown = struct('id', {}, 'plot', {});   % confirmed tracks shown: number and newest plot
end
if isfield(params, 'method') && strcmp(params.method, 'tomht')
    tr.tracks = struct('id', {}, 'confirmed', {}, 'x', {}, 'P', {}, 't', {}, 'hits', {}, 'misses', {}, ...
                       'r', {}, 'mu', {}, 'score', {});
    tr.T = struct('id', {}, 'confirmed', {}, 'h', {}, 'rMarg', {}, 'rTarget', {}, 'plot', {}, 't0', {});   % track trees
    tr.G = zeros(1, 0);        % global hypotheses: the leaf of every tree (one row each)
    tr.logw = 0;               % their scores (sum of the LLR of their tracks)
    tr.batch = 0;              % batches associated
    tr.tMaint = -Inf;          % last check of the missed looks
    tr.tNow = -Inf;            % newest plot or look time (N-scan pruning)
    tr.logged = 0;             % newest plot in tr.log
    tr.shown = struct('id', {}, 'plot', {});   % confirmed tracks shown: number and newest plot
end
end

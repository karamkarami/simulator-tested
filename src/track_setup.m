function p = track_setup(cfg, sys)
%TRACK_SETUP  Parameters of the tracker (track_update.m) from cfg.track.
%   p.method is 'ekf' (EKF, GNN) or 'pmbm' (PMBM filter, IMM-UKF: p.pmbm,
%   track_pmbm.m); the settings below are common to both.
%   p.dim is 3 (range, azimuth, elevation measured; needs the elevation of
%   the DBF) or 2 (range and azimuth in the horizontal plane). With
%   cfg.track.dim = 'auto' it is 3 with beamforming and 2 without.
%   p.model is 'cv' (state [position; velocity], p.nd = 2 values per axis)
%   or 'ca' ([position; velocity; acceleration], p.nd = 3); p.q holds the
%   white-acceleration ('cv') or white-jerk ('ca') spectral density per axis.
%   Radial velocity (cfg.track.useRadialVelocity): p.useVr, its variance
%   p.Rvr, the span of the Doppler axis p.vrSpan (lambda/(2*row PRI)),
%   p.vrUnambiguous (maxSpeed below half the span: the plot velocity is the
%   true one) and the chi-square gate with the radial velocity p.gateVr.
%   Speed gate: p.minSpeed_mps (cfg.track.minSpeed_mps with
%   cfg.track.minSpeedEnable, else 0 = off); a track with at least
%   p.minSpeedHits updates that is slower is deleted (track_update.m).
t = cfg.track;
p.enable = logical(t.enable);
p.method = 'ekf';
if isfield(t, 'method') && ~isempty(t.method), p.method = lower(t.method); end
if ~any(strcmp(p.method, {'ekf', 'pmbm'}))
    error('cfg.track.method must be ''ekf'' or ''pmbm''.');
end
switch lower(t.dim)
    case 'auto', p.dim = 2 + sys.dbf.enabled;
    case '2d',   p.dim = 2;
    case '3d'
        p.dim = 3;
        if ~sys.dbf.enabled
            warning('cfg.track.dim = ''3d'' needs the elevation of the DBF (1 element: no DBF): tracking in 2-D.');
            p.dim = 2;
        end
    otherwise
        error('cfg.track.dim must be ''auto'', ''2d'' or ''3d''.');
end
p.model = lower(t.model);
switch p.model
    case 'cv'
        p.nd = 2;
        q = [t.sigmaAccel_mps2, t.sigmaAccel_mps2, t.sigmaAccelZ_mps2].^2;   % [m^2/s^3]
    case 'ca'
        p.nd = 3;
        q = [t.sigmaJerk_mps3, t.sigmaJerk_mps3, t.sigmaJerkZ_mps3].^2;      % [m^2/s^5]
    otherwise
        error('cfg.track.model must be ''cv'' or ''ca''.');
end
sig = [t.sigmaRange_m, deg2rad(t.sigmaAz_deg), deg2rad(t.sigmaEl_deg)];
p.R = diag(sig(1:p.dim).^2);                       % measurement noise [range az (el)]
p.q = q(1:p.dim);                                  % process noise spectral density per axis
p.sigmaV0 = t.maxSpeed_mps / sqrt(3);              % velocity std of a new track (uniform in +-maxSpeed)
p.sigmaA0 = t.maxAccel_mps2 / sqrt(3);             % acceleration std of a new track ('ca')
p.gate = 2 * gammaincinv(t.gateProb, p.dim/2);     % chi-square gate with dim degrees of freedom
p.useVr = logical(t.useRadialVelocity);
sigmaVr = t.sigmaVr_mps;
if isempty(sigmaVr), sigmaVr = 0.4 * abs(diff(sys.dop.v_mps(1:2))); end   % 0.4 Doppler bins
p.Rvr = sigmaVr^2;
p.vrSpan = sys.lambda / (2 * sys.dop.rowPri_s);   % the plot velocity is known modulo this
p.vrUnambiguous = t.maxSpeed_mps < p.vrSpan / 2;
p.gateVr = 2 * gammaincinv(t.gateProb, (p.dim + 1)/2);
p.confirmHits = t.confirmHits;
p.maxMissTentative = t.maxMissTentative;
p.maxMiss = t.maxMiss;
p.minSpeed_mps = 0;                                % 0: no deletion of slow tracks
if isfield(t, 'minSpeedEnable') && t.minSpeedEnable
    if ~(isscalar(t.minSpeed_mps) && t.minSpeed_mps > 0 && isfinite(t.minSpeed_mps))
        error('cfg.track.minSpeed_mps must be a positive speed [m/s] when cfg.track.minSpeedEnable is true.');
    end
    p.minSpeed_mps = t.minSpeed_mps;
end
p.minSpeedHits = 2;                                % updates before the speed is tested
p.scanTime_s = 360 / sys.ant.rate_dps;
p.batchGap_s = t.gnnSector_deg / sys.ant.rate_dps;     % GNN batch: time without a new plot
p.batchMax_s = t.gnnMaxSector_deg / sys.ant.rate_dps;  % GNN batch: longest collection time
p.height_m = sys.ant.height_m;                     % the state is relative to the antenna
p.az0_deg = cfg.ant.az0_deg;                       % antenna azimuth at t = 0 (look times, PMBM)
p.rate_dps = sys.ant.rate_dps;
if strcmp(p.method, 'pmbm')
    p.pmbm = pmbm_setup(t, sys, p);
end
end

function m = pmbm_setup(t, sys, p)
% Parameters of the PMBM tracker with IMM-UKF state estimation (track_pmbm.m).
s = t.pmbm;
m = struct();
m.pd = s.pd;  m.ps = s.ps;
if ~(m.pd > 0 && m.pd < 1) || ~(m.ps > 0 && m.ps <= 1)
    error('cfg.track.pmbm: pd must lie in (0, 1) and ps in (0, 1].');
end
if ~(s.falsePlotsPerScan > 0) || ~(s.undetectedTargets > 0)
    error('cfg.track.pmbm: falsePlotsPerScan and undetectedTargets must be positive.');
end
% measurement space of a plot: range x azimuth (x elevation) (x radial velocity); the
% false plots and the undetected targets are uniform in it
rCover = [sys.wave.short.range_m(1), sys.wave.long.range_m(end)];
vol = diff(rCover) * 2*pi;                         % [m rad]
if p.dim == 3
    el = sys.dbf.beamEl_deg(:);
    de = 1;  if numel(el) > 1, de = abs(diff(el(end-1:end))); end
    vol = vol * deg2rad(max(el) - min(el) + de);   % [m rad^2]
end
m.logLambdaC = log(s.falsePlotsPerScan / vol);     % false plots per scan and unit volume
m.logLambdaU = log(s.undetectedTargets / vol);     % undetected targets per unit volume
m.logVrDensity = -log(p.vrSpan);                   % uniform radial velocity of the plot (when used)
m.maxHyp = s.maxHypotheses;
m.logHypPrune = log(s.hypPrune);
m.rPrune = s.rPrune;  m.rConfirm = s.rConfirm;  m.rDelete = s.rDelete;
if ~(m.maxHyp >= 1) || ~(m.rPrune < m.rDelete && m.rDelete < m.rConfirm && m.rConfirm <= 1)
    error('cfg.track.pmbm: maxHypotheses >= 1 and rPrune < rDelete < rConfirm <= 1 needed.');
end
% IMM: models and the Markov chain per scan (uniform switching)
names = lower(s.imm.models);
if ischar(names), names = {names}; end
if isempty(names) || ~all(ismember(names, {'cv', 'ca', 'ct'})) || numel(unique(names)) < numel(names)
    error('cfg.track.pmbm.imm.models must be distinct models of ''cv'', ''ca'', ''ct''.');
end
M = numel(names);
m.models = names;
mu0 = s.imm.mu0(:).';
if numel(mu0) ~= M || any(mu0 < 0) || sum(mu0) <= 0
    error('cfg.track.pmbm.imm.mu0 needs one non-negative probability per model.');
end
m.mu0 = mu0 / sum(mu0);
sw = s.imm.switchProb;
if M == 1, sw = 0; end
if ~(sw >= 0 && sw < (M - 1)/M)
    error('cfg.track.pmbm.imm.switchProb must lie in [0, (M-1)/M) for M models.');
end
m.switchProb = sw;
m.switchLambda = 1 - sw*M/(M - 1 + (M == 1));     % second eigenvalue of the transition matrix per scan
d = p.dim;
m.n = 3*d + 1;                                     % state [position; velocity; acceleration; turn rate]
qz = [t.sigmaAccelZ_mps2^2, t.sigmaJerkZ_mps3^2];  % vertical (3-D)
m.qCv = [s.imm.sigmaAccelCv_mps2^2 * [1 1], qz(1)];  m.qCv = m.qCv(1:d);   % [m^2/s^3]
m.qCa = [s.imm.sigmaJerkCa_mps3^2 * [1 1], qz(2)];   m.qCa = m.qCa(1:d);   % [m^2/s^5]
m.qCt = [s.imm.sigmaAccelCt_mps2^2 * [1 1], qz(1)];  m.qCt = m.qCt(1:d);   % [m^2/s^3]
m.qTurn = deg2rad(s.imm.sigmaTurn_dps)^2;          % [rad^2/s^3]
m.sigmaTurn0 = deg2rad(s.imm.maxTurn_dps) / sqrt(3);
m.sigmaClimb0 = min(t.maxSpeed_mps, s.maxClimb_mps) / sqrt(3);   % vertical velocity std of a new track (3-D)
% A model without acceleration (CV, CT) or turn rate (CV, CA) predicts it as 0 with the
% uncertainty of a new track: it does not change the motion of that model, but when the
% IMM mixes it into the CA or CT model that model can take up a new maneuver at once.
m.varUnusedA = p.sigmaA0^2;                        % [m^2/s^4]
m.varUnusedW = m.sigmaTurn0^2;                     % [rad^2/s^2]
% unscented transform: weights of the 2n+1 sigma points
n = m.n;
a = s.ukf.alpha;  b = s.ukf.beta;  k = s.ukf.kappa;
lam = a^2 * (n + k) - n;
if ~(n + lam > 0), error('cfg.track.pmbm.ukf: alpha^2*(n + kappa) must be positive (n = %d).', n); end
m.ukfScale = sqrt(n + lam);
m.Wm = [lam/(n + lam), repmat(1/(2*(n + lam)), 1, 2*n)];
m.Wc = m.Wm;  m.Wc(1) = m.Wc(1) + 1 - a^2 + b;
% chi-square gates of the measurement with and without the radial velocity (as the EKF)
m.gate = p.gate;  m.gateVr = p.gateVr;
m.latency_s = t.gnnMaxSector_deg / sys.ant.rate_dps + 2 * sys.cpi.time_s;   % a plot is in the tracker this
                                                   % long after the antenna passed it (extraction, batch)
m.maintStep_s = 1 / sys.ant.rate_dps;             % missed looks are checked every degree of rotation
end

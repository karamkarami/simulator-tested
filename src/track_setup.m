function p = track_setup(cfg, sys)
%TRACK_SETUP  Parameters of the EKF tracker (track_update.m) from cfg.track.
%   p.dim is 3 (range, azimuth, elevation measured; needs the elevation of
%   the DBF) or 2 (range and azimuth in the horizontal plane). With
%   cfg.track.dim = 'auto' it is 3 with beamforming and 2 without.
%   p.model is 'cv' (state [position; velocity], p.nd = 2 values per axis)
%   or 'ca' ([position; velocity; acceleration], p.nd = 3); p.q holds the
%   white-acceleration ('cv') or white-jerk ('ca') spectral density per axis.
t = cfg.track;
p.enable = logical(t.enable);
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
p.confirmHits = t.confirmHits;
p.maxMissTentative = t.maxMissTentative;
p.maxMiss = t.maxMiss;
p.scanTime_s = 360 / sys.ant.rate_dps;
p.batchGap_s = t.gnnSector_deg / sys.ant.rate_dps;     % GNN batch: time without a new plot
p.batchMax_s = t.gnnMaxSector_deg / sys.ant.rate_dps;  % GNN batch: longest collection time
p.height_m = sys.ant.height_m;                     % the state is relative to the antenna
end

function tr = track_update(tr, plots, t_s, flush)
%TRACK_UPDATE  Tracking block: the new plots of time t_s go to the tracker of cfg.track.method.
%   tr = TRACK_UPDATE(tr, plots, t_s) takes the new plots of time t_s (call
%   it every PRI, also without plots); tr = TRACK_UPDATE(tr, plots, t_s, true)
%   also associates the plots still waiting in the batch (end of the data,
%   skipped sector). tr comes from track_init(sys.track) (track_setup.m):
%     'ekf'   track_ekf.m   extended Kalman filter per track, GNN association
%     'pmbm'  track_pmbm.m  PMBM multi-target filter, IMM-UKF state estimation
%   Both keep tr.tracks (the current tracks, track_report.m) and tr.log (one
%   row per track update: [time id confirmed x y z vx vy vz hits plotIndex]).
if nargin < 4, flush = false; end
switch tr.p.method
    case 'pmbm', tr = track_pmbm(tr, plots, t_s, flush);
    otherwise,   tr = track_ekf(tr, plots, t_s, flush);
end
end

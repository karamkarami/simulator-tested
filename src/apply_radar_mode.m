function cfg = apply_radar_mode(cfg)
%APPLY_RADAR_MODE  Apply the radar mode cfg.modes.active and its sub-mode cfg.modes.subMode.
%
%   cfg = APPLY_RADAR_MODE(cfg) runs the mode file radar_config_mode_<k>.m
%   (folder modes/, k = cfg.modes.active), which sets the waveforms of both
%   sub-pulses, cfg.wave.singlePulse and the PRF of the mode, and then the
%   sub-mode: 'dual' keeps both sub-pulses, 'single' switches off the
%   sub-pulse other than cfg.wave.singlePulse. radar_config.m calls it; to run
%   another mode in a script:
%
%       cfg = radar_config();
%       cfg.modes.active  = 2;          % one of cfg.modes.available
%       cfg.modes.subMode = 'single';
%       cfg = apply_radar_mode(cfg);
%
%   setup_radar.m checks cfg.modes (at most 6 modes 0 ... 5, the active one
%   among cfg.modes.available, the sub-mode). This file is also evaluated by
%   radar_export_config (no if/for blocks: only assignments).
cfg = feval(sprintf('radar_config_mode_%d', cfg.modes.active), cfg);
isSingle = strcmpi(cfg.modes.subMode, 'single');
cfg.wave.short.enable = cfg.wave.short.enable && (~isSingle || strcmpi(cfg.wave.singlePulse, 'short'));
cfg.wave.long.enable  = cfg.wave.long.enable  && (~isSingle || strcmpi(cfg.wave.singlePulse, 'long'));
end

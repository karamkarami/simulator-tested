function cfg = apply_radar_mode(cfg)
%APPLY_RADAR_MODE  Apply the radar mode cfg.modes.active.
%
%   cfg = APPLY_RADAR_MODE(cfg) runs the mode file radar_config_mode_<k>.m
%   (folder modes/, k = cfg.modes.active), which sets the waveforms of both
%   sub-pulses (each one on or off with cfg.wave.<pulse>.enable: a mode with
%   only the short or only the long pulse switches the other one off), the PRF
%   and the whole processing of the mode. radar_config.m calls it; to run
%   another mode in a script:
%
%       cfg = radar_config();
%       cfg.modes.active = 2;           % one of cfg.modes.available
%       cfg = apply_radar_mode(cfg);
%
%   setup_radar.m checks cfg.modes (at most 8 modes 0 ... 7, the active one
%   among cfg.modes.available). This file is also evaluated by
%   radar_export_config (no if/for blocks: only assignments).
cfg = feval(sprintf('radar_config_mode_%d', cfg.modes.active), cfg);
end

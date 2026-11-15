function cfg = lfm_mode_config()
%LFM_MODE_CONFIG  radar_config() in radar mode 1 'lfm_mode': both sub-pulses LFM, fixed PRF.
%   cfg = LFM_MODE_CONFIG() is the configuration the tests were written and
%   checked with (the LFM waveforms of the radar before the radar modes; the
%   default mode 0 'medium_mode' has phase codes). Tests that compare MATLAB
%   with the C tools reading radar_config.m keep radar_config() instead.
cfg = radar_config();
cfg.modes.active = 1;
if ~any(cfg.modes.available == 1), cfg.modes.available = unique([cfg.modes.available(:).', 1]); end
cfg = apply_radar_mode(cfg);
end

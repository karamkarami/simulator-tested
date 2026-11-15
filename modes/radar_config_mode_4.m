function cfg = radar_config_mode_4(cfg)
%RADAR_CONFIG_MODE_4  Settings of radar mode 4 (see "Radar modes" in radar_config.m).
%
%   Mode 4 'mode_4': not yet defined for a program; it holds the
%   waveforms of mode 1 (LFM, fixed PRF). Edit the values below to define it.
%
%   cfg = RADAR_CONFIG_MODE_4(cfg) sets the waveforms of both sub-pulses, the
%   sub-pulse of the single-pulse sub-mode and the PRF of this mode in cfg.
%   radar_config.m calls it (through apply_radar_mode.m) when this mode is
%   cfg.modes.active. Every mode file sets the same fields (all of them), so
%   that selecting another mode with apply_radar_mode replaces them all. A mode
%   may set further fields of cfg (then set them in every mode file too).
%   Each mode has two sub-modes (cfg.modes.subMode): 'dual' (both sub-pulses,
%   as set here) and 'single' (only the sub-pulse cfg.wave.singlePulse).
%   The meaning of every field: section "Waveforms" of README.md and
%   radar_config.m (cfg.wave.*, cfg.prf.*).

cfg.modes.name = 'mode_4';                % name of the mode (run summary, receiver configuration)

% Sub-mode 'single' (one sub-pulse per PRI): the sub-pulse it keeps, 'short' or 'long'
cfg.wave.singlePulse = 'long';

% Waveforms: every PRI holds a short and a long sub-pulse (short first)
cfg.wave.short.enable       = true;     % short sub-pulse on/off (in the 'dual' sub-mode)
cfg.wave.short.fc_Hz        = [];       % its carrier [Hz] ([] = cfg.fc)
cfg.wave.short.type         = 'lfm';        % 'lfm' or 'phasecode' (code from cfg.files.phaseCode)
cfg.wave.short.width_s      = 25e-6;    % pulse width [s]
cfg.wave.short.bandwidth_Hz = 1e6;      % LFM sweep bandwidth [Hz]
cfg.wave.short.chip_s       = 1e-6;     % phase-code chip length [s] (width_s/chip_s chips)
cfg.wave.short.compression  = 'matched';   % phase code: 'matched' or 'mismatch' filter (LFM: 'matched')
cfg.wave.short.mismatchLength = [];     % coefficients (one per chip) of the mismatch filter in the default
                                        % phase-code file ([] = chips + 10, at most 2x the chips)
cfg.wave.short.tTx_s        = 0;        % transmit time inside the PRI [s]

cfg.wave.long.enable        = true;     % long sub-pulse on/off (in the 'dual' sub-mode)
cfg.wave.long.fc_Hz         = [];
cfg.wave.long.type          = 'lfm';
cfg.wave.long.width_s       = 150e-6;
cfg.wave.long.bandwidth_Hz  = 1e6;
cfg.wave.long.chip_s        = 1e-6;
cfg.wave.long.compression   = 'matched';
cfg.wave.long.mismatchLength = [];
cfg.wave.long.tTx_s         = 200e-6;

% PRF: 'fixed' (pri_s), 'staggered' (PRIs of cfg.files.stagger) or 'scan2scan' (one PRF per
% antenna scan from cfg.files.scan2scan, constant during the scan)
cfg.prf.mode  = 'fixed';
cfg.prf.pri_s = 1e-3;      % PRI in fixed mode [s]
end

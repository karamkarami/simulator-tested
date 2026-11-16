function cfg = radar_config_mode_4(cfg)
%RADAR_CONFIG_MODE_4  Settings of radar mode 4 (see "Radar modes" in radar_config.m).
%
%   Mode 4 'mode_4': not yet defined for a program; it holds the settings
%   of mode 1 (LFM, fixed PRF). Edit the values below to define it.
%
%   cfg = RADAR_CONFIG_MODE_4(cfg) sets everything of this mode in cfg: the
%   waveforms, transmit times and receive windows of both sub-pulses (each
%   one on or off), the PRF, the antenna rotation speed, the input files of the
%   processing (stagger, scan-to-scan PRFs, phase codes, FIR clutter filter),
%   the STC and the whole receiver processing chain (pulse compression, CPI
%   buffer, Doppler window and FFT, clutter filter, CFAR, pruning, reports,
%   plot extraction, tracking). Common to all modes (radar_config.m): the
%   antenna (except its rotation speed), transmitter, receiver gain, ADC, DDC, clutter, targets, recording
%   and display. radar_config.m calls this file (through apply_radar_mode.m)
%   when this mode is cfg.modes.active. Every mode file sets the same fields
%   (all of them), so that selecting another mode with apply_radar_mode
%   replaces them all. A mode with one sub-pulse per PRI switches the other
%   one off (cfg.wave.<pulse>.enable = false).
%   The meaning of every setting: radar_config_mode_0.m and README.md.
%   The file is also evaluated by radar_export_config (C receiver): only
%   assignments, no if/for blocks.

root = fileparts(fileparts(mfilename('fullpath')));   % folder of radar_config.m

cfg.modes.name = 'mode_4';              % name of the mode (run summary, receiver configuration)

%% Waveforms and timing: every PRI holds a short and a long sub-pulse (short first)
% Times in the PRI [s] from its start (the first ADC/DDC sample): transmission at tTx_s, receive
% window [rxStart_s rxEnd_s) ([] = from the end of the pulse; short pulse: to the transmission of
% the long pulse, long pulse: to the end of the echoes from cfg.rx.maxRange_m). The DDC output
% runs from the PRI start to the end of the last receive window; the receiver takes only the
% receive windows out of it (the samples while transmitting and between the windows are not used).
% Range cells of a sub-pulse: the echoes that lie whole inside its window, i.e. the delays
% rxStart_s - tTx_s ... rxEnd_s - tTx_s - width_s.
cfg.wave.short.enable       = true;     % short sub-pulse on/off (false: only the long pulse)
cfg.wave.short.fc_Hz        = [];       % its carrier [Hz] ([] = cfg.fc)
cfg.wave.short.type         = 'lfm';    % 'lfm' or 'phasecode' (code from cfg.files.phaseCode)
cfg.wave.short.width_s      = 25e-6;    % pulse width [s]
cfg.wave.short.bandwidth_Hz = 1e6;      % LFM sweep bandwidth [Hz]
cfg.wave.short.chip_s       = 1e-6;     % phase-code chip length [s] (width_s/chip_s chips)
cfg.wave.short.compression  = 'matched'; % phase code: 'matched' or 'mismatch' filter (LFM: 'matched')
cfg.wave.short.availableCodes = [1 2 3 4]; % phase code: pairs of cfg.files.phaseCode this mode may use
cfg.wave.short.codeIndex    = 1;        % phase code: code / mismatch pair of cfg.files.phaseCode (one of availableCodes)
cfg.wave.short.mismatchLength = [];     % coefficients (one per chip) of the mismatch filter in the default
                                        % phase-code file ([] = chips + 10, at most 2x the chips)
cfg.wave.short.tTx_s        = 0;        % transmit time in the PRI [s]
cfg.wave.short.rxStart_s    = [];       % start of its receive window in the PRI [s] ([] = end of the pulse)
cfg.wave.short.rxEnd_s      = [];       % end of its receive window [s] ([] = transmission of the long pulse)

cfg.wave.long.enable        = true;     % long sub-pulse on/off (false: only the short pulse)
cfg.wave.long.fc_Hz         = [];
cfg.wave.long.type          = 'lfm';
cfg.wave.long.width_s       = 150e-6;
cfg.wave.long.bandwidth_Hz  = 1e6;
cfg.wave.long.chip_s        = 1e-6;
cfg.wave.long.compression   = 'matched';
cfg.wave.long.availableCodes = [1 2 3 4];
cfg.wave.long.codeIndex     = 1;
cfg.wave.long.mismatchLength = [];
cfg.wave.long.tTx_s         = 200e-6;   % transmit time in the PRI [s]
cfg.wave.long.rxStart_s     = [];       % start of its receive window [s] ([] = end of the pulse)
cfg.wave.long.rxEnd_s       = [];       % end of its receive window = end of the DDC record [s]
                                        % ([] = end of the echoes from cfg.rx.maxRange_m)
cfg.rx.maxRange_m      = 80e3;   % instrumented range [m]: the end of the long-pulse receive window when its
                                 % rxEnd_s is [] (else set it to match the window); the PPI range
cfg.rx.manualDdc = false;        % true: the DDC record of this mode is set by hand (nDdc, ddcT0_s)
cfg.rx.nDdc      = [];           % manualDdc: DDC samples per PRI
cfg.rx.ddcT0_s   = [];           % manualDdc: time of the first DDC sample in the PRI [s] ([] = automatic)

%% PRF: 'fixed' (pri_s), 'staggered' (PRIs of cfg.files.stagger, with their MTI weights) or
% 'scan2scan' (one PRF per antenna scan from cfg.files.scan2scan, constant during the scan)
cfg.prf.mode  = 'fixed';
cfg.prf.pri_s = 1e-3;                   % PRI in fixed mode [s] (at least the end of the last receive window)
cfg.prf.availableStaggerSets = [1 2 3 4]; % staggered PRF: (T, w) pairs of cfg.files.stagger this mode may use
cfg.prf.staggerSet           = 1;         % staggered PRF: pair in use (one of availableStaggerSets)

%% Antenna rotation of this mode
cfg.ant.availableRpm = [6 12];          % rotation speeds this mode may use [rev/min]
cfg.ant.rpm          = 6;               % rotation speed in use [rev/min] (one of availableRpm; [] = the first)

%% Input files of this mode (each mode may have its own files)
cfg.files.stagger   = fullfile(root, 'inputs', 'stagger_prf.mat');    % staggered PRIs + MTI weights (cfg.prf.mode = 'staggered')
cfg.files.scan2scan = fullfile(root, 'inputs', 'scan2scan_prf.mat');  % PRFs of the scans (cfg.prf.mode = 'scan2scan')
cfg.files.phaseCode = fullfile(root, 'inputs', 'phase_code.mat');     % phase codes + mismatch filters of both sub-pulses
cfg.files.mtiFir    = fullfile(root, 'inputs', 'mti_fir.mat');        % coefficients of the FIR clutter filter (cfg.mti.type = 'fir')

%% STC (sensitivity time control)
cfg.rx.stc.enable = true;        % true: STC on; false: STC off (no attenuation). Without the STC lower
                                 % cfg.rx.gain_dB to 30 dB, otherwise the near clutter clips the ADC
                                 % (6 - 9 % of the receive samples at 60 dB: false plots; run_radar_sim warns)
cfg.rx.stc.curve  = [    0 40;  1000 40;  2000 32;  4000 25;  6000 22;  8000 16;  9000 12;  10000 5;  10500 0];
cfg.rx.stc.curveLong = [0 60;  3000 60;  6000 45;  9000 30;  12000 18;  15000 10;  18000 5;  22500 0];

%% Pulse compression and Doppler processing (CPI buffer)
cfg.proc.taylorNbar   = 4;     % Taylor window (LFM only): number of nearly constant sidelobes
cfg.proc.taylorSll_dB = -35;   % Taylor window peak sidelobe level [dB]
cfg.proc.cpiMode      = 'auto'; % pulses per CPI buffer: 'auto' from the radar parameters (the pulses
                               % arriving while a target is in the 3 dB azimuth beam:
                               % round(beamwidth / (rotation rate * PRI)) = 114), or 'manual':
                               % cfg.proc.cpiLength below
cfg.proc.cpiLength    = 114;   % pulses per CPI buffer with cpiMode = 'manual' (ignored with 'auto'; with a
                               % staggered PRF rounded down to whole stagger periods)
cfg.proc.nfft         = [];    % Doppler FFT length ([] = next power of 2: of the pulses in the buffer
                               % with a fixed PRF (114 -> 128), of the samples per row when staggered (28 -> 32))
cfg.proc.dopplerWindow = 'hamming'; % window of the Doppler FFT: 'hamming', 'hann', 'blackman' or 'rect'
                               % (mainlobe +-2, +-2, +-3, +-1 bins x nfft/samples: sets the MTI notch width)
cfg.proc.updateEvery  = 1;     % run MTI/Doppler/CFAR every N new pulses (1 = every pulse)

%% Clutter filter along slow time (MTI), in front of the Doppler FFT
cfg.mti.type = 'canceller';    % 'canceller' or 'fir'
cfg.mti.maxRange_m = Inf;      % the clutter filter is applied to the range cells up to this range [m]
                               % (both sub-pulses); the cells beyond it go to the Doppler FFT unfiltered
                               % (the newest L pulses' last one, so every cell has the same CPI), with
                               % their own noise-gain equalisation. Inf = every range cell (as before).
                               % Use cfg.cfar.notchMaxRange_m to test the zero-Doppler bins beyond it.

%% Detection: CFAR on every range-Doppler map
cfg.cfar.type   = 'ca';   % 'ca', 'go', 'so' or 'os'
cfg.cfar.method = '2d';   % '2d', 'range' or 'doppler'
cfg.cfar.maxPfa = 1e-4;   % largest Pfa allowed (per cell and CPI): the lowest threshold that can be set
                          % (pfa / threshold_dB of the window and the PPI threshold field, MATLAB and C)
cfg.cfar.osRank = 0.9;    % OS: rank of the statistic as a fraction of the training cells (0 ... 1]. The
                          % textbook 0.75 suits isolated interfering targets; here the eclipsed long-pulse
                          % echo of a near target (short-pulse range, below 22.5 km) smears over the next
                          % 15 km of the long window as a ripple of similar peaks, whose 75th percentile is
                          % below its mean: with 0.75 the OS detects tens of its peaks per scan
                          % (test_multi_scan: 62 false plots in 3 scans, 0.9: 1, as CA; 0.95: 0)

cfg.cfar.twoD.guardRange   = 3;     % guard cells on each side in range
cfg.cfar.twoD.trainRange   = 8;     % training cells on each side in range
cfg.cfar.twoD.guardDoppler = 3;     % guard cells on each side in Doppler
cfg.cfar.twoD.trainDoppler = 4;     % training cells on each side in Doppler
cfg.cfar.twoD.pfa          = 1e-8;  % design probability of false alarm (per cell and CPI)
cfg.cfar.twoD.threshold_dB = [];    % [] = from pfa; else the threshold factor in dB

cfg.cfar.range.guard        = 3;    % guard cells on each side
cfg.cfar.range.train        = 16;   % training cells on each side
cfg.cfar.range.pfa          = 1e-8;
cfg.cfar.range.threshold_dB = [];

cfg.cfar.doppler.guard        = 2;  % guard bins on each side
cfg.cfar.doppler.train        = 8;  % training bins on each side (2*(guard+train)+1 <= FFT length)
cfg.cfar.doppler.pfa          = 1e-8;
cfg.cfar.doppler.threshold_dB = [];

cfg.cfar.notchWidth   = 1;     % half-width of the zero-Doppler band that is not tested (MTI notch,
                               % clutter residue), in Doppler window mainlobes (1 -> +-2 bins)
cfg.cfar.notchMaxRange_m = Inf; % apply notchWidth only through this range [m]; beyond it, test every Doppler bin
cfg.cfar.short.minRange_m = 0;
cfg.cfar.short.maxRange_m = Inf;
cfg.cfar.long.minRange_m  = 0;
cfg.cfar.long.maxRange_m  = Inf;

%% Pruning of the CFAR detections
cfg.prune.sidelobeMargin_dB = 20;  % a peak this much weaker than a peak on the same range, beam
                                   % or Doppler line is treated as its sidelobe and removed

%% Detection reports: noise floor of the SNR = power / noise - 1 of every detection and plot
cfg.report.snrNoise = 'median';
cfg.report.dopplerInterp   = false;
cfg.report.elevationInterp = true;

%% Plot extraction (azimuth peak picking over successive CPIs)
cfg.extract.minExtent_deg = 2; % minimum azimuth extent of the detections of a plot [deg]
                               % (noise false alarms last only a few CPIs, targets about one beamwidth)
cfg.extract.maxMiss = 2;       % CPIs without a detection before a plot is closed

%% Target tracking
cfg.track.method          = 'ekf';   % 'ekf', 'pmbm' or 'tomht'

cfg.track.enable          = true;
cfg.track.dim             = 'auto';  % 'auto' (3d with DBF, 2d without), '2d' or '3d'
cfg.track.model           = 'cv';    % motion model: 'cv' constant velocity or 'ca' constant acceleration
cfg.track.sigmaRange_m    = 50;      % range error of a plot (cell size, range-Doppler coupling) [m]
cfg.track.sigmaAz_deg     = 0.3;     % azimuth error of a plot [deg]
cfg.track.sigmaEl_deg     = 0.5;     % elevation error of a plot (3d) [deg]
cfg.track.useRadialVelocity = true;  % the radial velocity of the plot (Doppler bin) is a measurement too.
                                     % It is unambiguous when maxSpeed_mps is below half the span of the
                                     % Doppler axis (fixed PRF at 170 MHz: +-441 m/s); otherwise (staggered
                                     % PRF: +-112 m/s) it is unfolded with the range rate predicted by the
                                     % track, once that is known to 1/6 of the span
cfg.track.sigmaVr_mps     = [];      % radial velocity error of a plot ([] = 0.4 Doppler bins: 2.8 m/s) [m/s]
cfg.track.sigmaAccel_mps2 = 5;       % 'cv': horizontal process noise (white acceleration) [m/s^2]
cfg.track.sigmaAccelZ_mps2 = 1;      % 'cv': vertical process noise (3d) [m/s^2]
cfg.track.sigmaJerk_mps3  = 0.5;     % 'ca': horizontal process noise (white jerk) [m/s^3]
cfg.track.sigmaJerkZ_mps3 = 0.1;     % 'ca': vertical process noise (3d) [m/s^3]
cfg.track.maxSpeed_mps    = 400;     % initial velocity uncertainty of a new track: std maxSpeed/sqrt(3)
cfg.track.maxAccel_mps2   = 5;       % 'ca': initial acceleration uncertainty: std maxAccel/sqrt(3); with
                                     % 3 unknowns per axis 3 plots always fit a track, so a large value
                                     % makes false plots confirm tracks more easily
cfg.track.gnnSector_deg   = 5;       % GNN batch: associated after this antenna turn without a new plot [deg]
cfg.track.gnnMaxSector_deg = 30;     % GNN batch: at the latest this antenna turn after its first plot [deg]
cfg.track.gateProb        = 0.999;   % probability that the gate holds the correct plot (chi-square)
cfg.track.confirmHits     = 3;       % updates (incl. the first plot) that confirm a tentative track
cfg.track.maxMissTentative = 1;      % scans without a plot before a tentative track is deleted
cfg.track.maxMiss         = 3;       % scans without a plot before a confirmed track is deleted
cfg.track.minSpeedEnable  = false;   % true: delete the tracks slower than minSpeed_mps (clutter residue,
                                     % birds, ground vehicles). Speed of the state: sqrt(vx^2+vy^2+vz^2) in
                                     % 3-D, sqrt(vx^2+vy^2) in 2-D. Tested on every track with at least 2
                                     % updates (a new track has no measured velocity yet), tentative or confirmed
cfg.track.minSpeed_mps    = 5;       % minimum speed of a track when minSpeedEnable is true [m/s]

%% Target tracking with cfg.track.method = 'pmbm': PMBM filter, IMM-UKF state estimation
cfg.track.pmbm.pd                = 0.9;    % detection probability of a target in the beam (per scan)
cfg.track.pmbm.ps                = 0.99;   % survival probability of a target per scan
cfg.track.pmbm.falsePlotsPerScan = 2;      % expected false plots per scan (uniform in the measurement space)
cfg.track.pmbm.undetectedTargets = 1;      % expected targets in the coverage not yet detected (Poisson part)
cfg.track.pmbm.maxClimb_mps      = 100;    % 3-D: vertical velocity uncertainty of a new track: std
                                           % maxClimb/sqrt(3) (horizontal: cfg.track.maxSpeed_mps)
cfg.track.pmbm.maxHypotheses     = 25;     % global hypotheses kept
cfg.track.pmbm.hypPrune          = 1e-3;   % global hypotheses below this weight (of the best one) are dropped
cfg.track.pmbm.rPrune            = 1e-5;   % single-target hypotheses below this existence probability: no target
cfg.track.pmbm.rConfirm          = 0.9;    % probability of a target that confirms a track (with confirmHits plots)
cfg.track.pmbm.rDelete           = 1e-4;   % a track whose existence probability falls below is deleted (a
                                           % low value keeps the evidence of a new track that was missed once
                                           % in dense false plots; maxMissTentative / maxMiss end it as well)
cfg.track.pmbm.imm.models        = {'cv', 'ca', 'ct'};   % models of the IMM (any of them, in any order)
cfg.track.pmbm.imm.mu0           = [0.8 0.1 0.1];        % model probabilities of a new track
cfg.track.pmbm.imm.switchProb    = 0.05;   % probability per scan of switching to another model
cfg.track.pmbm.imm.sigmaAccelCv_mps2 = 2;  % CV: horizontal white acceleration [m/s^2]
cfg.track.pmbm.imm.sigmaJerkCa_mps3  = 1;  % CA: horizontal white jerk [m/s^3]
cfg.track.pmbm.imm.sigmaAccelCt_mps2 = 2;  % CT: horizontal white acceleration [m/s^2]
cfg.track.pmbm.imm.sigmaTurn_dps     = 0.5;% CT: change of the turn rate (random walk) [deg/s per sqrt(s)]
cfg.track.pmbm.imm.maxTurn_dps       = 6;  % CT: turn rate uncertainty of a new track: std maxTurn/sqrt(3) [deg/s]
                                           % (vertical: cfg.track.sigmaAccelZ_mps2, all models)
cfg.track.pmbm.ukf.alpha = 1;              % unscented transform: spread of the sigma points,
cfg.track.pmbm.ukf.beta  = 2;              % prior knowledge (2: Gaussian) and
cfg.track.pmbm.ukf.kappa = 0;              % secondary scaling

%% Target tracking with cfg.track.method = 'tomht': track-oriented MHT, IMM-UKF filter
cfg.track.tomht.pd                = 0.9;   % detection probability of a target in the beam (per scan)
cfg.track.tomht.falsePlotsPerScan = 2;     % expected false plots per scan (uniform in the measurement space)
cfg.track.tomht.undetectedTargets = 1;     % expected new targets per scan in the coverage (track initiation)
cfg.track.tomht.maxClimb_mps      = 100;   % 3-D: vertical velocity uncertainty of a new track (std maxClimb/sqrt(3))
cfg.track.tomht.maxHypotheses     = 25;    % k-best: global hypotheses kept
cfg.track.tomht.kBest             = 10;    % k-best: at most this many assignments (Murty) per global hypothesis
cfg.track.tomht.hypPrune          = 1e-4;  % global hypotheses below this likelihood ratio to the best one are dropped
cfg.track.tomht.nScan             = 3;     % N-scan pruning: decisions older than this many scans (looks) are fixed
cfg.track.tomht.confirmScore      = 10;    % LLR that confirms a track (with confirmHits plots; 3 plots on a
                                           % target: 15 - 30, 3 false plots that fit a track: rarely above 5)
cfg.track.tomht.deleteScore       = 10;    % a track hypothesis whose LLR fell this much below its maximum is deleted
cfg.track.tomht.pDelete           = 1e-3;  % a track whose probability over the global hypotheses is below is deleted
cfg.track.tomht.imm.models        = {'cv', 'ca', 'ct'};   % IMM-UKF filter of every track hypothesis (as the PMBM)
cfg.track.tomht.imm.mu0           = [0.8 0.1 0.1];
cfg.track.tomht.imm.switchProb    = 0.05;
cfg.track.tomht.imm.sigmaAccelCv_mps2 = 2;
cfg.track.tomht.imm.sigmaJerkCa_mps3  = 1;
cfg.track.tomht.imm.sigmaAccelCt_mps2 = 2;
cfg.track.tomht.imm.sigmaTurn_dps     = 0.5;
cfg.track.tomht.imm.maxTurn_dps       = 6;
cfg.track.tomht.ukf.alpha = 1;
cfg.track.tomht.ukf.beta  = 2;
cfg.track.tomht.ukf.kappa = 0;
end

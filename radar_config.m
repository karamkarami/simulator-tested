function cfg = radar_config()
%RADAR_CONFIG  All user-editable parameters of the rotating search radar simulator.
%
%   cfg = RADAR_CONFIG() returns a struct holding every parameter of the
%   simulation. Edit the values below, or change fields of cfg after calling
%   this function (this is what the tests do), e.g.
%
%       cfg = radar_config();
%       cfg.prf.mode = 'staggered';
%       cfg.clutter.enable = false;
%       results = run_radar_sim(cfg);
%
%   Quantities that come from input .mat files (antenna patterns, DDC
%   filters, DBF matrix, staggered PRIs and MTI weights) are listed under
%   "Input files". Default files are created by make_input_files.m when they
%   are missing.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));

%% Physical constants
cfg.const.c  = 299792458;       % speed of light [m/s]
cfg.const.kB = 1.380649e-23;    % Boltzmann constant [J/K]
cfg.const.T0 = 290;             % reference noise temperature [K]

%% Input files
cfg.files.antenna = fullfile(root, 'inputs', 'antenna_pattern.mat'); % patterns, element spacing, gains
cfg.files.ddc     = fullfile(root, 'inputs', 'ddc_filters.mat');     % low-pass filters of the 2 DDC stages, 2 MHz output
cfg.files.ddc4MHz = fullfile(root, 'inputs', 'ddc_filters_4MHz.mat'); % the same for the 4 MHz output
cfg.files.dbf     = fullfile(root, 'inputs', 'dbf_matrix.mat');      % elevation beamforming matrix
cfg.files.stagger = fullfile(root, 'inputs', 'stagger_prf.mat');     % staggered PRIs + MTI weights
cfg.files.phaseCode = fullfile(root, 'inputs', 'phase_code.mat');    % phase codes + mismatch filters of both sub-pulses

%% Radar site and antenna
cfg.fc            = 170e6; % center frequency [Hz] (VHF; wavelength 1.76 m)
cfg.ant.nElements = 48;   % receive elements in elevation (1 = single antenna, beamforming disabled)
cfg.ant.elementSpacing_lambda = 0.66;  % element spacing in elevation [wavelengths]
% (azimuth: fixed 4.1 deg beam of every element; patterns come from the antenna file)
cfg.ant.height_m  = 6;    % antenna height above the ground [m]
cfg.ant.tilt_deg  = 11;   % array tilt above the horizon [deg]
cfg.ant.rpm       = 6;    % rotation speed [rev/min]
cfg.ant.az0_deg   = 0;    % antenna azimuth at t = 0, clockwise from North [deg]

%% Transmitter and receiver
cfg.tx.peakPower_W     = 20e3;   % peak transmit power [W]
cfg.rx.noiseFigure_dB  = 4;      % receiver noise figure [dB]
cfg.rx.losses_dB       = 3;      % system losses not modelled elsewhere (RF, radome, ...) [dB]
cfg.rx.gain_dB         = 60;     % gain of the receiver path from the antenna to the ADC [dB]: the
                                 % thermal noise is 180 LSB rms at the ADC (no loss by the ADC noise or
                                 % the 16-bit DDC output); the STC below keeps the near clutter
                                 % (+32 dB above full scale without it) within the ADC range
cfg.rx.txIsolation_dB  = 40;     % transmitter-to-receiver isolation (leakage while transmitting) [dB]
cfg.rx.maxRange_m      = 80e3;   % instrumented range of the long pulse: sets the receive window [m]
cfg.rx.addThermalNoise = true;   % false only for debugging / SNR measurements
cfg.rx.addLeakage      = true;   % false only for debugging / SNR measurements

%% STC (sensitivity time control)
cfg.rx.stc.enable = true;        % true: STC on; false: STC off (no attenuation). Without the STC lower
                                 % cfg.rx.gain_dB to 30 dB, otherwise the near clutter clips the ADC
                                 % (6 - 9 % of the receive samples at 60 dB: false plots; run_radar_sim warns)
% An attenuator in front of the receiver that attenuates the echoes (not the
% receiver noise) early in each receive window. The attenuation depends on
% the time dt since the end of the transmitted pulse, as the equivalent range
% c*dt/2: the nearest range from which echoes still arrive then, where the
% strongest clutter comes from (at the start of the long-pulse window the
% eclipsed clutter of 1 - 10 km). Rows [range_m attenuation_dB], linear in
% dB, none beyond the last row; see stc_demo.m (it also derives the
% attenuation the clutter needs, and compares with and without the STC).
% short-pulse window (and while transmitting): the clutter up to the horizon (10 km)
cfg.rx.stc.curve  = [    0 40;  1000 40;  2000 32;  4000 25;  6000 22;  8000 16;  9000 12;  10000 5;  10500 0];
% long-pulse window ([] = the curve above): its start also holds the eclipsed
% long-pulse echo of every target nearer than the long-pulse blind range
% (22.5 km), whose compressed response would give false plots up to 22.5 km
% behind the target; the curve goes on to there. It costs the long-pulse
% targets at 23 - 40 km up to 13 dB of SNR (1 m^2 at 23 km: 57 dB) and
% tapers their echo: range sidelobes -26 dB with the LFM pulse (pruned), but
% -11 ... -24 dB with phase codes (see README "170 MHz").
cfg.rx.stc.curveLong = [0 60;  3000 60;  6000 45;  9000 30;  12000 18;  15000 10;  18000 5;  22500 0];

%% ADC
cfg.adc.fs_Hz            = 120e6;  % sampling frequency [Hz]
cfg.adc.nBits            = 16;     % ADC word length [bit]
cfg.adc.enob             = 12.7;   % effective number of bits
cfg.adc.fullScale_dBm    = 10;     % power of a full-scale sine at the ADC input (2 Vpp into 50 ohm) [dBm]
cfg.adc.addAdcNoise      = true;   % false only for debugging / SNR measurements
cfg.adc.quantize         = true;   % false only for SNR measurements (keeps the analog values)
cfg.adc.fIF_Hz           = 30e6;   % IF carrier of the received signal (receiver LO = fc - 30 MHz)
cfg.adc.overloadWarn     = 0.01;   % warn when more than this fraction of the receive samples of a PRI
                                   % is at full scale (cfg.rx.gain_dB too high for the echoes)

%% Waveforms: every PRI holds a short and a long sub-pulse (short first)
% PRI timeline: short pulse at its tTx -> receive (near range) -> long pulse
% at its tTx -> receive (far range) until the end of the record. While a
% pulse is transmitted the receiver sees the transmitter leakage. Both
% sub-pulses are on the same 30 MHz IF and are separated in time. The long
% pulse is blind up to c*150us/2 = 22.5 km; with tTx = 200 us the short
% pulse covers 4.2 ... 23.7 km and hands over to the long pulse.
% Pulse compression: LFM with the Taylor-windowed matched filter; a phase
% code (chips from cfg.files.phaseCode, default P4) with its matched filter
% or with the mismatch filter of the same file (lower range sidelobes for a
% small SNR loss). The mismatch filter may have any length (not shorter than
% the code): it sees zeros outside the receive interval, so the range cells,
% the receive window and the DDC output are the same as with the matched filter.
cfg.wave.short.type         = 'lfm';    % 'lfm' or 'phasecode' (code from cfg.files.phaseCode)
cfg.wave.short.width_s      = 25e-6;    % pulse width [s]
cfg.wave.short.bandwidth_Hz = 1e6;      % LFM sweep bandwidth [Hz]
cfg.wave.short.chip_s       = 1e-6;     % phase-code chip length [s] (width_s/chip_s chips)
cfg.wave.short.compression  = 'matched'; % phase code: 'matched' or 'mismatch' filter (LFM: 'matched')
cfg.wave.short.mismatchLength = [];     % coefficients (one per chip) of the mismatch filter in the default
                                        % phase-code file, e.g. 3*25 ([] = chips + 10, at most 2x the chips)
cfg.wave.short.tTx_s        = 0;        % transmit time inside the PRI [s]

cfg.wave.long.type          = 'lfm';
cfg.wave.long.width_s       = 150e-6;
cfg.wave.long.bandwidth_Hz  = 1e6;
cfg.wave.long.chip_s        = 1e-6;
cfg.wave.long.compression   = 'matched';
cfg.wave.long.mismatchLength = [];
cfg.wave.long.tTx_s         = 200e-6;

%% PRF
cfg.prf.mode  = 'fixed';   % 'fixed' or 'staggered' (staggered PRIs come from cfg.files.stagger)
cfg.prf.pri_s = 1e-3;      % PRI in fixed mode [s]

%% DDC (filters come from cfg.files.ddc for 2 MHz, cfg.files.ddc4MHz for 4 MHz)
cfg.ddc.outputRate_Hz = 2e6;   % DDC output sample rate: 2e6 or 4e6 [Hz]. Stage 1 decimates by 10
                               % (120 -> 12 MHz), stage 2 by 6 (-> 2 MHz) or by 3 (-> 4 MHz). The
                               % 4 MHz output halves the range cell (37.5 m instead of 75 m) and doubles
                               % the matched filter taps and range cells; the CFAR windows are in cells.
cfg.ddc.stage1Decimation = 10; % decimation of stage 1 (120 -> 12 MHz)
cfg.ddc.coefFracBits = 13;     % filter coefficients in the file: signed 16-bit, 13 fractional bits
cfg.ddc.shift      = [13 13];  % right shift of the accumulator of each stage (13 = unit gain);
                               % every stage output is rounded and saturated to signed 16 bit
cfg.ddc.fixedPoint = true;     % false only for SNR measurements (no rounding / saturation)

%% Matched filter and Doppler processing
cfg.proc.taylorNbar   = 4;     % Taylor window (LFM only): number of nearly constant sidelobes
cfg.proc.taylorSll_dB = -35;   % Taylor window peak sidelobe level [dB]
cfg.proc.cpiLength    = [];    % pulses per CPI buffer ([] = pulses arriving while a target is in the
                               % 3 dB azimuth beam: round(beamwidth / (rotation rate * PRI)) = 114)
cfg.proc.nfft         = [];    % Doppler FFT length ([] = next power of 2: of the pulses in the buffer
                               % with a fixed PRF (114 -> 128), of the samples per row when staggered (28 -> 32))
cfg.proc.updateEvery  = 1;     % run MTI/Doppler/CFAR every N new pulses (1 = every pulse)

%% Detection: cell-averaging CFAR on every range-Doppler map
% method: '2d'      window in range and Doppler around the cell under test
%         'range'   1-D window along range (same Doppler bin)
%         'doppler' 1-D window along Doppler (same range cell; circular)
% Every method has its own window and threshold. The threshold factor alpha
% comes from pfa (and the training cells), unless threshold_dB is given:
% alpha = 10^(threshold_dB/10). The PPI threshold field changes the threshold
% of the selected method while the radar runs.
% Doppler is circular, so the window is always complete. At the first and
% last range cells of a sub-pulse the range window is shifted inwards (same
% number of training cells, the cell under test is no longer in the middle),
% so every range cell can be tested.
cfg.cfar.method = '2d';

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
% Range interval in which each sub-pulse is tested (all methods); the cells
% outside it are only used as training cells. The short pulse also stops
% where the long pulse starts (handover).
cfg.cfar.short.minRange_m = 0;
cfg.cfar.short.maxRange_m = Inf;
cfg.cfar.long.minRange_m  = 0;
cfg.cfar.long.maxRange_m  = Inf;

%% Pruning of the CFAR detections
cfg.prune.sidelobeMargin_dB = 20;  % a peak this much weaker than a peak on the same range, beam
                                   % or Doppler line is treated as its sidelobe and removed

%% Detection reports: noise floor of the SNR = power / noise - 1 of every detection and plot
% 'median' the median of the range-Doppler map of the beam, scaled to the mean: the receiver
%          noise, not raised by the range sidelobes of strong targets or by a second target
%          near the detection
% 'cfar'   the noise estimate of the CFAR window of the cell (mean of its training cells):
%          no extra computation (no median of the map), but the SNR of a target near a
%          strong one (range sidelobes, Doppler mainlobe) or near clutter residue is too low,
%          and with few training cells it scatters more (32 cells: +-0.8 dB)
cfg.report.snrNoise = 'median';

%% Plot extraction (azimuth peak picking over successive CPIs)
cfg.extract.minExtent_deg = 2; % minimum azimuth extent of the detections of a plot [deg]
                               % (noise false alarms last only a few CPIs, targets about one beamwidth)
cfg.extract.maxMiss = 2;       % CPIs without a detection before a plot is closed

%% Target tracking: extended Kalman filter (EKF) on the plots
% State relative to the antenna (x = East, y = North, z = up): position and
% velocity ('cv': nearly constant velocity, white acceleration) or position,
% velocity and acceleration ('ca': nearly constant acceleration, white jerk).
% Measurements: range and azimuth of every plot ('2d', in the horizontal
% plane, ground range), or range, azimuth and elevation ('3d', needs the
% elevation from the DBF).
% Association: global nearest neighbour (GNN). The plots are collected until
% the antenna has turned gnnSector_deg without a new plot (or for at most
% gnnMaxSector_deg after the first one); then every track gets at most one plot of the
% batch: the assignment with the most plot-track pairs inside their gates and,
% among those, the smallest sum of d^2 + ln|S| (Hungarian algorithm). A plot
% without a track starts a tentative track (velocity, acceleration unknown).
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

%% Ground clutter
cfg.clutter.enable           = true;
cfg.clutter.sigma0_dB        = -30;   % land clutter reflectivity [dB m^2/m^2]
cfg.clutter.minRange_m       = 1000;  % clutter starts here [m]
cfg.clutter.maxRange_m       = [];    % [] = radar horizon (4/3 earth model)
cfg.clutter.azStep_deg       = 0.25;  % azimuth size of a clutter cell [deg]
cfg.clutter.rangeStepSamples = 12;    % range size of a clutter cell in ADC samples (12 -> 15 m)
% Doppler spectrum of the internal clutter motion (wind-blown vegetation):
%   'gaussian'    each clutter cell moves with a Gaussian random radial velocity
%                 of std icmStd_mps (Gaussian spectrum)
%   'exponential' Billingsley's model of land clutter: a DC part (stationary
%                 scatterers) and an AC part with an exponential spectrum
%                 P(v) = beta/2 * exp(-beta*|v|); DC/AC power ratio r. From the
%                 wind speed: 1/beta = 0.1*log10(w) + 0.068 m/s and
%                 r = 489.8 * w^-1.55 * f^-1.21 (w in mph, f in GHz)
cfg.clutter.spectrum         = 'gaussian';
cfg.clutter.icmStd_mps       = 0.1;   % 'gaussian': std of internal clutter motion [m/s]
cfg.clutter.windSpeed_mps    = 4.5;   % 'exponential': wind speed [m/s] (4.5 m/s = 10 mph, breezy)
cfg.clutter.expBeta_spm      = [];    % 'exponential': beta [s/m] ([] = from the wind speed)
cfg.clutter.dcToAc           = [];    % 'exponential': DC/AC power ratio r ([] = from the wind speed)
cfg.clutter.elevation_deg    = 0;     % elevation of the clutter returns [deg]

%% Simulation length
cfg.sim.nScans      = 1;      % number of antenna revolutions
cfg.sim.sectors_deg = [];     % [] = full 360 deg; else rows [azStart azEnd]: only these sectors are simulated
cfg.sim.seed        = 1;      % random seed
cfg.sim.swerling    = 1;      % target RCS model: 1 = Swerling 1 (new RCS every scan), 0 = constant

%% Targets (see make_target.m)
% Test scenario:
%   'six'       six_targets_scenario.m: six targets around the radar (0.5 - 10 m^2, 8 - 46 km)
%   'challenge' challenge_targets_scenario.m: ten targets that are hard for the radar
%               (low SNR, slow in the clutter, tangential in the MTI notch, behind a strong
%               near target at the pulse handover, two targets 300 m apart, 420 m/s, high
%               elevation, low altitude near the horizon); some are expected to be missed
%   'random'    random_targets_scenario.m: cfg.scenario.random.n targets with random position,
%               RCS and velocity in the limits below (the same seed gives the same targets)
% cfg.targets is set from the scenario below; after changing cfg.scenario in a
% script, set cfg.targets = target_scenario(cfg) again (or give cfg.targets directly).
cfg.scenario.type = 'six';                 % 'six', 'challenge' or 'random'
cfg.scenario.random.n         = 10;        % number of random targets
cfg.scenario.random.seed      = [];        % [] = cfg.sim.seed
cfg.scenario.random.rcs_m2    = [0.5 10];  % mean RCS [min max], log-uniform [m^2]
cfg.scenario.random.range_m   = [5e3 70e3];% ground range at t = 0 [min max] [m] (at most 95% of cfg.rx.maxRange_m)
cfg.scenario.random.alt_m     = [300 10000];% altitude [min max] [m]
cfg.scenario.random.maxEl_deg = 40;        % at most this elevation at t = 0 (limits the altitude of near targets) [deg]
cfg.scenario.random.speed_mps = [30 300];  % horizontal speed [min max] [m/s], heading uniform 0 - 360 deg
cfg.targets = target_scenario(cfg);

%% Recording of the raw ADC data and of the DDC output, streaming of the DDC output
cfg.io.outFolder = fullfile(root, 'output');
cfg.io.saveRaw   = true;                                    % raw ADC data: one .mat file per PRI
cfg.io.rawFolder = fullfile(root, 'output', 'raw_pulses');  % -> pulse_<n>.mat (adcData, hdr)
cfg.io.saveDdc   = false;                                   % DDC output: one .mat file per PRI
cfg.io.ddcFolder = fullfile(root, 'output', 'ddc_pulses');  % -> ddc_<n>.mat (I, Q, hdr); replay_ddc_udp.m
                                                            %    sends these files with the PRI timing
cfg.io.matFormat = '-v6';                                   % '-v6' fast/uncompressed, '-v7' compressed

cfg.udp.enable       = false;        % stream the DDC output of every PRI over Ethernet/UDP
cfg.udp.host         = '127.0.0.1';
cfg.udp.port         = 50000;
cfg.udp.payloadBytes = 1400;         % DDC data bytes per datagram (int16 I/Q; + 64-byte header fits a 1500-byte
                                     % Ethernet MTU; use up to 65400 with jumbo frames / localhost)

%% Console output
% Printing every plot can take longer than the processing when a low CFAR
% threshold gives many plots: switch cfg.console.plots off then.
cfg.console.enable = true;   % MATLAB: summary, SNR budget, progress lines, CFAR threshold changes,
                             % warnings of the run and the lists below
cfg.console.plots  = true;   % every plot and every track: MATLAB lists them at the end of the run, the C
                             % receiver prints them while it runs (its default; --quiet / --print-plots)

%% Display
cfg.display.enable            = true;
cfg.display.visible           = true;     % false = draw off-screen (figures are still saved)
cfg.display.updateEvery       = 50;       % PPI refresh period [pulses]
cfg.display.maxRange_m        = [];       % [] = cfg.rx.maxRange_m
cfg.display.showTruth         = true;     % draw the true target positions
cfg.display.showMode          = 'plots';  % PPI shows 'plots' (after the plot extractor), 'detections'
                                          % (CFAR detections before it: arcs) or 'both' (also in the panel)
cfg.display.maxDetections     = 100000;   % CFAR detections kept for the PPI (the newest)
cfg.display.showTracks        = true;     % draw the tracks (hollow triangles; also a checkbox in the panel)
cfg.display.trackLeader_s     = 30;       % length of the velocity leader of a track [s of flight]
cfg.display.rangeProfilePulse = [];       % pulse number whose matched filter output is plotted
cfg.display.rangeProfileAxis  = 'range';  % 'range', 'delay' or 'index'
end

function cfg = radar_config_mode_0(cfg)
%RADAR_CONFIG_MODE_0  Settings of radar mode 0 (see "Radar modes" in radar_config.m).
%
%   Mode 0 'medium_mode': both sub-pulses are phase codes (P4 of
%   cfg.files.phaseCode: 25 chips for the short, 150 chips for the long pulse,
%   1 us chips) compressed with their mismatch filters (lower range sidelobes
%   than the matched filter for a small SNR loss), fixed PRF of 1 kHz.
%
%   cfg = RADAR_CONFIG_MODE_0(cfg) sets everything of this mode in cfg: the
%   waveforms, transmit times and receive windows of both sub-pulses (each
%   one on or off), the PRF, the input files of the
%   processing (stagger, scan-to-scan PRFs, phase codes, FIR clutter filter),
%   the STC and the whole receiver processing chain (pulse compression, CPI
%   buffer, Doppler window and FFT, clutter filter, CFAR, pruning, reports,
%   plot extraction, tracking). Common to all modes (radar_config.m): the
%   antenna, transmitter, receiver gain, ADC, DDC, clutter, targets, recording
%   and display. radar_config.m calls this file (through apply_radar_mode.m)
%   when this mode is cfg.modes.active. Every mode file sets the same fields
%   (all of them), so that selecting another mode with apply_radar_mode
%   replaces them all. A mode with one sub-pulse per PRI switches the other
%   one off (cfg.wave.<pulse>.enable = false).
%   The meaning of every setting is explained here and in README.md.
%   The file is also evaluated by radar_export_config (C receiver): only
%   assignments, no if/for blocks.

root = fileparts(fileparts(mfilename('fullpath')));   % folder of radar_config.m

cfg.modes.name = 'medium_mode';         % name of the mode (run summary, receiver configuration)

%% Waveforms and timing: every PRI holds a short and a long sub-pulse (short first)
% Times in the PRI [s] from its start (the first ADC/DDC sample): transmission at tTx_s, receive
% window [rxStart_s rxEnd_s) ([] = from the end of the pulse; short pulse: to the transmission of
% the long pulse, long pulse: to the end of the echoes from cfg.rx.maxRange_m). The DDC output
% runs from the PRI start to the end of the last receive window; the receiver takes only the
% receive windows out of it (the samples while transmitting and between the windows are not used).
% Range cells of a sub-pulse: the echoes that lie whole inside its window, i.e. the delays
% rxStart_s - tTx_s ... rxEnd_s - tTx_s - width_s.
% PRI timeline: short pulse at its tTx -> receive (near range) -> long pulse
% at its tTx -> receive (far range) until the end of its receive window. While a
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
% Each sub-pulse can be switched off (enable = false: not transmitted, not processed; with the
% long pulse off the short-pulse receive window lasts to cfg.rx.maxRange_m) and can have its own
% carrier (fc_Hz, [] = cfg.fc). The IF stays at cfg.adc.fIF_Hz (30 MHz at the DDC input for both):
% the receiver LO is tuned to fc_Hz - cfg.adc.fIF_Hz of each sub-pulse (LO of the short pulse, from
% the transmission of the long pulse on that of the long pulse), the DDC NCO stays at 30 MHz, and
% the wavelength of each sub-pulse sets its radial velocity
% axis, received power (~ lambda^2), element phases (beam directions) and clutter Doppler. Carriers
% more than the pulse bandwidth apart also keep the short-pulse echo of far targets out of the
% long window (the DDC filters reject it).
%
% Code agility: cfg.files.phaseCode holds up to cfg.wave.maxCodes code / mismatch-filter pairs per
% sub-pulse (<pulse>Code_<k>, <pulse>Mismatch_<k>, k = 1, 2, ...; the default file: maxCodes pairs,
% code 1 P4, the others polyphase codes of other roots); cfg.wave.<pulse>.codeIndex selects the pair
% of each sub-pulse. A file with more pairs than maxCodes, or a codeIndex above maxCodes or above
% the pairs of the file, stops with a message.
cfg.wave.maxCodes           = 4;        % phase codes: largest number of code / mismatch pairs per sub-pulse
cfg.wave.short.enable       = true;     % short sub-pulse on/off (false: only the long pulse)
cfg.wave.short.fc_Hz        = [];       % its carrier [Hz] ([] = cfg.fc)
cfg.wave.short.type         = 'phasecode'; % 'lfm' or 'phasecode' (code from cfg.files.phaseCode)
cfg.wave.short.width_s      = 25e-6;    % pulse width [s]
cfg.wave.short.bandwidth_Hz = 1e6;      % LFM sweep bandwidth [Hz]
cfg.wave.short.chip_s       = 1e-6;     % phase-code chip length [s] (width_s/chip_s chips)
cfg.wave.short.compression  = 'mismatch'; % phase code: 'matched' or 'mismatch' filter (LFM: 'matched')
cfg.wave.short.codeIndex    = 1;        % phase code: code / mismatch pair of cfg.files.phaseCode (1 ... maxCodes)
cfg.wave.short.mismatchLength = [];     % coefficients (one per chip) of the mismatch filter in the default
                                        % phase-code file ([] = chips + 10, at most 2x the chips)
cfg.wave.short.tTx_s        = 0;        % transmit time in the PRI [s]
cfg.wave.short.rxStart_s    = [];       % start of its receive window in the PRI [s] ([] = end of the pulse)
cfg.wave.short.rxEnd_s      = [];       % end of its receive window [s] ([] = transmission of the long pulse)

cfg.wave.long.enable        = true;     % long sub-pulse on/off (false: only the short pulse)
cfg.wave.long.fc_Hz         = [];
cfg.wave.long.type          = 'phasecode';
cfg.wave.long.width_s       = 150e-6;
cfg.wave.long.bandwidth_Hz  = 1e6;
cfg.wave.long.chip_s        = 1e-6;
cfg.wave.long.compression   = 'mismatch';
cfg.wave.long.codeIndex     = 1;
cfg.wave.long.mismatchLength = [];
cfg.wave.long.tTx_s         = 200e-6;   % transmit time in the PRI [s]
cfg.wave.long.rxStart_s     = [];       % start of its receive window [s] ([] = end of the pulse)
cfg.wave.long.rxEnd_s       = [];       % end of its receive window = end of the DDC record [s]
                                        % ([] = end of the echoes from cfg.rx.maxRange_m)
cfg.rx.maxRange_m      = 80e3;   % instrumented range [m]: the end of the long-pulse receive window when its
                                 % rxEnd_s is [] (else set it to match the window); the PPI range
% DDC record of every PRI: automatic (from the PRI start to the end of the last receive window) or,
% with manualDdc = true, nDdc DDC samples whose first one is at ddcT0_s in the PRI ([] = the
% automatic time: half the input span of a DDC sample after the PRI start). ddcT0_s is rounded to
% the ADC sample grid; it may be before the PRI start (e.g. 0). The receive windows are taken out
% of the record as with the automatic one (a set rxEnd_s must lie inside it); the record must end
% within the shortest PRI.
cfg.rx.manualDdc = false;        % true: the DDC record of this mode is set by hand (nDdc, ddcT0_s)
cfg.rx.nDdc      = [];           % manualDdc: DDC samples per PRI
cfg.rx.ddcT0_s   = [];           % manualDdc: time of the first DDC sample in the PRI [s] ([] = automatic)

%% PRF: 'fixed' (pri_s), 'staggered' (PRIs of cfg.files.stagger, with their MTI weights) or
% 'scan2scan' (one PRF per antenna scan from cfg.files.scan2scan, constant during the scan)
% PRF modes:
%   'fixed'     the same PRI cfg.prf.pri_s in every PRI
%   'staggered' the PRIs of cfg.files.stagger in turn (with their MTI weights)
%   'scan2scan' one PRF per antenna scan, constant during the scan: the PRFs
%               of cfg.files.scan2scan (variable prf_Hz, one per scan, used in
%               turn and repeated). Inside a scan the processing is that of the
%               fixed PRF (MTI canceller or FIR filter, Doppler FFT, CFAR, notch);
%               the CPI length ('auto'), the FFT length, the Doppler axis and the
%               blind speeds follow the PRF of the scan.
cfg.prf.mode  = 'fixed';
cfg.prf.pri_s = 1e-3;                   % PRI in fixed mode [s] (at least the end of the last receive window)
% Staggered PRF: cfg.files.stagger holds up to maxStaggerSets pairs of a PRI sequence T
% (priVector_s_<k>) and its MTI weights w (mtiCoeffs_<k>: one row per stagger position, applied
% time-varying: each MTI output uses the row of the stagger position of its oldest pulse);
% staggerSet selects the pair. More pairs in the file than maxStaggerSets, or a staggerSet above
% maxStaggerSets or above the pairs of the file, stops with a message.
cfg.prf.staggerSet     = 1;     % staggered PRF: pair (T, w) of cfg.files.stagger (1 ... maxStaggerSets)
cfg.prf.maxStaggerSets = 4;     % staggered PRF: largest number of (T, w) pairs

%% Input files of this mode (each mode may have its own files)
cfg.files.stagger   = fullfile(root, 'inputs', 'stagger_prf.mat');    % staggered PRIs + MTI weights (cfg.prf.mode = 'staggered')
cfg.files.scan2scan = fullfile(root, 'inputs', 'scan2scan_prf.mat');  % PRFs of the scans (cfg.prf.mode = 'scan2scan')
cfg.files.phaseCode = fullfile(root, 'inputs', 'phase_code.mat');     % phase codes + mismatch filters of both sub-pulses
cfg.files.mtiFir    = fullfile(root, 'inputs', 'mti_fir.mat');        % coefficients of the FIR clutter filter (cfg.mti.type = 'fir')
% A file made by make_input_files.m (default stagger, codes, FIR) is rebuilt
% when it is missing or was made for other settings; a file of your own (without
% the variable 'madeFor') is used as it is. Modes with different phase-code
% lengths need their own phaseCode file (e.g. 'phase_code_mode_2.mat').

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
% 'canceller' the three-pulse canceller [1 -2 1] with a fixed PRF, the
%             time-varying weights of cfg.files.stagger with a staggered PRF
% 'fir'       an FIR filter of your own: the variable 'coeffs' of cfg.files.mtiFir,
%             real, one row of L taps (any L >= 1, used for every pulse) or one row per
%             stagger position (staggered PRF: row k for the outputs whose oldest pulse
%             has stagger position k). Output m = sum_l coeffs(l) * x(m + l - 1), x
%             oldest first (the newest pulse meets the last tap). The default file
%             (make_input_files.m) holds a 5-tap high-pass with a deeper and wider notch
%             than the canceller (against moving clutter) and a narrower passband. The
%             CPI loses L - 1 pulses to the filter; the noise gain of every Doppler bin
%             and the CFAR looks follow from the coefficients.
cfg.mti.type = 'canceller';    % 'canceller' or 'fir'
cfg.mti.maxRange_m = Inf;      % the clutter filter is applied to the range cells up to this range [m]
                               % (both sub-pulses); the cells beyond it go to the Doppler FFT unfiltered
                               % (the newest L pulses' last one, so every cell has the same CPI), with
                               % their own noise-gain equalisation. Inf = every range cell (as before).
                               % Use cfg.cfar.notchMaxRange_m to test the zero-Doppler bins beyond it.

%% Detection: CFAR on every range-Doppler map
% type (detector; the noise statistic of the training cells of the window):
%         'ca'      cell averaging: the mean of the training cells (best in
%                   homogeneous noise)
%         'go'      greatest of: the larger of the means of the two halves of
%                   the window (fewer false alarms at clutter edges, more
%                   masking of close targets)
%         'so'      smallest of: the smaller of them (resolves close targets,
%                   more false alarms at clutter edges)
%         'os'      order statistic: the k-th smallest training cell,
%                   k = round(osRank * training cells) (robust against other
%                   targets in the window; more work)
%         Halves (GO, SO): with the 2-D and the range window the cells at a
%         smaller and at a larger range than the cell under test, with the
%         Doppler window the bins below and above it.
% method (window):
%         '2d'      window in range and Doppler around the cell under test
%         'range'   1-D window along range (same Doppler bin)
%         'doppler' 1-D window along Doppler (same range cell; circular)
% Every type works with every window. Every window has its own size and
% threshold. The threshold factor alpha multiplies the statistic; it comes
% from pfa (and the training cells, the type and the looks: cfar_alpha.m),
% unless threshold_dB is given: alpha = 10^(threshold_dB/10). The PPI
% threshold field changes the threshold of the selected type and window
% while the radar runs, not below the threshold of cfg.cfar.maxPfa.
% Doppler is circular, so the window is always complete. At the first and
% last range cells of a sub-pulse the range window is shifted inwards (same
% number of training cells, the cell under test is no longer in the middle),
% so every range cell can be tested.
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
% Interpolation of the detected peaks (block "Peak interp." of the timing panel), after the pruning:
% dopplerInterp   the radial velocity from a parabola through the powers [dB] of the Doppler bins
%                 d-1, d, d+1 of the peak (circular), instead of the centre of bin d
% elevationInterp the elevation from a parabola through the powers [dB] of the beams b-1, b, b+1
%                 of the peak, when it has both neighbours (not the lowest or highest beam) and
%                 the parabola has a maximum; otherwise (and when false) the centre of beam b
cfg.report.dopplerInterp   = false;
cfg.report.elevationInterp = true;

%% Plot extraction (azimuth peak picking over successive CPIs)
cfg.extract.minExtent_deg = 2; % minimum azimuth extent of the detections of a plot [deg]
                               % (noise false alarms last only a few CPIs, targets about one beamwidth)
cfg.extract.maxMiss = 2;       % CPIs without a detection before a plot is closed

%% Target tracking
% Method of the tracking block:
%   'ekf'   extended Kalman filter (EKF) per track, global nearest neighbour (GNN)
%           association, M-of-N confirmation (this section below)
%   'pmbm'  multi-target Bayesian tracker: Poisson multi-Bernoulli mixture (PMBM)
%           filter with an IMM-UKF estimator of every target state (interacting
%           multiple models CV / CA / CT, each an unscented Kalman filter);
%           settings in cfg.track.pmbm (next section) and the common ones here
%           (dim, measurement errors, radial velocity, maxSpeed, maxAccel, batch,
%           gateProb, confirmHits, maxMiss*, speed gate)
%   'tomht' track-oriented multiple hypothesis tracker (TOMHT): track trees scored
%           by their log-likelihood ratio, global hypotheses by the k-best
%           assignments (Murty) and N-scan pruning, IMM-UKF filter of every track
%           hypothesis (as the PMBM); settings in cfg.track.tomht (section after
%           the PMBM) and the common ones here
% Both work in 2-D (one element: range and azimuth, as the 2-D EKF) and in
% 3-D (DBF: range, azimuth and elevation), cfg.track.dim. The C receiver
% runs the same trackers (tracker_ekf.c, tracker_pmbm.c).
cfg.track.method          = 'ekf';   % 'ekf', 'pmbm' or 'tomht'

% EKF: extended Kalman filter on the plots
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

%% Target tracking with cfg.track.method = 'pmbm': PMBM filter, IMM-UKF state estimation
% Poisson multi-Bernoulli mixture (PMBM, the conjugate prior of the standard
% point-target model): the targets not yet detected are a Poisson process
% (uniform intensity: undetectedTargets in the coverage), every detected
% target a Bernoulli component (existence probability r and state density),
% and the data association hypotheses a mixture: global hypotheses, each one
% single-target hypothesis per Bernoulli, with its weight. The plots of a GNN
% batch (gnnSector_deg / gnnMaxSector_deg above) update the mixture: every
% global hypothesis gives its best assignments (Murty's k-best, as many as its
% weight earns of maxHypotheses); a plot may also be a new target (a new
% Bernoulli from the Poisson part) or a false plot. The antenna looks at a
% Bernoulli once per scan: when it passes its predicted azimuth without a
% plot the Bernoulli is missed (r falls with pd). Track = Bernoulli; the
% tracks shown are those of the most likely global hypothesis (with their
% state). A track is confirmed when the probability that a target made its
% newest plot (over all global hypotheses, whatever the earlier plots of
% that target) is at least rConfirm and it has confirmHits plots in the
% most likely hypothesis; when the most likely hypothesis changes, the track
% that holds the newest plot of a confirmed track takes its number and
% confirmation (the same target with another history). It is deleted when its existence probability over all
% global hypotheses falls below rDelete (or after maxMiss / maxMissTentative
% missed scans).
% IMM-UKF: per Bernoulli the models CV (nearly constant velocity), CA
% (nearly constant acceleration) and CT (coordinated turn, turn rate in the
% state) on the common state [position; velocity; acceleration; turn rate]
% (2-D: 7, 3-D: 10 values), mixed with the Markov switching probability
% switchProb per scan. The unscented transform propagates the CT model and
% maps every model to the measurement (range, azimuth (, elevation) and the
% radial velocity).
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
% Every track is a tree of track hypotheses (the plot or the miss of every antenna
% look); its score is the log-likelihood ratio (LLR) target / false plots:
% log(pd*lambdaU/lambdaC) for a new track, + log(pd*l(z)/lambdaC) per plot (l the
% IMM-UKF likelihood), + log(1 - pd) per missed look (lambdaC, lambdaU: false plots
% and new targets per scan and unit measurement volume). The plots of a GNN batch
% (gnnSector_deg / gnnMaxSector_deg) update every global hypothesis (a compatible
% set of track hypotheses, score = the sum of theirs) with its k best assignments
% (Murty: plot to track, new track or false plot). Hypothesis management:
% k-best (at most maxHypotheses global hypotheses, none below hypPrune of the best
% one), N-scan pruning (the decisions older than nScan looks are those of the best
% global hypothesis: the other branches of each tree are cut), deletion of a track
% hypothesis whose score fell deleteScore below its maximum (or after maxMiss /
% maxMissTentative missed scans) and of a tree whose probability is below pDelete.
% A track is confirmed when its hypothesis in the best global hypothesis has a
% score of at least confirmScore and confirmHits plots. The tracks shown are those
% of the best global hypothesis. Filter: the IMM-UKF of the PMBM (imm, ukf below).
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

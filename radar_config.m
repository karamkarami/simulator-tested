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
%   The waveforms, timing (transmit times, receive windows), PRF, STC and the
%   whole receiver processing chain (CFAR, tracker, windows, CPI buffer, ...)
%   belong to the radar modes: the files modes/radar_config_mode_<k>.m
%   (section "Radar modes" below); this file holds what all modes share.
%
%   Quantities that come from input .mat files are listed under "Input
%   files" (antenna patterns, DDC filters, DBF matrix) and in the mode files
%   (staggered PRIs and MTI weights, PRFs per scan, phase codes, FIR clutter
%   filter). Default files are created by make_input_files.m when they are
%   missing.

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
% (the files of the processing (staggered PRIs, PRFs per scan, phase codes, FIR clutter filter)
%  belong to the radar modes: cfg.files.* in modes/radar_config_mode_<k>.m)

%% Radar site and antenna
cfg.fc            = 170e6; % center frequency [Hz] (VHF; wavelength 1.76 m)
cfg.fcMin_Hz      = 150e6; % carrier range of the VBDR board [Hz]: the commands F1 (short pulse) and F2
cfg.fcMax_Hz      = 190e6; % (long pulse) set the carrier fcMin_Hz + 0 ... 255 MHz in steps of 1 MHz;
                           % with cfg.udp.format = 'vbdr' the carriers of the modes must be whole MHz in this range
cfg.ant.nElements = 48;   % receive elements in elevation (1 = single antenna, beamforming disabled)
cfg.rx.nChannelsUsed = [];  % channels of the DDC output the receiver processes: the first nChannelsUsed
                            % of the cfg.ant.nElements in every PRI ([] = all; 1 = only the first element,
                            % no DBF). The DDC output (simulator, files, UDP stream, logs) keeps all of
                            % them; the DBF uses the first nChannelsUsed columns of its matrix
cfg.ant.elementSpacing_lambda = 0.66;  % element spacing in elevation [wavelengths]
% (azimuth: fixed 4.1 deg beam of every element; patterns come from the antenna file)
cfg.ant.height_m  = 6;    % antenna height above the ground [m]
cfg.ant.tilt_deg  = 11;   % array tilt above the horizon [deg]
% (rotation speed: cfg.ant.availableRpm / cfg.ant.rpm of the radar mode, modes/radar_config_mode_<k>.m)
cfg.ant.az0_deg   = 0;    % antenna azimuth at t = 0, clockwise from North [deg]

%% Transmitter and receiver
cfg.tx.peakPower_W     = 4e3;    % peak transmit power [W]
cfg.rx.noiseFigure_dB  = 4;      % receiver noise figure [dB]
cfg.rx.losses_dB       = 3;      % system losses not modelled elsewhere (RF, radome, ...) [dB]
cfg.rx.gain_dB         = 60;     % gain of the receiver path from the antenna to the ADC [dB]: the
                                 % thermal noise is 180 LSB rms at the ADC (no loss by the ADC noise or
                                 % the 16-bit DDC output); the STC below keeps the near clutter
                                 % (+32 dB above full scale without it) within the ADC range
cfg.rx.txIsolation_dB  = 40;     % transmitter-to-receiver isolation (leakage while transmitting) [dB]
cfg.rx.addThermalNoise = true;   % false only for debugging / SNR measurements
cfg.rx.addLeakage      = true;   % false only for debugging / SNR measurements

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

%% Radar modes: waveforms, timing, PRF, STC and the whole receiver processing chain
% The radar has at most 8 modes (numbers 0 ... 7). Each mode has its own
% settings in its mode file modes/radar_config_mode_<k>.m: the waveforms,
% transmit times and receive windows of both sub-pulses (cfg.wave.short,
% cfg.wave.long, each one on or off), the DDC record, the instrumented range (cfg.rx.maxRange_m), the PRF
% (cfg.prf) with its files (cfg.files.stagger, cfg.files.scan2scan), the
% phase codes (cfg.files.phaseCode), the STC (cfg.rx.stc) and the whole
% processing chain: pulse compression and CPI buffer (cfg.proc, incl. the
% Doppler window), clutter filter (cfg.mti, cfg.files.mtiFir), CFAR
% (cfg.cfar), pruning (cfg.prune), detection reports (cfg.report), plot
% extraction (cfg.extract) and tracking (cfg.track: method, model, ...).
% Common to all modes (this file): the antenna, transmitter, receiver gain,
% ADC, DDC, clutter, targets, recording and display. A program uses some of
% the modes (cfg.modes.available, e.g. [0 1 2]: modes 0, 1 and 2) and runs one
% of them (cfg.modes.active), at most 8 modes 0 ... 7. A mode with one sub-pulse per PRI (only
% the short or only the long pulse) switches the other one off in its file
% (cfg.wave.<pulse>.enable = false).
%   mode 0 'medium_mode'      both sub-pulses phase codes (P4, mismatch filters), fixed PRF 1 kHz
%   mode 1 'lfm_mode'         both sub-pulses LFM (Taylor-windowed matched filter), fixed PRF 1 kHz
%   mode 2 'long_range_mode'  example of set transmit times and receive windows: short pulse at
%                             11 us (received 40 ... 248 us), long pulse of 125 us at 321 us
%                             (received 533 ... 1780 us, to 200 km), LFM, PRI 1.8 ms
%   modes 3 ... 7             not yet defined (the settings of mode 1): edit their files
% radar_config_mode_0.m explains every setting. apply_radar_mode.m (called
% below) runs the mode file of cfg.modes.active.
% Another mode in a script: set cfg.modes.active and call
% cfg = apply_radar_mode(cfg) again: it sets every setting of that mode
% (C receiver: radar_export_config --set "cfg.modes.active = 1; cfg = apply_radar_mode(cfg);").
% Changes of these settings after that (as the tests do) are kept. A new
% setting of the processing goes into every mode file.
cfg.modes.available = [0 1 2];  % modes of this program (at most 8 of the numbers 0 ... 7)
cfg.modes.active    = 0;        % mode the radar runs in (one of cfg.modes.available)
addpath(fullfile(root, 'modes'));
cfg = apply_radar_mode(cfg);    % every setting of the active mode (modes/radar_config_mode_<k>.m)

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

%% Ground clutter
cfg.clutter.enable           = true;
% Clutter power: CNR (clutter-to-noise ratio) of land at the reference range:
% the mean clutter power over the noise power at the matched filter output of
% one element (short pulse, one pulse, no STC, before the MTI), averaged over
% the antenna azimuth (clutter_cnr.m; the pulse compression gives clutter and
% noise the same gain, so it is also the CNR of a range cell). It sets the
% reflectivity sigma0 of land (printed in the run summary); every other
% terrain type is sigma0Rel_dB above or below it, so cnr_dB scales all the
% clutter. 72.7 dB at 5 km = sigma0 -30 dB. A higher CNR may need more STC
% attenuation (run_radar_sim warns about an ADC overload).
cfg.clutter.cnr_dB           = 72.7;  % CNR of land at cnrRange_m [dB] ([] = sigma0_dB below)
cfg.clutter.cnrRange_m       = 5000;  % reference range of the CNR [m]
cfg.clutter.sigma0_dB        = -30;   % land reflectivity [dB m^2/m^2], used only when cnr_dB = []
cfg.clutter.minRange_m       = 1000;  % clutter starts here [m]
cfg.clutter.maxRange_m       = [];    % [] = radar horizon (4/3 earth model) of the flat ground; elevated
                                      % terrain (hills, mountains) is seen beyond it
% Terrain around the radar (clutter_terrain.m):
%   'uniform'  land everywhere, to the horizon (the model before the terrain types)
%   'regions'  land, and the regions below on top of it (later rows win). Rows:
%              {type, [az1 az2] sector clockwise from az1 to az2 [deg, from North],
%               [r1 r2] ground range [m], height_m (mountain: of the ridge, hill: of the crests)}.
%              A cell is seen only when no nearer terrain on its azimuth is higher
%              in elevation (4/3 earth): back slopes are shadowed, sea and desert
%              end at the radar horizon (10 km), mountains are seen far beyond it.
cfg.clutter.terrain = 'regions';
cfg.clutter.regions = {
    'mountain', [ 30  75], [42e3 62e3], 2500   % mountain range in the North-East (strong, seen beyond the horizon)
    'hill',     [100 160], [ 3e3 12e3],  200   % hills in the East-South-East (patchy: shadowed valleys)
    'desert',   [180 250], [ 1e3 20e3],    0   % desert in the South-West (weak, nearly no motion)
    'sea',      [280 340], [ 3e3 30e3],    0   % sea in the West-North-West (waves: moving clutter)
};
% Terrain types: reflectivity relative to land [dB], internal motion (icmStd_mps = []: the spectrum
% of cfg.clutter.spectrum below, else Gaussian with this std around the mean radial velocity
% meanVel_mps*cos(az - meanDir_deg) of a drift toward meanDir_deg), shape: 'flat', 'ridge', 'rolling'.
cfg.clutter.types.land     = struct('sigma0Rel_dB',   0, 'icmStd_mps', [],   'meanVel_mps', 0, 'meanDir_deg',  0, 'shape', 'flat');
cfg.clutter.types.hill     = struct('sigma0Rel_dB',   3, 'icmStd_mps', [],   'meanVel_mps', 0, 'meanDir_deg',  0, 'shape', 'rolling');
cfg.clutter.types.mountain = struct('sigma0Rel_dB',   4, 'icmStd_mps', 0.05, 'meanVel_mps', 0, 'meanDir_deg',  0, 'shape', 'ridge');
cfg.clutter.types.desert   = struct('sigma0Rel_dB', -15, 'icmStd_mps', 0.02, 'meanVel_mps', 0, 'meanDir_deg',  0, 'shape', 'flat');
cfg.clutter.types.sea      = struct('sigma0Rel_dB', -15, 'icmStd_mps', 1.0,  'meanVel_mps', 2, 'meanDir_deg', 90, 'shape', 'flat');
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
cfg.clutter.elevation_deg    = 0;     % elevation of the returns of flat terrain (land, desert, sea) [deg];
                                      % hills and mountains return from their own elevation

%% Simulation length
cfg.sim.nScans      = 1;      % number of antenna revolutions
cfg.sim.sectors_deg = [];     % [] = full 360 deg; else rows [azStart azEnd]: only these sectors are simulated
cfg.sim.liveFidelity = 'full'; % radar_live_sim (real-time simulator): 'full' = as radar_generate_ddc (ADC at
                               % 120 MHz and the fixed-point DDC), 'fast' = the signal at the DDC rate (no ADC
                               % and DDC filters: much less computation for real time with many elements)
cfg.sim.seed        = 1;      % random seed
cfg.sim.swerling    = 1;      % target RCS model: 1 = Swerling 1 (new RCS every scan), 0 = constant

%% Targets (see make_target.m: constant velocity, or maneuvers: coordinated turns and constant
% accelerations, target_state.m)
% Test scenario:
%   'six'       six_targets_scenario.m: six targets around the radar (0.5 - 10 m^2, 8 - 46 km)
%   'challenge' challenge_targets_scenario.m: thirteen targets that are hard for the radar
%               (low SNR over the mountains, slow over the hills, tangential in the MTI notch,
%               behind a strong near target at the pulse handover, two targets 300 m apart,
%               420 m/s, high elevation over the sea, low altitude near the horizon) and for
%               the tracker (a 5 deg/s weave, a constant acceleration of 8 m/s^2, a turn
%               while slowing down and climbing); some are expected to be missed
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

%% Recording of the raw ADC data, the DDC, matched filter and Doppler outputs; streaming of the DDC output
% Every recording has its own switch (saveRaw, saveDdc, saveMf, saveDoppler) and folder. The matched
% filter and Doppler files are large (Doppler, 48 elements: about 12 MB per PRI and sub-pulse).
cfg.io.outFolder = fullfile(root, 'output');
cfg.io.saveRaw   = true;                                    % raw ADC data: one .mat file per PRI
cfg.io.rawFolder = fullfile(root, 'output', 'raw_pulses');  % -> pulse_<n>.mat (adcData, hdr)
cfg.io.saveDdc   = false;                                   % DDC output: one .mat file per PRI
cfg.io.ddcFolder = fullfile(root, 'output', 'ddc_pulses');  % -> ddc_<n>.mat (I, Q, hdr); replay_ddc_udp.m
                                                            %    sends these files with the PRI timing
cfg.io.saveMf    = false;                                   % matched filter output: one .mat file per PRI
cfg.io.mfFolder  = fullfile(root, 'output', 'mf_pulses');   % -> mf_<n>.mat (short, long: complex single
                                                            %    [nBeams x nRange], their range axes, hdr)
cfg.io.saveDoppler = false;                                 % Doppler processing output: one .mat file per PRI
cfg.io.dopplerFolder = fullfile(root, 'output', 'doppler_pulses');  % with a CPI update -> doppler_<n>.mat
                                                            %    (short, long: power single [nBeams x nRange x
                                                            %    nfft] after the noise-gain equalisation, [] for a
                                                            %    sub-pulse without an update; axes, hdr)
cfg.io.matFormat = '-v6';                                   % '-v6' fast/uncompressed, '-v7' compressed

cfg.udp.enable       = false;        % stream the DDC output of every PRI over Ethernet/UDP
cfg.udp.format       = 'rddc';       % packet format of the DDC stream: 'rddc' (64-byte header in every datagram,
                                     % all elements; README "Raw data, DDC files and UDP format") or 'vbdr' (the
                                     % VBDR board: one stream of N packets with a 28-byte label per PRI, the 4
                                     % channels CH1 (main) and CH2 ... CH4 (J1 ... J3); payload and N per mode:
                                     % cfg.udp.vbdr of the mode file; README "VBDR stream and board control")
cfg.udp.host         = '127.0.0.1';  % receiver PC (destination of the DDC stream)
cfg.udp.port         = 50000;        % its data port
cfg.udp.payloadBytes = 1400;         % 'rddc': DDC data bytes per datagram (int16 I/Q; + 64-byte header fits a 1500-byte
                                     % Ethernet MTU; use up to 65400 with jumbo frames / localhost)

%% Board control (VBDR board): ASCII commands over UDP, e.g. "CNTL W 0x02 0x1C" (write) and
% "CNTL R 0x02 0x00" (read); the board answers "CNTL W 0x02 0x1C Returned Value = 0x00000000" and
% "CNTL R 0x02 0x00 Returned Value = 0x0000001C" to the receiver PC. The PPI (tab Controls, Radar
% Settings) sends every setting and reads it back; the simulator (run_radar_sim, replay_ddc_udp,
% radar_stream_send, radar_live_sim) answers like the board. Addresses: README.
cfg.board.host      = '127.0.0.1';   % board (real board: 192.168.100.17)
cfg.board.port      = 5000;          % its command port
cfg.board.pcHost    = '127.0.0.1';   % receiver PC the board answers to (real board: 192.168.100.18)
cfg.board.replyPort = 5100;          % port of the answers on the receiver PC
cfg.board.timeout_s = 0.5;           % wait this long for an answer ...
cfg.board.retries   = 1;             % ... and send a command this many times more without one
cfg.board.sendRpm   = false;         % true: the rotation speed of the PPI also goes to the board (address
                                     % 0x32; only the simulator knows it), false: only the receiver uses it
cfg.board.nOperationModes = 4;       % options of the Operation Mode (address 0x28): 1 ... nOperationModes
cfg.board.operationMode   = 1;       % Operation Mode at the start (one value for all radar modes)

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
                                          % (CFAR detections before it: arcs), 'both' or 'video' (the
                                          % DDC output of the elements averaged non-coherently, receive
                                          % intervals of both sub-pulses joined, colour = power; also in
                                          % the panel), 'dbf' (DBF output of the selected beam) or 'mf'
                                          % (pulse compression output of the selected beam)
cfg.display.beam              = 1;        % beam after the DBF shown by 'dbf', 'mf' and the range-Doppler
                                          % tab (also in the panel: Beam)
cfg.display.tab               = 'ppi';    % tab at the start: 'ppi', or 'rd' (range-Doppler map of the
                                          % selected beam next to a small PPI). The map and the plot below
                                          % it are always computed, whichever tab is shown
cfg.display.rdLowerPlot       = 'spectrum'; % range-Doppler tab, the plot below the map: 'spectrum' (Doppler
                                          % spectrum of a range bin) or 'range' (range profile of a Doppler
                                          % bin); also in the panel (C receiver: --ppi-rd-plot)
cfg.display.rdRate_Hz         = 10;       % range-Doppler tab: maps handed to the display per second (C
                                          % receiver, radar_config.bin; MATLAB: at most once per PPI
                                          % refresh). The display redraws at most 25 times a second
cfg.display.maxDetections     = 100000;   % CFAR detections kept for the PPI (the newest)
cfg.display.showTracks        = true;     % draw the tracks (hollow triangles; also a checkbox in the panel)
cfg.display.trackLeader_s     = 30;       % length of the velocity leader of a track [s of flight]
cfg.display.rangeProfilePulse = [];       % pulse number whose matched filter output is plotted
cfg.display.rangeProfileAxis  = 'range';  % 'range', 'delay' or 'index'
end

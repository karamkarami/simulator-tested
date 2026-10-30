function info = generate_ddc_files(cfg, nWorkers)
%GENERATE_DDC_FILES  Simulate every PRI up to the DDC output and save it, in parallel (parfor).
%
%   info = GENERATE_DDC_FILES() uses radar_config().
%   info = GENERATE_DDC_FILES(cfg) uses a modified configuration.
%   info = GENERATE_DDC_FILES(cfg, nWorkers) runs on nWorkers parallel workers
%   (Parallel Computing Toolbox: a pool of this size is started if needed).
%   nWorkers = 0 runs every pulse in this MATLAB session. Default: the pool
%   that is open, else a pool of the default cluster profile. Without the
%   toolbox, and in GNU Octave, parfor runs as a plain loop.
%
%   Every active PRI of the schedule (cfg.sim.nScans, cfg.sim.sectors_deg)
%   goes through steps 1-3 of run_radar_sim only, independently of the others:
%     1. targets, clutter, transmitter leakage and thermal noise on the IF,
%        receiver gain, ADC (16 bit)
%     2. raw ADC data -> <cfg.io.rawFolder>/pulse_<n>.mat (only with cfg.io.saveRaw)
%     3. DDC (2 or 4 MHz output, cfg.ddc.outputRate_Hz) -> <cfg.io.ddcFolder>/ddc_<n>.mat
%   The DDC files are the same as those of run_radar_sim with cfg.io.saveDdc
%   (I, Q int16 [nDdc x nElements] and hdr), so replay_ddc_udp.m or the C
%   sender (radar_receiver/build/radar_stream_send) can send them to the
%   receiver. No processing after the DDC is done.
%
%   Random numbers: the target RCS of every scan and the clutter map are drawn
%   first with cfg.sim.seed, as in run_radar_sim (same targets and clutter).
%   Every pulse then has its own random stream for the noise (substream
%   <pulse number> of an mrg32k3a stream seeded with cfg.sim.seed; in Octave
%   rand/randn seeded with [seed; pulse number]), so the files do not depend
%   on the number of workers or on the order in which the pulses are done.
%   The noise is not the same random numbers as in run_radar_sim (which draws
%   from one stream in pulse order), only the same statistics.
%
%   info holds the number of files, their folder, the run time and the
%   pulse numbers.
%
%   Example:
%     cfg = radar_config;  cfg.io.saveRaw = false;
%     cfg.sim.sectors_deg = [20 45];
%     info = generate_ddc_files(cfg, 8);          % 8 workers
%     replay_ddc_udp(info.folder);                % stream them to the receiver

if nargin < 1 || isempty(cfg), cfg = radar_config(); end
if nargin < 2, nWorkers = []; end
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'));
tStart = tic;

%% Setup (as run_radar_sim: same random draws for the RCS and the clutter map)
rng(cfg.sim.seed);
sys     = setup_radar(cfg);
sched   = pulse_schedule(sys);
targets = draw_target_rcs(cfg.targets, max(sched.scan) + 1, cfg.sim.swerling);
clut = [];
if cfg.clutter.enable
    clut = clutter_init(sys);
end
ddcFolder = cfg.io.ddcFolder;
rawFolder = cfg.io.rawFolder;
saveRaw   = cfg.io.saveRaw;
matFormat = cfg.io.matFormat;
seed      = cfg.sim.seed;
if ~exist(ddcFolder, 'dir'), mkdir(ddcFolder); end
if saveRaw && ~exist(rawFolder, 'dir'), mkdir(rawFolder); end

pulses = find(sched.active);
nPulses = numel(pulses);
maxWorkers = open_pool(nWorkers);
if cfg.sim.verbose
    fprintf('Generating the DDC output of %d PRIs (%d element(s), DDC %g MHz) -> %s\n', nPulses, ...
            sys.nCh, sys.fsOut/1e6, ddcFolder);
end

%% One PRI per iteration, in any order
parfor (j = 1:nPulses, maxWorkers)
    k = pulses(j);
    seed_pulse(seed, k);
    t  = sched.time_s(k);
    az = sched.az_deg(k);
    clutterIf = [];
    if ~isempty(clut)
        clutterIf = clutter_signal(clut, sys, az, t);
    end
    adcData = adc_sample(generate_rx_signal(sys, targets, t, az, clutterIf), sys);
    hdr = pulse_header(sys, sched, k);
    if saveRaw
        save_pulse_file(rawFolder, matFormat, hdr, adcData);
    end
    [I, Q] = ddc(adcData, sys);
    save_ddc_file(ddcFolder, matFormat, hdr, I, Q);
end

info.nFiles    = nPulses;
info.folder    = ddcFolder;
info.pulses    = pulses;
info.runTime_s = toc(tStart);
info.sys       = sys;
if cfg.sim.verbose
    fprintf('%d DDC files written in %.1f s (%.1f ms per PRI)\n', nPulses, info.runTime_s, ...
            1e3*info.runTime_s/max(nPulses, 1));
end
end

function M = open_pool(nWorkers)
% Number of parfor workers (0 = this session); starts a pool when needed.
M = 0;
if exist('OCTAVE_VERSION', 'builtin') || exist('gcp', 'file') ~= 2 || isequal(nWorkers, 0)
    return                                           % Octave, no Parallel Computing Toolbox, or none asked
end
pool = gcp('nocreate');
if ~isempty(nWorkers) && (isempty(pool) || pool.NumWorkers ~= nWorkers)
    delete(pool);
    pool = parpool(nWorkers);
elseif isempty(pool)
    pool = parpool();                                % size of the default cluster profile
end
M = pool.NumWorkers;
end

function seed_pulse(seed, k)
% Random stream of pulse k: independent of the worker and of the order of the pulses.
if exist('OCTAVE_VERSION', 'builtin')
    rand('state', [seed; k]);   %#ok<RAND>
    randn('state', [seed; k]);  %#ok<RAND>
else
    s = RandStream('mrg32k3a', 'Seed', seed);
    s.Substream = k;
    RandStream.setGlobalStream(s);
end
end

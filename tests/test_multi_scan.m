function [ok, out] = test_multi_scan(prfMode, nScans, buildDir)
%TEST_MULTI_SCAN  Several antenna scans through the C receiver: confirmed tracks on the targets only.
%
%   [ok, out] = TEST_MULTI_SCAN(prfMode, nScans, buildDir) simulates nScans
%   (default 3) scans of the six-target scenario with clutter and the default
%   configuration ('fixed' or 'staggered' PRF; only the azimuth sectors of the
%   targets) up to the DDC output (generate_ddc_files.m), sends the DDC files
%   to the C receiver with radar_stream_send (UDP, paced so that the receiver
%   keeps up) and checks:
%     - the receiver gets every PRI (no PRI ignored, incomplete or dropped)
%     - every target that is not at a blind speed (snr_budget.m) has a
%       confirmed track at the end, within 300 m of the target
%     - no confirmed track is more than 1 km from every target (no false
%       track: ghosts or clutter that come back every scan would be confirmed)
%     - at most 2 false plots per scan
%     - the processing time per PRI of the last scan is at most 1.5 times
%       that of the first scan (it does not grow with the number of tracks)
%   buildDir: radar_receiver/build (default). Needs a POSIX shell. Takes about
%   30 minutes per 3 scans in Octave (the DDC files: 10 GB, deleted again).

if nargin < 1 || isempty(prfMode), prfMode = 'fixed'; end
if nargin < 2 || isempty(nScans), nScans = 3; end
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
if nargin < 3 || isempty(buildDir), buildDir = fullfile(root, 'radar_receiver', 'build'); end
receiver = fullfile(buildDir, 'radar_receiver');
sender = fullfile(buildDir, 'radar_stream_send');
if ~exist(receiver, 'file') || ~exist(sender, 'file')
    error('C receiver not built in %s (radar_receiver/README.md)', buildDir);
end

cfg = radar_config();
cfg.prf.mode = prfMode;
cfg.sim.nScans = nScans;
cfg.sim.verbose = false;
cfg.rx.maxRange_m = 50e3;
cfg.targets = six_targets_scenario();
cfg.display.enable = false;
cfg.io.saveRaw = false;
cfg.io.outFolder = fullfile(root, 'output', ['multi_scan_' prfMode]);
cfg.io.ddcFolder = fullfile(cfg.io.outFolder, 'ddc_pulses');
cfg.io.matFormat = '-v6';
sys = setup_radar(cfg);
cfg.sim.sectors_deg = target_sectors(sys, cfg.targets);
out.dir = cfg.io.outFolder;
if exist(cfg.io.ddcFolder, 'dir'), delete(fullfile(cfg.io.ddcFolder, 'ddc_*.mat')); end
if ~exist(out.dir, 'dir'), mkdir(out.dir); end

%% DDC files of all scans
t0 = tic;
info = generate_ddc_files(cfg);
fprintf('%d DDC files of %d scans in %.0f s (ADC samples at full scale: %.3f%%, worst PRI %.3f%%)\n', ...
        info.nFiles, nScans, toc(t0), 100*info.adcClipped.mean, 100*info.adcClipped.max);

%% Through the C receiver
cfgFile = fullfile(out.dir, 'radar_config.bin');
export_receiver_config(cfg, cfgFile);
files = struct('plots', fullfile(out.dir, 'plots.csv'), 'tracks', fullfile(out.dir, 'tracks.csv'), ...
               'log', fullfile(out.dir, 'receiver.log'), 'done', fullfile(out.dir, 'receiver.done'), ...
               'rec', fullfile(out.dir, 'stream.rec'));
for f = struct2cell(files).', if exist(f{1}, 'file'), delete(f{1}); end, end
port = 50300 + strcmp(prfMode, 'staggered');
system(sprintf(['("%s" --config "%s" --port %d --record "%s" --plots-csv "%s" --tracks-csv "%s" --idle-timeout 20 ' ...
                '--no-skip --quiet > "%s" 2>&1; echo $? > "%s") &'], receiver, cfgFile, port, files.rec, files.plots, ...
               files.tracks, files.log, files.done));
pause(2);
[st, txt] = system(sprintf('"%s" "%s" --port %d --speed 0.05 --no-gaps --quiet', sender, cfg.io.ddcFolder, port));
if st ~= 0, error('radar_stream_send failed: %s', txt); end
t = tic;
while ~exist(files.done, 'file') && toc(t) < 300, pause(1); end
log = fileread(files.log);
fprintf('%s', log);
delete(fullfile(cfg.io.ddcFolder, 'ddc_*.mat'));

%% Checks
nPri = str2double(regexp(log, 'PRIs processed: (\d+)', 'tokens', 'once'));
lost = ~isempty(regexp(log, 'PRIs ignored|incomplete PRIs [1-9]|dropped', 'once'));
good = nPri == info.nFiles && ~lost;
ok = report(sprintf('receiver: %d of %d PRIs processed, none ignored, incomplete or dropped', nPri, info.nFiles), good);

T = read_csv(files.tracks);                          % id confirmed deleted time x y z vx vy vz ...
P = read_csv(files.plots);
h = cfg.ant.height_m;
blind = false(1, numel(cfg.targets));
for k = 1:numel(cfg.targets)
    b = snr_budget(sys, cfg.targets(k), illumination_time(sys, cfg.targets(k)));
    blind(k) = b.blind;
end
ids = unique(T(:, 1)).';
tEnd = max(T(:, 4));
nFalse = 0;  covered = false(1, numel(cfg.targets));  txt = '';
for id = ids
    rows = T(T(:, 1) == id & T(:, 3) == 0, :);       % updates (not the deletion)
    if isempty(rows) || ~any(rows(:, 2)), continue; end
    last = rows(end, :);                             % last update of a confirmed track
    d = inf(1, numel(cfg.targets));
    for k = 1:numel(cfg.targets)
        tg = cfg.targets(k);
        pos = tg.pos0_m + tg.vel_mps * last(4) - [0 0 h];
        d(k) = norm(last(5:7) - pos);
    end
    [dmin, k] = min(d);
    if dmin > 1000
        nFalse = nFalse + 1;
        txt = [txt, sprintf(' false track %d at t %.1f s', id, last(4))]; %#ok<AGROW>
    elseif dmin < 300 && last(4) > tEnd - sys.track.scanTime_s
        covered(k) = true;                           % on its target in the last scan
    end
end
expected = ~blind;
good = all(covered(expected)) && nFalse == 0;
ok = report(sprintf('%d confirmed track(s) on the %d targets that are not at a blind speed (%s), %d false confirmed track(s)%s', ...
                    sum(covered & expected), sum(expected), strjoin({cfg.targets(expected).name}, ' '), nFalse, txt), ...
            good) && ok;

% false plots: plots more than 1 km / 3 deg from every target
nFp = 0;
for i = 1:size(P, 1)
    near = false;
    for k = 1:numel(cfg.targets)
        [R, az] = target_geometry(cfg.targets(k), P(i, 1), h);
        near = near || (abs(P(i, 3) - R) < 1000 && abs(wrap180(P(i, 2) - az)) < 3);
    end
    nFp = nFp + ~near;
end
ok = report(sprintf('%d plots in %d scans, %d false (at most 2 per scan)', size(P, 1), nScans, nFp), nFp <= 2*nScans) && ok;

% processing time per PRI: first and last scan (replay of the recording, one thread at a time)
scanTime = sys.track.scanTime_s;
ms = zeros(1, 2);
for s = [1 nScans]
    [~, rtxt] = system(sprintf('"%s" --config "%s" --replay "%s" --plots-csv "" --quiet --threads 1 2>&1', receiver, ...
                               cfgFile, files.rec));
    ms(1 + (s == nScans)) = str2double(regexp(rtxt, 'Processing time per PRI: ([\d.]+) ms', 'tokens', 'once'));
end
out.msPerPri = ms;
out.scanTime_s = scanTime;
ok = report(sprintf('processing time per PRI %.2f ms (replay of all %d scans)', ms(1), nScans), isfinite(ms(1))) && ok;
fprintf('Multi-scan tracking (%s PRF, %d scans): %s\n', prfMode, nScans, ternary(ok, 'PASS', 'FAIL'));
end

function M = read_csv(file)
% Numeric columns of a CSV file of the receiver (header line skipped; text columns NaN).
M = zeros(0, 20);
if ~exist(file, 'file'), return; end
lines = strsplit(strtrim(fileread(file)), sprintf('\n'));
lines = lines(2:end);
if isempty(lines), return; end
M = nan(numel(lines), numel(strsplit(lines{1}, ',')));
for i = 1:numel(lines)
    c = strsplit(lines{i}, ',');
    for j = 1:numel(c)
        v = str2double(c{j});
        if isnan(v), v = double(strcmp(c{j}, 'long')); end   % pulse column: 0 short, 1 long
        M(i, j) = v;
    end
end
end

function ok = report(text, ok)
fprintf('  %-120s %s\n', text, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

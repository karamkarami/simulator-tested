function stats = replay_ddc_udp(folder, varargin)
%REPLAY_DDC_UDP  Send recorded DDC output files over UDP with the original PRI timing.
%
%   REPLAY_DDC_UDP() sends the files ddc_<n>.mat of cfg.io.ddcFolder (written
%   by run_radar_sim with cfg.io.saveDdc = true) to cfg.udp.host:cfg.udp.port
%   of radar_config.m.
%   REPLAY_DDC_UDP(folder, name, value, ...) sends the files of another folder
%   and/or changes options:
%     'host', 'port', 'payloadBytes'  destination and DDC bytes per datagram
%                       (default: cfg.udp of radar_config.m)
%     'speed'    1      replay speed (2 = twice as fast, 0.5 = half as fast)
%     'spread'   0      the datagrams of a PRI are sent evenly over this fraction
%                       of the PRI (0 = all at the PRI start, as the simulator does;
%                       e.g. 0.8 gives a smooth data rate for small receive buffers)
%     'keepGaps' true   a skipped sector (gap in the pulse numbers) keeps its real
%                       duration; false = the next PRI follows after one PRI
%     'catchUp'  false  after a PRI that went out late, the following PRIs keep the
%                       PRI spacing (the replay gets longer); true = they are sent
%                       at once until the original schedule is reached again
%     'preloadMB' 2048  files read and split into datagrams before sending starts,
%                       up to this total file size [MB] (inf = all files); the
%                       others are read in the idle time between PRIs and datagrams
%     'verbose'  true
%   stats = REPLAY_DDC_UDP(...) returns the number of PRIs, datagrams and
%   bytes, the run time and the timing: stats.sent_s is the time each PRI went
%   out (relative to the first one), stats.late_s its delay against its schedule.
%
%   Nothing is processed: the stored int16 I/Q samples and header go out in
%   exactly the same datagrams as the live stream of run_radar_sim
%   (udp_make_packets.m, header format in README.md), so the receiver sees
%   no difference between a replay and the simulator.
%
%   Timing: the send time of every PRI comes from the PRI start times in the
%   file headers, on one clock (no accumulated drift). The loop sleeps until
%   20 ms before a send time and then busy-waits, so a PRI leaves within a few
%   microseconds of its time when the computer is fast enough. Reading the
%   files never delays a PRI: it is done only while the next datagram is not
%   yet due. Keep in mind:
%   * One PRI of the default radar (48 elements) is 342 kB, i.e. 2.7 Gbit/s
%     at a 1 ms PRI: this needs localhost or 10 Gbit Ethernet.
%   * Every datagram costs a fixed time in MATLAB (tens of microseconds). With
%     1400-byte payloads a PRI is 245 datagrams, too many for 1 ms; use
%     'payloadBytes' 65400 (6 datagrams; the network stack fragments them;
%     the C receiver accepts any payload size).
%   * Reading a file can take longer than one PRI; preload the files
%     ('preloadMB') when the idle time is too short to read them.
%   stats.late_s shows whether the timing was kept. For exact real-time
%   timing use the C sender, which reads the same files and sends the same
%   datagrams: radar_receiver/build/radar_stream_send <folder> (see
%   radar_receiver/README.md).
%
%   Example:
%     cfg = radar_config;  cfg.io.saveDdc = true;  cfg.io.saveRaw = false;
%     cfg.sim.sectors_deg = [20 45];
%     run_radar_sim(cfg);                          % records output/ddc_pulses/ddc_<n>.mat
%     stats = replay_ddc_udp([], 'host', '192.168.1.20', 'port', 50000);

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'));
cfg = radar_config();
if nargin < 1 || isempty(folder), folder = cfg.io.ddcFolder; end
opt = struct('host', cfg.udp.host, 'port', cfg.udp.port, 'payloadBytes', cfg.udp.payloadBytes, ...
             'speed', 1, 'spread', 0, 'keepGaps', true, 'catchUp', false, 'preloadMB', 2048, 'verbose', true);
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i}), error('replay_ddc_udp: unknown option ''%s''', varargin{i}); end
    opt.(varargin{i}) = varargin{i+1};
end
spin_s    = 20e-3;                                   % busy-wait before every send time (pause() can
                                                     % wake up late, on Windows by several ms) [s]
lateTol_s = 1e-4;                                    % a PRI sent later than this counts as late

files = dir(fullfile(folder, 'ddc_*.mat'));
if isempty(files), error('replay_ddc_udp: no ddc_*.mat files in %s', folder); end
[names, order] = sort({files.name});                 % zero-padded pulse numbers: file order = pulse order
names = cellfun(@(n) fullfile(folder, n), names, 'UniformOutput', false);
fileMB = cumsum([files(order).bytes]) / 2^20;
nPri = numel(names);

u = udp_sender_open(struct('host', opt.host, 'port', opt.port, 'payloadBytes', opt.payloadBytes));
closer = onCleanup(@() udp_sender_close(u));         % also on Ctrl+C or error

% Read and split the first PRIs before the clock starts
pri = cell(1, nPri);
nReady = max(1, sum(fileMB <= opt.preloadMB));
prepTime = 0;                                        % (pessimistic) time to prepare one PRI [s]
for k = 1:nReady
    t1 = tic;
    pri{k} = prepare_pri(names{k}, opt.payloadBytes);
    prepTime = max(prepTime, toc(t1));
end
if opt.verbose
    fprintf('Replaying %d PRIs from %s to %s:%d (speed %g, %d preloaded)\n', nPri, folder, opt.host, ...
            opt.port, opt.speed, nReady);
end

late = zeros(1, nPri);                               % delay of every PRI against its schedule [s]
sent = zeros(1, nPri);                               % time the first datagram of every PRI went out [s]
nDatagrams = 0;
nBytes = 0;
nRead = 0;                                           % PRIs that had to be read when they were due
clk = tic;
tStart = toc(clk) + 0.01;                            % send time of the first PRI
tNextPrint = 5;
prev = [];
for k = 1:nPri
    if k > nReady                                    % not prepared in time: read it now
        pri{k} = prepare_pri(names{k}, opt.payloadBytes);
        nReady = k;
        nRead = nRead + 1;
    end
    p = pri{k};
    pri{k} = [];                                     % free the memory
    if isempty(prev)
        tPri = tStart;
    else
        dt = p.time_s - prev.time_s;
        if ~opt.keepGaps && p.pulseNumber ~= prev.pulseNumber + 1
            dt = prev.pri_s;                         % skipped sector -> one PRI
        end
        tPri = tPri + dt/opt.speed;
    end
    prev = p;
    nPkt = numel(p.packets);
    dtPkt = opt.spread * p.pri_s / opt.speed / nPkt; % spacing of the datagrams of this PRI

    for j = 1:nPkt
        tSend = tPri + (j-1)*dtPkt;
        % use the time until tSend to read the next files
        while nReady < nPri && tSend - toc(clk) > prepTime + 1e-3
            t1 = toc(clk);
            pri{nReady+1} = prepare_pri(names{nReady+1}, opt.payloadBytes);
            nReady = nReady + 1;
            prepTime = max(0.95*prepTime, toc(clk) - t1);
        end
        wait_until(clk, tSend, spin_s);
        if j == 1
            sent(k) = toc(clk);
            late(k) = sent(k) - tPri;
            if ~opt.catchUp && late(k) > lateTol_s
                tPri = tPri + late(k);               % keep the PRI spacing from here on
            end
        end
        udp_send_datagram(u, p.packets{j});
        nBytes = nBytes + numel(p.packets{j});
    end
    nDatagrams = nDatagrams + nPkt;

    if opt.verbose && toc(clk) > tNextPrint
        fprintf('  PRI %d of %d (pulse %d), max delay so far %.3f ms\n', k, nPri, p.pulseNumber, 1e3*max(late(1:k)));
        tNextPrint = tNextPrint + 5;
    end
end
runTime = toc(clk) - tStart;
clear closer                                         % closes the socket

stats.nPri        = nPri;
stats.nDatagrams  = nDatagrams;
stats.nBytes      = nBytes;
stats.runTime_s   = runTime;
stats.late_s      = late;
stats.sent_s      = sent - sent(1);
stats.maxLate_s   = max(late);
stats.nLate       = sum(late > lateTol_s);
stats.nReadLate   = nRead;
if opt.verbose
    fprintf(['Sent %d PRIs (%d datagrams, %.1f MB) in %.2f s: %d PRIs more than 0.1 ms late, ', ...
             'max delay %.3f ms, median %.1f us\n'], nPri, nDatagrams, nBytes/2^20, runTime, stats.nLate, ...
            1e3*stats.maxLate_s, 1e6*median(late));
    if stats.nLate > 0
        fprintf(['  To keep the PRI timing: fewer datagrams per PRI (''payloadBytes''), a larger ''preloadMB'', ', ...
                 'or a lower ''speed''\n']);
    end
end
end

function p = prepare_pri(fileName, payloadBytes)
% Read one DDC file and split it into the datagrams of the live stream.
s = load(fileName);
p.packets     = udp_make_packets(s.hdr, s.I, s.Q, payloadBytes);
p.time_s      = s.hdr.time_s;
p.pri_s       = s.hdr.pri_s;
p.pulseNumber = s.hdr.pulseNumber;
end

function wait_until(clk, t, spin_s)
% Wait until toc(clk) = t: sleep for the coarse part, busy-wait for the last spin_s.
dt = t - toc(clk);
if dt > spin_s
    pause(dt - spin_s);
end
while toc(clk) < t
end
end

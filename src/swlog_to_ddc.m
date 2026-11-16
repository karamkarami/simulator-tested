function info = swlog_to_ddc(logPath, outFolder, varargin)
%SWLOG_TO_DDC  Convert a DDC log (.swlog) of the real receiver into DDC files of this project.
%
%   SWLOG_TO_DDC(logPath) reads every record (one PRI) of the log and saves
%   it as <folder>/ddc_<n>.mat in the format of src/save_ddc_file.m, the same
%   files that run_radar_sim (cfg.io.saveDdc = true) and generate_ddc_files
%   write:
%     I, Q   int16 [nDdc x nChannels]  DDC output of the PRI (column = channel)
%     hdr    header of the PRI (fields of src/pulse_header.m, see below)
%   Default folder: output/ddc_pulses_<log name>. replay_ddc_udp.m and
%   radar_receiver_<os>/build/radar_stream_send send these files to the
%   receiver.
%
%   The log is read directly (no SwlogReader object needed), with the
%   format of SwlogReader.m:
%     record = 84-byte header + float32 [n x 2 x 4]:
%       header bytes 0-3  uint32 PRF counter, 4-5 uint16 azimuth
%       (0 ... 16383 = 0 ... 360 deg), byte 10 log mode (0 ... 3)
%       samples: channel after channel, first the n I values, then the n Q values
%     log mode           0     1     2     3
%     samples n in log   7200  6660  2910  2660
%     DDC record nDdc    9317  9387  4937  3059
%   The log stores only the receive windows; they are put at their place in
%   the DDC record (SwlogReader.windowsFor), the rest stays 0:
%     log mode 0: log samples 1-552 -> record 9-560 (short pulse),
%                 log samples 553-7140 -> record 1298-7885 (long pulse)
%     log mode 1: log samples 1-6660 -> record 7-6666 (one window)
%     log mode 2: log samples 1-416 -> record 7-422 (short pulse),
%                 log samples 417-2910 -> record 993-3486 (long pulse)
%     log mode 3: log samples 1-2660 -> record 7-2666 (one window)
%   The values are rounded (and saturated) to int16 as in SwlogReader
%   (int32(round(x))), after the factor 'scale'.
%
%   SWLOG_TO_DDC(logPath, outFolder, name, value, ...) changes the folder
%   ([] = default) and the options:
%     'pri_s'        1e-3    PRI [s]; scalar or [1 x 4] for the log modes 0 ... 3
%     'ddcFs_Hz'     2e6     DDC output rate [Hz]; scalar or [1 x 4]
%     'ddcT0_s'      0       time of the first DDC sample in the PRI [s]; scalar or [1 x 4]
%     'oneWindow'    'long'  a log mode with one window (1, 3): it belongs to the
%                            'long' or the 'short' sub-pulse (hdr.longDdcIdx / shortDdcIdx)
%     'channels'     1:4     channels of the log that are saved (columns of I, Q)
%     'scale'        1       factor applied before the conversion to int16
%     'staggerIndex' 1       hdr.staggerIndex of every PRI
%     'skipZeroPrf'  true    skip records with PRF counter 0 (as SwlogReader.sample)
%     'firstRecord'  1       first record of the log that is converted
%     'maxPulses'    inf     stop after this many saved PRIs
%     'matFormat'    '-v6'   save format (cfg.io.matFormat; the C programs read -v6)
%     'utcOffset_min' 210    time zone of the time in the log name (SwlogReader.startUsFromName)
%     'verbose'      true
%   info = SWLOG_TO_DDC(...) returns the folder, the number of records read
%   and saved, the log modes found, the number of scans, of saturated samples
%   and of PRF counter jumps backwards.
%
%   Header of every PRI:
%     pulseNumber    PRF counter - first PRF counter + 1 (gaps in the counter stay
%                    gaps; a jump backwards of the counter is bridged, so the
%                    numbers always grow). Also the file number: ddc_<n>.mat
%     scanNumber     1 + number of azimuth wraps (the azimuth drops by more than 180 deg)
%     time_s         PRI start time: the counter steps times the PRI, summed
%     antennaAz_deg  azimuth * 360 / 16384
%     pri_s, staggerIndex, ddcFs_Hz, ddcT0_s   from the options
%     shortDdcIdx,   [first last] DDC sample (1-based) of the receive windows
%     longDdcIdx     ([0 0] = no window)
%     fs_Hz, fIF_Hz, fc_Hz, rxStart_s, shortTx_s, longTx_s, shortFc_Hz, longFc_Hz
%                    NaN (not in the log; the receiver takes them from its
%                    configuration, radar_config.bin)
%     logMode        log mode (0 ... 3)
%     logRecord      record number in the log (1-based)
%     logPrfCounter  PRF counter of the log
%     logStartUs     start time of the log from its name [us since 1970 UTC], 0 if unknown
%   The receiver accepts a PRI only when its configuration matches: the
%   number of channels, ddcFs_Hz, at least nDdc samples and the receive
%   windows shortDdcIdx / longDdcIdx.
%
%   Example:
%     info = swlog_to_ddc('H:\log\SoftDcdLog20250930_1718.swlog', [], ...
%                         'pri_s', 1e-3, 'ddcFs_Hz', 2e6, 'maxPulses', 5000);
%     d = load(fullfile(info.folder, 'ddc_0000001.mat'));   % d.I, d.Q, d.hdr
%     replay_ddc_udp(info.folder);

root = fileparts(fileparts(mfilename('fullpath')));    % project folder (this file is in src/)
addpath(root, fullfile(root, 'src'));
opt = struct('pri_s', 1e-3, 'ddcFs_Hz', 2e6, 'ddcT0_s', 0, 'oneWindow', 'long', 'channels', 1:4, ...
             'scale', 1, 'staggerIndex', 1, 'skipZeroPrf', true, 'firstRecord', 1, 'maxPulses', inf, ...
             'matFormat', '-v6', 'utcOffset_min', 210, 'verbose', true);
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i}), error('swlog_to_ddc: unknown option ''%s''', varargin{i}); end
    opt.(varargin{i}) = varargin{i+1};
end
if ~any(strcmp(opt.oneWindow, {'long', 'short'}))
    error('swlog_to_ddc: ''oneWindow'' must be ''long'' or ''short''.');
end

% log format (SwlogReader.m)
headerBytes = 84;
modeOffset  = 10;
nChLog      = 4;
logSamples  = [7200, 6660, 2910, 2660];      % samples per channel in the log, log mode 0 ... 3
recSamples  = [9317, 9387, 4937, 3059];      % DDC record (SwlogReader RawSize)
winTable = {[0, 552, 8; 552, 6588, 1297], ... % [from count to] (0-based): log samples from+1 ... from+count
            [0, 6660, 6; 0, 0, 0], ...        % -> record samples to+1 ... to+count
            [0, 416, 6; 416, 2494, 992], ...
            [0, 2660, 6; 0, 0, 0]};
ch = opt.channels(:).';
if isempty(ch) || any(ch < 1 | ch > nChLog | ch ~= round(ch))
    error('swlog_to_ddc: ''channels'' must be channel numbers 1 ... %d.', nChLog);
end

logPath = char(logPath);
d = dir(logPath);
if isempty(d), error('swlog_to_ddc: cannot find %s', logPath); end
fileBytes = d(1).bytes;
if nargin < 2 || isempty(outFolder)
    [~, stem] = fileparts(logPath);
    outFolder = fullfile(root, 'output', ['ddc_pulses_' regexprep(stem, '[^\w]', '_')]);
end
if ~exist(outFolder, 'dir'), mkdir(outFolder); end
fid = fopen(logPath, 'rb', 'ieee-le');
if fid < 0, error('swlog_to_ddc: cannot open %s', logPath); end
closer = onCleanup(@() fclose(fid));
startUs = start_us_from_name(logPath, opt.utcOffset_min);

nRead = 0; nSaved = 0; nClipped = 0; nBack = 0;
modes = [];
scan = 1; prevAz = [];
prevCount = []; pulseNo = 0; t = 0;
truncated = false; badMode = false;
lastPrint = tic;
while nSaved < opt.maxPulses
    head = fread(fid, headerBytes, 'uint8=>uint8');
    if isempty(head), break; end
    if numel(head) ~= headerBytes, truncated = true; break; end
    m = double(head(modeOffset + 1));
    if m >= numel(logSamples)
        badMode = true;                      % SwlogReader stops here too
        break;
    end
    n = logSamples(m + 1);
    nRead = nRead + 1;
    if nRead < opt.firstRecord                % skip without reading the samples
        if ftell(fid) + n*2*nChLog*4 > fileBytes || fseek(fid, n*2*nChLog*4, 'cof') ~= 0
            truncated = true; nRead = nRead - 1; break;
        end
        continue;
    end
    x = fread(fid, n*2*nChLog, 'float32=>double');
    if numel(x) ~= n*2*nChLog, truncated = true; nRead = nRead - 1; break; end
    count = double(typecast(head(1:4), 'uint32'));
    if opt.skipZeroPrf && count == 0, continue; end
    modes = union(modes, m);

    % pulse number, time, scan, azimuth
    pri = per_mode(opt.pri_s, m, 'pri_s');
    if isempty(prevCount)
        step = 1;
    else
        step = count - prevCount;
        if step <= 0                          % counter jumped back (reset): next number
            step = 1; nBack = nBack + 1;
        end
        t = t + step * pri;
    end
    prevCount = count;
    pulseNo = pulseNo + step;
    az = mod(double(typecast(head(5:6), 'uint16')), 16384) * 360 / 16384;
    if ~isempty(prevAz) && az < prevAz - 180, scan = scan + 1; end
    prevAz = az;

    % samples: [n x (I,Q) x channel] -> record [nDdc x channel], windows at their place
    x = reshape(x, n, 2, nChLog);
    x = x(:, :, ch) * opt.scale;
    w = winTable{m + 1};
    w = w(w(:, 2) > 0, :);
    nDdc = recSamples(m + 1);
    Id = zeros(nDdc, numel(ch));
    Qd = zeros(nDdc, numel(ch));
    for r = 1:size(w, 1)
        src = w(r, 1) + (1:w(r, 2));
        dst = w(r, 3) + (1:w(r, 2));
        Id(dst, :) = reshape(x(src, 1, :), [], numel(ch));
        Qd(dst, :) = reshape(x(src, 2, :), [], numel(ch));
    end
    Id = round(Id); Qd = round(Qd);
    nClipped = nClipped + sum(abs(Id(:)) > 32767 | Id(:) < -32768) + sum(abs(Qd(:)) > 32767 | Qd(:) < -32768);
    I = int16(Id);                            % saturates
    Q = int16(Qd);
    [sIdx, lIdx] = windows(w, opt.oneWindow);

    hdr = struct();
    hdr.pulseNumber   = pulseNo;
    hdr.scanNumber    = scan;
    hdr.time_s        = t;
    hdr.antennaAz_deg = az;
    hdr.pri_s         = pri;
    hdr.staggerIndex  = opt.staggerIndex;
    hdr.fs_Hz         = NaN;                  % not in the log
    hdr.fIF_Hz        = NaN;
    hdr.fc_Hz         = NaN;
    hdr.rxStart_s     = NaN;
    hdr.shortTx_s     = [NaN NaN];
    hdr.longTx_s      = [NaN NaN];
    hdr.ddcFs_Hz      = per_mode(opt.ddcFs_Hz, m, 'ddcFs_Hz');
    hdr.ddcT0_s       = per_mode(opt.ddcT0_s, m, 'ddcT0_s');
    hdr.shortDdcIdx   = sIdx;
    hdr.longDdcIdx    = lIdx;
    hdr.shortFc_Hz    = NaN;
    hdr.longFc_Hz     = NaN;
    hdr.logMode       = m;
    hdr.logRecord     = nRead;
    hdr.logPrfCounter = count;
    hdr.logStartUs    = startUs;
    save_ddc_file(outFolder, opt.matFormat, hdr, I, Q);
    nSaved = nSaved + 1;

    if opt.verbose && toc(lastPrint) > 2
        fprintf('  record %d of ~%d, %d PRIs saved, scan %d, az %.1f deg\n', ...
                nRead, round(fileBytes / record_bytes(logSamples(m + 1))), nSaved, scan, az);
        lastPrint = tic;
    end
end
if truncated
    warning('swlog_to_ddc: the last record of %s is incomplete and was skipped.', logPath);
end
if badMode
    warning('swlog_to_ddc: record %d of %s has log mode %d (only 0 ... 3 exist): stopped there.', ...
            nRead + 1, logPath, m);
end
if nClipped > 0
    warning('swlog_to_ddc: %d samples were outside int16 and saturated (use ''scale'').', nClipped);
end
if nBack > 0
    warning('swlog_to_ddc: the PRF counter jumped backwards %d times (bridged in pulseNumber).', nBack);
end
info = struct('folder', outFolder, 'nRecords', nRead, 'nSaved', nSaved, 'logModes', modes, ...
              'scans', scan, 'nClipped', nClipped, 'nCounterBack', nBack, 'logStartUs', startUs);
if opt.verbose
    fprintf('swlog_to_ddc: %d of %d records saved to %s (log modes %s, %d scans)\n', ...
            nSaved, nRead, outFolder, mat2str(modes), scan);
end
end

function b = record_bytes(n)
b = 84 + n * 2 * 4 * 4;
end

function v = per_mode(x, m, name)
if isscalar(x), v = x; return; end
if numel(x) < m + 1
    error('swlog_to_ddc: ''%s'' has %d values, log mode %d needs %d.', name, numel(x), m, m + 1);
end
v = x(m + 1);
end

function [sIdx, lIdx] = windows(w, oneWindow)
% rows [from count to] of the windows: record samples to+1 ... to+count (1-based)
iv = @(r) [w(r, 3) + 1, w(r, 3) + w(r, 2)];
sIdx = [0 0]; lIdx = [0 0];
if size(w, 1) >= 2
    sIdx = iv(1); lIdx = iv(2);
elseif size(w, 1) == 1
    if strcmp(oneWindow, 'short'), sIdx = iv(1); else, lIdx = iv(1); end
end
end

function us = start_us_from_name(path, utcOffset_min)
% SwlogReader.startUsFromName: ...SoftDcdLogYYYYMMDD_HHMM... (local time) -> us since 1970 UTC
us = 0;
[~, stem] = fileparts(path);
tok = regexp(stem, 'SoftDcdLog(\d{8})_(\d{4})', 'tokens', 'once');
if isempty(tok), return; end
v = str2double({tok{1}(1:4), tok{1}(5:6), tok{1}(7:8), tok{2}(1:2), tok{2}(3:4)});
if any(isnan(v)) || v(2) < 1 || v(3) < 1, return; end
days = datenum(v(1), v(2), v(3), v(4), v(5), 0) - datenum(1970, 1, 1);
us = (round(days * 86400) - round(utcOffset_min * 60)) * 1e6;
end

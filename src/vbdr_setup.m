function v = vbdr_setup(cfg, nDdc, nChUsed)
%VBDR_SETUP  Layout of the VBDR stream of the mode and the board settings, with their checks.
%   v = VBDR_SETUP(cfg, nDdc, nChUsed) (setup_radar.m: sys.vbdr) for the DDC record of nDdc
%   samples per PRI. Fields:
%     format        'rddc' or 'vbdr' (cfg.udp.format)
%     payloadBytes  bytes of the half-stream per packet (cfg.udp.vbdr.payloadBytes)
%     nPackets      packets per stream (both halves; cfg.udp.vbdr.nPackets or computed)
%     halfBytes     bytes of a half-stream: 64 (32 uint16 parameters) + 8*nDdc (2 channels x I, Q int16)
%     fcMin_Hz, fcMax_Hz, f1, f2   carrier range and the F1 / F2 register values (carrier - fcMin
%                   in MHz; NaN when the carrier is not a whole MHz in the range or the sub-pulse is off)
%     stcOptions, stc, gainOptions, gain   board settings of the mode (1 x 4: Main, J1, J2, J3)
%     nOperationModes, operationMode       (cfg.board)
%   With cfg.udp.format = 'vbdr' the stream carries 4 channels (the receiver processes at most
%   these) and the carriers of the enabled sub-pulses must be whole MHz in [fcMin_Hz, fcMax_Hz]
%   (at most 255 MHz above fcMin_Hz). README "VBDR stream and board control".
v.format = 'rddc';
if isfield(cfg, 'udp') && isfield(cfg.udp, 'format') && ~isempty(cfg.udp.format)
    v.format = lower(cfg.udp.format);
end
if ~ischar(v.format) || ~any(strcmp(v.format, {'rddc', 'vbdr'}))
    error('cfg.udp.format must be ''rddc'' or ''vbdr''.');
end
isVbdr = strcmp(v.format, 'vbdr');

% layout of the stream
v.payloadBytes = 1400;
nPk = [];
if isfield(cfg, 'udp') && isfield(cfg.udp, 'vbdr')
    if isfield(cfg.udp.vbdr, 'payloadBytes') && ~isempty(cfg.udp.vbdr.payloadBytes)
        v.payloadBytes = cfg.udp.vbdr.payloadBytes;
    end
    if isfield(cfg.udp.vbdr, 'nPackets'), nPk = cfg.udp.vbdr.nPackets; end
end
if ~isscalar(v.payloadBytes) || ~isreal(v.payloadBytes) || v.payloadBytes ~= round(v.payloadBytes) || ...
        v.payloadBytes < 2 || v.payloadBytes > 65479
    error('cfg.udp.vbdr.payloadBytes must be a whole number 2 ... 65479 (bytes per packet after the 28-byte label).');
end
v.halfBytes = 64 + 8*nDdc;
v.nPackets = 2*ceil(v.halfBytes / v.payloadBytes);
if ~isempty(nPk) && isVbdr && ~isequal(nPk, v.nPackets)
    error(['cfg.udp.vbdr.nPackets = %s does not fit the DDC record of the mode: %d samples per PRI give ', ...
           '2 x %d bytes, %d packets of %d bytes (2*ceil((64 + 8*nDdc)/payloadBytes)); set nPackets = %d ', ...
           '(or [] = computed).'], mat2str(nPk), nDdc, v.halfBytes, v.nPackets, v.payloadBytes, v.nPackets);
end
if v.nPackets > 65535
    error('cfg.udp.vbdr.payloadBytes = %d gives %d packets per stream: use a larger payload.', v.payloadBytes, v.nPackets);
end

% carriers (F1 = short pulse, F2 = long pulse)
v.fcMin_Hz = field_or(cfg, 'fcMin_Hz', 150e6);
v.fcMax_Hz = field_or(cfg, 'fcMax_Hz', 190e6);
if ~isscalar(v.fcMin_Hz) || ~isscalar(v.fcMax_Hz) || ~(v.fcMin_Hz > 0) || ~(v.fcMax_Hz >= v.fcMin_Hz) || ...
        mod(v.fcMin_Hz, 1e6) ~= 0 || mod(v.fcMax_Hz, 1e6) ~= 0 || v.fcMax_Hz - v.fcMin_Hz > 255e6
    error(['cfg.fcMin_Hz and cfg.fcMax_Hz must be whole MHz, fcMin_Hz <= fcMax_Hz <= fcMin_Hz + 255 MHz ', ...
           '(the F1 / F2 registers hold carrier - fcMin in MHz, one byte).']);
end
names = {'short', 'long'};
regs = {'F1', 'F2'};
f = [NaN NaN];
for i = 1:2
    wv = cfg.wave.(names{i});
    if isfield(wv, 'enable') && ~wv.enable, continue; end
    fc = cfg.fc;
    if isfield(wv, 'fc_Hz') && ~isempty(wv.fc_Hz), fc = wv.fc_Hz; end
    ok = mod(fc, 1e6) == 0 && fc >= v.fcMin_Hz && fc <= v.fcMax_Hz;
    if ok
        f(i) = (fc - v.fcMin_Hz) / 1e6;
    elseif isVbdr
        error(['The %s-pulse carrier %.6g MHz (%s) is not a whole MHz in cfg.fcMin_Hz ... cfg.fcMax_Hz = ', ...
               '%g ... %g MHz (the VBDR board sets it with %s).'], names{i}, fc/1e6, ...
              ['cfg.wave.' names{i} '.fc_Hz or cfg.fc'], v.fcMin_Hz/1e6, v.fcMax_Hz/1e6, regs{i});
    end
end
v.f1 = f(1);
v.f2 = f(2);

% board settings of the mode and the operation mode
b = struct();
if isfield(cfg, 'board'), b = cfg.board; end
v.stcOptions  = options(b, 'stcOptions',  'cfg.board.stcOptions');
v.stc         = values(b, 'stc', v.stcOptions, 'cfg.board.stc', 'stcOptions');
v.gainOptions = options(b, 'gainOptions', 'cfg.board.gainOptions');
v.gain        = values(b, 'gain', v.gainOptions, 'cfg.board.gain', 'gainOptions');
v.nOperationModes = field_or(b, 'nOperationModes', 1);
v.operationMode   = field_or(b, 'operationMode', 1);
if ~isscalar(v.nOperationModes) || v.nOperationModes ~= round(v.nOperationModes) || v.nOperationModes < 1 || ...
        v.nOperationModes > 255
    error('cfg.board.nOperationModes must be a whole number 1 ... 255.');
end
if ~isscalar(v.operationMode) || v.operationMode ~= round(v.operationMode) || v.operationMode < 1 || ...
        v.operationMode > v.nOperationModes
    error('cfg.board.operationMode must be a whole number 1 ... cfg.board.nOperationModes = %d.', v.nOperationModes);
end

if isVbdr && nChUsed > 4
    error(['cfg.udp.format = ''vbdr'' carries 4 channels (CH1 main, CH2 ... CH4 = J1 ... J3): the receiver ', ...
           'can process at most 4 (cfg.rx.nChannelsUsed = %d; set cfg.rx.nChannelsUsed or cfg.ant.nElements ', ...
           '<= 4).'], nChUsed);
end
end

function x = field_or(s, f, def)
x = def;
if isfield(s, f) && ~isempty(s.(f)), x = s.(f); end
end

function n = options(b, f, path)
n = [4 4 4 4];
if isfield(b, f) && ~isempty(b.(f)), n = b.(f); end
if ~isnumeric(n) || numel(n) ~= 4 || any(n ~= round(n)) || any(n < 1) || any(n > 255)
    error('%s must hold 4 whole numbers 1 ... 255 (options of Main, J1, J2, J3).', path);
end
n = reshape(double(n), 1, 4);
end

function x = values(b, f, n, path, listName)
x = [1 1 1 1];
if isfield(b, f) && ~isempty(b.(f)), x = b.(f); end
if ~isnumeric(x) || numel(x) ~= 4 || any(x ~= round(x)) || any(x < 1) || any(reshape(x, 1, 4) > n)
    error('%s must hold 4 whole numbers 1 ... cfg.board.%s (Main, J1, J2, J3).', path, listName);
end
x = reshape(double(x), 1, 4);
end

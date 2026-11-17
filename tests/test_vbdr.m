function ok = test_vbdr()
%TEST_VBDR  VBDR stream (cfg.udp.format = 'vbdr') and the board commands.
%   1) udp_make_vbdr: label, packet sizes, system parameters and the interleaved channels
%   2) udp_send_pulse / udp_receive_pulse(rx, 'vbdr') over localhost: the 4 channels intact
%   3) setup_radar checks: nPackets, carriers in whole MHz in fcMin ... fcMax, at most 4
%      channels processed, the board options
%   4) board_command: writes are stored, reads return them, the answers of the slides
%   5) board_sim_open / board_sim_poll over localhost: a command gets its answer
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;

%% 1) packets
nS = 1770;
I = int16(randi([-3000 3000], nS, 3));            % 3 channels: CH4 is sent as zeros
Q = int16(randi([-3000 3000], nS, 3));
hdr = struct('pulseNumber', 70000, 'antennaAz_deg', 90, 'pri_s', 1e-3, 'prevPri_s', 1.2e-3);
vb = struct('format', 'vbdr', 'payloadBytes', 1400, 'nPackets', 22);
p = udp_make_vbdr(hdr, I, Q, vb, 3);
lens = cellfun(@numel, p);
ok = check(numel(p) == 22 && all(lens([1:10 12:21]) == 1428) && lens(11) == 28 + 14224 - 10*1400 && ...
           lens(22) == lens(11), 'udp_make_vbdr: 22 packets of 1428 bytes, the last of each half shorter') && ok;
lab = double(typecast(p{12}(1:12), 'uint32'));
ok = check(isequal(lab, [70000 12 22]) && strcmp(char(p{12}(13:18)), 'Loss=3') && all(p{12}(19:28) == 0), ...
           'label: stream number, packet number 1-based, packets per stream, "Loss=3" zero padded') && ok;
half1 = [p{1:11}];  half1 = half1(~ismember(1:numel(half1), reshape((0:10)'*1428 + (1:28), 1, [])));
par = double(typecast(half1(1:64), 'uint16'));
ok = check(par(1) == mod(70000, 65536) && par(2) == 16384 && par(4) == hex2dec('5555') && par(5) == 1200 && ...
           par(6) == 1000 && all(par([3 7:32]) == 0), ...
           'system parameters: PRF number, azimuth 90 deg = 16384, 0x5555, PRT 1200 / 1000 us') && ok;
d = typecast(half1(65:end), 'int16');
ok = check(isequal(d(1:4), [I(1,1) Q(1,1) I(1,2) Q(1,2)]) && isequal(d(5:8), [I(2,1) Q(2,1) I(2,2) Q(2,2)]), ...
           'half 1: CH1_I, CH1_Q, CH2_I, CH2_Q of every sample') && ok;
try
    udp_make_vbdr(hdr, I, Q, struct('payloadBytes', 1400, 'nPackets', 20), 0);
    ok = check(false, 'udp_make_vbdr refuses a wrong nPackets') && ok;
catch
    ok = check(true, 'udp_make_vbdr refuses a wrong nPackets') && ok;
end

%% 2) localhost
port = 50124;
try
    rx = udp_receiver_open(port, 2);
    tx = udp_sender_open(struct('host', '127.0.0.1', 'port', port, 'payloadBytes', 1400, 'format', 'vbdr'), ...
                         struct('format', 'vbdr', 'payloadBytes', 1000, 'nPackets', []));
    I2 = I(1:200, :);  Q2 = Q(1:200, :);
    tx = udp_send_pulse(tx, hdr, I2, Q2);
    [h, Ir, Qr] = udp_receive_pulse(rx, 'vbdr');
    udp_sender_close(tx);
    clear rx
    ok = check(isequal(Ir(:, 1:3), I2) && isequal(Qr(:, 1:3), Q2) && all(Ir(:, 4) == 0) && h.pulseNumber == 70000 && ...
               abs(h.antennaAz_deg - 90) < 1e-9 && abs(h.pri_s - 1e-3) < 1e-12 && abs(h.prevPri_s - 1.2e-3) < 1e-12 && ...
               h.loss == 0 && h.marker == hex2dec('5555') && h.nPackets == 2*ceil((64 + 8*200)/1000), ...
               'udp_send_pulse / udp_receive_pulse ''vbdr'' over localhost') && ok;
catch err
    fprintf('  localhost VBDR stream: SKIPPED (%s)\n', err.message);
end

%% 3) checks of setup_radar
cfg = radar_config();
cfg.udp.format = 'vbdr';
cfg.rx.nChannelsUsed = 4;                          % (48 elements: the input files stay as they are)
cfg.clutter.enable = false;
sys = setup_radar(cfg);
ok = check(sys.vbdr.nPackets == 22 && sys.vbdr.f1 == 20 && sys.vbdr.f2 == 20 && isequal(sys.vbdr.stcOptions, [4 4 4 4]), ...
           'setup_radar: 22 packets, F1 = F2 = 20 MHz above fcMin') && ok;
ok = expect_error(setfield_path(cfg, 'udp.vbdr.nPackets', 20), 'nPackets', 'a wrong nPackets is refused') && ok;
c = cfg;  c.wave.long.fc_Hz = 170.5e6;
ok = expect_error(c, 'whole MHz', 'a carrier of 170.5 MHz is refused (vbdr)') && ok;
c = cfg;  c.wave.short.fc_Hz = 210e6;
ok = expect_error(c, 'whole MHz', 'a carrier above fcMax is refused (vbdr)') && ok;
c = cfg;  c.udp.format = 'rddc';  c.wave.long.fc_Hz = 170.5e6;
try
    s = setup_radar(c);
    ok = check(isnan(s.vbdr.f2), 'rddc: any carrier (F2 NaN)') && ok;
catch err
    ok = check(false, ['rddc: any carrier: ' err.message]) && ok;
end
c = cfg;  c.rx.nChannelsUsed = 8;
ok = expect_error(c, '4 channels', 'vbdr with 8 channels processed is refused') && ok;
c.rx.nChannelsUsed = 1;
try
    setup_radar(c);
    ok = check(true, 'vbdr with 48 elements, 1 processed') && ok;
catch err
    ok = check(false, ['vbdr with 48 elements, 1 processed: ' err.message]) && ok;
end
ok = expect_error(setfield_path(cfg, 'board.stc', [1 5 1 1]), 'cfg.board.stc', 'an STC value above its options is refused') && ok;
ok = expect_error(setfield_path(cfg, 'board.operationMode', 9), 'operationMode', 'an operation mode above the options is refused') && ok;
ok = expect_error(setfield_path(cfg, 'udp.format', 'xyz'), 'udp.format', 'an unknown format is refused') && ok;

%% 4) board commands
regs = zeros(256, 1);
[r, regs] = board_command('CNTL W 0x02 0x1C', regs);
ok = check(strcmp(r, 'CNTL W 0x02 0x1C Returned Value = 0x00000000') && regs(3) == 28, 'write: answer of the slides') && ok;
r = board_command('CNTL R 0x02 0x00', regs);
ok = check(strcmp(r, 'CNTL R 0x02 0x00 Returned Value = 0x0000001C'), 'read: answer of the slides') && ok;
r = board_command('hello', regs);
ok = check(isempty(r), 'not a command: no answer') && ok;

%% 5) the simulator as the board over localhost
try
    cfgB = cfg;
    cfgB.board.port = 50125;  cfgB.board.replyPort = 50126;
    b = board_sim_open(cfgB, sys, false);
    rxB = udp_receiver_open(50126, 2);
    txB = udp_sender_open(struct('host', '127.0.0.1', 'port', 50125, 'payloadBytes', 100));
    udp_send_datagram(txB, uint8('CNTL R 0x00 0x00'));
    t0 = tic;
    while b.nCommands == 0 && toc(t0) < 2
        b = board_sim_poll(b);
        pause(0.01);
    end
    rep = char(udp_read_any(rxB));
    udp_sender_close(txB);
    board_sim_close(b);
    clear rxB
    ok = check(strcmp(rep, 'CNTL R 0x00 0x00 Returned Value = 0x00000000'), ...
               'board_sim_poll answers a read of Radar_Mode (mode 0)') && ok;
catch err
    fprintf('  board commands over localhost: SKIPPED (%s)\n', err.message);
end
fprintf('VBDR stream and board commands: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function pkt = udp_read_any(rx)
switch rx.kind
    case 'java'
        dp = javaObject('java.net.DatagramPacket', zeros(1, 2048, 'int8'), 2048);
        rx.obj.receive(dp);
        data = typecast(int8(dp.getData()), 'uint8');
        pkt = data(1:dp.getLength());
    case 'matlab'
        dg = read(rx.obj, 1, 'uint8');  pkt = uint8(dg.Data);
    otherwise
        pkt = uint8(read(rx.obj, rx.obj.NumBytesAvailable, 'uint8'));
end
pkt = pkt(:).';
end

function c = setfield_path(c, path, v)
f = strsplit(path, '.');
c = setfield(c, f{:}, v);
end

function ok = expect_error(cfg, pattern, what)
try
    setup_radar(cfg);
    ok = check(false, [what ' (no error)']);
catch err
    ok = check(~isempty(strfind(err.message, pattern)), [what ': ' err.message]);
end
end

function ok = check(cond, what)
ok = logical(cond);
fprintf('  %-4s %s\n', ternary(ok, 'ok', 'FAIL'), what);
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

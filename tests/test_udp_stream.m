function ok = test_udp_stream()
%TEST_UDP_STREAM  Send the DDC output of one PRI over UDP (localhost) and check it is received intact.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
port = 50123;
cfgUdp = struct('host', '127.0.0.1', 'port', port, 'payloadBytes', 1400);
try
    rx = udp_receiver_open(port, 2);
    tx = udp_sender_open(cfgUdp);
catch err
    fprintf('UDP stream test: SKIPPED (no UDP support: %s)\n', err.message);
    ok = true;
    return
end

% Small signed 16-bit DDC output (fits the receive socket buffer) and its header
I = int16(randi([-3000 3000], 100, 48));
Q = int16(randi([-3000 3000], 100, 48));
hdr = struct('pulseNumber', 1234, 'scanNumber', 2, 'antennaAz_deg', 123.5, 'staggerIndex', 3, ...
             'time_s', 12.345, 'ddcFs_Hz', 2e6, 'ddcT0_s', 28e-6, 'shortDdcIdx', [1 40], 'longDdcIdx', [70 100]);
udp_send_pulse(tx, hdr, I, Q);
[hdrRx, Irx, Qrx] = udp_receive_pulse(rx);
udp_sender_close(tx);
clear rx                                              % closes the udpport receiver

ok = isequal(Irx, I) && isequal(Qrx, Q) && hdrRx.pulseNumber == 1234 && hdrRx.scanNumber == 2 ...
     && abs(hdrRx.antennaAz_deg - 123.5) < 1e-4 && hdrRx.staggerIndex == 3 && hdrRx.time_s == 12.345 ...
     && hdrRx.ddcFs_Hz == 2e6 && isequal(hdrRx.shortDdcIdx, [1 40]) && isequal(hdrRx.longDdcIdx, [70 100]);
fprintf('UDP stream test (int16 DDC output, %d bytes in %d datagrams): %s\n', 4*numel(I), ...
        ceil(4*numel(I)/cfgUdp.payloadBytes), ternary(ok, 'PASS', 'FAIL'));
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

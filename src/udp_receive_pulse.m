function [hdr, I, Q] = udp_receive_pulse(rx)
%UDP_RECEIVE_PULSE  Receive and reassemble the DDC output of one PRI sent by udp_send_pulse.m.
%   [hdr, I, Q] = UDP_RECEIVE_PULSE(rx) reads datagrams from the receiver
%   rx (see udp_receiver_open.m) until all packets of one PRI have arrived.
%   Packets of another PRI restart the reassembly. I and Q are int16
%   [nSamples x nChannels]. This function is a reference implementation of
%   the stream format for a receiver.
bytes = [];
received = [];
nPackets = inf;
pulse = -1;
while numel(received) < nPackets
    pkt = read_datagram(rx);
    if isempty(pkt), error('udp_receive_pulse: timeout'); end
    if typecast(pkt(1:4), 'uint32') ~= hex2dec('52444443'), continue; end   % not ours
    w32a = double(typecast(pkt(5:12),  'uint32'));     % pulse, scan
    if w32a(1) ~= pulse                                % first packet of a new PRI
        pulse = w32a(1);
        bytes = [];
        received = [];
    end
    w16  = double(typecast(pkt(17:24), 'uint16'));     % packet index, count, channels, stagger
    w32b = double(typecast(pkt(25:32), 'uint32'));     % samples per channel, first sample
    wS   = double(typecast(pkt(41:48), 'single'));     % DDC rate, time of the first DDC sample
    w32c = double(typecast(pkt(49:64), 'uint32'));     % segments of the two sub-pulses
    hdr.pulseNumber   = w32a(1);
    hdr.scanNumber    = w32a(2);
    hdr.antennaAz_deg = double(typecast(pkt(13:16), 'single'));
    hdr.staggerIndex  = w16(4);
    hdr.time_s        = typecast(pkt(33:40), 'double');
    hdr.ddcFs_Hz      = wS(1);
    hdr.ddcT0_s       = wS(2);
    hdr.shortDdcIdx   = reshape(w32c(1:2), 1, 2);
    hdr.longDdcIdx    = reshape(w32c(3:4), 1, 2);
    nPackets  = w16(2);
    nChannels = w16(3);
    nSamples  = w32b(1);
    if isempty(bytes)
        bytes = zeros(1, 4*nSamples*nChannels, 'uint8');
    end
    payload = pkt(65:end);
    bytes(4*w32b(2) + (1:numel(payload))) = payload;
    received = unique([received, w16(1)]);
end
iq = reshape(typecast(bytes, 'int16'), 2, []);
I = reshape(iq(1, :), nSamples, nChannels);
Q = reshape(iq(2, :), nSamples, nChannels);
end

function pkt = read_datagram(rx)
% One datagram as a uint8 row vector (empty on timeout).
switch rx.kind
    case 'matlab'
        dg = read(rx.obj, 1, 'uint8');
        if isempty(dg), pkt = []; return; end
        pkt = uint8(dg.Data);
    case 'octave'
        t0 = tic;
        while rx.obj.NumBytesAvailable == 0 && toc(t0) < rx.timeout_s
            pause(0.001);
        end
        if rx.obj.NumBytesAvailable == 0, pkt = []; return; end
        pkt = uint8(read(rx.obj, rx.obj.NumBytesAvailable, 'uint8'));
    case 'java'
        dp = javaObject('java.net.DatagramPacket', zeros(1, 65536, 'int8'), 65536);
        try
            rx.obj.receive(dp);
        catch
            pkt = []; return;                          % timeout
        end
        data = typecast(int8(dp.getData()), 'uint8');
        pkt = data(1:dp.getLength());
end
pkt = pkt(:).';
end

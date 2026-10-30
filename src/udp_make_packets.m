function packets = udp_make_packets(hdr, I, Q, payloadBytes)
%UDP_MAKE_PACKETS  Split the DDC output of one PRI into UDP datagrams.
%   packets = UDP_MAKE_PACKETS(hdr, I, Q, payloadBytes) returns a cell array
%   of uint8 row vectors, one per datagram: a 64-byte header + up to
%   payloadBytes bytes of the signed 16-bit DDC output I, Q [nSamples x
%   nChannels] of the PRI described by hdr (see
%   pulse_header.m). Header (little-endian):
%     bytes  0-3   uint32  magic 0x52444443 ('RDDC')
%     bytes  4-7   uint32  pulse number
%     bytes  8-11  uint32  scan number
%     bytes 12-15  single  antenna azimuth [deg]
%     bytes 16-19  uint16  packet index (0-based), packet count
%     bytes 20-23  uint16  number of channels, stagger index
%     bytes 24-31  uint32  samples per channel, first sample in this packet (0-based)
%     bytes 32-39  double  PRI start time [s]
%     bytes 40-47  single  DDC sample rate [Hz], time of the first DDC sample in the PRI [s]
%     bytes 48-63  uint32  short pulse first/last DDC sample, long pulse first/last DDC sample (1-based)
%   Payload: int16 pairs (I, Q), channel after channel (column-major order).
iq = [reshape(int16(I), 1, []); reshape(int16(Q), 1, [])];
bytes = typecast(iq(:).', 'uint8');
perPacket = 4*floor(payloadBytes/4);                 % whole complex samples
nPackets = ceil(numel(bytes)/perPacket);
if perPacket < 4 || nPackets > 65535
    error(['udp_make_packets: %d payload bytes per datagram give %d datagrams per PRI (at most 65535): ' ...
           'use a larger payload.'], payloadBytes, nPackets);
end
% header fields that are the same in every packet of the PRI
headA = [typecast(uint32(hex2dec('52444443')), 'uint8'), ...
         typecast(uint32([hdr.pulseNumber, hdr.scanNumber]), 'uint8'), ...
         typecast(single(hdr.antennaAz_deg), 'uint8')];
headB = [typecast(double(hdr.time_s), 'uint8'), ...
         typecast(single([hdr.ddcFs_Hz, hdr.ddcT0_s]), 'uint8'), ...
         typecast(uint32([hdr.shortDdcIdx, hdr.longDdcIdx]), 'uint8')];
packets = cell(1, nPackets);
for i = 1:nPackets
    first = (i-1)*perPacket + 1;
    last  = min(i*perPacket, numel(bytes));
    packets{i} = [headA, ...
                  typecast(uint16([i-1, nPackets, size(I, 2), hdr.staggerIndex]), 'uint8'), ...
                  typecast(uint32([size(I, 1), (first-1)/4]), 'uint8'), ...
                  headB, bytes(first:last)];
end
end

function udp_send_pulse(u, hdr, I, Q)
%UDP_SEND_PULSE  Send the DDC output of one PRI as a sequence of UDP datagrams.
%   UDP_SEND_PULSE(u, hdr, I, Q) sends the signed 16-bit DDC output I, Q
%   [nSamples x nChannels] of the PRI described by hdr (see
%   pulse_header.m). Every datagram = 64-byte header + up to u.payloadBytes
%   bytes of samples; the format is described in udp_make_packets.m.
udp_send_packets(u, udp_make_packets(hdr, I, Q, u.payloadBytes));
end

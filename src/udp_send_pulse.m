function u = udp_send_pulse(u, hdr, I, Q)
%UDP_SEND_PULSE  Send the DDC output of one PRI as a sequence of UDP datagrams.
%   u = UDP_SEND_PULSE(u, hdr, I, Q) sends the signed 16-bit DDC output I, Q
%   [nSamples x nChannels] of the PRI described by hdr (see pulse_header.m) with the
%   sender u (udp_sender_open.m), in its format: 'rddc' = every datagram a 64-byte header +
%   up to u.payloadBytes bytes of samples (udp_make_packets.m); 'vbdr' = the stream of the
%   VBDR board (udp_make_vbdr.m), the first 4 channels. u.loss counts the packets that
%   could not be sent (VBDR: the "Loss=" of the next streams).
if isfield(u, 'format') && strcmp(u.format, 'vbdr')
    packets = udp_make_vbdr(hdr, I, Q, u.vbdr, u.loss);
    for i = 1:numel(packets)
        try
            udp_send_datagram(u, packets{i});
        catch
            u.loss = u.loss + 1;
        end
    end
else
    udp_send_packets(u, udp_make_packets(hdr, I, Q, u.payloadBytes));
end
end

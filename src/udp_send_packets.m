function udp_send_packets(u, packets)
%UDP_SEND_PACKETS  Send datagrams made by udp_make_packets.m with the sender u (udp_sender_open.m).
for i = 1:numel(packets)
    udp_send_datagram(u, packets{i});
end
end

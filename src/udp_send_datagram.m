function udp_send_datagram(u, packet)
%UDP_SEND_DATAGRAM  Send one datagram (uint8 row vector) with the sender u (udp_sender_open.m).
if strcmp(u.kind, 'udpport')
    write(u.obj, packet, 'uint8', u.host, u.port);
else
    dp = javaObject('java.net.DatagramPacket', typecast(packet, 'int8'), numel(packet), u.addr, u.port);
    u.obj.send(dp);
end
end

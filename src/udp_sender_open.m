function u = udp_sender_open(cfgUdp)
%UDP_SENDER_OPEN  Open a UDP sender for streaming the DDC output.
%   Uses udpport (MATLAB R2020b+, or the Octave instrument-control package)
%   and falls back to a Java DatagramSocket.
u = cfgUdp;
datagramSize = 64 + u.payloadBytes;                  % header + payload of one datagram
try
    u.obj  = udpport('datagram');                    % MATLAB
    % every write must go out as ONE datagram (MATLAB's default limit is 512 bytes)
    u.obj.OutputDatagramSize = datagramSize;
    u.kind = 'udpport';
catch
    try
        pkg('load', 'instrument-control');           % Octave
        u.obj  = udpport();
        u.kind = 'udpport';
    catch
        u.obj  = javaObject('java.net.DatagramSocket');
        u.addr = javaMethod('getByName', 'java.net.InetAddress', u.host);
        u.kind = 'java';
    end
end
end

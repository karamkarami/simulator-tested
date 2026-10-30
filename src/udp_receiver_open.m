function rx = udp_receiver_open(port, timeout_s)
%UDP_RECEIVER_OPEN  Open a UDP receiver (used by udp_receive_pulse.m).
%   Uses udpport (MATLAB R2020b+, or the Octave instrument-control package)
%   and falls back to a Java DatagramSocket.
try
    rx.obj  = udpport('datagram', 'LocalPort', port, 'Timeout', timeout_s);   % MATLAB
    rx.kind = 'matlab';
catch
    try
        pkg('load', 'instrument-control');                                     % Octave
        rx.obj  = udpport('LocalPort', port, 'Timeout', timeout_s);
        rx.kind = 'octave';
    catch
        rx.obj  = javaObject('java.net.DatagramSocket', port);
        rx.obj.setSoTimeout(round(1000*timeout_s));
        rx.kind = 'java';
    end
end
rx.timeout_s = timeout_s;
end

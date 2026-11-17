function u = udp_sender_open(cfgUdp, vbdr)
%UDP_SENDER_OPEN  Open a UDP sender for streaming the DDC output.
%   u = UDP_SENDER_OPEN(cfg.udp) sends in the format cfg.udp.format ('rddc' when missing);
%   u = UDP_SENDER_OPEN(cfg.udp, sys.vbdr) is needed for the format 'vbdr' (payload and
%   packets of the stream of the mode, setup_radar.m). u.loss counts the packets that could
%   not be sent (VBDR label "Loss=").
%   Uses udpport (MATLAB R2020b+, or the Octave instrument-control package)
%   and falls back to a Java DatagramSocket.
u = cfgUdp;
if ~isfield(u, 'format') || isempty(u.format), u.format = 'rddc'; end
u.format = lower(u.format);
if strcmp(u.format, 'vbdr')
    if nargin < 2 || isempty(vbdr)
        error('udp_sender_open: the format ''vbdr'' needs the layout of the stream (sys.vbdr of setup_radar.m).');
    end
    u.vbdr = vbdr;
    datagramSize = 28 + vbdr.payloadBytes;           % label + payload of one packet
else
    datagramSize = 64 + u.payloadBytes;              % header + payload of one datagram
end
u.loss = 0;
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

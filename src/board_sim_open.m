function b = board_sim_open(cfg, sys, verbose)
%BOARD_SIM_OPEN  The simulator answers the commands of the receiver like the VBDR board.
%   b = BOARD_SIM_OPEN(cfg, sys) listens on UDP port cfg.board.port; board_sim_poll.m reads
%   the commands that arrived, prints them, stores the written values and answers to
%   cfg.board.pcHost:cfg.board.replyPort (board_command.m). The registers start with the
%   settings of the simulation (board_registers.m). The simulation itself does not change
%   (rerun it in the new mode). b = [] (with a warning) when the port cannot be opened.
%   Close with board_sim_close(b).
if nargin < 3, verbose = true; end
b = [];
if ~isfield(cfg, 'board'), return; end
b.regs = board_registers(cfg, sys);
b.verbose = verbose;
b.replyHost = cfg.board.pcHost;
b.replyPort = cfg.board.replyPort;
b.nCommands = 0;
try
    b.obj = udpport('datagram', 'LocalPort', cfg.board.port);              % MATLAB
    b.kind = 'matlab';
catch
    try
        ch = javaMethod('open', 'java.nio.channels.DatagramChannel');      % Octave: Java, not blocking
        ch.configureBlocking(false);
        ch.socket().bind(javaObject('java.net.InetSocketAddress', cfg.board.port));
        b.obj = ch;
        b.buf = javaMethod('allocate', 'java.nio.ByteBuffer', 4096);
        b.addr = javaMethod('getByName', 'java.net.InetAddress', b.replyHost);
        b.kind = 'java';
    catch err
        warning('board_sim_open: cannot listen for the board commands on UDP port %d: %s', cfg.board.port, err.message);
        b = [];
        return
    end
end
if verbose
    fprintf('Board commands: listening on UDP port %d, answers to %s:%d\n', cfg.board.port, b.replyHost, b.replyPort);
end
end

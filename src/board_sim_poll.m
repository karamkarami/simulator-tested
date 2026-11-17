function b = board_sim_poll(b)
%BOARD_SIM_POLL  Read and answer the board commands that arrived (board_sim_open.m); never waits.
if isempty(b), return; end
for guard = 1:64
    line = '';
    switch b.kind
        case 'matlab'
            if b.obj.NumDatagramsAvailable == 0, break; end
            dg = read(b.obj, 1, 'uint8');
            line = char(uint8(dg.Data(:).'));
        case 'java'
            b.buf.clear();
            from = b.obj.receive(b.buf);
            if isempty(from), break; end
            n = b.buf.position();
            arr = b.buf.array();
            raw = typecast(int8(arr(1:n)), 'uint8');
            line = char(raw(:).');
    end
    [reply, b.regs, info] = board_command(line, b.regs);
    b.nCommands = b.nCommands + 1;
    if b.verbose, fprintf('board command: %s\n', info); end
    if isempty(reply), continue; end
    switch b.kind
        case 'matlab'
            write(b.obj, uint8(reply), 'uint8', b.replyHost, b.replyPort);
        case 'java'
            msg = typecast(uint8(reply), 'int8');
            bb = javaMethod('wrap', 'java.nio.ByteBuffer', msg);
            b.obj.send(bb, javaObject('java.net.InetSocketAddress', b.replyHost, b.replyPort));
    end
end
end

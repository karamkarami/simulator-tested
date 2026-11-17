function [reply, regs, info] = board_command(line, regs)
%BOARD_COMMAND  Answer of the VBDR board to one ASCII command (the simulator acts as the board).
%   [reply, regs, info] = BOARD_COMMAND(line, regs): line = "CNTL W 0xAA 0xVV" (write value VV
%   to address AA) or "CNTL R 0xAA 0x00" (read address AA); regs (256 x 1) holds the register
%   values (regs(AA+1)). reply = the command followed by " Returned Value = 0x%08X" (write: 0,
%   read: the value), '' for a line that is not a command; info = a text for the console, e.g.
%   "write Radar_Mode (0x00) = 1". README "VBDR stream and board control".
reply = '';
info = '';
tok = regexp(strtrim(line), '^CNTL\s+([WRwr])\s+0[xX]([0-9A-Fa-f]{1,2})\s+0[xX]([0-9A-Fa-f]{1,8})$', 'tokens', 'once');
if isempty(tok)
    info = sprintf('unknown command "%s"', strtrim(line));
    return
end
op = upper(tok{1});
addr = hex2dec(tok{2});
val = hex2dec(tok{3});
cmd = sprintf('CNTL %s 0x%02X 0x%02X', op, addr, mod(val, 256));
if op == 'W'
    regs(addr+1) = mod(val, 256);
    reply = sprintf('%s Returned Value = 0x%08X', cmd, 0);
    info = sprintf('write %s (0x%02X) = %d', board_register_name(addr), addr, mod(val, 256));
else
    reply = sprintf('%s Returned Value = 0x%08X', cmd, regs(addr+1));
    info = sprintf('read  %s (0x%02X) -> %d', board_register_name(addr), addr, regs(addr+1));
end
end

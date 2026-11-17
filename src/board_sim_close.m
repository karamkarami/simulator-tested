function board_sim_close(b)
%BOARD_SIM_CLOSE  Stop answering the board commands (board_sim_open.m).
if isempty(b), return; end
try
    if strcmp(b.kind, 'java'), b.obj.close(); else, delete(b.obj); end
catch
end
end

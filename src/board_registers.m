function regs = board_registers(cfg, sys)
%BOARD_REGISTERS  Register values of the VBDR board for the settings of the simulator (setup_radar.m).
%   regs (256 x 1): regs(addr+1) = value of address addr. The simulator answers reads with
%   them until the receiver writes new values (board_sim_poll.m).
regs = zeros(256, 1);
regs(0+1) = cfg.modes.active;
code = max([sys.wave.short.codeIndex, sys.wave.long.codeIndex, 0]);
regs(1+1) = code;
if ~isnan(sys.vbdr.f1), regs(2+1) = sys.vbdr.f1; end
if ~isnan(sys.vbdr.f2), regs(3+1) = sys.vbdr.f2; end
regs(5+1) = sys.staggerSet;                      % 0: no staggered PRF
regs(40+1) = sys.vbdr.operationMode;
regs((41:44)+1) = sys.vbdr.stc;
regs((45:48)+1) = sys.vbdr.gain;
regs(50+1) = mod(round(sys.ant.rpm), 256);
end

function name = board_register_name(addr)
%BOARD_REGISTER_NAME  Name of a register of the VBDR board (README "VBDR stream and board control").
names = containers.Map('KeyType', 'double', 'ValueType', 'char');
names(0) = 'Radar_Mode';   names(1) = 'Code_index';  names(2) = 'F1_Freq';    names(3) = 'F2_Freq';
names(5) = 'PRF_index';    names(40) = 'Operation_Mode';
names(41) = 'STC_Main';    names(42) = 'STC_J1';     names(43) = 'STC_J2';    names(44) = 'STC_J3';
names(45) = 'Gain_Main';   names(46) = 'Gain_J1';    names(47) = 'Gain_J2';   names(48) = 'Gain_J3';
names(49) = 'Ext_Trig';    names(50) = 'RPM (simulator)';
if isKey(names, addr), name = names(addr); else, name = 'unknown'; end
end

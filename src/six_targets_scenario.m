function targets = six_targets_scenario()
%SIX_TARGETS_SCENARIO  Six targets around the radar with different RCS and kinematics.
%                    name  RCS[m2] range[m] az[deg] alt[m] vRadial[m/s] vCross[m/s]
targets(1) = make_target('T1', 1.0,  12e3,   30,   1000,  -237,    0);
targets(2) = make_target('T2', 0.5,  20e3,   95,   3000,   130,   40);
targets(3) = make_target('T3', 5.0,  35e3,  160,   6000,  -290,    0);
targets(4) = make_target('T4', 10,   42e3,  215,   8000,    80,  150);
targets(5) = make_target('T5', 2.0,   8e3,  270,    300,   -40,   20);
targets(6) = make_target('T6', 3.0,  46e3,  330,   2000,   210,  -50);
end

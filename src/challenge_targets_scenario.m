function targets = challenge_targets_scenario()
%CHALLENGE_TARGETS_SCENARIO  Ten targets that are hard for the 170 MHz radar.
%   Every target stands for one of the limits of the radar (see README
%   "170 MHz" and tests/test_vhf_challenges.m); some are expected to be missed
%   or detected only in some scans. With the terrain clutter of radar_config
%   (cfg.clutter.terrain = 'regions') some of them fly over its regions:
%   C1      0.3 m^2 at 48 km over the mountain range (North-East, 42 - 62 km):
%           low SNR inside the strongest clutter (MTI residue, sidelobes)
%   C2      slow (25 m/s radial) and low (300 m) at 8 km over the hills
%           (East-South-East), just outside the MTI notch (+-2 bins of 6.9 m/s with
%           the fixed PRF)
%   C3      tangential crosser: radial velocity 4 m/s, inside the MTI notch
%           (blind while it crosses), 220 m/s across the line of sight
%   C4, C5  a strong near target (10 m^2 at 5 km over the hills, receding,
%           short pulse, deep STC) and
%           a 1 m^2 target behind it on the same azimuth at 23 km: just behind
%           the long-pulse blind range (22.5 km), at the short/long-pulse
%           handover and where the long-window STC ends
%   C6, C7  two targets 300 m apart in range on the same azimuth (other
%           Doppler: 40 m/s apart; 525 m in slant range at the first look):
%           C7 lies in the 2-D CFAR window of the stronger C6 and is masked
%           at first (CA-CFAR), then resolved: two plots, two tracks
%   C8      fast: 420 m/s radial, near the unambiguous limit of the fixed PRF
%           (+-441 m/s; with the staggered PRF it folds)
%   C9      high elevation: 9 km altitude at 12 km (37 deg, receding: below
%           the highest DBF beam, 40 deg), over the sea
%   C10     low altitude far away: 100 m at 45 km (0.1 deg elevation) over the
%           sea, near the radar horizon (51 km for 100 m), at the lower edge
%           of the beams
%                    name   RCS[m2] range[m] az[deg] alt[m] vRadial[m/s] vCross[m/s]
targets(1)  = make_target('C1',  0.3,  48e3,    50,   4500,  -150,    50);
targets(2)  = make_target('C2',  2.0,   8e3,   125,    300,   -25,     4);
targets(3)  = make_target('C3',  5.0,  28e3,    90,   6000,     4,  -220);
targets(4)  = make_target('C4', 10,     5e3,   140,    400,   120,     0);
targets(5)  = make_target('C5',  1.0,  23e3,   140,   3000,   150,     0);
targets(6)  = make_target('C6',  2.0,  33e3,   185,   5000,  -200,     0);
targets(7)  = make_target('C7',  2.0,  33.3e3, 185,   5200,  -160,     0);
targets(8)  = make_target('C8',  3.0,  50e3,   235,   9000,  -420,     0);
targets(9)  = make_target('C9',  2.0,  12e3,   280,   9000,   150,   100);
targets(10) = make_target('C10', 5.0,  45e3,   325,    100,  -250,    30);
end

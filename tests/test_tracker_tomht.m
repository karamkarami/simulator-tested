function ok = test_tracker_tomht()
%TEST_TRACKER_TOMHT  TOMHT tracker with IMM-UKF (cfg.track.method = 'tomht') on synthetic plots.
%   The cases of test_tracker_pmbm with the track-oriented MHT (track trees
%   scored by their log-likelihood ratio, k-best global hypotheses by Murty's
%   algorithm, N-scan pruning) instead of the PMBM: six targets with false
%   plots (fixed and staggered PRF), a 4 deg/s turn and a 6 m/s^2 acceleration
%   (smaller errors than the EKF, the CT / CA model dominant), Pd 0.7 with 30
%   false plots per scan, two targets 200 m apart, two crossing targets (no
%   swap), track termination, speed gate, the 13 challenge targets. When the C
%   receiver is built, the C TOMHT tracker (tracker_tomht.c) runs on the same
%   calls of every case: every track update must match MATLAB.
ok = test_tracker_pmbm('tomht');
end

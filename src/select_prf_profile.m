function sys = select_prf_profile(sys, j)
%SELECT_PRF_PROFILE  Slow-time settings of the PRF j of the scan-to-scan PRF (cfg.prf.mode = 'scan2scan').
%   sys = SELECT_PRF_PROFILE(sys, j) puts the settings of sys.prfProfiles(j)
%   (setup_radar.m: PRI, CPI length, Doppler rows, FFT length, window, axes,
%   noise gains, notch) into sys, and the spans of the Doppler axis of the
%   plots (sys.track.vrSpanPulse, vrSpan). run_radar_sim.m calls it at the
%   first PRI of every scan whose PRF differs from that of the previous one.
q = sys.prfProfiles(j);
sys.pri = q.pri;
sys.cpi = q.cpi;
sys.dop = q.dop;
sys.wave.short.v_mps = q.v_mps{1};
sys.wave.long.v_mps  = q.v_mps{2};
sys.track.vrSpanPulse.short = sys.wave.short.lambda_m / (2 * q.dop.rowPri_s);
sys.track.vrSpanPulse.long  = sys.wave.long.lambda_m / (2 * q.dop.rowPri_s);
spans = [sys.track.vrSpanPulse.short, sys.track.vrSpanPulse.long];
sys.track.vrSpan = min(spans(sys.pulseOn));
sys.prfProfile = j;
end

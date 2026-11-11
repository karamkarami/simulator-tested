function hdr = pulse_header(sys, sched, k)
%PULSE_HEADER  Header of PRI k: stored with its raw data and DDC output, sent with its DDC output.
%   hdr = PULSE_HEADER(sys, sched, k) with the schedule of pulse_schedule.m.
hdr.pulseNumber    = k;
hdr.scanNumber     = sched.scan(k);
hdr.time_s         = sched.time_s(k);
hdr.antennaAz_deg  = sched.az_deg(k);
hdr.pri_s          = sched.pri_s(k);
hdr.staggerIndex   = sched.staggerIdx(k);
hdr.fs_Hz          = sys.fs;
hdr.fIF_Hz         = sys.fIF;
hdr.fc_Hz          = sys.cfg.fc;
hdr.rxStart_s      = sys.rxStart_s;                  % time of the first ADC sample in the PRI
hdr.shortTx_s      = [sys.wave.short.tTx_s, sys.wave.short.width_s];   % [start width]
hdr.longTx_s       = [sys.wave.long.tTx_s,  sys.wave.long.width_s];
hdr.ddcFs_Hz       = sys.fsOut;                  % DDC output: signed 16-bit I/Q
hdr.ddcT0_s        = sys.tDdc(1);                    % time of the first DDC sample in the PRI
hdr.shortDdcIdx    = interval(sys.wave.short);       % receive interval of the short pulse [first last]
hdr.longDdcIdx     = interval(sys.wave.long);        % (DDC samples, independent of CFAR and filters;
                                                     %  [0 0]: the sub-pulse is switched off)
hdr.shortFc_Hz     = sys.wave.short.fc_Hz;           % carriers of the sub-pulses (cfg.wave.<pulse>.fc_Hz)
hdr.longFc_Hz      = sys.wave.long.fc_Hz;
end

function iv = interval(wv)
iv = [0 0];
if wv.enable && ~isempty(wv.rxIdx), iv = wv.rxIdx([1 end]); end
end

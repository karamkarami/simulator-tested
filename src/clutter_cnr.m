function cnr = clutter_cnr(sys, sigma0, R_m, el_deg)
%CLUTTER_CNR  Clutter-to-noise ratio of uniform clutter at the matched filter output.
%   cnr = CLUTTER_CNR(sys, sigma0, R_m, el_deg) returns the linear CNR of
%   clutter of reflectivity sigma0 [m^2/m^2, linear] that covers the ground
%   around range R_m, arriving from elevation el_deg: the mean clutter power
%   over the noise power at the matched filter output of one element (short
%   pulse, one pulse, no STC, before the MTI), averaged over the azimuth of
%   the antenna. This is the definition of cfg.clutter.cnr_dB.
%   The clutter of the cells within the range response of the matched filter
%   adds up in power: per element
%     CNR = radarConst*sigma0*gEl/R^3 * int(gAz^2 dtheta) * c/(2 fs) * S
%           / (N0ddc * noiseGain * fs)
%   with gAz the one-way azimuth power pattern (two-way amplitude), gEl the
%   two-way elevation element pattern, S the energy of the response of the
%   DDC filters and the matched filter to a unit echo, summed over the
%   delays of one ADC sample (sys.fs) and noiseGain the noise gain of the
%   same filters (as snr_budget.m). The pulse compression gives the clutter
%   and the noise the same gain, so this is the CNR of the range cell.
c  = sys.cfg.const.c;
fs = sys.fs;
wv = sys.wave.short;
[~, u] = steering_vector(sys, 0, el_deg);
gEl = pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, asind(u)) ...
    * pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, asind(u));
% azimuth: integral of the two-way power pattern over the angle [rad]
az = sys.ant.az_deg(:);
g  = pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, az);
sel = abs(az) <= sys.ant.azSpan_deg;
azInt = trapz(deg2rad(az(sel)), g(sel).^2);
% range: energy of the response to a unit echo, summed over the ADC-sample delays
D = prod(sys.ddc.decim);
nPad = numel(sys.ddc.hEq);
tapsUp = zeros(D*(numel(wv.mfTaps)-1) + 1, 1);
tapsUp(1:D:end) = wv.mfTaps;
noiseGain = sum(abs(conv(sys.ddc.hEq, tapsUp)).^2);
S = 0;
for shift = 0:D-1
    t = ((0:round(wv.width_s*fs) + 2*nPad - 1)' - nPad - shift) / fs;
    s = waveform_envelope(wv, t);
    y = fir_decimate_valid(fir_decimate_valid(s, sys.ddc.h1real, sys.ddc.decim(1)), sys.ddc.h2real, sys.ddc.decim(2));
    S = S + sum(abs(conv(y, wv.mfTaps.')).^2);
end
cnr = sys.radarConst * sigma0 * gEl / R_m^3 * azInt * c/(2*fs) * S / (sys.N0ddc * noiseGain * fs);
end

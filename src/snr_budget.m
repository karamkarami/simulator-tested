function b = snr_budget(sys, tgt, t_s)
%SNR_BUDGET  Expected SNR of a target at the output of every processing block.
%   b = SNR_BUDGET(sys, tgt, t_s) computes the SNR (dB) of target tgt at time
%   t_s, assuming the CPI is centred on the target azimuth:
%     b.adc      ADC output, one element, one real sample (noise in fs/2, incl. ADC noise)
%     b.ddc      DDC output, one element, one complex sample
%     b.dbf      beamformer output (best beam b.beam; = b.ddc without beamforming)
%     b.mf       matched filter output (peak), incl. the window, DDC filtering
%                and intra-pulse Doppler losses (b.mfLoss)
%     b.mti      MTI output (gain at the target Doppler)
%     b.doppler  Doppler filter output incl. Hamming window, Doppler straddle and
%                the two-way azimuth beam shape over the CPI. With several rows
%                this is the mean S/N per row before non-coherent integration.
%     b.required SNR per row needed for Pd = 0.9 at the CFAR Pfa with the
%                effective number of independent rows (sys.dop.nLooks)
%                integrated non-coherently: Swerling 1 (Shnidman) when
%                cfg.sim.swerling = 1, else a steady target (Albersheim);
%                b.margin = doppler - required
%   The slow-time gains are obtained by running the actual MTI and Doppler
%   functions on a noise-free target and on unit impulses (white noise).

[R, ~, el, vr] = target_geometry(tgt, t_s, sys.ant.height_m);
b.name    = tgt.name;
b.rcs_m2  = tgt.rcs_m2;
b.range_m = R;
b.el_deg  = el;
b.vr_mps  = vr;
if R < sys.handover_m
    b.pulse = 'short';                     % the short pulse is tested up to the handover range
else
    b.pulse = 'long';
end
wv = sys.wave.(b.pulse);

%% Fast time: received power with the beam pointing at the target
[v, u] = steering_vector(sys, 0, el);
gEl = pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, asind(u)) ...
    * pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, asind(u));
Pr  = sys.radarConst * tgt.rcs_m2 / R^4 * gEl;                     % per element [W]
snrAdc = Pr / (sys.N0eff * sys.fs/2);
snrDdc = Pr / (sys.N0eff * sys.ddc.noiseBw_Hz);
if sys.dbf.enabled
    arrayGain = abs(sys.dbf.W * v).^2 ./ sum(abs(sys.dbf.W).^2, 2);
    [gDbf, b.beam] = max(arrayGain);
else
    gDbf = 1;                                                       % single element: no array gain
    b.beam = 1;
end
snrDbf = snrDdc * gDbf;
% Matched filter: the noise-free pulse (with its Doppler shift) goes through the DDC
% filters and the matched filter. tEff = |peak|^2 / (noise gain * fs) is the effective
% pulse length: equal to the pulse width for an ideal matched filter, smaller with the
% window, the DDC filtering of the waveform and the Doppler mismatch. The echo is
% shifted by fractions of an output sample and the best one is kept (target on a
% range cell, no straddle loss).
fd    = -2*vr/wv.lambda_m;                                          % Doppler [Hz]
D     = prod(sys.ddc.decim);
nPad  = numel(sys.ddc.hEq);
tapsUp = zeros(D*(numel(wv.mfTaps)-1) + 1, 1);
tapsUp(1:D:end) = wv.mfTaps;
noiseGain = sum(abs(conv(sys.ddc.hEq, tapsUp)).^2);                % ADC input -> MF output
peak = 0;
for shift = 0:D/6:D-1
    t = ((0:round(wv.width_s*sys.fs) + 2*nPad - 1)' - nPad - shift) / sys.fs;
    s = waveform_envelope(wv, t) .* exp(1j*2*pi*fd*t);              % baseband echo at the ADC rate
    y = fir_decimate_valid(fir_decimate_valid(s, sys.ddc.h1real, sys.ddc.decim(1)), sys.ddc.h2real, sys.ddc.decim(2));
    peak = max(peak, max(abs(conv(y, wv.mfTaps.')).^2));
end
tEff  = peak / (noiseGain * sys.fs);
b.mfLoss = 10*log10(tEff / wv.width_s);
snrMf = Pr * gDbf * tEff / sys.N0eff;                               % energy / noise density

%% Slow time: one CPI of pulses centred on the target
n   = sys.cpi.length;
st  = mod(0:n-1, sys.nStagger) + 1;                                 % stagger positions
pri = sys.pri(st(1:end-1));
tp  = [0, cumsum(pri(:).')];                                        % pulse times [s]
gAz = pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, sys.ant.rate_dps*(tp - mean(tp)));
sigConst = reshape(exp(1j*2*pi*fd*tp), 1, 1, n);                    % constant amplitude
sigScan  = reshape(gAz .* exp(1j*2*pi*fd*tp), 1, 1, n);             % with beam shape
impulses = reshape(eye(n), 1, n, n);                                % white noise, one impulse per cell

yN = mti_filter(impulses, sys.mtiCoeffs, st);
gMti = mean(abs(mti_filter(sigConst, sys.mtiCoeffs, st)).^2) / mean(sum(abs(yN).^2, 2));
pS = doppler_process(mti_filter(sigScan, sys.mtiCoeffs, st), sys);
pN = sum(doppler_process(yN, sys), 2);
[~, k] = max(pS);
gSlow = pS(k) / pN(k);

b.nRows    = sys.cpi.nRows;
b.adc      = 10*log10(snrAdc);
b.ddc      = 10*log10(snrDdc);
b.dbf      = 10*log10(snrDbf);
b.mf       = 10*log10(snrMf);
b.mti      = 10*log10(snrMf * gMti);
b.doppler  = 10*log10(snrMf * gSlow);
if isfield(sys.cfg, 'sim') && isfield(sys.cfg.sim, 'swerling') && sys.cfg.sim.swerling == 1
    b.required = shnidman_swerling1_snr(0.9, sys.cfar.pfa, sys.dop.nLooks);
else
    b.required = albersheim_snr(0.9, sys.cfar.pfa, sys.dop.nLooks);
end
b.margin   = b.doppler - b.required;
end

function snr_dB = shnidman_swerling1_snr(pd, pfa, n)
% SNR per sample for a Swerling 1 target (constant during the dwell) with n
% samples integrated non-coherently: exact for one sample (Pd = Pfa^(1/(1+SNR))),
% else Shnidman's approximation (+-0.5 dB, 0.1 <= pd <= 0.99, n < 40).
if abs(n - 1) < 1e-9
    snr_dB = 10*log10(log(pfa)/log(pd) - 1);
    return
end
eta  = sqrt(-0.8*log(4*pfa*(1 - pfa))) + sign(pd - 0.5)*sqrt(-0.8*log(4*pd*(1 - pd)));
xInf = eta * (eta + 2*sqrt(n/2 - 0.25));
c1 = ((17.7006*pd - 18.4496)*pd + 14.5339)*pd - 3.525;
c2 = exp(27.31*pd - 25.14) + (pd - 0.8)*(0.7*log(1e-5/pfa) + (2*n - 20)/80);
if pd <= 0.872
    cdB = c1;
else
    cdB = c1 + c2;
end
snr_dB = 10*log10(10^(cdB/10) * xInf / n);
end

function snr_dB = albersheim_snr(pd, pfa, n)
% Albersheim: SNR per sample for a steady target with n samples integrated non-coherently.
A = log(0.62/pfa);
B = log(pd/(1 - pd));
snr_dB = -5*log10(n) + (6.2 + 4.54/sqrt(n + 0.44)) * log10(A + 0.12*A*B + 1.7*B);
end

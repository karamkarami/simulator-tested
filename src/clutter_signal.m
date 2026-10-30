function xc = clutter_signal(clut, sys, az_deg, t_s)
%CLUTTER_SIGNAL  Complex IF clutter echo of one PRI (common to all elements).
%   xc = CLUTTER_SIGNAL(clut, sys, az_deg, t_s) returns the clutter echo of
%   both sub-pulses during the receive window [nAdc x 1] for the antenna
%   azimuth az_deg at time t_s. The two-way azimuth pattern weights the
%   clutter cells (this also creates the scanning modulation); every cell
%   rotates with its internal-motion Doppler (plus a stationary part with
%   the exponential spectrum model, see clutter_init.m). Element i receives
%   real(xc * sys.clutterSteer(i)).

azRel = wrap180(clut.az_deg - az_deg);
sel   = abs(azRel) <= sys.ant.azSpan_deg;
gAz   = pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, azRel(sel));     % one-way power = two-way amplitude
refl  = gAz.' * (clut.rho(sel, :) .* exp(1j*2*pi*clut.fd_Hz(sel, :)*t_s)); % [1 x nRange]
if ~isempty(clut.dc)                                                     % stationary part (exponential model)
    refl = refl + gAz.' * clut.dc(sel, :);
end

xc = zeros(sys.nAdc, 1);
names = {'short', 'long'};
for i = 1:2
    w = clut.(names{i});
    % Echo of each range cell starts at ADC sample n0 (0-based, may be negative)
    n0  = w.start0 + clut.delaySamples;
    use = n0 > -w.nTx & n0 < sys.nAdc;
    % Impulses shifted by nTx-1 so that partially received echoes are kept
    imp = zeros(w.nFft, 1);
    imp(n0(use) + w.nTx) = refl(use);
    y = ifft(fft(imp) .* w.txFft);                         % convolution with the IF pulse
    xc = xc + y(w.nTx : w.nTx + sys.nAdc - 1);
end
end

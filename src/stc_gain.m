function g = stc_gain(sys, t)
%STC_GAIN  Amplitude gain of the STC (sensitivity time control) at times t of the PRI.
%   g = STC_GAIN(sys, t) returns the gain (<= 1) of the attenuator in front
%   of the receiver at the times t [s] (from the PRI start; any shape). The
%   attenuation depends on the time dt since the end of the pulse whose
%   receive window t lies in, given as the equivalent range R = c*dt/2 (the
%   nearest range from which echoes still arrive at that time; this is where
%   the strongest clutter comes from): cfg.rx.stc.curve holds rows
%   [R_m attenuation_dB], interpolated linearly in dB, no attenuation beyond
%   the last row; cfg.rx.stc.curveLong, when not empty, is the curve of the
%   long-pulse window (its start also holds the eclipsed echoes of the
%   targets nearer than the long-pulse blind range). While a pulse is
%   transmitted the attenuation of the first row applies. With
%   cfg.rx.stc.enable = false the gain is 1.
%   The receiver noise is added after the attenuator (generate_rx_signal.m):
%   the STC lowers the sensitivity at short range on purpose and keeps the
%   strong near clutter below the full scale of the ADC.
g = ones(size(t));
stc = sys.cfg.rx.stc;
if ~stc.enable, return; end
c = sys.cfg.const.c;
ws = sys.wave.short;  wl = sys.wave.long;
curveLong = stc.curve;
if isfield(stc, 'curveLong') && ~isempty(stc.curveLong), curveLong = stc.curveLong; end
A = zeros(size(t)) + stc.curve(1, 2);                          % while transmitting
inShort = t >= ws.tTx_s + ws.width_s & t < wl.tTx_s;          % receive windows
inLong  = t >= wl.tTx_s + wl.width_s;
A(inShort) = attenuation(stc.curve, c * (t(inShort) - ws.tTx_s - ws.width_s) / 2);
A(inLong)  = attenuation(curveLong, c * (t(inLong) - wl.tTx_s - wl.width_s) / 2);
g = 10.^(-A/20);
end

function A = attenuation(curve, R)
% Attenuation [dB] at the equivalent ranges R: linear in dB between the rows, 0 beyond the last.
A = interp1(curve(:, 1), curve(:, 2), min(max(R, curve(1, 1)), curve(end, 1)));
A(R > curve(end, 1)) = 0;
end

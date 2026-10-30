function x = generate_rx_signal(sys, targets, t_s, azAnt_deg, clutterIf)
%GENERATE_RX_SIGNAL  Received IF signal of all elements at the ADC input for one PRI.
%   x = GENERATE_RX_SIGNAL(sys, targets, t_s, azAnt_deg, clutterIf)
%   returns the real signal [nAdc x nElements] at the ADC input in sqrt(W)
%   (a sine of amplitude A has power A^2/2), for the PRI starting at t_s with
%   the antenna at azAnt_deg:
%     1. every target: position from its motion, received power from the
%        radar equation (RCS of this antenna pass, Swerling 1), echo of both
%        sub-pulses on the 30 MHz IF with the element phases
%     2. ground clutter (clutterIf from clutter_signal.m, [] = none)
%     3. transmitter leakage while a pulse is transmitted (isolation)
%     4. thermal noise of the receiver (k*T0*F over fs/2)
%     5. receiver path gain
%   adc_sample.m then converts the result to 16-bit ADC samples.

c   = sys.cfg.const.c;
fs  = sys.fs;
% Single precision is plenty for a 16-bit ADC and about twice as fast
x   = zeros(sys.nAdc, sys.nCh, 'single');
names = {'short', 'long'};

%% 1) Target echoes at the receiver input
for i = 1:numel(targets)
    [R, az, el, vr] = target_geometry(targets(i), t_s, sys.ant.height_m);
    azRel = wrap180(az - azAnt_deg);
    if abs(azRel) > sys.ant.azSpan_deg
        continue                                   % far outside the azimuth beam
    end
    [v, u] = steering_vector(sys, azRel, el);
    elRel  = asind(u);
    % Two-way pattern factor: product of the one-way power patterns
    gain = pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, azRel)^2 ...
         * pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, elRel) ...
         * pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, elRel);
    rcs   = target_rcs(targets(i), target_look(sys, az, t_s));        % Swerling 1: one value per pass
    power = sys.radarConst * rcs / R^4 * gain;                         % per element [W]
    amp   = sqrt(2*power);                                                % [sqrt(W)]
    for p = 1:2
        wv   = sys.wave.(names{p});
        tau0 = 2*R/c;
        % ADC samples that contain the echo
        first = max(1, ceil((wv.tTx_s + tau0 - sys.rxStart_s)*fs) + 1);
        last  = min(sys.nAdc, floor((wv.tTx_s + tau0 + wv.width_s - sys.rxStart_s)*fs) + 1);
        if first > last, continue; end
        n  = (first:last)';
        t  = sys.tAdc(n);                           % time inside the PRI
        tau = 2*(R + vr*t)/c;                       % delay of every sample (intra-pulse Doppler)
        % RF echo mixed down by the LO: IF carrier plus the carrier phase -2*pi*fLO*tau
        echo = amp * waveform_envelope(wv, t - wv.tTx_s - tau0) ...
             .* exp(1j*2*pi*(sys.fIF*(t - wv.tTx_s - tau) - sys.fLO*tau));
        x(n, :) = x(n, :) + real(single(echo) * single(v).');
    end
end

%% 2) Clutter (same echo on all elements, different phase)
if ~isempty(clutterIf)
    x = x + real(single(clutterIf) * single(sys.clutterSteer).');
end

%% 3) Transmitter leakage into every receive channel while a pulse is transmitted
if sys.cfg.rx.addLeakage
    for p = 1:2
        wv = sys.wave.(names{p});
        n  = find(sys.tAdc >= wv.tTx_s & sys.tAdc < wv.tTx_s + wv.width_s);
        t  = sys.tAdc(n) - wv.tTx_s;
        leak = sys.rx.leakAmp * real(waveform_envelope(wv, t) .* exp(1j*2*pi*sys.fIF*t));
        x(n, :) = x(n, :) + single(leak);
    end
end

%% 4) Thermal noise of the receiver
if sys.cfg.rx.addThermalNoise
    x = x + sys.rx.sigmaTh * randn(sys.nAdc, sys.nCh, 'single');
end

%% 5) Receiver path gain
x = x * sys.rx.voltGain;
end

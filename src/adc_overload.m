function frac = adc_overload(adcData, sys)
%ADC_OVERLOAD  Fraction of the ADC samples of the receive intervals that are at full scale.
%   frac = ADC_OVERLOAD(adcData, sys) counts the samples of adcData
%   [nAdc x nElements] (adc_sample.m) inside the receive intervals of both
%   sub-pulses (sys.adc.rxSamples; the transmitter leakage always saturates
%   the ADC while a pulse is sent) that are clipped to the ADC range. Clipped
%   clutter is no longer cancelled by the MTI, and clipped targets give
%   intermodulation products at other Doppler bins and beams (false plots):
%   a fraction above cfg.adc.overloadWarn means that cfg.rx.gain_dB is too
%   high for the echoes (e.g. after a change of cfg.fc: the echo power grows
%   with the wavelength squared).
x = adcData(sys.adc.rxSamples, :);
fs = sys.adc.fullScale;
frac = nnz(x >= fs - 1 | x <= -fs) / numel(x);
end

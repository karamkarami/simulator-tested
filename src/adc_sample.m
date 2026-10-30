function adcData = adc_sample(x, sys)
%ADC_SAMPLE  16-bit ADC with 12.7 effective bits.
%   adcData = ADC_SAMPLE(x, sys) converts the ADC input x (sqrt(W), see
%   generate_rx_signal.m) to LSB and quantizes it to 16 bits (rounding and
%   clipping to int16). Then the ADC noise that sets the effective number of
%   bits is added, and the result is rounded and clipped to int16 again, so
%   the samples stay signed 16-bit integers. With cfg.adc.quantize = false
%   (SNR tests only) the analog values in LSB plus the ADC noise are
%   returned instead.
x = x * sys.adc.lsbPerUnit;                           % [sqrt(W)] -> LSB
if ~sys.cfg.adc.quantize
    if sys.cfg.adc.addAdcNoise
        x = x + sys.adc.sigmaAdc * randn(size(x), 'single');
    end
    adcData = x;
    return
end
adcData = quantize(x, sys);                           % 16-bit quantization
if sys.cfg.adc.addAdcNoise                            % ADC noise (ENOB 12.7)
    adcData = quantize(single(adcData) + sys.adc.sigmaAdc * randn(size(x), 'single'), sys);
end
end

function q = quantize(x, sys)
% Round and clip to the ADC word length (stored as int16).
q = int16(x);                                         % rounds and clips to 16 bits
if sys.cfg.adc.nBits < 16
    fs = sys.adc.fullScale;
    q = min(max(q, -fs), fs - 1);
end
end

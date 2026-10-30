function [I, Q] = ddc(adcData, sys)
%DDC  Two-stage fixed-point digital down-converter of all channels.
%   [I, Q] = DDC(adcData, sys) converts the 16-bit ADC samples [nAdc x nElements]
%   on the 30 MHz IF to complex baseband, returned as signed 16-bit I and Q
%   [nDdc x nElements]:
%     mixer   : I line x .* cos(2*pi*fIF*t), Q line x .* (-sin(2*pi*fIF*t))
%     stage 1 : low-pass filter h1 (valid time-domain convolution), decimation by 10
%     stage 2 : low-pass filter h2 (valid time-domain convolution), decimation by 6
%               (2 MHz output) or 3 (4 MHz output, cfg.ddc.outputRate_Hz)
%   The filter coefficients are signed 16-bit integers with 13 fractional bits.
%   The accumulator of every stage is shifted right by cfg.ddc.shift bits
%   (13 = unit gain), rounded and saturated to signed 16 bit.
%   The time of every output sample is sys.tDdc; the samples of the short and
%   of the long sub-pulse are sys.wave.short.ddcIdx and sys.wave.long.ddcIdx.
fixedPoint = sys.cfg.ddc.fixedPoint;
x  = double(adcData);                              % integer arithmetic is exact in double
xI = x .* sys.ddc.cos;                             % real (I) line
xQ = x .* sys.ddc.minusSin;                        % imaginary (Q) line
D  = sys.ddc.decim;
s  = sys.ddc.shift;
I  = stage(stage(xI, sys.ddc.h1, D(1), s(1), fixedPoint), sys.ddc.h2, D(2), s(2), fixedPoint);
Q  = stage(stage(xQ, sys.ddc.h1, D(1), s(1), fixedPoint), sys.ddc.h2, D(2), s(2), fixedPoint);
if fixedPoint
    I = int16(I);
    Q = int16(Q);
end
end

function y = stage(x, h, D, shift, fixedPoint)
% One DDC stage: filter + decimation, then scaling to 16 bit.
y = fir_decimate_valid(x, h, D) / 2^shift;
if fixedPoint
    y = min(max(round(y), -32768), 32767);         % round and saturate to signed 16 bit
end
end

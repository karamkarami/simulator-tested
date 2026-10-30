function Y = ddc_pulse(sys, wv, phases)
%DDC_PULSE  A sub-pulse after the DDC filters, sampled at the DDC output rate.
%   Y = DDC_PULSE(sys, wv, phases) returns the noise-free baseband pulse wv
%   (no Doppler) after the two DDC filters, one column per sampling phase:
%   phases are delays of the pulse in DDC output samples (e.g. 0:0.25:0.75;
%   default 0). The first and last samples hold the filter transients.
if nargin < 3, phases = 0; end
nPad = numel(sys.ddc.hEq);
D = prod(sys.ddc.decim);
Y = [];
for ph = phases(:).'
    t = ((0:round(wv.width_s*sys.fs) + 2*nPad - 1)' - nPad - ph*D) / sys.fs;
    y = fir_decimate_valid(fir_decimate_valid(waveform_envelope(wv, t), sys.ddc.h1real, sys.ddc.decim(1)), ...
                           sys.ddc.h2real, sys.ddc.decim(2));
    Y(:, end+1) = y; %#ok<AGROW>
end
end

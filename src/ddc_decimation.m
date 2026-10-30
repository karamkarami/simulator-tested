function [D, fileKey] = ddc_decimation(cfg, outputRate_Hz)
%DDC_DECIMATION  Decimation of the two DDC stages and the filter file of the output rate.
%   [D, fileKey] = DDC_DECIMATION(cfg) returns D = [D1 D2] for the output
%   rate cfg.ddc.outputRate_Hz (2 or 4 MHz): stage 1 decimates the ADC rate
%   by cfg.ddc.stage1Decimation, stage 2 the stage-1 rate down to the output
%   rate. fileKey names the filter file in cfg.files: 'ddc' (2 MHz) or
%   'ddc4MHz' (4 MHz).
%   DDC_DECIMATION(cfg, outputRate_Hz) does the same for another rate.
if nargin < 2, outputRate_Hz = cfg.ddc.outputRate_Hz; end
switch outputRate_Hz
    case 2e6, fileKey = 'ddc';
    case 4e6, fileKey = 'ddc4MHz';
    otherwise
        error('cfg.ddc.outputRate_Hz must be 2e6 or 4e6 (is %g).', outputRate_Hz);
end
D1 = cfg.ddc.stage1Decimation;
D2 = cfg.adc.fs_Hz / D1 / outputRate_Hz;
if D2 < 1 || abs(D2 - round(D2)) > 1e-9
    error('The stage-1 rate %g MHz is not a whole multiple of the DDC output rate %g MHz.', ...
          cfg.adc.fs_Hz/D1/1e6, outputRate_Hz/1e6);
end
D = [D1 round(D2)];
end

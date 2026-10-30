function save_pulse_file(folder, matFormat, hdr, adcData) %#ok<INUSD>
%SAVE_PULSE_FILE  Save the raw ADC data of one PRI to <folder>/pulse_<number>.mat.
%   The file holds adcData (int16 [nSamples x nElements]) and hdr (pulse
%   number, scan number, antenna azimuth, time, PRI, sampling parameters ...).
fileName = fullfile(folder, sprintf('pulse_%07d.mat', hdr.pulseNumber));
save(fileName, 'adcData', 'hdr', matFormat);
end

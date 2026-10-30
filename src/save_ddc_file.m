function save_ddc_file(folder, matFormat, hdr, I, Q) %#ok<INUSD>
%SAVE_DDC_FILE  Save the DDC output of one PRI to <folder>/ddc_<number>.mat.
%   The file holds I and Q (int16 [nDdc x nElements], both sub-pulses) and
%   hdr (the same header as the raw data file). replay_ddc_udp.m sends these
%   files over UDP with the original PRI timing.
fileName = fullfile(folder, sprintf('ddc_%07d.mat', hdr.pulseNumber));
save(fileName, 'I', 'Q', 'hdr', matFormat);
end

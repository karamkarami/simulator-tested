function save_block_file(folder, prefix, matFormat, hdr, data) %#ok<INUSL>
%SAVE_BLOCK_FILE  Save the output of a processing block of one PRI to <folder>/<prefix>_<number>.mat.
%   SAVE_BLOCK_FILE(folder, prefix, matFormat, hdr, data) writes every field
%   of the struct data as a variable of its own, and hdr (the header of the
%   PRI, pulse_header.m). Used by run_radar_sim for the matched filter output
%   (prefix 'mf', cfg.io.saveMf) and the Doppler processing output (prefix
%   'doppler', cfg.io.saveDoppler).
data.hdr = hdr;
fileName = fullfile(folder, sprintf('%s_%07d.mat', prefix, hdr.pulseNumber));
save(fileName, '-struct', 'data', matFormat);
end

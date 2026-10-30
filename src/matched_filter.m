function rp = matched_filter(beams, wv)
%MATCHED_FILTER  Pulse compression of every beam (valid convolution).
%   rp = MATCHED_FILTER(beams, wv) filters every row of beams [nBeams x nSamples]
%   with the matched filter of sub-pulse wv and returns the range profiles
%   rp [nBeams x nRange] on the range axis wv.range_m (delay axis wv.delay_s,
%   index axis 1:wv.nRange).
rp = conv2(beams, wv.mfTaps, 'valid');
end

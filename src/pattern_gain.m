function g = pattern_gain(angles_deg, gain_dB, query_deg)
%PATTERN_GAIN  Linear (power) gain of a tabulated antenna pattern.
%   g = PATTERN_GAIN(angles_deg, gain_dB, query_deg) interpolates the pattern
%   (given in dB) at the angles query_deg and returns the linear power gain.
%   Angles outside the table get the lowest gain of the table.
g = 10.^(interp1(angles_deg, gain_dB, query_deg, 'linear', min(gain_dB)) / 10);
end

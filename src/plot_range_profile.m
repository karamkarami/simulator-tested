function plot_range_profile(rp, wv, axisType, beams)
%PLOT_RANGE_PROFILE  Plot matched-filter outputs against range, delay or sample index.
%   PLOT_RANGE_PROFILE(rp, wv, axisType, beams) plots |rp|^2 in dB for the
%   selected beams (default: all). axisType is 'range', 'delay' or 'index'.
if nargin < 3, axisType = 'range'; end
if nargin < 4, beams = 1:size(rp, 1); end
switch lower(axisType)
    case 'range', ax = wv.range_m / 1e3;   label = 'Range [km]';
    case 'delay', ax = wv.delay_s * 1e6;   label = 'Delay [\mus]';
    case 'index', ax = 1:size(rp, 2);      label = 'Range cell index';
    otherwise, error('axisType must be ''range'', ''delay'' or ''index''.');
end
figure('Name', sprintf('Matched filter output (%s pulse)', wv.name));
plot(ax, 10*log10(abs(rp(beams, :)).^2 + eps).');
grid on;
xlabel(label);
ylabel('Power [dB]');
title(sprintf('Matched filter output, %s pulse', wv.name));
end

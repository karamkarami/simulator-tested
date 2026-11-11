function q = plot_params(p, b)
%PLOT_PARAMS  Tracker settings for one plot of a batch (b = [time az range el plotIndex vr vrSpan]).
%   q = PLOT_PARAMS(p, b) is p with the span of the plot's Doppler axis
%   (column 7, the carrier of its sub-pulse: plot_vr_span.m) as q.vrSpan.
q = p;
if numel(b) >= 7 && isfinite(b(7)) && b(7) > 0, q.vrSpan = b(7); end
end

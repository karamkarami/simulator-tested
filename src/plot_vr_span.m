function s = plot_vr_span(p, plot)
%PLOT_VR_SPAN  Span of the Doppler axis of a plot [m/s] (its radial velocity is known modulo it).
%   s = PLOT_VR_SPAN(p, plot) with the tracker settings p (track_setup.m):
%   the span on the carrier of the plot's sub-pulse (plot.pulse 'short' or
%   'long', cfg.wave.<pulse>.fc_Hz), or p.vrSpan.
s = p.vrSpan;
if isfield(p, 'vrSpanPulse') && isfield(plot, 'pulse') && ischar(plot.pulse) && isfield(p.vrSpanPulse, plot.pulse)
    s = p.vrSpanPulse.(plot.pulse);
end
end

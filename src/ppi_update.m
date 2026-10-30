function sys = ppi_update(ppi, sys, az_deg, t_s, scan, plots, dets, tracks, timing_us)
%PPI_UPDATE  Refresh the PPI: antenna, true targets, plots, detections, tracks, block times.
%   sys = PPI_UPDATE(ppi, sys, az_deg, t_s, scan, plots, dets, tracks, timing_us)
%   stores the data in the figure, so the interactive options of
%   ppi_callbacks.m (filters, zoom, data tips) can redraw it at any time.
%   dets are the CFAR detections before the plot extractor (rows of
%   detection_log_init.m), tracks the current tracks (track_report.m) and
%   timing_us the mean time per PRI of DBF, matched filter, MTI, Doppler,
%   detection, plot extractor, tracking, console output and display [us]
%   (all optional).
%   While the display is paused (space key or the Pause button) the data is
%   still stored but not drawn.
%   When the user applied a new CFAR threshold in the panel, the returned sys
%   holds it (sys.cfar.alpha), so the processing uses it from now on.
if nargin < 7, dets = zeros(0, 10); end
if nargin < 8, tracks = []; end
if nargin < 9, timing_us = []; end
if ~ishandle(ppi.fig), return; end                       % the window was closed
dB = ppi_callbacks('takeCfar', ppi.fig);
if ~isempty(dB)
    sys.cfar.alpha = 10^(dB/10);
    if sys.cfg.sim.verbose
        fprintf('  CFAR threshold set to %.2f dB (Pfa %.1e) at t = %.3f s\n', dB, ...
                cfar_pfa(sys.cfar.alpha, sys.cfar.nTrain, sys.dop.nLooks), t_s);
    end
end
S = getappdata(ppi.fig, 'ppi');
S.sys = sys;
S.data = struct('az_deg', az_deg, 't_s', t_s, 'scan', scan, 'plots', plots, 'dets', dets, ...
                'tracks', tracks, 'timing_us', timing_us);
setappdata(ppi.fig, 'ppi', S);
if ~isempty(dB)
    ppi_callbacks('cfarApplied', ppi.fig, dB);
end
if ~S.opt.paused
    ppi_callbacks('redraw', ppi.fig);
else
    ppi_callbacks('timing', ppi.fig);                     % the block times are shown also when paused
end
drawnow;
end

function [sys, want] = ppi_update(ppi, sys, az_deg, t_s, scan, plots, dets, tracks, timing_us, timingPris, video, rd)
%PPI_UPDATE  Refresh the PPI: antenna, true targets, plots, detections, tracks, block times.
%   [sys, want] = PPI_UPDATE(ppi, sys, az_deg, t_s, scan, plots, dets, tracks, timing_us,
%                            timingPris, video, rd)
%   stores the data in the figure, so the interactive options of
%   ppi_callbacks.m (filters, zoom, data tips) can redraw it at any time.
%   dets are the CFAR detections before the plot extractor (rows of
%   detection_log_init.m), tracks the current tracks (track_report.m) and
%   timing_us the mean time per PRI of DBF, matched filter, MTI, Doppler,
%   detection, peak interpolation, plot extractor, tracking, console output
%   and display [us]
%   over the last timingPris PRIs ([] = the whole run), video the raw data map
%   (ppi_video_add.m, ppi_video_add_row.m; [] = the one stored), rd the
%   range-Doppler map of the selected beam (ppi_rd_map.m; [] = the one stored)
%   (all optional). want says what the display needs from the caller:
%   want.raw = 'video', 'dbf', 'mf' (fill the map every PRI with that data of
%   beam want.beam) or '' (none), want.rd = true while the range-Doppler tab
%   is shown (pass the map of beam want.beam).
%   While the display is paused (space key or the Pause button) the data is
%   kept aside and not drawn (the display, its data tips and the range-Doppler
%   map stay as they were); Resume shows the newest data.
%   When the user applied a new CFAR threshold in the panel, the returned sys
%   holds it (sys.cfar.alpha), so the processing uses it from now on.
if nargin < 7, dets = zeros(0, 10); end
if nargin < 8, tracks = []; end
if nargin < 9, timing_us = []; end
if nargin < 10, timingPris = []; end
if nargin < 11, video = []; end
if nargin < 12, rd = []; end
want = struct('raw', '', 'beam', 1, 'rd', false);
if ~ishandle(ppi.fig), return; end                       % the window was closed
dB = ppi_callbacks('takeCfar', ppi.fig);
if ~isempty(dB)
    if isfield(sys.cfar, 'alphaMin'), dB = max(dB, 10*log10(sys.cfar.alphaMin)); end   % cfg.cfar.maxPfa
    sys.cfar.alpha = 10^(dB/10);
    if sys.cfg.console.enable
        fprintf('  CFAR threshold set to %.2f dB (Pfa %.1e) at t = %.3f s\n', dB, ...
                cfar_pfa(sys.cfar.alpha, sys.cfar, sys.dop.nLooks), t_s);
    end
end
S = getappdata(ppi.fig, 'ppi');
S.sys = sys;
base = S.data;
if S.opt.paused && ~isempty(S.pending), base = S.pending; end
if isempty(video), video = base.video; end
if isempty(rd), rd = base.rd; end
new = struct('az_deg', az_deg, 't_s', t_s, 'scan', scan, 'plots', plots, 'dets', dets, ...
             'tracks', tracks, 'timing_us', timing_us, 'timingPris', timingPris, 'video', video, 'rd', rd);
if S.opt.paused
    S.pending = new;                                      % shown on resume
else
    S.data = new;
    S.pending = [];
end
if any(strcmp(S.opt.showMode, {'video', 'dbf', 'mf'})), want.raw = S.opt.showMode; end
want.beam = S.opt.beam;
want.rd = S.opt.tab == 2;
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

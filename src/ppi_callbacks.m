function out = ppi_callbacks(action, fig, varargin)
%PPI_CALLBACKS  Interaction with the PPI created by ppi_init.m.
%   out = PPI_CALLBACKS(action, fig, ...) is called by the figure and its
%   controls; it can also be called directly (tests, scripts):
%     'redraw'                    draw the stored data with the current options
%     'scroll', evt               mouse wheel: zoom around the cursor
%     'zoom', factor, [x y]       zoom by factor (> 1 = in) around [x y] km
%     'pan', [dx dy]              move the view by [dx dy] km
%     'down' / 'move' / 'up'      mouse buttons and motion (drag = pan, click = pin)
%     'hover', [x y]              data tip of the object nearest to [x y] km
%                                 (out = its text, '' if nothing is close)
%     'pin', [x y]                pin (or remove) the data tip of the nearest object
%     'key', evt                  + - arrows r c d space
%     'option', name, value       showMode ('plots', 'detections', 'both', 'video',
%                                 'dbf' or 'mf'), beam (1 ... beams),
%                                 showShort, showLong, showTruth, showRings,
%                                 showLabels, showAntenna, showTracks,
%                                 minSnr_dB, history
%     'timing'                    write the block processing times into the panel
%     'tab', k                    tab 1 (PPI) or 2 (range-Doppler map of the
%                                 selected beam, small PPI next to it)
%     'rdhover', [range_km v]     readout and spectrum of the range-Doppler map at
%                                 that point (out = the readout text)
%     'reset' 'clear' 'cleartips' 'pause'   buttons of the panel
%     'save' [, fileName]         save a screenshot of the PPI as a PNG file
%                                 (without fileName: file dialog; out = file)
%     'cfar', dB                  request a new CFAR threshold (Apply button);
%                                 the processing takes it at its next display update
%     'takeCfar'                  out = the requested threshold [dB] ([] if none),
%                                 and clears the request (used by ppi_update.m)
%     'cfarApplied', dB           show the threshold the processing uses
out = [];
if ~ishandle(fig), return; end
S = getappdata(fig, 'ppi');
switch action
    case 'redraw'
        S = redraw(S);
    case 'timing'
        show_timing(S);
    case 'tab'
        S.opt.tab = varargin{1};
        S = layout_tab(S);
        S = redraw(S);
    case 'rdhover'
        [S, out] = rd_hover(S, varargin{1});
    case 'scroll'
        evt = varargin{1};
        p = cursor_km(S);
        if ~isempty(p)
            S = zoom_at(S, 1.25^(-evt.VerticalScrollCount), p);
        end
    case 'zoom'
        S = zoom_at(S, varargin{1}, varargin{2});
    case 'pan'
        S = set_limits(S, get(S.ax, 'XLim') + varargin{1}(1), get(S.ax, 'YLim') + varargin{1}(2));
    case 'down'
        p = cursor_km(S);
        if isempty(p), return; end                          % outside the axes
        if strcmp(get(fig, 'SelectionType'), 'open')        % double click
            S = set_limits(S, S.home(1:2), S.home(3:4));
            S.drag = [];
        else
            S.drag = struct('px', figure_px(fig), 'xl', get(S.ax, 'XLim'), 'yl', get(S.ax, 'YLim'), ...
                            'moved', false, 'p', p);
        end
    case 'move'
        if ~isempty(S.drag)                                 % pan with the mouse
            d = figure_px(fig) - S.drag.px;
            if any(abs(d) > 3), S.drag.moved = true; end
            if S.drag.moved
                k = km_per_pixel(S, S.drag.xl, S.drag.yl);
                S = set_limits(S, S.drag.xl - d(1)*k, S.drag.yl - d(2)*k);
            end
        else
            if S.opt.tab == 2                               % over the range-Doppler map: readout, spectrum
                q = get(S.h.rdAx, 'CurrentPoint');
                xl = get(S.h.rdAx, 'XLim');  yl = get(S.h.rdAx, 'YLim');
                if q(1, 1) >= xl(1) && q(1, 1) <= xl(2) && q(1, 2) >= yl(1) && q(1, 2) <= yl(2)
                    [S, out] = rd_hover(S, q(1, 1:2));
                    setappdata(fig, 'ppi', S);
                    return
                end
            end
            p = cursor_km(S);
            S = show_cursor(S, p);
            [S, out] = hover(S, p);
        end
    case 'up'
        if ~isempty(S.drag) && ~S.drag.moved                % a click: pin the data tip
            S = pin(S, S.drag.p);
        end
        S.drag = [];
    case 'hover'
        [S, out] = hover(S, varargin{1});
    case 'pin'
        S = pin(S, varargin{1});
        out = numel(S.h.pins);
    case 'key'
        S = key(S, varargin{1});
    case 'option'
        S.opt.(varargin{1}) = varargin{2};
        if strcmp(varargin{1}, 'showMode')
            set(S.h.ui.showMode, 'Value', find(strcmp(varargin{2}, {'plots', 'detections', 'both', 'video', 'dbf', 'mf'})));
        elseif strcmp(varargin{1}, 'beam')
            set(S.h.ui.beam, 'Value', varargin{2});
            S.data.rd = [];                                 % the map of the new beam comes with the next update
            S.rdCursorBin = [];
        end
        S = redraw(S);
    case 'cfar'
        dB = varargin{1};
        if isscalar(dB) && isfinite(dB) && dB >= 0 && dB <= 50
            S.opt.cfarRequest_dB = dB;
            set(S.h.ui.cfar, 'String', sprintf('%.2f', dB));
            set(S.h.ui.cfarInfo, 'String', sprintf('%.2f dB requested ...', dB));
        else
            set(S.h.ui.cfarInfo, 'String', 'enter 0 ... 50 dB');
        end
    case 'takeCfar'
        out = S.opt.cfarRequest_dB;
        S.opt.cfarRequest_dB = [];
    case 'cfarApplied'
        dB = varargin{1};
        set(S.h.ui.cfar, 'String', sprintf('%.2f', dB));
        set(S.h.ui.cfarInfo, 'String', sprintf('%s: %.2f dB (Pfa %.1e)', S.sys.cfar.label, dB, ...
            cfar_pfa(10^(dB/10), S.sys.cfar, S.sys.dop.nLooks)));
    case 'reset'
        S = set_limits(S, S.home(1:2), S.home(3:4));
    case 'clear'
        S.opt.clearTime_s = S.data.t_s;
        S = redraw(S);
    case 'cleartips'
        S = clear_pins(S);
    case 'pause'
        S.opt.paused = ~S.opt.paused;
        labels = {'Pause display', 'Resume display'};
        set(S.h.ui.pause, 'String', labels{1 + S.opt.paused});
        if ~S.opt.paused, S = redraw(S); end
    case 'save'
        if ~isempty(varargin)
            file = varargin{1};
        else
            [name, folder] = uiputfile('*.png', 'Save a screenshot of the PPI', ...
                                       sprintf('ppi_%s.png', datestr(now, 'yyyymmdd_HHMMSS')));
            if isequal(name, 0), return; end
            file = fullfile(folder, name);
        end
        setappdata(fig, 'ppi', S);
        print(fig, file, '-dpng');
        out = file;
        return
    otherwise
        error('ppi_callbacks: unknown action ''%s''', action);
end
setappdata(fig, 'ppi', S);
end

%% ------------------------------------------------------------------------
function S = redraw(S)
% Antenna, true targets and the plots that pass the filters.
h = S.h;  o = S.opt;  d = S.data;  sys = S.sys;
onOff = {'off', 'on'};
set(h.ant, 'XData', [0 S.rMax*sind(d.az_deg)], 'YData', [0 S.rMax*cosd(d.az_deg)], ...
           'Visible', onOff{1 + o.showAntenna});
set(h.rings, 'Visible', onOff{1 + o.showRings});

targets = sys.cfg.targets;
if o.showTruth && ~isempty(targets)
    n = numel(targets);
    xy = zeros(n, 2);
    for i = 1:n
        p = target_state(targets(i), d.t_s);
        xy(i, :) = p(1:2) / 1e3;
        tt = [0, d.t_s];                                    % true path since t = 0
        if isfield(targets(i), 'maneuvers') && ~isempty(targets(i).maneuvers)
            tt = linspace(0, d.t_s, 41);                     % maneuvers: curved path
        end
        path = zeros(numel(tt), 2);
        for k = 1:numel(tt)
            q = target_state(targets(i), tt(k));
            path(k, :) = q(1:2) / 1e3;
        end
        set(h.names(i), 'Position', [xy(i, :) 0], 'Visible', 'on');
        set(h.tracks(i), 'XData', path(:, 1).', 'YData', path(:, 2).', 'Visible', 'on');
    end
    set(h.truth, 'XData', xy(:, 1), 'YData', xy(:, 2), 'Visible', 'on');
else
    set(h.truth, 'Visible', 'off');
    set(h.names, 'Visible', 'off');
    set(h.tracks, 'Visible', 'off');
end

% Plots that pass the filters (shown unless only the detections are selected)
showVideo = any(strcmp(o.showMode, {'video', 'dbf', 'mf'}));   % raw data instead of plots and detections
showPlots = ~strcmp(o.showMode, 'detections') && ~showVideo;
showDets  = ~strcmp(o.showMode, 'plots') && ~showVideo;
S = draw_video(S, showVideo);
plots = d.plots;
nAll = numel(plots);
x = zeros(1, nAll);  y = x;  show = false(1, nAll);  isShort = show;
if nAll > 0
    isShort = strcmp({plots.pulse}, 'short');
    el = [plots.el_deg];
    el(isnan(el)) = 0;                                     % no elevation without beamforming
    groundRange = [plots.range_m] .* cosd(el) / 1e3;
    x = groundRange .* sind([plots.az_deg]);
    y = groundRange .* cosd([plots.az_deg]);
    t = [plots.time_s];
    show = ((isShort & o.showShort) | (~isShort & o.showLong)) ...
         & [plots.snr_dB] >= o.minSnr_dB & t > o.clearTime_s;
    if o.history > 0
        show = show & (d.t_s - t) < o.history * 360 / sys.ant.rate_dps;
    end
    show = show & showPlots;
end
s = show & isShort;  l = show & ~isShort;
set(h.short, 'XData', x(s), 'YData', y(s));
set(h.long,  'XData', x(l), 'YData', y(l));
set_datatip_rows(h.short, plots(s));
set_datatip_rows(h.long,  plots(l));
S.vis = struct('x', x(show), 'y', y(show), 'idx', find(show));

% CFAR detections (before the plot extractor): an arc over the azimuth step of the CPI update
D = d.dets;
xd = zeros(1, 0);  yd = xd;  showD = false(0, 1);
if showDets && ~isempty(D)
    el = D(:, 4);
    el(isnan(el)) = 0;
    g = D(:, 3) .* cosd(el) / 1e3;
    xd = (g .* sind(D(:, 2))).';
    yd = (g .* cosd(D(:, 2))).';
    showD = ((D(:, 7) == 1 & o.showShort) | (D(:, 7) == 2 & o.showLong)) & D(:, 6) >= o.minSnr_dB ...
          & D(:, 1) > o.clearTime_s;
    if o.history > 0
        showD = showD & (d.t_s - D(:, 1)) < o.history * 360 / sys.ant.rate_dps;
    end
    half = max(S.detArc_deg, 1e-3) / 2;
    hDet = [h.detShort h.detLong];
    for pul = 1:2
        k = find(showD & D(:, 7) == pul);
        az = D(k, 2).';
        gk = g(k).';
        X = [gk .* sind(az - half); gk .* sind(az + half); NaN(1, numel(k))];
        Y = [gk .* cosd(az - half); gk .* cosd(az + half); NaN(1, numel(k))];
        set(hDet(pul), 'XData', X(:).', 'YData', Y(:).', 'Visible', 'on');
    end
else
    set([h.detShort h.detLong], 'XData', NaN, 'YData', NaN, 'Visible', 'off');
end
S.visDet = struct('x', xd(showD), 'y', yd(showD), 'idx', find(showD).');

% Tracks: hollow triangle at the last update (around its plot), velocity leader, number
S = draw_tracks(S);
S = draw_rd(S);
show_timing(S);

% Labels next to the plots
delete(h.labels(ishandle(h.labels)));
S.h.labels = [];
if o.showLabels
    idx = find(show);
    for k = idx(1:min(end, 300))
        S.h.labels(end+1) = text(x(k), y(k), sprintf('  %.1f km, %+.1f m/s', plots(k).range_m/1e3, ...
            plots(k).vr_mps), 'Parent', S.ax, 'Color', [0.8 0.8 0.8], 'FontSize', 8, 'HitTest', 'off', 'Clipping', 'on');
    end
end
counts = '';
if showPlots, counts = sprintf('plots %d / %d   ', sum(show), nAll); end
if showDets, counts = [counts sprintf('detections %d   ', sum(showD))]; end
if showVideo
    what = 'video';
    if ~strcmp(o.showMode, 'video'), what = sprintf('%s out, beam %d', upper(o.showMode), o.beam); end
    counts = sprintf('%s (noise %.1f dB)   ', what, d.video.noise_dB);
end
if o.showTracks && ~isempty(d.tracks)
    counts = [counts sprintf('tracks %d   ', sum([d.tracks.confirmed]))];
end
set(h.title, 'String', sprintf('Scan %d   t = %.2f s   az = %5.1f deg   %szoom x%.1f', ...
    d.scan, d.t_s, d.az_deg, counts, zoom_factor(S)));
end

function S = draw_video(S, on)
% Raw data (video, DBF out, MF out): the map as an image (transparent where there is no data), and
% its colour scale. The map must hold the data of the selected mode (and beam).
h = S.h;
vis = {'off', 'on'};
set([h.video h.videoBar h.videoBarImg h.videoBarLabel], 'Visible', vis{1 + on});
if ~on || ~isfield(S.data, 'video') || isempty(S.data.video), return; end
V = S.data.video;
names = struct('video', 'video', 'dbf', 'DBF out', 'mf', 'MF out');
if strcmp(S.opt.showMode, 'video')
    set(h.videoBarLabel, 'String', 'video [dB above noise]');
else
    set(h.videoBarLabel, 'String', sprintf('%s, beam %d [dB above noise]', names.(S.opt.showMode), S.opt.beam));
end
if ~isfield(V, 'kind') || ~strcmp(V.kind, S.opt.showMode) || (~strcmp(V.kind, 'video') && V.beam ~= S.opt.beam)
    set(h.video, 'AlphaData', zeros(size(get(h.video, 'AlphaData'))));   % other data: not shown
    return
end
tMin = S.opt.clearTime_s;
if S.opt.history > 0
    tMin = max(tMin, S.data.t_s - S.opt.history * 360 / S.sys.ant.rate_dps);
end
[rgb, alpha] = ppi_video_image(S.data.video, tMin);
set(h.video, 'CData', rgb, 'AlphaData', alpha, 'AlphaDataMapping', 'none');
end

function S = layout_tab(S)
% Tab 1: the PPI over the whole display area. Tab 2: the range-Doppler map of the selected beam
% with its colour scale and the Doppler spectrum on the left, the PPI small on the right.
h = S.h;
two = S.opt.tab == 2;
vis = {'off', 'on'};
set(h.tabBtn(1), 'Value', ~two, 'BackgroundColor', [0.12 0.12 0.12] + 0.13*~two);
set(h.tabBtn(2), 'Value', two, 'BackgroundColor', [0.12 0.12 0.12] + 0.13*two);
set([h.rdAx h.rdBar h.specAx], 'Visible', vis{1 + two});
set([h.rdImg h.rdNotch h.rdDet h.rdCut h.rdInfo h.rdBarImg h.specLine h.specDet], 'Visible', vis{1 + two});
labels = [get(h.rdAx, 'XLabel') get(h.rdAx, 'YLabel') h.rdTitle get(h.specAx, 'XLabel') get(h.specAx, 'YLabel') h.specTitle];
set(labels, 'Visible', vis{1 + two});                      % (Octave shows them also on hidden axes)
if two
    set(S.ax, 'Position', [0.56 0.22 0.19 0.62]);
    set(h.videoBar, 'Position', [0.59 0.87 0.15 0.012]);
    set(h.legend, 'Visible', 'off');
    set(h.title, 'FontSize', 8);
else
    set(S.ax, 'Position', h.axPos);
    set(h.videoBar, 'Position', [0.57 0.872 0.15 0.012]);
    set(h.legend, 'Visible', 'on');
    set(h.title, 'FontSize', get(0, 'DefaultAxesFontSize') * 1.1);
end
end

function S = draw_rd(S)
% Range-Doppler map of the selected beam (tab 2): dB above the noise floor of the map (its
% median), -3 dB (black) ... the strongest cell (at least 20 dB); the beam's CFAR detections
% (circles), the MTI notch (dashed), and the Doppler spectrum of the range bin under the cursor
% or else of the strongest detection (strongest cell).
h = S.h;
if S.opt.tab ~= 2, return; end
rd = S.data.rd;
if isempty(rd) || rd.beam ~= S.opt.beam || all(isnan(rd.P(:)))
    set(h.rdTitle, 'String', sprintf('Range-Doppler map, beam %d: waiting for a CPI update ...', S.opt.beam));
    set([h.rdImg h.specLine h.specDet h.rdDet h.rdNotch h.rdCut], 'Visible', 'off');
    return
end
cols = find(any(~isnan(rd.P), 2));
b0 = cols(1);  b1 = cols(end);
dB = 10*log10(rd.P(b0:b1, :).' / rd.noise);              % [nfft x nBins]
hi = max(20, ceil(max(dB(:)) / 5) * 5);  lo = -3;
nB = b1 - b0 + 1;
W = min(nB, 600);                                          % columns: the largest bin of each (peaks stay)
if W < nB
    e = floor((0:W) * nB / W);
    D = NaN(size(dB, 1), W);
    for c = 1:W, D(:, c) = max(dB(:, e(c)+1:max(e(c+1), e(c)+1)), [], 2); end
else
    D = dB;
end
cmap = S.data.video.cmap;
k = round((D - lo) / (hi - lo) * 255) + 1;
k = min(max(k, 1), 256);
rgb = reshape(cmap(k(:), :), [size(D) 3]);
rgb(repmat(isnan(D), [1 1 3])) = 0;
x = [b0 - 0.5, b1 - 0.5] * rd.dR_m / 1e3;  y = rd.vel_mps([1 end]);
set(h.rdImg, 'CData', rgb, 'XData', x, 'YData', y, 'Visible', 'on');
dv = abs(diff(rd.vel_mps(1:2)));
set(h.rdAx, 'XLim', [b0 - 1, b1] * rd.dR_m / 1e3, 'YLim', sort(y) + [-dv dv]/2);
z = floor(numel(rd.vel_mps)/2) + 1;                        % zero Doppler
yn = rd.vel_mps(z) + [1 -1] * (rd.notchBins + 0.5) * dv;
xl = get(h.rdAx, 'XLim');
set(h.rdNotch, 'XData', [xl NaN xl], 'YData', [yn(1) yn(1) NaN yn(2) yn(2)], 'Visible', 'on');
if ~isempty(rd.dets)
    set(h.rdDet, 'XData', rd.dets(:, 1).' / 1e3, 'YData', rd.vel_mps(rd.dets(:, 2)), 'Visible', 'on');
else
    set(h.rdDet, 'XData', NaN, 'YData', NaN);
end
set(h.rdBarImg, 'XData', [0 1], 'YData', [lo hi], 'Visible', 'on');
set(h.rdBar, 'YLim', [lo hi], 'YTick', 0:10:hi);
set(h.rdTitle, 'FontSize', 9, 'String', sprintf('Range-Doppler, beam %d   t = %.2f s   az %.1f deg   peak %.1f dB', ...
    rd.beam, rd.t_s, rd.az_deg, max(dB(:))));
% Doppler spectrum: the bin under the cursor, else the strongest detection, else the strongest cell
bin = S.rdCursorBin;
if isempty(bin) || bin < b0 || bin > b1
    if ~isempty(rd.dets)
        [~, i] = max(rd.dets(:, 3));
        bin = floor(rd.dets(i, 1) / rd.dR_m) + 1;
    else
        [~, i] = max(rd.P(:));
        [bin, ~] = ind2sub(size(rd.P), i);
    end
end
s = 10*log10(rd.P(bin, :) / rd.noise);
set(h.specLine, 'XData', rd.vel_mps, 'YData', s, 'Visible', 'on');
m = [];
if ~isempty(rd.dets), m = rd.dets(floor(rd.dets(:, 1) / rd.dR_m) + 1 == bin, 2); end
set(h.specDet, 'XData', rd.vel_mps(m), 'YData', s(m), 'Visible', 'on');
set(h.specAx, 'XLim', sort(rd.vel_mps([1 end])), 'YLim', [-10 max(hi, 10)]);
set(h.rdCut, 'XData', [1 1] * (bin - 0.5) * rd.dR_m / 1e3, 'YData', get(h.rdAx, 'YLim'), 'Visible', 'on');
src = 'strongest peak';  if ~isempty(S.rdCursorBin), src = 'cursor'; end
set(h.specTitle, 'String', sprintf('Doppler spectrum at %.2f km (%s), circles: detections', ...
    (bin - 0.5) * rd.dR_m / 1e3, src));
end

function [S, txt] = rd_hover(S, p)
% Readout of the range-Doppler cell at p = [range_km velocity]; the spectrum follows the cursor.
txt = '';
rd = S.data.rd;
if isempty(rd) || S.opt.tab ~= 2, return; end
bin = floor(p(1) * 1e3 / rd.dR_m) + 1;
[~, d] = min(abs(rd.vel_mps - p(2)));
if bin < 1 || bin > size(rd.P, 1), return; end
v = rd.P(bin, d);
if isnan(v)
    txt = sprintf('R %.2f km   v %.1f m/s   no cell', (bin - 0.5) * rd.dR_m / 1e3, rd.vel_mps(d));
else
    txt = sprintf('R %.2f km   v %.1f m/s   %.1f dB above noise', (bin - 0.5) * rd.dR_m / 1e3, rd.vel_mps(d), ...
                  10*log10(v / rd.noise));
end
set(S.h.rdInfo, 'String', txt);
S.rdCursorBin = bin;
S = draw_rd(S);
end

function set_datatip_rows(hLine, plots)
% MATLAB data tips (axes toolbar) with the plot values (not in Octave).
if isempty(plots) || ~isprop(hLine, 'DataTipTemplate'), return; end
try
    tpl = hLine.DataTipTemplate;
    el = [plots.el_deg];
    tpl.DataTipRows = [dataTipTextRow('Range [km]', [plots.range_m]/1e3, '%.3f'), ...
                       dataTipTextRow('Azimuth [deg]', [plots.az_deg], '%.2f'), ...
                       dataTipTextRow('Elevation [deg]', el, '%.2f'), ...
                       dataTipTextRow('Radial vel. [m/s]', [plots.vr_mps], '%.2f'), ...
                       dataTipTextRow('SNR [dB]', [plots.snr_dB], '%.1f'), ...
                       dataTipTextRow('Time [s]', [plots.time_s], '%.3f'), ...
                       dataTipTextRow('CPI hits', [plots.hits], '%d')];
catch
end
end

%% ------------------------------------------------------------------------
function [S, txt] = hover(S, p)
% Data tip of the plot or true target nearest to p (within 8 pixels).
txt = '';
[txt, pos] = nearest_text(S, p);
if isempty(txt)
    set(S.h.tip, 'Visible', 'off');
    set(S.h.sel, 'XData', NaN, 'YData', NaN);
    return
end
k = km_per_pixel(S, get(S.ax, 'XLim'), get(S.ax, 'YLim'));
xl = get(S.ax, 'XLim');
if pos(1) > mean(xl)                                       % keep the tip inside the axes
    set(S.h.tip, 'Position', [pos(1) - 10*k, pos(2) + 6*k, 0], 'HorizontalAlignment', 'right');
else
    set(S.h.tip, 'Position', [pos(1) + 10*k, pos(2) + 6*k, 0], 'HorizontalAlignment', 'left');
end
set(S.h.tip, 'String', txt, 'Visible', 'on');
set(S.h.sel, 'XData', pos(1), 'YData', pos(2));
end

function [txt, pos] = nearest_text(S, p)
% Text and position of the nearest track, plot, detection or true target, '' if none is close.
% A track (its triangle surrounds its plot) wins unless another object is 3 pixels closer.
txt = '';  pos = [];
if isempty(p), return; end
kpx = km_per_pixel(S, get(S.ax, 'XLim'), get(S.ax, 'YLim'));
best = 8 * kpx;
[dist, j] = min(hypot(S.visTrk.x - p(1), S.visTrk.y - p(2)));   % tracks
if ~isempty(dist) && dist <= best
    best = max(dist - 3*kpx, 0);
    pos = [S.visTrk.x(j) S.visTrk.y(j)];
    txt = track_text(S.data.tracks(S.visTrk.idx(j)), S.sys.ant.height_m);
end
[dist, j] = min(hypot(S.vis.x - p(1), S.vis.y - p(2)));     % plots
if ~isempty(dist) && dist <= best
    best = dist;
    pos = [S.vis.x(j) S.vis.y(j)];
    txt = plot_text(S.data.plots(S.vis.idx(j)));
end
[dist, j] = min(hypot(S.visDet.x - p(1), S.visDet.y - p(2)));   % detections (plots first)
if ~isempty(dist) && dist < best
    best = dist;
    pos = [S.visDet.x(j) S.visDet.y(j)];
    txt = detection_text(S.data.dets(S.visDet.idx(j), :));
end
if S.opt.showTruth                                         % true targets
    xt = get(S.h.truth, 'XData');  yt = get(S.h.truth, 'YData');
    for i = 1:numel(xt)
        dist = hypot(xt(i) - p(1), yt(i) - p(2));
        if dist < best
            best = dist;
            pos = [xt(i) yt(i)];
            txt = truth_text(S, i);
        end
    end
end
end

function txt = plot_text(p)
if isnan(p.el_deg), el = '- (no DBF)'; else, el = sprintf('%.2f deg', p.el_deg); end
txt = sprintf(['Plot (%s pulse)\nRange       %.3f km\nAzimuth     %.2f deg\nElevation   %s\n' ...
               'Radial vel. %.2f m/s\nSNR         %.1f dB\nTime        %.3f s\nCPI hits    %d\n' ...
               'Beam %d, range cell %d, Doppler bin %d'], p.pulse, p.range_m/1e3, p.az_deg, el, ...
               p.vr_mps, p.snr_dB, p.time_s, p.hits, p.beam, p.rIdx, p.dIdx);
end

function txt = detection_text(r)
% r = [time az range el vr snr pulse beam rIdx dIdx] (detection_log_init.m)
names = {'short', 'long'};
if isnan(r(4)), el = '- (no DBF)'; else, el = sprintf('%.2f deg', r(4)); end
txt = sprintf(['CFAR detection (%s pulse)\nRange       %.3f km\nAzimuth     %.2f deg\nElevation   %s\n' ...
               'Radial vel. %.2f m/s\nSNR         %.1f dB\nTime        %.3f s\n' ...
               'Beam %d, range cell %d, Doppler bin %d'], names{r(7)}, r(3)/1e3, r(2), el, r(5), r(6), r(1), ...
               r(8), r(9), r(10));
end

function txt = track_text(tk, height_m)
% Data tip of a track.
states = {'tentative', 'confirmed'};
p = tk.pos_m;  v = tk.vel_mps;
R = norm(p(1:tk.dim));
method = 'EKF';
if isfield(tk, 'method') && ~isempty(tk.method), method = tk.method; end
txt = sprintf(['Track %d (%s, %d-D %s)\nRange       %.3f km\nAzimuth     %.2f deg\n'], tk.id, ...
              states{1 + tk.confirmed}, tk.dim, method, R/1e3, mod(atan2(p(1), p(2))*180/pi, 360));
if tk.dim == 3
    txt = [txt sprintf('Height      %.0f m\nClimb       %.1f m/s\n', p(3) + height_m, v(3))];
end
txt = [txt sprintf('Speed       %.1f m/s\nHeading     %.1f deg\n', norm(v(1:2)), mod(atan2(v(1), v(2))*180/pi, 360))];
txt = [txt sprintf('Radial vel. %.1f m/s\n', p(1:tk.dim)*v(1:tk.dim).'/R)];   % from the track: unambiguous
if isfield(tk, 'acc_mps2') && any(tk.acc_mps2 ~= 0)                     % constant-acceleration model
    txt = [txt sprintf('Accel.      %.2f m/s^2\n', norm(tk.acc_mps2(1:2)))];
end
if isfield(tk, 'existence') && isfinite(tk.existence)                 % PMBM: existence, IMM models
    txt = [txt sprintf('Existence   %.3f\n', tk.existence)];
    if isfield(tk, 'modeProb') && ~isempty(tk.modeProb)
        txt = [txt sprintf('Models      %s\n', tk.modeProb)];
    end
end
txt = [txt sprintf('Updates     %d (missed scans %d)\nLast update %.3f s', tk.hits, tk.misses, tk.t_s)];
end

function S = draw_tracks(S)
% Confirmed (magenta) and tentative (gray) tracks at their last update, with leaders and numbers.
h = S.h;  d = S.data;
delete(h.trackIds(ishandle(h.trackIds)));
S.h.trackIds = [];
tk = d.tracks;
if ~S.opt.showTracks || isempty(tk)
    set([h.track h.trackTent h.trackLead], 'XData', NaN, 'YData', NaN);
    S.visTrk = struct('x', [], 'y', [], 'idx', []);
    return
end
P = reshape([tk.pos_m], 3, []) / 1e3;                     % [km]
V = reshape([tk.vel_mps], 3, []) / 1e3;
conf = logical([tk.confirmed]);
set(h.track, 'XData', P(1, conf), 'YData', P(2, conf));
set(h.trackTent, 'XData', P(1, ~conf), 'YData', P(2, ~conf));
L = S.sys.cfg.display.trackLeader_s;
X = [P(1, conf); P(1, conf) + L*V(1, conf); NaN(1, sum(conf))];
Y = [P(2, conf); P(2, conf) + L*V(2, conf); NaN(1, sum(conf))];
set(h.trackLead, 'XData', X(:).', 'YData', Y(:).');
for i = find(conf)
    S.h.trackIds(end+1) = text(P(1, i), P(2, i), sprintf('   T%d', tk(i).id), 'Parent', S.ax, ...
        'Color', [1 0.35 1], 'FontSize', 8, 'VerticalAlignment', 'top', 'HitTest', 'off', 'Clipping', 'on');
end
S.visTrk = struct('x', P(1, :), 'y', P(2, :), 'idx', 1:numel(tk));
end

function show_timing(S)
% Mean time per PRI of every block, of the console output and of the display in the panel.
% The labels are at most 15 characters, so the value column stays inside the panel.
us = S.data.timing_us;
if isempty(us), return; end
names = {'DBF', 'Matched filter', 'MTI (clutter)', 'Doppler FFT+NCI', compact_label(['CFAR ' S.sys.cfar.label]), ...
         'Peak interp.', 'Plot extractor', tracker_label(S.sys.track), ...
         'Console output', 'Display (PPI)'};
if numel(us) == 9, names(6) = []; end                       % (block times without the peak interpolation)
n = min(numel(us), numel(names));
rows = cell(1, n + 3);
rows{1} = 'Mean time per PRI [us]';
for i = 1:n
    rows{i+1} = sprintf('%-15s%9.1f', names{i}(1:min(end, 15)), us(i));
end
rows{n+2} = sprintf('%-15s%9.1f', 'Total', sum(us(1:n)));       % incl. console and display
if isfield(S.data, 'timingPris') && ~isempty(S.data.timingPris)
    rows{n+3} = sprintf('over the last %d PRIs', S.data.timingPris);
else
    rows{n+3} = 'over the whole run';
end
set(S.h.timing, 'String', rows);
end

function txt = truth_text(S, i)
tgt = S.sys.cfg.targets(i);
[R, az, el, vr] = target_geometry(tgt, S.data.t_s, S.sys.ant.height_m);
[pos, vel, acc] = target_state(tgt, S.data.t_s);
txt = sprintf(['True target %s\nRange       %.3f km\nAzimuth     %.2f deg\nElevation   %.2f deg\n' ...
               'Height      %.0f m\nRadial vel. %.1f m/s\nSpeed       %.1f m/s\nMean RCS    %.2f m^2'], ...
               tgt.name, R/1e3, az, el, pos(3), vr, norm(vel), tgt.rcs_m2);
if any(acc ~= 0)                                            % maneuvering now
    txt = [txt sprintf('\nAccel.      %.1f m/s^2', norm(acc))];
end
end

function S = pin(S, p)
% Pin the data tip of the nearest object, or remove a pinned tip at p.
k = km_per_pixel(S, get(S.ax, 'XLim'), get(S.ax, 'YLim'));
for i = 1:numel(S.h.pins)                                  % click on a pinned tip: remove it
    if hypot(S.h.pins{i}.x - p(1), S.h.pins{i}.y - p(2)) <= 8*k
        delete([S.h.pins{i}.txt S.h.pins{i}.mk]);
        S.h.pins(i) = [];
        return
    end
end
[txt, pos] = nearest_text(S, p);
if isempty(txt), return; end
pinTxt = text(pos(1) + 10*k, pos(2) - 6*k, txt, 'Parent', S.ax, 'Color', 'w', ...
              'BackgroundColor', [0.05 0.2 0.05], 'EdgeColor', [0.3 0.9 0.3], 'FontSize', 9, ...
              'VerticalAlignment', 'top', 'Margin', 4, 'HitTest', 'off');
mk = plot(S.ax, pos(1), pos(2), 's', 'Color', [0.3 0.9 0.3], 'MarkerSize', 13, 'LineWidth', 1.5, 'HitTest', 'off', 'Clipping', 'on');
S.h.pins{end+1} = struct('txt', pinTxt, 'mk', mk, 'x', pos(1), 'y', pos(2));
end

function S = clear_pins(S)
for i = 1:numel(S.h.pins)
    delete([S.h.pins{i}.txt S.h.pins{i}.mk]);
end
S.h.pins = {};
end

function S = key(S, evt)
xl = get(S.ax, 'XLim');  yl = get(S.ax, 'YLim');
step = 0.1 * diff(xl);
c = [mean(xl) mean(yl)];
switch evt.Key
    case {'add', 'equal', 'plus'},       S = zoom_at(S, 1.5, c);
    case {'subtract', 'hyphen', 'minus'}, S = zoom_at(S, 1/1.5, c);
    case 'leftarrow',  S = set_limits(S, xl - step, yl);
    case 'rightarrow', S = set_limits(S, xl + step, yl);
    case 'uparrow',    S = set_limits(S, xl, yl + step);
    case 'downarrow',  S = set_limits(S, xl, yl - step);
    case {'r', 'home'}, S = set_limits(S, S.home(1:2), S.home(3:4));
    case 'c',           S = clear_pins(S);
    case 'd'                                                % next display mode
        modes = {'plots', 'detections', 'both', 'video', 'dbf', 'mf'};
        S.opt.showMode = modes{mod(find(strcmp(S.opt.showMode, modes)), 6) + 1};
        set(S.h.ui.showMode, 'Value', find(strcmp(S.opt.showMode, modes)));
        S = redraw(S);
    case 'space'
        S.opt.paused = ~S.opt.paused;
        labels = {'Pause display', 'Resume display'};
        set(S.h.ui.pause, 'String', labels{1 + S.opt.paused});
        if ~S.opt.paused, S = redraw(S); end
end
end

%% ------------------------------------------------------------------------
function S = zoom_at(S, factor, c)
% Zoom by factor around the point c [km] (the point stays under the cursor).
xl = get(S.ax, 'XLim');  yl = get(S.ax, 'YLim');
span = diff(S.home(1:2));
f = max(min(factor, diff(xl)/(span/500)), diff(xl)/span); % between full view and x500
S = set_limits(S, c(1) + (xl - c(1))/f, c(2) + (yl - c(2))/f);
end

function S = set_limits(S, xl, yl)
set(S.ax, 'XLim', xl, 'YLim', yl);
set(S.h.title, 'String', regexprep(get(S.h.title, 'String'), 'zoom x[\d.]+', ...
    sprintf('zoom x%.1f', zoom_factor(S))));
end

function z = zoom_factor(S)
z = diff(S.home(1:2)) / diff(get(S.ax, 'XLim'));
end

function k = km_per_pixel(S, xl, yl)
% Data units per screen pixel (axis equal: the larger of the two scales).
pos = getpixelposition(S.ax);
k = max(diff(xl)/pos(3), diff(yl)/pos(4));
end

function p = cursor_km(S)
% Cursor position in the axes [km], [] when it is outside the axes.
cp = get(S.ax, 'CurrentPoint');
p = cp(1, 1:2);
xl = get(S.ax, 'XLim');  yl = get(S.ax, 'YLim');
if p(1) < xl(1) || p(1) > xl(2) || p(2) < yl(1) || p(2) > yl(2), p = []; end
end

function S = show_cursor(S, p)
if isempty(p)
    set(S.h.cursor, 'String', '');
    return
end
az = mod(atan2(p(1), p(2)) * 180/pi, 360);
set(S.h.cursor, 'String', {sprintf('Cursor: %.2f km, %.1f deg', hypot(p(1), p(2)), az), ...
                           sprintf('E %.2f km, N %.2f km', p(1), p(2))});
end

function px = figure_px(fig)
% Mouse position in the figure in pixels.
units = get(fig, 'Units');
set(fig, 'Units', 'pixels');
px = get(fig, 'CurrentPoint');
set(fig, 'Units', units);
px = px(1:2);
end

function s = tracker_label(p)
% Name of the tracking block in the timing panel (at most 15 characters).
if ~p.enable
    s = 'Track (off)';
elseif isfield(p, 'method') && strcmp(p.method, 'pmbm')
    s = sprintf('Track %dD PMBM', p.dim);
else
    s = sprintf('Track %dD %s', p.dim, upper(p.model));
end
end

function s = compact_label(s)
% A label of the timing panel in at most 15 characters: '1-D' -> '1D', then
% 'range' -> 'rng', 'Doppler' -> 'Dop', then cut.
s = strrep(s, '-D', 'D');
if numel(s) > 15, s = strrep(strrep(s, ' range', ' rng'), ' Doppler', ' Dop'); end
s = s(1:min(end, 15));
end

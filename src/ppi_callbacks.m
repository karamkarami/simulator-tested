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
%     'option', name, value       showMode ('plots', 'detections' or 'both'),
%                                 showShort, showLong, showTruth, showRings,
%                                 showLabels, showAntenna, showTracks,
%                                 minSnr_dB, history
%     'timing'                    write the block processing times into the panel
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
            set(S.h.ui.showMode, 'Value', find(strcmp(varargin{2}, {'plots', 'detections', 'both'})));
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
            cfar_pfa(10^(dB/10), S.sys.cfar.nTrain, S.sys.dop.nLooks)));
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
        p = targets(i).pos0_m + targets(i).vel_mps*d.t_s;
        xy(i, :) = p(1:2) / 1e3;
        p0 = targets(i).pos0_m / 1e3;
        set(h.names(i), 'Position', [xy(i, :) 0], 'Visible', 'on');
        set(h.tracks(i), 'XData', [p0(1) xy(i, 1)], 'YData', [p0(2) xy(i, 2)], 'Visible', 'on');
    end
    set(h.truth, 'XData', xy(:, 1), 'YData', xy(:, 2), 'Visible', 'on');
else
    set(h.truth, 'Visible', 'off');
    set(h.names, 'Visible', 'off');
    set(h.tracks, 'Visible', 'off');
end

% Plots that pass the filters (shown unless only the detections are selected)
showPlots = ~strcmp(o.showMode, 'detections');
showDets  = ~strcmp(o.showMode, 'plots');
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
if o.showTracks && ~isempty(d.tracks)
    counts = [counts sprintf('tracks %d   ', sum([d.tracks.confirmed]))];
end
set(h.title, 'String', sprintf('Scan %d   t = %.2f s   az = %5.1f deg   %szoom x%.1f', ...
    d.scan, d.t_s, d.az_deg, counts, zoom_factor(S)));
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
txt = sprintf(['Track %d (%s, %d-D EKF)\nRange       %.3f km\nAzimuth     %.2f deg\n'], tk.id, ...
              states{1 + tk.confirmed}, tk.dim, R/1e3, mod(atan2(p(1), p(2))*180/pi, 360));
if tk.dim == 3
    txt = [txt sprintf('Height      %.0f m\nClimb       %.1f m/s\n', p(3) + height_m, v(3))];
end
txt = [txt sprintf('Speed       %.1f m/s\nHeading     %.1f deg\n', norm(v(1:2)), mod(atan2(v(1), v(2))*180/pi, 360))];
if isfield(tk, 'acc_mps2') && any(tk.acc_mps2 ~= 0)                     % constant-acceleration model
    txt = [txt sprintf('Accel.      %.2f m/s^2\n', norm(tk.acc_mps2(1:2)))];
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
us = S.data.timing_us;
if isempty(us), return; end
names = {'DBF', 'Matched filter', 'MTI (clutter)', 'Doppler FFT+NCI', ['CFAR ' S.sys.cfar.label], ...
         'Plot extractor', sprintf('Tracking %d-D %s', S.sys.track.dim, upper(S.sys.track.model)), ...
         'Console output', 'Display (PPI)'};
if ~S.sys.track.enable, names{7} = 'Tracking (off)'; end
n = min(numel(us), numel(names));
rows = cell(1, n + 2);
rows{1} = 'Mean time per PRI [us]';
for i = 1:n
    rows{i+1} = sprintf('%-17s%9.1f', names{i}, us(i));
end
rows{n+2} = sprintf('%-17s%9.1f', 'Total', sum(us(1:n)));
set(S.h.timing, 'String', rows);
end

function txt = truth_text(S, i)
tgt = S.sys.cfg.targets(i);
[R, az, el, vr] = target_geometry(tgt, S.data.t_s, S.sys.ant.height_m);
pos = tgt.pos0_m + tgt.vel_mps * S.data.t_s;
txt = sprintf(['True target %s\nRange       %.3f km\nAzimuth     %.2f deg\nElevation   %.2f deg\n' ...
               'Height      %.0f m\nRadial vel. %.1f m/s\nSpeed       %.1f m/s\nMean RCS    %.2f m^2'], ...
               tgt.name, R/1e3, az, el, pos(3), vr, norm(tgt.vel_mps), tgt.rcs_m2);
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
        modes = {'plots', 'detections', 'both'};
        S.opt.showMode = modes{mod(find(strcmp(S.opt.showMode, modes)), 3) + 1};
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

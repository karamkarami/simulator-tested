function ppi = ppi_init(sys)
%PPI_INIT  Create the interactive PPI (plan position indicator) display.
%   ppi = PPI_INIT(sys) creates the figure. ppi_update.m draws the antenna,
%   the true targets and the plots; ppi_callbacks.m handles the interaction:
%     mouse wheel          zoom in/out around the cursor
%     left drag            pan
%     double click         reset the view (whole instrumented range)
%     hover over a plot    data tip with range, azimuth, elevation, radial
%                          velocity, SNR, time, sub-pulse, hits, beam
%     click on a plot      pin its data tip (click again to remove it)
%     hover over a track   data tip with number, state, position, height (3-D),
%                          speed, heading, updates
%     keys                 + / - zoom, arrows pan, r reset, c clear the
%                          pinned data tips, space pause/resume the display,
%                          d next display mode
%   The panel on the right selects what is shown: the plots (after the plot
%   extractor), the CFAR detections before the extractor (each detection is
%   an arc over the azimuth step of its CPI update, so a target appears as an
%   arc about one beamwidth long) or both; shows or hides the short/long-pulse
%   data, the true targets, range rings, plot labels and the antenna; filters
%   by SNR and age; shows or hides the tracks (hollow triangles: confirmed
%   magenta, tentative gray, with a velocity leader); shows the mean
%   time per PRI of every block, of the console output and of the display
%   [us]; sets the CFAR threshold [dB]
%   of the selected CFAR method of the running processing
%   (Apply; Default = the value from the method's Pfa); clears the display; and
%   saves a screenshot (PNG). The MATLAB zoom/pan/data tip tools of the axes
%   toolbar also work.
cfg  = sys.cfg;
rMax = cfg.display.maxRange_m;
if isempty(rMax), rMax = cfg.rx.maxRange_m; end
rMax = rMax / 1e3;                                         % km
if cfg.display.visible, vis = 'on'; else, vis = 'off'; end

fig = figure('Name', 'PPI', 'Color', 'k', 'Visible', vis, 'InvertHardcopy', 'off', ...
             'NumberTitle', 'off', 'Position', [100 60 1040 880]);
ax  = axes('Parent', fig, 'Color', 'k', 'XColor', [0.6 0.6 0.6], 'YColor', [0.6 0.6 0.6], ...
           'Units', 'normalized', 'Position', [0.07 0.12 0.66 0.8]);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, rMax*[-1.05 1.05 -1.05 1.05]);
xlabel(ax, 'East [km]');
ylabel(ax, 'North [km]');
if exist('disableDefaultInteractivity', 'file')
    disableDefaultInteractivity(ax);                       % our own wheel zoom / drag pan
end

%% Range rings and azimuth spokes
ringStep = 10^floor(log10(rMax/2));
if rMax/ringStep > 6, ringStep = 2*ringStep; end
phi  = linspace(0, 2*pi, 361);
gray = [0.25 0.25 0.25];
h.rings = [];
for r = ringStep:ringStep:rMax
    h.rings(end+1) = plot(ax, r*sin(phi), r*cos(phi), '-', 'Color', gray, 'HitTest', 'off', 'Clipping', 'on');
    h.rings(end+1) = text(r*sind(45), r*cosd(45), sprintf('%g km', r), 'Color', [0.5 0.5 0.5], ...
                          'Parent', ax, 'HitTest', 'off', 'Clipping', 'on');
end
for a = 0:30:330
    h.rings(end+1) = plot(ax, [0 rMax*sind(a)], [0 rMax*cosd(a)], '-', 'Color', gray, 'HitTest', 'off', 'Clipping', 'on');
    h.rings(end+1) = text(1.03*rMax*sind(a), 1.03*rMax*cosd(a), sprintf('%d', a), 'Color', [0.5 0.5 0.5], ...
                          'Parent', ax, 'HorizontalAlignment', 'center', 'HitTest', 'off', 'Clipping', 'on');
end

%% Antenna, true targets, plots
h.ant    = plot(ax, [0 0], [0 rMax], '-', 'Color', [0 0.8 0], 'LineWidth', 1.5, 'HitTest', 'off', 'Clipping', 'on');
h.truth  = plot(ax, NaN, NaN, 's', 'Color', [1 0.3 0.3], 'MarkerSize', 10);
h.names  = gobjects_compat(numel(cfg.targets));
h.tracks = gobjects_compat(numel(cfg.targets));
for i = 1:numel(cfg.targets)
    h.names(i)  = text(NaN, NaN, ['  ' cfg.targets(i).name], 'Color', [1 0.3 0.3], 'Parent', ax, 'HitTest', 'off', 'Clipping', 'on');
    h.tracks(i) = plot(ax, NaN, NaN, ':', 'Color', [1 0.3 0.3], 'HitTest', 'off', 'Clipping', 'on');   % true track since t = 0
end
% CFAR detections before the plot extractor: one arc segment per detection
h.detShort = plot(ax, NaN, NaN, '-', 'Color', [0 0.55 0.65], 'LineWidth', 2.5, 'HitTest', 'off', 'Clipping', 'on');
h.detLong  = plot(ax, NaN, NaN, '-', 'Color', [0.75 0.55 0], 'LineWidth', 2.5, 'HitTest', 'off', 'Clipping', 'on');
h.short = plot(ax, NaN, NaN, 'o', 'Color', [0 0.9 1], 'MarkerFaceColor', [0 0.9 1], 'MarkerSize', 5);
h.long  = plot(ax, NaN, NaN, 'o', 'Color', [1 0.9 0], 'MarkerFaceColor', [1 0.9 0], 'MarkerSize', 5);
% Tracks: hollow triangles around the plot (confirmed: magenta, tentative: gray), velocity leaders
h.trackLead = plot(ax, NaN, NaN, '-', 'Color', [0.9 0.3 0.9], 'LineWidth', 1, 'HitTest', 'off', 'Clipping', 'on');
h.trackTent = plot(ax, NaN, NaN, '^', 'Color', [0.55 0.55 0.55], 'MarkerSize', 14, 'LineWidth', 1, ...
                   'HitTest', 'off', 'Clipping', 'on');
h.track = plot(ax, NaN, NaN, '^', 'Color', [1 0.35 1], 'MarkerSize', 14, 'LineWidth', 1.5, ...
               'HitTest', 'off', 'Clipping', 'on');
h.trackIds = [];                                           % track numbers next to confirmed tracks
h.labels = [];                                             % plot labels (option)
hLeg = legend([h.truth h.short h.long h.detShort h.detLong h.track], {'true target', 'plot (short)', ...
       'plot (long)', 'detection (short)', 'detection (long)', 'track'}, ...
       'TextColor', 'w', 'Color', 'k', 'Location', 'southoutside', 'Orientation', 'horizontal');
try
    set(hLeg, 'AutoUpdate', 'off');                        % data tips and labels stay out of the legend
catch
end
h.title = title(ax, '', 'Color', 'w');

%% Data tip (hover) and selection marker
h.tip = text(0, 0, '', 'Parent', ax, 'Color', 'w', 'BackgroundColor', [0.12 0.12 0.12], ...
             'EdgeColor', [0.8 0.8 0.8], 'FontSize', 9, 'VerticalAlignment', 'bottom', ...
             'Visible', 'off', 'HitTest', 'off', 'Margin', 4);
h.sel = plot(ax, NaN, NaN, 'o', 'Color', 'w', 'MarkerSize', 13, 'LineWidth', 1.5, 'HitTest', 'off', 'Clipping', 'on');
h.pins = {};                                               % pinned data tips: {text, marker}

%% Control panel (right)
bg = [0.15 0.15 0.15];  fg = [0.9 0.9 0.9];
x0 = 0.76;  w = 0.22;  rowH = 0.03;
y = 0.95;
label = @(str, pos) uicontrol(fig, 'Style', 'text', 'String', str, 'Units', 'normalized', 'Position', pos, ...
                              'BackgroundColor', 'k', 'ForegroundColor', fg, 'HorizontalAlignment', 'left');
uicontrol(fig, 'Style', 'text', 'String', 'Display', 'Units', 'normalized', 'Position', [x0 y w rowH], ...
          'BackgroundColor', 'k', 'ForegroundColor', 'w', 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
y = y - rowH;
label('Show', [x0 y 0.05 rowH]);
modes = {'plots', 'detections', 'both'};
h.ui.showMode = uicontrol(fig, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [x0+0.05 y w-0.05 rowH], ...
    'String', {'Plots (after extractor)', 'Detections (before extractor)', 'Plots and detections'}, ...
    'Value', find(strcmp(cfg.display.showMode, modes)), ...
    'BackgroundColor', bg, 'ForegroundColor', fg, ...
    'Callback', @(src, evt) ppi_callbacks('option', fig, 'showMode', modes{get(src, 'Value')}));
opt.showMode = cfg.display.showMode;
checks = {'showShort', 'Short pulse', true;  'showLong', 'Long pulse', true;
          'showTruth', 'True targets', logical(cfg.display.showTruth);  'showRings', 'Range rings', true;
          'showLabels', 'Plot labels', false;  'showAntenna', 'Antenna', true;
          'showTracks', 'Tracks', logical(cfg.display.showTracks)};
for i = 1:size(checks, 1)
    y = y - rowH;
    h.ui.(checks{i, 1}) = uicontrol(fig, 'Style', 'checkbox', 'String', checks{i, 2}, 'Value', checks{i, 3}, ...
        'Units', 'normalized', 'Position', [x0 y w rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
        'Callback', @(src, evt) ppi_callbacks('option', fig, checks{i, 1}, logical(get(src, 'Value'))));
    opt.(checks{i, 1}) = checks{i, 3};
end
y = y - 1.3*rowH;
label('Min. SNR [dB]', [x0 y 0.12 rowH]);
h.ui.minSnr = uicontrol(fig, 'Style', 'edit', 'String', '0', 'Units', 'normalized', ...
    'Position', [x0+0.13 y 0.09 rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
    'Callback', @(src, evt) ppi_callbacks('option', fig, 'minSnr_dB', str2double(get(src, 'String'))));
opt.minSnr_dB = 0;
y = y - rowH;
label('Shown', [x0 y 0.12 rowH]);
h.ui.history = uicontrol(fig, 'Style', 'popupmenu', 'String', {'all', 'last scan', 'last 2 scans'}, ...
    'Value', 1, 'Units', 'normalized', 'Position', [x0+0.13 y 0.09 rowH], 'BackgroundColor', bg, ...
    'ForegroundColor', fg, 'Callback', @(src, evt) ppi_callbacks('option', fig, 'history', get(src, 'Value') - 1));
opt.history = 0;                                           % 0 all, 1 last scan, 2 last two scans

% CFAR threshold of the running processing
cfar_dB = 10*log10(sys.cfar.alpha);
y = y - 1.4*rowH;
label('CFAR thr. [dB]', [x0 y 0.12 rowH]);
h.ui.cfar = uicontrol(fig, 'Style', 'edit', 'String', sprintf('%.2f', cfar_dB), 'Units', 'normalized', ...
    'Position', [x0+0.13 y 0.09 rowH], 'BackgroundColor', bg, 'ForegroundColor', fg);
y = y - 1.1*rowH;
h.ui.cfarApply = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Apply', 'Units', 'normalized', ...
    'Position', [x0 y 0.105 rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
    'Callback', @(src, evt) ppi_callbacks('cfar', fig, str2double(get(h.ui.cfar, 'String'))));
h.ui.cfarDefault = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Default', 'Units', 'normalized', ...
    'Position', [x0+0.115 y 0.105 rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
    'Callback', @(src, evt) ppi_callbacks('cfar', fig, cfar_dB));
y = y - rowH;
h.ui.cfarInfo = label('', [x0 y w rowH]);                 % set by ppi_callbacks('cfarApplied')
set(h.ui.cfarInfo, 'ForegroundColor', [0.7 0.9 0.7]);

buttons = {'Reset view', 'reset';  'Clear display', 'clear';  'Remove data tips', 'cleartips';
           'Pause display', 'pause';  'Save screenshot (PNG) ...', 'save'};
y = y - 0.3*rowH;
for i = 1:size(buttons, 1)
    y = y - 1.15*rowH;
    h.ui.(buttons{i, 2}) = uicontrol(fig, 'Style', 'pushbutton', 'String', buttons{i, 1}, 'Units', 'normalized', ...
        'Position', [x0 y w rowH], 'BackgroundColor', bg, 'ForegroundColor', fg, ...
        'Callback', @(src, evt) ppi_callbacks(buttons{i, 2}, fig));
end
y = y - 2.1*rowH;
h.cursor = uicontrol(fig, 'Style', 'text', 'String', '', 'Units', 'normalized', 'Position', [x0 y w 2*rowH], ...
                     'BackgroundColor', 'k', 'ForegroundColor', [0.7 0.9 0.7], 'HorizontalAlignment', 'left');
% Mean processing time of every block (set by ppi_callbacks('redraw') from ppi_update),
% bottom right above the help
helpH = 0.085;
mono = 'FixedWidth';                                       % MATLAB's monospaced font
if exist('OCTAVE_VERSION', 'builtin'), mono = 'Monospace'; end
h.timing = uicontrol(fig, 'Style', 'text', 'String', '', 'Units', 'normalized', ...
    'Position', [x0 0.015+helpH w 0.19], 'BackgroundColor', [0.08 0.08 0.08], ...
    'ForegroundColor', [0.85 0.85 0.85], 'HorizontalAlignment', 'left', 'FontName', mono, 'FontSize', 8);
uicontrol(fig, 'Style', 'text', 'Units', 'normalized', 'Position', [x0 0.01 w helpH], ...
    'BackgroundColor', 'k', 'ForegroundColor', [0.6 0.6 0.6], 'HorizontalAlignment', 'left', ...
    'String', {'Wheel: zoom   Drag: pan', 'Double click: reset', 'Hover: tip   Click: pin', ...
               'Keys: + - arrows r c d space'});

%% State (used by ppi_update and ppi_callbacks)
opt.clearTime_s = -Inf;                                    % plots up to this time are cleared
opt.paused = false;
opt.cfarDefault_dB = cfar_dB;
opt.cfarRequest_dB = [];                                   % threshold to apply (ppi_update hands it over)
S.h = h;
S.ax = ax;
S.rMax = rMax;
S.home = rMax*[-1.05 1.05 -1.05 1.05];
S.opt = opt;
S.sys = sys;
S.data = struct('az_deg', 0, 't_s', 0, 'scan', 1, 'plots', [], 'dets', zeros(0, 10), ...
                'tracks', [], 'timing_us', []);
S.vis = struct('x', [], 'y', [], 'idx', []);               % visible plots (for picking)
S.visDet = struct('x', [], 'y', [], 'idx', []);            % visible detections
S.visTrk = struct('x', [], 'y', [], 'idx', []);            % visible tracks
S.detArc_deg = sys.ant.rate_dps * mean(sys.pri) * cfg.proc.updateEvery;   % azimuth step of a CPI update
S.drag = [];
setappdata(fig, 'ppi', S);

set(fig, 'WindowScrollWheelFcn',  @(src, evt) ppi_callbacks('scroll', fig, evt), ...
         'WindowButtonDownFcn',   @(src, evt) ppi_callbacks('down', fig), ...
         'WindowButtonMotionFcn', @(src, evt) ppi_callbacks('move', fig), ...
         'WindowButtonUpFcn',     @(src, evt) ppi_callbacks('up', fig), ...
         'KeyPressFcn',           @(src, evt) ppi_callbacks('key', fig, evt));
ppi_callbacks('cfarApplied', fig, cfar_dB);

ppi.fig = fig;
ppi.ax = ax;
ppi.rMax = rMax;
end

function g = gobjects_compat(n)
% Empty graphics handle array (gobjects in MATLAB, doubles in Octave).
if exist('gobjects', 'builtin') || exist('gobjects', 'file')
    g = gobjects(1, n);
else
    g = zeros(1, n);
end
end

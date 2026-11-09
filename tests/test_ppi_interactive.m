function ok = test_ppi_interactive()
%TEST_PPI_INTERACTIVE  Interaction with the PPI: data tips, zoom, pan, filters.
%   Creates the PPI off-screen with a few synthetic plots and drives it
%   through ppi_callbacks.m like the mouse, keyboard and panel would:
%   hover data tips (plot and true target), wheel zoom around the cursor,
%   pan, keys, show/hide of short- and long-pulse plots, SNR and age
%   filters, clear, pinned data tips, labels, reset and PNG export; the
%   display of the CFAR detections before the plot extractor (arcs, data
%   tips, filters, key d) and the CFAR threshold field (request, hand-over
%   to the processing through ppi_update, invalid input, default), and the
%   video mode (ppi_video_add.m, image of the sector the antenna passed),
%   DBF out / MF out of the selected beam (beam selector, key d) and the
%   range-Doppler tab (map, detections, spectrum, hover, small PPI).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
if exist('OCTAVE_VERSION', 'builtin') && isempty(available_graphics_toolkits())
    fprintf('Interactive PPI test: SKIPPED (no graphics toolkit: run octave with a display, e.g. xvfb-run octave)\n');
    ok = true;
    return
end
cfg = radar_config();
cfg.display.visible = false;
cfg.targets = make_target('T1', 1, 20e3, 60, 1000, 100, 0);
sys = setup_radar(cfg);
ppi = ppi_init(sys);
fig = ppi.fig;

% Three plots: short at 10 km / 30 deg, long at 40 km / 200 deg (low SNR, old)
p = struct('time_s', {2.0, 0.5, 2.1}, 'az_deg', {30, 200, 95}, 'range_m', {10e3, 40e3, 30e3}, ...
           'el_deg', {NaN, 3, 5}, 'vr_mps', {12.5, -7.25, 3}, 'snr_dB', {25, 14, 30}, 'power', 1, ...
           'beam', {3, 4, 5}, 'rIdx', {100, 200, 150}, 'dIdx', {10, 20, 30}, ...
           'pulse', {'short', 'long', 'long'}, 'hits', {80, 60, 70});
ppi_update(ppi, sys, 100, 12.0, 2, p);
pos = getpixelposition(ppi.ax);
xl0 = get(ppi.ax, 'XLim');
tol = 8 * max(diff(xl0)/pos(3), diff(get(ppi.ax, 'YLim'))/pos(4));
txtNear = ppi_callbacks('hover', fig, plot_xy(p(3)) + [0.75*tol 0]);
xy = @(k) plot_xy(p(k));
S = getappdata(fig, 'ppi');
ok = true;
check = @(cond, msg) report(cond, msg);

%% Data tips
txt = ppi_callbacks('hover', fig, xy(1) + [0.05 0.05]);
ok = check(~isempty(strfind(txt, 'short pulse')) && ~isempty(strfind(txt, '10.000 km')) ...
           && ~isempty(strfind(txt, '12.50 m/s')) && ~isempty(strfind(txt, 'no DBF')), ...
           'hover over the short-pulse plot shows its values') && ok;
ok = check(strcmp(get(S.h.tip, 'Visible'), 'on'), 'data tip is visible') && ok;
txt = ppi_callbacks('hover', fig, [-60 -60]);
ok = check(isempty(txt) && strcmp(get(S.h.tip, 'Visible'), 'off'), 'no data tip far from all plots') && ok;
tgtXY = (cfg.targets.pos0_m(1:2) + cfg.targets.vel_mps(1:2)*12) / 1e3;
txt = ppi_callbacks('hover', fig, tgtXY);
ok = check(~isempty(strfind(txt, 'True target T1')), 'hover over a true target') && ok;

ok = check(~isempty(strfind(txtNear, '30.000 km')), 'data tip within 8 px at the full view') && ok;

%% Zoom (wheel and keys) and pan
xl0 = get(ppi.ax, 'XLim');
ppi_callbacks('zoom', fig, 2, xy(3));
xl = get(ppi.ax, 'XLim');
c = xy(3);
ok = check(abs(diff(xl) - diff(xl0)/2) < 1e-9 && abs((c(1) - xl(1))/diff(xl) - (c(1) - xl0(1))/diff(xl0)) < 1e-9, ...
           'zoom x2 around a point keeps the point in place') && ok;
txt = ppi_callbacks('hover', fig, xy(3) + [0.75*tol 0]);   % 6 px at full view, 12 px zoomed
ok = check(isempty(txt), 'data tip tolerance (8 px) follows the zoom') && ok;
ppi_callbacks('pan', fig, [5 -3]);
ok = check(all(abs(get(ppi.ax, 'XLim') - (xl + 5)) < 1e-9), 'pan moves the view') && ok;
ppi_callbacks('key', fig, struct('Key', 'add', 'Character', '+'));
ok = check(abs(diff(get(ppi.ax, 'XLim')) - diff(xl)/1.5) < 1e-9, 'key + zooms in') && ok;
ppi_callbacks('key', fig, struct('Key', 'r', 'Character', 'r'));
ok = check(all(abs(get(ppi.ax, 'XLim') - xl0) < 1e-9), 'key r resets the view') && ok;
ppi_callbacks('zoom', fig, 1e6, [0 0]);
ok = check(abs(diff(xl0)/diff(get(ppi.ax, 'XLim')) - 500) < 1e-6, 'zoom is limited to x500') && ok;
ppi_callbacks('zoom', fig, 1e-6, [0 0]);
ok = check(all(abs(get(ppi.ax, 'XLim') - xl0) < 1e-6), 'zoom out stops at the full view') && ok;

%% Filters
ppi_callbacks('option', fig, 'showShort', false);
ok = check(isempty(get(S.h.short, 'XData')) && isempty(ppi_callbacks('hover', fig, xy(1))), ...
           'hiding the short-pulse plots removes them and their data tips') && ok;
ppi_callbacks('option', fig, 'showShort', true);
ppi_callbacks('option', fig, 'minSnr_dB', 20);
ok = check(numel(get(S.h.long, 'XData')) == 1, 'SNR filter hides the 14 dB plot') && ok;
ppi_callbacks('option', fig, 'minSnr_dB', 0);
ppi_callbacks('option', fig, 'history', 1);                % last scan: 10 s at 6 rpm
ok = check(numel(get(S.h.long, 'XData')) == 1, 'age filter hides the plot of the previous scan') && ok;
ppi_callbacks('option', fig, 'history', 0);
ok = check(numel(get(S.h.long, 'XData')) == 2, 'all plots shown again') && ok;
ppi_callbacks('option', fig, 'showLabels', true);
S = getappdata(fig, 'ppi');
ok = check(numel(S.h.labels) == 3, 'one label per plot') && ok;
ppi_callbacks('option', fig, 'showLabels', false);
ppi_callbacks('option', fig, 'showRings', false);
ok = check(all(strcmp(get(S.h.rings, 'Visible'), 'off')), 'range rings hidden') && ok;
ppi_callbacks('option', fig, 'showRings', true);

%% Pinned data tips
n = ppi_callbacks('pin', fig, xy(2));
ok = check(n == 1, 'click on a plot pins its data tip') && ok;
n = ppi_callbacks('pin', fig, xy(3));
ok = check(n == 2, 'second pinned data tip') && ok;
n = ppi_callbacks('pin', fig, xy(2));
ok = check(n == 1, 'click on a pinned plot removes its data tip') && ok;
ppi_callbacks('cleartips', fig);
S = getappdata(fig, 'ppi');
ok = check(isempty(S.h.pins), 'remove all data tips') && ok;

%% CFAR detections before the plot extractor (display modes)
% two detections of one target in successive CPI updates (long pulse) and a short-pulse one
dets = [11.0  95.0 30e3 5 3   28 2 5 150 30;
        11.1  95.2 30e3 5 3   27 2 5 150 30;
        11.5 300.0 15e3 2 -4  13 1 2  80 12];
ppi_update(ppi, sys, 100, 12.0, 2, p, dets);
S = getappdata(fig, 'ppi');
ok = check(strcmp(get(S.h.detLong, 'Visible'), 'off') && numel(get(S.h.long, 'XData')) == 2, ...
           'default mode: plots only, no detections') && ok;
ppi_callbacks('option', fig, 'showMode', 'detections');
xd = get(S.h.detLong, 'XData');
ok = check(isempty(get(S.h.long, 'XData')) && isempty(get(S.h.short, 'XData')) && numel(xd) == 6 ...
           && numel(get(S.h.detShort, 'XData')) == 3, 'detections mode: arcs only (3 points per arc)') && ok;
yd = get(S.h.detLong, 'YData');
arcLen = hypot(xd(2) - xd(1), yd(2) - yd(1));
ok = check(abs(arcLen - 30*cosd(5)*S.detArc_deg*pi/180) < 1e-6, ...
           'an arc spans the azimuth step of a CPI update') && ok;
dxy = 15*cosd(2)*[sind(300) cosd(300)];
txt = ppi_callbacks('hover', fig, dxy);
ok = check(~isempty(strfind(txt, 'CFAR detection (short pulse)')) && ~isempty(strfind(txt, '15.000 km')) ...
           && ~isempty(strfind(txt, 'Doppler bin 12')), 'hover over a detection shows its values') && ok;
ok = check(isempty(strfind(ppi_callbacks('hover', fig, xy(1)), 'Plot')), 'plots have no data tip in detections mode') && ok;
ppi_callbacks('option', fig, 'minSnr_dB', 20);
ok = check(all(isnan(get(S.h.detShort, 'XData'))) || isempty(get(S.h.detShort, 'XData')), ...
           'SNR filter also hides detections') && ok;
ppi_callbacks('option', fig, 'minSnr_dB', 0);
ppi_callbacks('key', fig, struct('Key', 'd', 'Character', 'd'));
S = getappdata(fig, 'ppi');
ok = check(strcmp(S.opt.showMode, 'both') && numel(get(S.h.long, 'XData')) == 2 && numel(get(S.h.detLong, 'XData')) == 6 ...
           && get(S.h.ui.showMode, 'Value') == 3, 'key d: plots and detections') && ok;
ok = check(~isempty(strfind(get(S.h.title, 'String'), 'detections 3')), 'title counts the detections') && ok;
n = ppi_callbacks('pin', fig, dxy);
ok = check(n == 1, 'a detection data tip can be pinned') && ok;
ppi_callbacks('cleartips', fig);
ppi_callbacks('option', fig, 'showMode', 'plots');

%% CFAR threshold
alpha0 = sys.cfar.alpha;
dB0 = 10*log10(alpha0);
ok = check(~isempty(strfind(get(S.h.ui.cfarInfo, 'String'), sprintf('2-D: %.2f dB', dB0))) ...
           && ~isempty(strfind(get(S.h.ui.cfarInfo, 'String'), 'Pfa 1.0e-08')), ...
           'the panel shows the threshold of cfg.cfar.pfa (1e-8)') && ok;
ppi_callbacks('cfar', fig, 18);
sys2 = ppi_update(ppi, sys, 101, 12.1, 2, p, dets);
ok = check(abs(sys2.cfar.alpha - 10^1.8) < 1e-9 && ~isempty(strfind(get(S.h.ui.cfarInfo, 'String'), '2-D: 18.00 dB')), ...
           'Apply: ppi_update hands the new threshold to the processing') && ok;
sys3 = ppi_update(ppi, sys2, 102, 12.2, 2, p, dets);
ok = check(sys3.cfar.alpha == sys2.cfar.alpha, 'the request is applied once') && ok;
ppi_callbacks('cfar', fig, NaN);
sys4 = ppi_update(ppi, sys3, 103, 12.3, 2, p, dets);
ok = check(sys4.cfar.alpha == sys3.cfar.alpha && ~isempty(strfind(get(S.h.ui.cfarInfo, 'String'), 'enter 0 ... 50 dB')), ...
           'an invalid value is refused') && ok;
ppi_callbacks('cfar', fig, dB0);                           % Default button
sys5 = ppi_update(ppi, sys4, 104, 12.4, 2, p, dets);
ok = check(abs(sys5.cfar.alpha - alpha0) < 1e-9*alpha0, 'Default restores the threshold of cfg.cfar.pfa') && ok;
ok = check(abs(cfar_pfa(alpha0, sys.cfar, sys.dop.nLooks) - cfg.cfar.twoD.pfa) < 1e-3*cfg.cfar.twoD.pfa, ...
           'cfar_pfa is the inverse of cfar_alpha') && ok;

%% Tracks and block processing times
tk = struct('id', {7, 8}, 'confirmed', {true, false}, 't_s', 11.9, 'pos_m', {[30e3*sind(95) 30e3*cosd(95) 2600], ...
            [5e3 5e3 0]}, 'vel_mps', {[-100 50 0], [0 0 0]}, 'alt_m', {2606, NaN}, 'hits', {5, 1}, ...
            'misses', 0, 'dim', 3);
us = [120 340 56 980 760 7 15 4 30 500];             % blocks, console output, display
ppi_update(ppi, sys, 105, 12.5, 2, p, dets, tk, us);
S = getappdata(fig, 'ppi');
ok = check(numel(get(S.h.track, 'XData')) == 1 && numel(get(S.h.trackTent, 'XData')) == 1, ...
           'confirmed and tentative tracks drawn as triangles') && ok;
c3 = xy(3);
ok = check(abs(get(S.h.track, 'XData') - c3(1)) < 0.3, 'the triangle of track 7 surrounds the plot at 30 km') && ok;
txt = ppi_callbacks('hover', fig, xy(3));
ok = check(~isempty(strfind(txt, 'Track 7 (confirmed, 3-D EKF)')) && ~isempty(strfind(txt, 'Speed       111.8 m/s')) ...
           && ~isempty(strfind(txt, 'Height      2606 m')), 'hover over a track (wins over its plot)') && ok;
str = get(S.h.timing, 'String');
ok = check(numel(str) == 13 && ~isempty(strfind(str{5}, '980.0')) && ~isempty(strfind(str{6}, 'CFAR CA 2D')) ...
           && ~isempty(strfind(str{7}, 'Peak interp.')) && ~isempty(strfind(str{9}, 'Track 3D CV')) ...
           && ~isempty(strfind(str{10}, 'Console output')) && ~isempty(strfind(str{11}, 'Display (PPI)')) ...
           && ~isempty(strfind(str{11}, '500.0')) && ~isempty(strfind(str{12}, '2812.0')) ...
           && ~isempty(strfind(str{13}, 'whole run')) && max(cellfun(@numel, str(2:12))) <= 24, ...
           'block, console and display times in the panel (Total includes the display)') && ok;
ppi_update(ppi, sys, 105, 12.5, 2, p, dets, tk, us, 50);           % means over the last 50 PRIs
str = get(S.h.timing, 'String');
ok = check(~isempty(strfind(str{end}, 'over the last 50 PRIs')), 'the panel says over how many PRIs') && ok;
ppi_callbacks('option', fig, 'showTracks', false);
ok = check(all(isnan(get(S.h.track, 'XData'))) && isempty(strfind(ppi_callbacks('hover', fig, xy(3)), 'Track')), ...
           'Tracks checkbox hides the tracks') && ok;
ppi_callbacks('option', fig, 'showTracks', true);

%% Pause, clear, save
ppi_callbacks('pause', fig);
ppi_update(ppi, sys, 120, 13.0, 2, p(1));
ok = check(numel(get(S.h.long, 'XData')) == 2, 'paused display is not redrawn') && ok;
ppi_callbacks('pause', fig);
ok = check(isempty(get(S.h.long, 'XData')), 'resumed display shows the new data') && ok;
ppi_update(ppi, sys, 120, 13.0, 2, p);
ppi_callbacks('clear', fig);
ok = check(isempty(get(S.h.short, 'XData')) && isempty(get(S.h.long, 'XData')), 'clear removes the plots') && ok;

%% Video: DDC output averaged over the elements, both receive intervals joined
ppi_callbacks('cleartips', fig);
V = ppi.video;
rng(3);
tgtR = 30e3;  dR = V.dR_m;
kTgt = find(sys.wave.long.rxIdx >= 1, 1) + round((tgtR - sys.wave.long.range_m(1)) / dR);   % rx sample of 30 km
for az = 30:0.036:70
    I = round(20*randn(sys.nDdc, sys.nCh));  Q = round(20*randn(sys.nDdc, sys.nCh));
    if az > 45 && az < 55                                 % echo 20 dB above the noise at 30 km
        k = sys.wave.long.rxIdx(kTgt) + (0:3);
        I(k, :) = I(k, :) + 200;
    end
    V = ppi_video_add(V, I, Q, az, 14 + az/36);
end
ppi_callbacks('option', fig, 'showMode', 'video');
[~, want] = ppi_update(ppi, sys, 70, 16.0, 2, p, dets, tk, [], [], V);
S = getappdata(fig, 'ppi');
ok = check(strcmp(want.raw, 'video') && ~want.rd && strcmp(get(S.h.video, 'Visible'), 'on') && isempty(get(S.h.long, 'XData')) ...
           && get(S.h.ui.showMode, 'Value') == 4, 'video mode: the video image instead of the plots') && ok;
A = get(S.h.video, 'AlphaData');  C = get(S.h.video, 'CData');
pix = @(r, az) round(([r*cosd(az) r*sind(az)] / (V.rMax_km*1e3) + 1) / 2 * (V.nPix - 1)) + 1;   % [row col]
pt = pix(tgtR + 2*dR, 50);  pn = pix(20e3, 40);  pe = pix(20e3, 150);
bright = @(q) max(max(sum(C(q(1)+(-1:1), q(2)+(-1:1), :), 3)));   % (a pixel is 267 m)
ok = check(A(pt(1), pt(2)) == 1 && A(pn(1), pn(2)) == 1 && A(pe(1), pe(2)) == 0, ...
           'video drawn where the antenna passed (30 - 70 deg) only') && ok;
ok = check(bright(pt) > bright(pn) + 0.5 && V.noise_dB < 10*log10(2*400) && V.noise_dB > 10*log10(2*400) - 1.5, ...
           'the echo is brighter than the noise; noise floor of the video') && ok;
ppi_callbacks('hover', fig, [0 0]);                        % (no data tip in the screenshot)
ppi_callbacks('save', fig, fullfile(tempdir, 'ppi_video_test.png'));
keyD = struct('Key', 'd', 'Character', 'd', 'Modifier', {{}});
ppi_callbacks('key', fig, keyD);
S = getappdata(fig, 'ppi');
ok = check(strcmp(S.opt.showMode, 'dbf') && get(S.h.ui.showMode, 'Value') == 5, 'key d: from the video to DBF out') && ok;

%% DBF out / MF out of the selected beam: the profile rows of ppi_profile_row.m
ok = check(numel(get(S.h.ui.beam, 'String')) == sys.nBeams && get(S.h.ui.beam, 'Value') == 1, ...
           'beam selector: every beam after the DBF, beam 1 at the start') && ok;
ppi_callbacks('option', fig, 'beam', 3);
[~, want] = ppi_update(ppi, sys, 70, 16.1, 2, p, dets, tk, [], [], V);   % (V still holds the video)
ok = check(strcmp(want.raw, 'dbf') && want.beam == 3, 'DBF out of beam 3 wanted') && ok;
A = get(S.h.video, 'AlphaData');
ok = check(~any(A(:)), 'the video is not shown as DBF out') && ok;
V = ppi_video_reset(ppi.video, 'dbf', 3);
ws = sys.wave.short;  wl = sys.wave.long;
for az = 100:0.036:120
    x = {complex(randn(1, numel(ws.ddcIdx)), randn(1, numel(ws.ddcIdx))), ...
         complex(randn(1, numel(wl.ddcIdx)), randn(1, numel(wl.ddcIdx)))};
    j = find(wl.ddcIdx == wl.rxIdx(1)) + round((tgtR - wl.range_m(1)) / dR);   % 30 km in the long pulse
    x{2}(j + (0:2)) = 30;
    V = ppi_video_add_row(V, ppi_profile_row(V, 'dbf', x), az, 17 + az/36);
end
ppi_update(ppi, sys, 120, 20.4, 2, p, dets, tk, [], [], V);
A = get(S.h.video, 'AlphaData');  C = get(S.h.video, 'CData');
pt = pix(tgtR + dR, 110);  pn = pix(15e3, 110);  pe = pix(20e3, 60);
bright = @(q) max(max(sum(C(q(1)+(-1:1), q(2)+(-1:1), :), 3)));   % (of this image)
ok = check(A(pt(1), pt(2)) == 1 && A(pe(1), pe(2)) == 0 && bright(pt) > bright(pn) + 0.5 && abs(V.noise_dB - 10*log10(2)) < 1, ...
           'DBF out: drawn at 100 - 120 deg only, the echo at 30 km bright, noise floor') && ok;
ok = check(~isempty(strfind(get(S.h.videoBarLabel, 'String'), 'DBF out, beam 3')), 'colour scale names DBF out and the beam') && ok;
ppi_callbacks('save', fig, fullfile(tempdir, 'ppi_dbf_test.png'));
ppi_callbacks('key', fig, keyD);
[~, want] = ppi_update(ppi, sys, 120, 20.5, 2, p, dets, tk, [], [], V);
A = get(S.h.video, 'AlphaData');
ok = check(strcmp(want.raw, 'mf') && ~any(A(:)), 'key d: MF out wanted, the DBF data is not shown as MF out') && ok;
rowMf = ppi_profile_row(V, 'mf', {zeros(1, ws.nRange), [zeros(1, 100) 5 zeros(1, wl.nRange - 101)]});
b = find(rowMf > 0);
ok = check(numel(b) == 1 && abs((b - 0.5) * dR - wl.range_m(101)) <= dR, 'MF out: a range cell in its range bin') && ok;
ppi_callbacks('key', fig, keyD);
S = getappdata(fig, 'ppi');
ok = check(strcmp(S.opt.showMode, 'plots') && strcmp(get(S.h.video, 'Visible'), 'off'), ...
           'key d: from MF out back to the plots') && ok;

%% Range-Doppler tab: the map of the selected beam next to a small PPI
pos1 = get(ppi.ax, 'Position');
ppi_callbacks('tab', fig, 2);
S = getappdata(fig, 'ppi');
[~, want] = ppi_update(ppi, sys, 120, 20.6, 2, p, dets, tk);
pos2 = get(ppi.ax, 'Position');
ok = check(want.rd && want.beam == 3 && strcmp(get(S.h.rdAx, 'Visible'), 'on') && pos2(3) < pos1(3) ...
           && get(S.h.tabBtn(2), 'Value') == 1 && ~isempty(strfind(get(S.h.rdTitle, 'String'), 'waiting')), ...
           'tab 2: the map wanted, the PPI small, the map waits for a CPI update') && ok;
nfft = sys.dop.nfft;
slices = {rand(ws.nRange, nfft) + 0.5, rand(wl.nRange, nfft) + 0.5};
slices{2}(101, 20 + (-2:2)) = 1000 * exp(-0.5 * (-2:2).^2);
rdm = ppi_rd_map(V, sys, slices, [wl.range_m(101), 20, 1000], 3, 20.6, 120);
ok = check(abs(rdm.noise - 1) < 0.05 && sum(~isnan(rdm.P(:, 1))) >= ws.nRange / 2 + wl.nRange - 2, ...
           'the map joins both sub-pulses in range bins; noise floor = median') && ok;
ppi_update(ppi, sys, 120, 20.7, 2, p, dets, tk, [], [], V, rdm);
img = get(S.h.rdImg, 'CData');
ok = check(strcmp(get(S.h.rdImg, 'Visible'), 'on') && max(img(:)) > 0.9 && ~isempty(strfind(get(S.h.rdTitle, 'String'), 'peak 30.0 dB')) ...
           && abs(get(S.h.rdDet, 'XData') - wl.range_m(101)/1e3) < 1e-9, 'map drawn: the 30 dB peak, its detection circled') && ok;
sp = get(S.h.specLine, 'YData');
ok = check(abs(max(sp) - 30) < 0.01, 'the spectrum at the strongest detection shows its peak') && ok;
txt = ppi_callbacks('rdhover', fig, [wl.range_m(101)/1e3, rdm.vel_mps(20)]);
ok = check(~isempty(strfind(txt, '30.0 dB above noise')), 'hover over the map: the cell in dB above noise') && ok;
S = getappdata(fig, 'ppi');
ok = check(~isempty(strfind(get(S.h.rdTip, 'String'), 'CFAR detection')) && strcmp(get(S.h.rdTip, 'Visible'), 'on'), ...
           'map: data tip of the detection under the cursor') && ok;
ok = check(strcmp(get(get(S.h.specAx, 'XLabel'), 'String'), 'Radial velocity [m/s]') && ...
           ~isempty(strfind(get(get(S.h.specAx, 'YLabel'), 'String'), 'dB above noise')), 'spectrum: axis labels') && ok;
% zoom, pan, pinned data tips, reset of the map
full = [get(S.h.rdAx, 'XLim') get(S.h.rdAx, 'YLim')];
c = [wl.range_m(101)/1e3, rdm.vel_mps(20)];
pkBin = floor(wl.range_m(101) / rdm.dR_m) + 1;
ppi_callbacks('rdzoom', fig, 'map', 4, c);
v = [get(S.h.rdAx, 'XLim') get(S.h.rdAx, 'YLim')];
ok = check(abs(diff(v(1:2)) - diff(full(1:2))/4) < 1e-6 && abs(diff(v(3:4)) - diff(full(3:4))/4) < 1e-6 && ...
           v(1) < c(1) && v(2) > c(1) && v(3) < c(2) && v(4) > c(2), 'map: zoom x4 around the peak (both axes)') && ok;
ppi_callbacks('rdzoom', fig, 'map', 2, c, [2 1]);
v1 = [get(S.h.rdAx, 'XLim') get(S.h.rdAx, 'YLim')];
ok = check(abs(diff(v1(1:2)) - diff(v(1:2))/2) < 1e-6 && abs(diff(v1(3:4)) - diff(v(3:4))) < 1e-6, ...
           'map: zoom of the range axis only (Shift + wheel)') && ok;
ppi_callbacks('rdpan', fig, 'map', [1 0]);
v2 = [get(S.h.rdAx, 'XLim') get(S.h.rdAx, 'YLim')];
ok = check(abs(v2(1) - v1(1) - 1) < 1e-6, 'map: pan by 1 km') && ok;
n = ppi_callbacks('rdpin', fig, 'map', c);
S = getappdata(fig, 'ppi');
ok = check(n == 1 && isequal(S.rdSel, [pkBin 20]) && numel(S.h.rdPinH) >= 2, 'map: click pins the data tip, selects the cell') && ok;
ppi_update(ppi, sys, 120, 20.75, 2, p, dets, tk, [], [], V, rdm);
S = getappdata(fig, 'ppi');
ok = check(max(abs([get(S.h.rdAx, 'XLim') get(S.h.rdAx, 'YLim')] - v2)) < 1e-9 && size(S.rdPins, 1) == 1 && numel(S.h.rdPinH) >= 2, ...
           'a new map keeps the zoom and the pinned tip') && ok;
ppi_callbacks('save', fig, fullfile(tempdir, 'ppi_rd_zoom_test.png'));
ppi_callbacks('rdreset', fig, 'map');
ok = check(max(abs([get(S.h.rdAx, 'XLim') get(S.h.rdAx, 'YLim')] - full)) < 1e-9, 'map: reset (double click)') && ok;
% the plot below the map: data tip, zoom, the range profile
txt = ppi_callbacks('lowhover', fig, [rdm.vel_mps(20), 25]);
ok = check(~isempty(strfind(txt, 'CFAR detection')) && strcmp(get(S.h.specTip, 'Visible'), 'on'), ...
           'spectrum: data tip of the nearest point') && ok;
lf = [get(S.h.specAx, 'XLim') get(S.h.specAx, 'YLim')];
ok = check(lf(4) > 30, 'spectrum: room above the peak') && ok;
ppi_callbacks('rdzoom', fig, 'lower', 2, [rdm.vel_mps(20), 20]);
lz = [get(S.h.specAx, 'XLim') get(S.h.specAx, 'YLim')];
ok = check(abs(diff(lz(1:2)) - diff(lf(1:2))/2) < 1e-6, 'spectrum: zoom around the cursor') && ok;
ppi_callbacks('option', fig, 'rdLower', 'range');
S = getappdata(fig, 'ppi');
xs = get(S.h.specLine, 'XData');  ys = get(S.h.specLine, 'YData');
[pk, ip] = max(ys);
ok = check(get(S.h.ui.rdLower, 'Value') == 2 && strcmp(get(get(S.h.specAx, 'XLabel'), 'String'), 'Range [km]') && ...
           abs(pk - 30) < 0.01 && abs(xs(ip) - (pkBin - 0.5) * rdm.dR_m / 1e3) < 1e-9 && ...
           ~isempty(strfind(get(S.h.specTitle, 'String'), 'Range profile')), ...
           'Below map box: the range profile through the Doppler bin of the peak') && ok;
txt = ppi_callbacks('lowhover', fig, [(pkBin - 0.5) * rdm.dR_m / 1e3, 30]);
ok = check(~isempty(strfind(txt, 'Doppler bin 20')), 'range profile: data tip of the peak') && ok;
n = ppi_callbacks('rdpin', fig, 'lower', [(pkBin - 0.5) * rdm.dR_m / 1e3, 30]);
ok = check(n == 0, 'a click on the pinned point removes the tip') && ok;
ppi_callbacks('save', fig, fullfile(tempdir, 'ppi_rd_range_test.png'));
ppi_callbacks('option', fig, 'rdLower', 'spectrum');
% pause: the map stays (also the data tips use it) until resumed
ppi_callbacks('pause', fig);
rdm2 = rdm;  rdm2.t_s = 21.5;
ppi_update(ppi, sys, 120, 21.5, 2, p, dets, tk, [], [], V, rdm2);
ppi_callbacks('rdhover', fig, c);
S = getappdata(fig, 'ppi');
ok = check(S.data.rd.t_s == rdm.t_s && ~isempty(strfind(get(S.h.rdTitle, 'String'), 'PAUSED')) && ...
           strcmp(get(S.h.ui.pause, 'String'), 'Resume display'), 'paused: the map is not replaced') && ok;
ppi_callbacks('pause', fig);
S = getappdata(fig, 'ppi');
ok = check(S.data.rd.t_s == 21.5 && isempty(S.pending) && isempty(strfind(get(S.h.rdTitle, 'String'), 'PAUSED')), ...
           'resume: the newest map') && ok;
ppi_callbacks('save', fig, fullfile(tempdir, 'ppi_rd_test.png'));
ppi_callbacks('option', fig, 'beam', 2);
S = getappdata(fig, 'ppi');
ok = check(isempty(S.data.rd) && ~isempty(strfind(get(S.h.rdTitle, 'String'), 'beam 2: waiting')), ...
           'another beam: the old map is dropped') && ok;
ppi_callbacks('tab', fig, 1);
[~, want] = ppi_update(ppi, sys, 120, 20.8, 2, p, dets, tk);
ok = check(~want.rd && max(abs(get(ppi.ax, 'Position') - pos1)) < 1e-6 && strcmp(get(S.h.rdAx, 'Visible'), 'off') ...
           && strcmp(get(S.h.rdTitle, 'Visible'), 'off') && strcmp(get(S.h.specTitle, 'Visible'), 'off'), ...
           'back to tab 1: the PPI full size, no map wanted') && ok;

file = [tempname() '.png'];
ppi_callbacks('save', fig, file);
ok = check(exist(file, 'file') == 2, 'PPI saved as PNG') && ok;
if exist(file, 'file'), delete(file); end
close(fig);
fprintf('Interactive PPI: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function v = plot_xy(p)
% Ground position of a plot on the PPI [km].
el = p.el_deg;
if isnan(el), el = 0; end
v = p.range_m * cosd(el) * [sind(p.az_deg) cosd(p.az_deg)] / 1e3;
end

function ok = report(cond, msg)
ok = logical(cond);
fprintf('  %-62s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

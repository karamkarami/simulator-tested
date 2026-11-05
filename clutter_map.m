function [img, info] = clutter_map(cfg, fileName)
%CLUTTER_MAP  Picture of the clutter terrain around the radar (PNG, no graphics toolkit needed).
%   [img, info] = CLUTTER_MAP(cfg, fileName) draws the terrain map of
%   cfg.clutter (clutter_terrain.m) seen from above, North up, out to the
%   farthest clutter cell: every visible cell in the colour of its terrain
%   type (land green, hill olive, mountain brown, desert yellow, sea blue),
%   brighter for a higher CNR (per range cell, from the CNR of land at
%   cfg.clutter.cnrRange_m and the R^3 law), shadowed and beyond-horizon
%   cells dark grey, range rings every 10 km. It is written to fileName
%   (default output/clutter_map.png). info holds the CNR range per type.
%   Example: cfg = radar_config(); clutter_map(cfg);
if nargin < 1 || isempty(cfg), cfg = radar_config(); end
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'));
if nargin < 2, fileName = fullfile(cfg.io.outFolder, 'clutter_map.png'); end
sys = setup_radar(cfg);
T = clutter_terrain(cfg, sys);
colors = struct('land', [0.25 0.75 0.25], 'hill', [0.6 0.65 0.15], 'mountain', [0.65 0.4 0.2], ...
                'desert', [0.95 0.85 0.35], 'sea', [0.2 0.45 0.95]);
% CNR of every cell's range cell [dB]
sig0 = 10^(sys.clutter.sigma0_dB/10);
cnr = sys.clutter.cnr_dB + 10*log10(max(T.sigma0, realmin) / sig0) + 30*log10(sys.clutter.cnrRange_m ./ T.range_m);
info = struct('type', {}, 'cnrMin_dB', {}, 'cnrMax_dB', {}, 'visibleCells', {});
for i = 1:numel(T.typeNames)
    m = T.type == i & T.visible;
    if any(m(:))
        info(end+1) = struct('type', T.typeNames{i}, 'cnrMin_dB', min(cnr(m)), 'cnrMax_dB', max(cnr(m)), ...
                             'visibleCells', nnz(m)); %#ok<AGROW>
    end
end
lo = 20;  hi = max(cnr(T.visible));                 % brightness scale [dB]
N = 801;
Rmax = T.range_m(end);
[X, Y] = meshgrid(linspace(-Rmax, Rmax, N), linspace(Rmax, -Rmax, N));
R = hypot(X, Y);
az = mod(atan2(X, Y) * 180/pi, 360);
ia = mod(round(az / cfg.clutter.azStep_deg), numel(T.az_deg)) + 1;
ir = round((R - T.range_m(1)) / T.dR_m) + 1;
inGrid = ir >= 1 & ir <= numel(T.range_m);
img = repmat(reshape([0.05 0.05 0.08], 1, 1, 3), N, N);
lin = sub2ind(size(T.type), ia(inGrid), ir(inGrid));
vis = T.visible(lin);
typ = T.type(lin);
b = min(max((cnr(lin) - lo) / (hi - lo), 0), 1);
pix = find(inGrid);
for c = 1:3
    ch = img(:, :, c);
    col = cellfun(@(n) colors.(n)(c), T.typeNames);
    v = 0.18 * ones(numel(pix), 1);                   % shadowed / beyond the horizon: grey
    v(vis) = reshape(col(typ(vis)), [], 1) .* (0.25 + 0.75 * b(vis));
    ch(pix) = v;
    ring = abs(mod(R + 500, 10e3) - 500) < Rmax / N & R <= Rmax;  % range rings every 10 km
    ch(ring) = 0.85;
    ch(R < Rmax / 100) = 1;                           % the radar
    img(:, :, c) = ch;
end
folder = fileparts(fileName);
if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
imwrite(img, fileName);
fprintf('Clutter map (%s terrain, CNR of land %.1f dB at %.1f km) -> %s\n', cfg.clutter.terrain, ...
        sys.clutter.cnr_dB, sys.clutter.cnrRange_m/1e3, fileName);
for i = 1:numel(info)
    fprintf('  %-8s %7d visible cells, CNR %5.1f ... %5.1f dB\n', info(i).type, info(i).visibleCells, ...
            info(i).cnrMin_dB, info(i).cnrMax_dB);
end
end

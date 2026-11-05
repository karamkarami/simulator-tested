function T = clutter_terrain(cfg, sys)
%CLUTTER_TERRAIN  Terrain map of the ground clutter (deterministic, no random numbers).
%   T = CLUTTER_TERRAIN(cfg, sys) returns the clutter grid around the radar
%   (azimuth x range cells of cfg.clutter.azStep_deg x
%   cfg.clutter.rangeStepSamples ADC samples) and for every cell its terrain
%   type, reflectivity, internal motion, elevation, visibility and
%   elevation group:
%     T.range_m [1 x nR], T.az_deg [nAz x 1], T.dR_m (range size of a cell)
%     T.typeNames  names of the terrain types (index 1 = 'land')
%     T.type       [nAz x nR] type index of every cell
%     T.visible    [nAz x nR] the radar sees the cell (not behind the horizon
%                  or behind higher terrain on the same azimuth)
%     T.sigma0     [nAz x nR] reflectivity (linear, 0 where not visible)
%     T.el_deg     [nAz x nR] elevation of the cell seen from the antenna
%     T.icmStd     [nAz x nR] std of the internal motion [m/s] (NaN: the
%                  spectrum of cfg.clutter: 'gaussian' icmStd_mps or 'exponential')
%     T.meanVel    [nAz x nR] mean radial velocity [m/s] (sea: wave drift)
%     T.group      [nAz x nR] elevation group (one steering vector each)
%     T.groupEl_deg [G x 1] elevation of every group
%   cfg.clutter.terrain:
%     'uniform'  land (reflectivity sigma0 of sys.clutter) everywhere from
%                cfg.clutter.minRange_m to the radar horizon of the flat
%                ground (4/3 earth) or cfg.clutter.maxRange_m, from the
%                elevation cfg.clutter.elevation_deg: the model before the
%                terrain types
%     'regions'  land everywhere, then every row {type, [az1 az2], [r1 r2],
%                height_m} of cfg.clutter.regions on top (later rows win;
%                the sector runs clockwise from az1 to az2 and may cross
%                North). Terrain shapes (cfg.clutter.types.<type>.shape):
%                'flat'    at height_m (0 for sea, desert, plain land)
%                'ridge'   a mountain range: height_m * sin(pi*u) in range
%                          (u = 0 ... 1 from r1 to r2), tapered in azimuth
%                          with sqrt(sin(pi*v)) (v = 0 ... 1 from az1 to az2),
%                          with a ripple of +-15%
%                'rolling' hills on ground rising from 0 to height_m from r1
%                          to r2: crests 1.5 km apart (sin^2), shifted and
%                          20% higher or lower with azimuth; the front slopes
%                          are seen, the valleys behind the crests shadowed
%                A cell is visible when its elevation angle (4/3 earth
%                curvature included) is not below that of any nearer cell
%                of its azimuth: the back slopes of ridges and hills are
%                shadowed, flat ground and sea end at the radar horizon,
%                elevated terrain is seen beyond it. Flat terrain at height 0
%                returns from cfg.clutter.elevation_deg, other terrain from
%                its own elevation. The reflectivity of a type is sigma0
%                (land) + sigma0Rel_dB.
%   The cells of height 0 form group 1 (cfg.clutter.elevation_deg); every
%   region of other height its own group, at the elevation of its visible
%   cells weighted by their power (sigma0 / R^3).

cc = cfg.clutter;
c = cfg.const.c;
hAnt = cfg.ant.height_m;
kae = 4/3 * 6371e3;                                       % 4/3 earth radius [m]
T.dR_m = cc.rangeStepSamples * c / (2*sys.fs);
sigma0Land = 10^(sys.clutter.sigma0_dB/10);
uniform = strcmpi(terrain_mode(cc), 'uniform');

regions = {};
if ~uniform, regions = check_regions(cc); end
rMax = cc.maxRange_m;
if isempty(rMax)
    rMax = sqrt(2 * kae * hAnt);                          % radar horizon of the flat ground
    for k = 1:size(regions, 1)
        if regions{k, 4} > 0, rMax = max(rMax, regions{k, 3}(2)); end   % elevated terrain beyond it
    end
end
T.range_m = (ceil(cc.minRange_m/T.dR_m):floor(rMax/T.dR_m)) * T.dR_m;
T.az_deg = (0:cc.azStep_deg:360-cc.azStep_deg)';
nAz = numel(T.az_deg);
nR = numel(T.range_m);

T.typeNames = {'land'};
if ~uniform
    T.typeNames = unique([{'land'}, regions(:, 1).'], 'stable');
end
nTypes = numel(T.typeNames);
types = cell(1, nTypes);
for i = 1:nTypes
    types{i} = terrain_type(cc, T.typeNames{i});
end

T.type = ones(nAz, nR, 'uint8');
height = zeros(nAz, nR);
regionOf = zeros(nAz, nR, 'uint8');                       % region of every cell (0 = plain land)
RR = repmat(T.range_m, nAz, 1);
for k = 1:size(regions, 1)
    [inAz, v] = sector(T.az_deg, regions{k, 2});
    rr = regions{k, 3};
    inR = T.range_m >= rr(1) & T.range_m <= rr(2);
    if ~any(inAz) || ~any(inR), continue; end
    ti = find(strcmp(T.typeNames, regions{k, 1}));
    H = regions{k, 4};
    u = (T.range_m(inR) - rr(1)) / max(rr(2) - rr(1), eps);
    vv = v(inAz);
    switch types{ti}.shape
        case 'flat'
            h = H * ones(numel(vv), numel(u));
        case 'ridge'
            h = H * sqrt(sin(pi*vv)) * sin(pi*u) .* (1 + 0.15 * cos(2*pi*4*vv) * sin(2*pi*3*u));
        case 'rolling'
            Rk = T.range_m(inR) - rr(1);
            crest = sin(pi * (ones(numel(vv), 1) * Rk / 1500 + 0.5 * cos(2*pi*2*vv) * ones(1, numel(u)))).^2;
            h = H * (0.8 + 0.2*cos(2*pi*3*vv)) * u .* (0.4 + 0.6*crest);
        otherwise
            error('cfg.clutter.types.%s.shape must be ''flat'', ''ridge'' or ''rolling''.', regions{k, 1});
    end
    T.type(inAz, inR) = ti;
    height(inAz, inR) = h;
    regionOf(inAz, inR) = k;
end

% elevation of every cell and the shadowing along every azimuth
tanEl = (height - hAnt - RR.^2 / (2*kae)) ./ RR;          % tangent of the elevation angle (curved earth)
if uniform
    T.visible = true(nAz, nR);                            % to the horizon or maxRange_m, as before
else
    T.visible = tanEl >= cummax(tanEl, 2) - 1e-12;
end
T.el_deg = atand(tanEl);
flat0 = height == 0;
T.el_deg(flat0) = cc.elevation_deg;

rel = cellfun(@(t) t.sigma0Rel_dB, types);
icm = cellfun(@(t) nan_if_empty(t.icmStd_mps), types);
mv  = cellfun(@(t) t.meanVel_mps, types);
md  = cellfun(@(t) t.meanDir_deg, types);
T.sigma0 = sigma0Land * 10.^(rel(T.type)/10) .* T.visible;
T.icmStd = icm(T.type);
AZ = repmat(T.az_deg, 1, nR);
T.meanVel = mv(T.type) .* cosd(AZ - md(T.type));         % radial component of the drift

% elevation groups: height 0, then every region of other height
T.group = ones(nAz, nR, 'uint8');
T.groupEl_deg = cc.elevation_deg;
for k = 1:size(regions, 1)
    m = regionOf == k & ~flat0 & T.visible;
    if ~any(m(:)), continue; end
    w = T.sigma0(m) ./ RR(m).^3;
    T.groupEl_deg(end+1, 1) = sum(w .* T.el_deg(m)) / sum(w);
    T.group(m) = numel(T.groupEl_deg);
end
end

function m = terrain_mode(cc)
m = 'uniform';
if isfield(cc, 'terrain') && ~isempty(cc.terrain), m = cc.terrain; end
if ~any(strcmpi(m, {'uniform', 'regions'}))
    error('cfg.clutter.terrain must be ''uniform'' or ''regions''.');
end
end

function regions = check_regions(cc)
regions = {};
if isfield(cc, 'regions'), regions = cc.regions; end
if isempty(regions), regions = cell(0, 4); end
if ~iscell(regions) || size(regions, 2) ~= 4
    error('cfg.clutter.regions must have rows {type, [az1 az2], [r1 r2], height_m}.');
end
for k = 1:size(regions, 1)
    if ~ischar(regions{k, 1}) || numel(regions{k, 2}) ~= 2 || numel(regions{k, 3}) ~= 2 || ...
            regions{k, 3}(1) >= regions{k, 3}(2) || ~isscalar(regions{k, 4}) || regions{k, 4} < 0
        error('cfg.clutter.regions row %d: {type, [az1 az2] deg, [r1 r2] m with r1 < r2, height_m >= 0}.', k);
    end
end
end

function t = terrain_type(cc, name)
% Parameters of a terrain type: cfg.clutter.types.<name> (land: defaults).
t = struct('sigma0Rel_dB', 0, 'icmStd_mps', [], 'meanVel_mps', 0, 'meanDir_deg', 0, 'shape', 'flat');
if isfield(cc, 'types') && isfield(cc.types, name)
    f = fieldnames(cc.types.(name));
    for i = 1:numel(f), t.(f{i}) = cc.types.(name).(f{i}); end
elseif ~strcmp(name, 'land')
    error('Terrain type ''%s'' of cfg.clutter.regions is not in cfg.clutter.types.', name);
end
end

function [in, v] = sector(az, s)
% Cells inside the clockwise sector s = [az1 az2] and their position v = 0 ... 1 in it.
width = mod(s(2) - s(1), 360);
if width == 0, width = 360; end
d = mod(az - s(1), 360);
in = d <= width;
v = d / width;
end

function x = nan_if_empty(x)
if isempty(x), x = NaN; end
end

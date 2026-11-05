function clut = clutter_init(sys)
%CLUTTER_INIT  Static ground clutter map around the radar.
%   clut = CLUTTER_INIT(sys) builds the clutter cells (azimuth x range) of
%   the terrain map (clutter_terrain.m: uniform land, or land with regions
%   of mountains, hills, desert and sea) with complex Gaussian (Rayleigh
%   amplitude) reflectivity of the sigma0 of their terrain type, a random
%   internal-motion Doppler per cell, and precomputes what clutter_signal.m
%   needs to synthesise the clutter echo of every PRI. The power of every
%   cell follows the radar equation with the elevation pattern at the
%   elevation of the cell; the reflectivity of land is set by
%   cfg.clutter.cnr_dB (sys.clutter.sigma0_dB, clutter_cnr.m). Only the
%   visible cells are kept, sorted by azimuth (clut.rowStart).
%   Internal clutter motion of a cell (type with icmStd_mps = []: the
%   spectrum of cfg.clutter.spectrum):
%     'gaussian'    every cell moves with a radial velocity drawn from
%                   N(0, icmStd_mps^2): a Gaussian clutter spectrum
%     'exponential' Billingsley's land clutter model: every cell has a
%                   stationary (DC) part with r/(1+r) of its power and a
%                   moving (AC) part with 1/(1+r), whose radial velocity is
%                   drawn from the exponential (Laplace) density
%                   beta/2*exp(-beta*|v|). beta and r follow from the wind
%                   speed unless cfg.clutter.expBeta_spm / dcToAc are given.
%   A type with its own icmStd_mps (mountain, desert, sea) moves with
%   N(meanVel, icmStd_mps^2), meanVel the radial component of its drift
%   (sea: the waves).

cfg = sys.cfg;
cc  = cfg.clutter;
c   = cfg.const.c;
fs  = sys.fs;

T = clutter_terrain(cfg, sys);
clut.range_m = T.range_m;
clut.delaySamples = round(2*clut.range_m/c*fs);           % two-way delay [ADC samples]
clut.az_deg = T.az_deg;
clut.nGroups = numel(T.groupEl_deg);
nAz = numel(clut.az_deg);
nR  = numel(clut.range_m);

% Mean power of a cell at the ADC (antenna boresight in azimuth) and its rms IF amplitude
tilt   = sys.ant.tilt_deg;                                % as steering_vector.m at the boresight azimuth
elRel  = asind(-sind(tilt)*cosd(T.el_deg) + cosd(tilt)*sind(T.el_deg));
gEl    = pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, elRel) ...
       .* pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, elRel);
area   = clut.range_m * deg2rad(cc.azStep_deg) * T.dR_m;  % cell area [m^2]
power  = sys.radarConst * T.sigma0 .* (area ./ clut.range_m.^4) .* gEl;  % [W]
ampRms = sqrt(2*power);                                   % rms IF amplitude at the receiver input [sqrt(W)]
rho = ampRms .* (randn(nAz, nR) + 1j*randn(nAz, nR)) / sqrt(2);
own = ~isnan(T.icmStd);                                   % types with their own Gaussian motion
cfgGauss = strcmpi(cc.spectrum, 'gaussian');
if ~cfgGauss && ~strcmpi(cc.spectrum, 'exponential')
    error('cfg.clutter.spectrum must be ''gaussian'' or ''exponential''.');
end
v = zeros(nAz, nR);
dc = [];
if any(own(:)) || (cfgGauss && ~all(own(:)))
    z = randn(nAz, nR);
    sd = T.icmStd;
    if cfgGauss, sd(~own) = cc.icmStd_mps; end
    g = own | cfgGauss;
    v(g) = T.meanVel(g) + sd(g) .* z(g);
end
if ~cfgGauss && ~all(own(:))                              % Billingsley for the other types
    [beta, r] = billingsley(cc, cfg.fc);
    uu = rand(nAz, nR) - 0.5;
    e = ~own;
    v(e) = T.meanVel(e) - sign(uu(e)) .* log(1 - 2*abs(uu(e))) / beta;   % Laplace: beta/2*exp(-beta*|v|)
    dc = zeros(nAz, nR);
    dc(e)  = rho(e) * sqrt(r/(1 + r));                    % stationary part
    rho(e) = rho(e) * sqrt(1/(1 + r));                    % moving part
    clut.beta_spm = beta;
    clut.dcToAc = r;
end

% keep the visible cells, sorted by azimuth row
keep = (T.sigma0 > 0).';                                  % [nR x nAz]: column-major order = by azimuth
idx = find(keep);
[ir, ia] = ind2sub([nR nAz], idx);
lin = sub2ind([nAz nR], ia, ir);
clut.cellRow   = int32(ia);
clut.cellRange = ir;
clut.cellGroup = double(T.group(lin));
clut.rho   = rho(lin);
clut.fd_Hz = 2*v(lin)/sys.lambda;
clut.dc    = [];
if ~isempty(dc), clut.dc = dc(lin); end
clut.rowStart = [0; cumsum(accumarray(ia, 1, [nAz 1]))] + 1;   % cells of row a: rowStart(a) ... rowStart(a+1)-1
clut.cellType = T.type(lin);                             % index into clut.typeNames
clut.typeNames = T.typeNames;
clut.nCells = numel(idx);

% IF transmit waveforms at the ADC rate and their FFTs (for the fast convolution)
names = {'short', 'long'};
for i = 1:2
    wv = sys.wave.(names{i});
    nTx = round(wv.width_s * fs);
    t   = (0:nTx-1)' / fs;
    tx  = waveform_envelope(wv, t) .* exp(1j*2*pi*sys.fIF*t);
    nFft = 2^nextpow2(sys.nAdc + nTx - 1);
    clut.(names{i}).nTx    = nTx;
    clut.(names{i}).nFft   = nFft;
    clut.(names{i}).txFft  = fft(tx, nFft);
    % ADC sample (0-based, relative to the window start) where an echo of zero delay starts
    clut.(names{i}).start0 = round((wv.tTx_s - sys.rxStart_s) * fs);
end
end

function [beta, r] = billingsley(cc, fc)
% Exponential spectrum parameters of Billingsley's model (wind speed in mph, f in GHz).
w = max(cc.windSpeed_mps / 0.44704, 0.1);
beta = cc.expBeta_spm;
if isempty(beta), beta = 1 / (0.1*log10(w) + 0.068); end
r = cc.dcToAc;
if isempty(r), r = 489.8 * w^-1.55 * (fc/1e9)^-1.21; end
end

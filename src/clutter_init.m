function clut = clutter_init(sys)
%CLUTTER_INIT  Static ground clutter map around the radar.
%   clut = CLUTTER_INIT(sys) builds a grid of clutter cells (azimuth x range)
%   with complex Gaussian (Rayleigh amplitude) reflectivity of constant sigma0,
%   a random internal-motion Doppler per cell, and precomputes what
%   clutter_signal.m needs to synthesise the clutter echo of every PRI.
%   Internal clutter motion (cfg.clutter.spectrum):
%     'gaussian'    every cell moves with a radial velocity drawn from
%                   N(0, icmStd_mps^2): a Gaussian clutter spectrum
%     'exponential' Billingsley's land clutter model: every cell has a
%                   stationary (DC) part with r/(1+r) of its power and a
%                   moving (AC) part with 1/(1+r), whose radial velocity is
%                   drawn from the exponential (Laplace) density
%                   beta/2*exp(-beta*|v|). beta and r follow from the wind
%                   speed unless cfg.clutter.expBeta_spm / dcToAc are given.

cfg = sys.cfg;
cc  = cfg.clutter;
c   = cfg.const.c;
fs  = sys.fs;

% Range grid: cell size is a whole number of ADC samples so delays fall on samples
dR = cc.rangeStepSamples * c / (2*fs);
rMax = cc.maxRange_m;
if isempty(rMax)
    rMax = sqrt(2 * 4/3 * 6371e3 * cfg.ant.height_m);       % radar horizon (4/3 earth)
end
clut.range_m = (ceil(cc.minRange_m/dR):floor(rMax/dR)) * dR;
clut.delaySamples = round(2*clut.range_m/c*fs);           % two-way delay [ADC samples]
clut.az_deg = (0:cc.azStep_deg:360-cc.azStep_deg)';
nAz = numel(clut.az_deg);
nR  = numel(clut.range_m);

% Mean power of a cell at the ADC (antenna boresight in azimuth) and its rms IF amplitude
[~, u] = steering_vector(sys, 0, cc.elevation_deg);
elRel  = asind(u);
gEl    = pattern_gain(sys.ant.elRel_deg, sys.ant.txElGain_dB, elRel) ...
       * pattern_gain(sys.ant.elRel_deg, sys.ant.elElemGain_dB, elRel);
area   = clut.range_m * deg2rad(cc.azStep_deg) * dR;      % cell area [m^2]
power  = sys.radarConst * 10^(cc.sigma0_dB/10) * area ./ clut.range_m.^4 * gEl;  % [W]
ampRms = sqrt(2*power);                                   % rms IF amplitude at the receiver input [sqrt(W)]
rho = ampRms .* (randn(nAz, nR) + 1j*randn(nAz, nR)) / sqrt(2);
switch lower(cc.spectrum)                                 % internal clutter motion
    case 'gaussian'
        v = cc.icmStd_mps * randn(nAz, nR);
        clut.dc  = [];                                    % no separate stationary part
        clut.rho = rho;
    case 'exponential'
        [beta, r] = billingsley(cc, cfg.fc);
        u = rand(nAz, nR) - 0.5;
        v = -sign(u) .* log(1 - 2*abs(u)) / beta;         % Laplace: density beta/2*exp(-beta*|v|)
        clut.dc  = rho * sqrt(r/(1 + r));                 % stationary part
        clut.rho = rho * sqrt(1/(1 + r));                 % moving part
        clut.beta_spm = beta;
        clut.dcToAc = r;
    otherwise
        error('cfg.clutter.spectrum must be ''gaussian'' or ''exponential''.');
end
clut.fd_Hz = 2*v/sys.lambda;

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

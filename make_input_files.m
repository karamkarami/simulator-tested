function make_input_files(cfg, which)
%MAKE_INPUT_FILES  Create the default input .mat files of the simulator.
%
%   MAKE_INPUT_FILES(cfg) writes the six input files named in radar_config.m:
%     antenna_pattern.mat : element spacing, azimuth pattern, elevation element
%                           pattern, transmit elevation pattern, peak gains
%     ddc_filters.mat     : low-pass FIR filters of the two DDC stages (h1, h2) for
%                           the 2 MHz output, signed 16-bit integers with 13
%                           fractional bits
%     ddc_filters_4MHz.mat: the same for the 4 MHz output (other h2)
%     dbf_matrix.mat      : elevation DBF matrix W [nBeams x cfg.ant.nElements]
%                           and beam elevations
%     stagger_prf.mat     : staggered PRI vector and the matching MTI weights
%     phase_code.mat      : phase codes of the short and the long sub-pulse
%                           (P4, width_s/chip_s chips) and their mismatch
%                           filters (one coefficient per chip)
%   MAKE_INPUT_FILES(cfg, which) writes only one of them: which = 'antenna',
%   'ddc', 'ddc4MHz', 'dbf', 'stagger' or 'phaseCode'.
%
%   Every default file stores in 'madeFor' the configuration values it was
%   built for. setup_radar.m rebuilds a default file when they no longer match
%   cfg (e.g. after changing cfg.ant.nElements or cfg.ant.tilt_deg). Replace
%   any file with your own data (same variable names, without 'madeFor') to
%   simulate a different radar.

if nargin < 1, cfg = radar_config(); end
if nargin < 2, which = 'all'; end
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'));
folder = fileparts(cfg.files.antenna);
if ~exist(folder, 'dir'), mkdir(folder); end
lambda = cfg.const.c / cfg.fc;
elementSpacing_m = cfg.ant.elementSpacing_lambda * lambda;   % element spacing in elevation [m]
doAll = strcmpi(which, 'all');

%% 1) Antenna
if doAll || strcmpi(which, 'antenna')
    azBeamwidth_deg = 4.1;       % 3 dB azimuth beamwidth (one way)
    elBeamwidth_deg = 70;        % 3 dB elevation beamwidth of a single element (row)

    % Azimuth pattern: Taylor-weighted horizontal line array scaled to the requested beamwidth
    az_deg  = (-180:0.02:180)';
    nAz     = 32;
    wAz     = taylor_window(nAz, 4, -30);
    afPower = @(d, az) abs(exp(1j*2*pi*d*(0:nAz-1)'*sind(az(:)')).' * wAz).^2 / sum(wAz)^2;
    halfBw  = @(d) 10*log10(afPower(d, azBeamwidth_deg/2)) + 3.0103;   % zero at the right spacing
    dAz     = fzero(halfBw, [0.3 0.7]);                                % element spacing [wavelengths]
    azGain_dB = max(10*log10(afPower(dAz, az_deg)), -50);
    azGain_dB(abs(az_deg) > 90) = -50;                                 % no back lobe

    % Elevation element pattern cos^q(theta), theta measured from the array broadside
    elRel_deg     = (-90:0.1:90)';
    q             = log(0.5) / log(cosd(elBeamwidth_deg/2));            % 3 dB at +-35 deg
    elElemGain_dB = max(10*log10(cosd(elRel_deg).^q + eps), -30);
    txElGain_dB   = elElemGain_dB;  % the transmitter uses the same wide elevation beam

    % Peak gains (rule of thumb G = 26000 / (az_bw * el_bw) with beamwidths in degrees)
    rxChannelGain_dBi = 10*log10(26000 / (azBeamwidth_deg * elBeamwidth_deg));
    txGain_dBi        = rxChannelGain_dBi;

    madeFor = struct('fc', cfg.fc, 'spacing', cfg.ant.elementSpacing_lambda);   %#ok<NASGU>
    save(cfg.files.antenna, 'elementSpacing_m', 'az_deg', 'azGain_dB', 'elRel_deg', ...
        'elElemGain_dB', 'txElGain_dB', 'rxChannelGain_dBi', 'txGain_dBi', 'madeFor', '-v7');
end

%% 2) DDC low-pass filters (windowed sinc, Blackman window, unit DC gain)
% After the 30 MHz mixer the signal occupies +-0.5 MHz around 0 Hz. Stage 1
% is the same for both output rates; stage 2 is a wider filter for 4 MHz.
fs = cfg.adc.fs_Hz;
q  = 2^cfg.ddc.coefFracBits;              % fixed point: signed 16 bit, 13 fractional bits
for rate = [2e6 4e6]
    [D, key] = ddc_decimation(cfg, rate);
    if ~(doAll || strcmpi(which, key)), continue; end
    fs1 = fs / D(1);
    h1  = int16(round(q * lowpass_fir(70, 5.8e6 / fs)));     % stage 1: pass +-1 MHz, stop above 10.5 MHz
    if rate == 2e6
        h2 = int16(round(q * lowpass_fir(66, 1.0e6 / fs1)));  % stage 2: pass +-0.5 MHz, stop above 1.5 MHz
    else
        h2 = int16(round(q * lowpass_fir(66, 1.5e6 / fs1)));  % stage 2: pass +-1 MHz, stop above 2 MHz
    end
    madeFor = struct('fs_Hz', fs, 'decimation', D, 'coefFracBits', cfg.ddc.coefFracBits); %#ok<NASGU>
    save(cfg.files.(key), 'h1', 'h2', 'madeFor', '-v7');
end

%% 3) Elevation DBF matrix: 21 beams equally spaced in u = sin(angle from broadside)
if doAll || strcmpi(which, 'dbf')
    nElements  = cfg.ant.nElements;
    nBeams     = 21;
    uFirst     = sind(0.5 - cfg.ant.tilt_deg);     % lowest beam 0.5 deg above the horizon
    uStep      = 1.2 * 0.886 / (nElements * cfg.ant.elementSpacing_lambda);   % 3 dB beamwidth in u
    uLast      = min(uFirst + (nBeams-1)*uStep, sind(60 - cfg.ant.tilt_deg)); % at most 60 deg elevation
    beamU      = linspace(uFirst, uLast, nBeams)';
    beamEl_deg = cfg.ant.tilt_deg + asind(beamU);  % beam elevations (at the azimuth boresight)
    taper      = taylor_window(nElements, 4, -30);
    n          = 0:nElements-1;
    W          = (taper(:).' .* exp(-1j*2*pi*elementSpacing_m/lambda * beamU * n));
    W          = W ./ sqrt(sum(abs(W).^2, 2));     % unit noise gain for every beam
    madeFor = struct('fc', cfg.fc, 'nElements', nElements, 'tilt_deg', cfg.ant.tilt_deg, ...
                     'spacing', cfg.ant.elementSpacing_lambda);     %#ok<NASGU>
    save(cfg.files.dbf, 'W', 'beamEl_deg', 'madeFor', '-v7');
end

%% 4) Staggered PRF: PRI vector and time-varying 3-pulse MTI weights
if doAll || strcmpi(which, 'stagger')
    priVector_s = [920; 1000; 960; 1040] * 1e-6;
    K = numel(priVector_s);
    mtiCoeffs = zeros(K, 3);
    for k = 1:K
        % Pulses k, k+1, k+2 are separated by T1 and T2. The weights null DC and
        % the first-order (linear phase) clutter term; they equal [1 -2 1] when T1 = T2.
        T1 = priVector_s(k);
        T2 = priVector_s(mod(k, K) + 1);
        mtiCoeffs(k, :) = 2 * [T2, -(T1 + T2), T1] / (T1 + T2);
    end
    save(cfg.files.stagger, 'priVector_s', 'mtiCoeffs', '-v7');
end

%% 5) Phase codes (P4) and their mismatch filters
% shortCode, longCode [nChips x 1]: complex chip values, the first chip is
% sent first. shortMismatch, longMismatch [nChips + 2*pad x 1]: mismatch
% filter with one coefficient per chip, as a convolution filter
% (output = conv(code, filter)); its main peak is in the middle of the full
% convolution. pad = min(5, nChips/2) chips on each side keeps the range
% handover between the sub-pulses (with a longer filter the short pulse
% ends before the long pulse starts) and the receive window of the
% staggered PRF; more padding gives little here (25 chips: PSL -25.9 dB with
% 5, -26.0 dB with 10). The filter minimises the peak sidelobe for a unit
% main peak (mismatch_filter.m). A filter with
% another sample rate (e.g. designed at the DDC output rate on the received
% pulse) can be given with <pulse>MismatchRate_Hz (setup_radar.m).
if doAll || strcmpi(which, 'phaseCode')
    names = {'short', 'long'};
    madeFor = struct();
    for i = 1:2
        wv = cfg.wave.(names{i});
        N = round(wv.width_s / wv.chip_s);
        m = (0:N-1)';
        code = exp(1j*(pi*m.^2/N - pi*m));            % P4 polyphase code
        pad = min(5, floor(N/2));
        out.([names{i} 'Code']) = code;
        out.([names{i} 'Mismatch']) = mismatch_filter(code, pad);
        madeFor.([names{i} 'Chips']) = N;
    end
    out.madeFor = madeFor;
    save(cfg.files.phaseCode, '-struct', 'out', '-v7');
end

fprintf('Default input file(s) "%s" written to %s\n', which, folder);
end

function h = lowpass_fir(nTaps, fCutNorm)
% Windowed-sinc low-pass FIR. fCutNorm = cut-off frequency / sampling frequency.
n = (0:nTaps-1)' - (nTaps-1)/2;
x = 2*fCutNorm*n;
s = ones(size(x));
s(x ~= 0) = sin(pi*x(x ~= 0)) ./ (pi*x(x ~= 0));
k = (0:nTaps-1)';
blackman = 0.42 - 0.5*cos(2*pi*k/(nTaps-1)) + 0.08*cos(4*pi*k/(nTaps-1));
h = s .* blackman;
h = h / sum(h);
end

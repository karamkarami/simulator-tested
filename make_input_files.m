function make_input_files(cfg, which)
%MAKE_INPUT_FILES  Create the default input .mat files of the simulator.
%
%   MAKE_INPUT_FILES(cfg) writes the input files named in radar_config.m:
%     antenna_pattern.mat : element spacing, azimuth pattern, elevation element
%                           pattern, transmit elevation pattern, peak gains
%     ddc_filters.mat     : low-pass FIR filters of the two DDC stages (h1, h2) for
%                           the 2 MHz output, signed 16-bit integers with 13
%                           fractional bits
%     ddc_filters_4MHz.mat: the same for the 4 MHz output (other h2)
%     dbf_matrix.mat      : elevation DBF matrix W [nBeams x cfg.ant.nElements]
%                           and beam elevations
%     stagger_prf.mat     : pairs of a staggered PRI vector and its MTI weights
%     scan2scan_prf.mat   : PRFs of the antenna scans (prf_Hz, used in turn) for
%                           cfg.prf.mode = 'scan2scan'
%     phase_code.mat      : phase codes of the short and the long sub-pulse
%                           (cfg.wave.maxCodes per sub-pulse, the first P4,
%                           width_s/chip_s chips) and their mismatch filters
%                           (one coefficient per chip)
%     mti_fir.mat         : FIR clutter filter along slow time (coeffs), used with
%                           cfg.mti.type = 'fir'
%   MAKE_INPUT_FILES(cfg, which) writes only one of them: which = 'antenna',
%   'ddc', 'ddc4MHz', 'dbf', 'stagger', 'scan2scan', 'phaseCode' or 'mtiFir'.
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

%% 4) Staggered PRF: pairs of a PRI sequence T and its time-varying 3-pulse MTI weights w
% priVector_s_<k> [K x 1] and mtiCoeffs_<k> [K x 3] (one row per stagger position: the weights of
% the MTI outputs whose oldest pulse has that position), k = 1 ... cfg.prf.maxStaggerSets (default
% 1): the mode selects pair cfg.prf.staggerSet. Pair 1 is the stagger of the radar before the
% pairs; the others are further sequences (the same mean PRI of about 1 ms, other blind speeds).
% A file of your own may hold fewer pairs (at least 1, numbered from 1) of any lengths, or only
% priVector_s and mtiCoeffs (one pair).
if doAll || strcmpi(which, 'stagger')
    maxSets = 1;
    if isfield(cfg.prf, 'maxStaggerSets') && ~isempty(cfg.prf.maxStaggerSets), maxSets = cfg.prf.maxStaggerSets; end
    seqs = {[920; 1000; 960; 1040], [940; 1000; 920; 1060], [960; 1020; 930; 1050; 990], [920; 980; 1040; 1100], ...
            [930; 1010; 970; 1050; 950; 1030]};
    out = struct();
    for j = 1:maxSets
        T = seqs{mod(j - 1, numel(seqs)) + 1} * 1e-6;
        if j > numel(seqs), T = circshift(T, floor((j - 1) / numel(seqs))); end   % (more pairs: rotated)
        K = numel(T);
        w = zeros(K, 3);
        for k = 1:K
            % Pulses k, k+1, k+2 are separated by T1 and T2. The weights null DC and the
            % first-order (linear phase) clutter term; they equal [1 -2 1] when T1 = T2.
            T1 = T(k);
            T2 = T(mod(k, K) + 1);
            w(k, :) = 2 * [T2, -(T1 + T2), T1] / (T1 + T2);
        end
        out.(sprintf('priVector_s_%d', j)) = T;
        out.(sprintf('mtiCoeffs_%d', j)) = w;
    end
    out.madeFor = struct('nSets', maxSets);
    save(cfg.files.stagger, '-struct', 'out', '-v7');
end

%% 4b) Scan-to-scan PRF: one PRF per antenna scan (used in turn, cfg.prf.mode = 'scan2scan')
if (doAll || strcmpi(which, 'scan2scan')) && isfield(cfg.files, 'scan2scan')
    prf_Hz = [1000 1040 960];                     %#ok<NASGU> PRFs of the scans [Hz] (PRIs 1000, 962, 1042 us)
    save(cfg.files.scan2scan, 'prf_Hz', '-v7');
end

%% 5) Phase codes and their mismatch filters (code agility: several pairs per sub-pulse)
% shortCode_<k>, longCode_<k> [nChips x 1]: complex chip values, the first chip is sent first;
% shortMismatch_<k>, longMismatch_<k> [K x 1]: the mismatch filter of code k with one
% coefficient per chip, as a convolution filter (output = conv(code, filter)); its main peak is
% in the middle of the full convolution. k = 1 ... cfg.wave.maxCodes (default 1); the mode selects
% the pair cfg.wave.<pulse>.codeIndex. Code 1 is P4; the others are polyphase codes with the phases
% pi*r*m^2/N - pi*m (r coprime to N; r = 1: P4), the roots r chosen one after the other among the
% first 16 so that the echo of any other code of the set through the mismatch filter of a code
% stays as far as possible below its own peak (agile_codes below). A file of your own may hold fewer pairs (numbered from 1)
% or only <pulse>Code and <pulse>Mismatch (one pair). Filter length K =
% cfg.wave.<pulse>.mismatchLength (any K >= nChips, e.g. 3*nChips; the extra coefficients are
% split between both ends), by default nChips + 2*pad with pad = min(5, nChips/2): more padding
% gives little here (25 chips: PSL -25.9 dB with 5, -26.0 dB with 10). Any length works with the
% rest of the chain: the range cells and the receive window do not depend on it (setup_radar.m).
% The filter minimises the peak sidelobe for a unit main peak (mismatch_filter.m). A filter with
% another sample rate (e.g. designed at the DDC output rate on the received pulse) can be given
% with <pulse>MismatchRate_Hz (for all codes of the sub-pulse; setup_radar.m).
if doAll || strcmpi(which, 'phaseCode')
    names = {'short', 'long'};
    madeFor = struct();
    out = struct();
    maxCodes = 1;
    if isfield(cfg.wave, 'maxCodes') && ~isempty(cfg.wave.maxCodes), maxCodes = cfg.wave.maxCodes; end
    for i = 1:2
        wv = cfg.wave.(names{i});
        N = round(wv.width_s / wv.chip_s);
        m = (0:N-1)';
        K = mismatch_length(wv);
        if isempty(K), K = N + 2*min(5, floor(N/2)); end
        if K < N || K ~= round(K)
            error('cfg.wave.%s.mismatchLength = %g: a whole number of coefficients, at least the %d chips of the code.', ...
                  names{i}, K, N);
        end
        pad = [floor((K - N)/2), K - N - floor((K - N)/2)];
        [codes, filters] = agile_codes(N, pad, maxCodes);
        for k = 1:maxCodes
            out.(sprintf('%sCode_%d', names{i}, k)) = codes{k};
            out.(sprintf('%sMismatch_%d', names{i}, k)) = filters{k};
        end
        madeFor.([names{i} 'Chips']) = N;
        madeFor.([names{i} 'MismatchLength']) = mismatch_length(wv);
    end
    madeFor.nCodes = maxCodes;
    out.madeFor = madeFor;
    save(cfg.files.phaseCode, '-struct', 'out', '-v7');
end

%% 6) FIR clutter filter along slow time (cfg.mti.type = 'fir')
% coeffs [1 x L] (or one row per stagger position): output m = sum_l coeffs(l) * x(m+l-1),
% oldest pulse first. Default: 5 taps with a double zero at zero Doppler (as the
% three-pulse canceller) and a zero pair at +-3% of the PRF: a deeper and wider notch
% (at 2% of the PRF -88 dB instead of -44 dB, at 5% -63 dB instead of -28 dB, relative
% to white noise) against moving clutter (sea, vegetation), for a narrower passband
% (-6.5 dB instead of -1.8 dB at a quarter of the PRF). Any real FIR filter of your own
% can replace it.
if (doAll || strcmpi(which, 'mtiFir')) && isfield(cfg.files, 'mtiFir')
    theta = 2*pi*0.03;
    coeffs = conv([1 -2 1], [1 -2*cos(theta) 1]);
    coeffs = coeffs / norm(coeffs) * sqrt(6);       % the noise gain of [1 -2 1]
    save(cfg.files.mtiFir, 'coeffs', '-v7');
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

function [codes, filters] = agile_codes(N, pad, n)
% n polyphase codes of N chips (phases pi*r*m^2/N - pi*m, root r coprime to N) with their mismatch
% filters: r = 1 (P4) first, then greedily the root (among the first 16 coprime to N) whose worst
% cross-talk with the codes already chosen (peak of a code through the filter of another, relative
% to the peak of the filter's own code, both ways) is the lowest.
m = (0:N-1)';
cand = find(gcd(1:max(N-1, 1), N) == 1);
cand = cand(1:min(end, max(16, n)));
C = cell(1, numel(cand));  H = C;  P = zeros(1, numel(cand));
for i = 1:numel(cand)
    C{i} = exp(1j*(pi*cand(i)*m.^2/N - pi*m));
    H{i} = mismatch_filter(C{i}, pad);
    P(i) = max(abs(conv(C{i}, H{i})));
end
X = zeros(numel(cand));                                   % X(i, j): code j through filter i [dB]
for i = 1:numel(cand)
    for j = 1:numel(cand)
        if i ~= j, X(i, j) = 20*log10(max(abs(conv(C{j}, H{i}))) / P(i)); end
    end
end
chosen = 1;
while numel(chosen) < min(n, numel(cand))
    best = -inf;  pick = 0;
    for j = setdiff(1:numel(cand), chosen)
        worst = max([X(chosen, j); X(j, chosen).']);
        if -worst > best, best = -worst;  pick = j; end
    end
    chosen(end+1) = pick; %#ok<AGROW>
end
while numel(chosen) < n, chosen(end+1) = chosen(mod(numel(chosen), numel(cand)) + 1); end %#ok<AGROW>
codes = C(chosen);
filters = H(chosen);
end

function K = mismatch_length(wv)
% cfg.wave.<pulse>.mismatchLength ([] = default length; also for a cfg without the field).
K = [];
if isfield(wv, 'mismatchLength'), K = wv.mismatchLength; end
end

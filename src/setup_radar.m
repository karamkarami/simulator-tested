function sys = setup_radar(cfg)
%SETUP_RADAR  Load the input files and derive every dependent radar quantity.
%
%   sys = SETUP_RADAR(cfg) returns a struct with the antenna, noise, timing,
%   waveform, filter, axis, CPI and CFAR parameters used by the processing
%   blocks. The raw configuration is kept in sys.cfg.

c  = cfg.const.c;
fs = cfg.adc.fs_Hz;
sys.cfg    = cfg;
sys.lambda = c / cfg.fc;
sys.fs     = fs;
sys.fIF    = cfg.adc.fIF_Hz;                                  % IF carrier at the ADC [Hz]
sys.fLO    = cfg.fc - sys.fIF;                                % receiver local oscillator [Hz]

sys.cfar = cfar_window(cfg);                                  % window of the selected CFAR method

%% Input files (defaults are built when missing or when made for other cfg values)
ant  = load_input(cfg, 'antenna', struct('fc', cfg.fc, 'spacing', cfg.ant.elementSpacing_lambda));
[ddcDecim, ddcFile] = ddc_decimation(cfg);                    % 2 or 4 MHz output: its filter file
ddcF = load_input(cfg, ddcFile, struct('fs_Hz', fs, 'decimation', ddcDecim, ...
                                       'coefFracBits', cfg.ddc.coefFracBits));

%% Antenna
sys.ant = ant;
sys.ant.height_m = cfg.ant.height_m;
sys.ant.tilt_deg = cfg.ant.tilt_deg;
sys.ant.rate_dps = cfg.ant.rpm * 360 / 60;                    % rotation rate [deg/s]
sys.ant.azBeamwidth_deg = beamwidth_3dB(ant.az_deg, ant.azGain_dB);
sys.ant.azSpan_deg = 3 * sys.ant.azBeamwidth_deg;             % echoes outside +-azSpan are ignored
sys.nCh = cfg.ant.nElements;

% Received power per element = radarConst * rcs / R^4 * (pattern factors)
sys.radarConst = cfg.tx.peakPower_W * 10^(ant.txGain_dBi/10) * 10^(ant.rxChannelGain_dBi/10) ...
               * sys.lambda^2 / ((4*pi)^3 * 10^(cfg.rx.losses_dB/10));

%% Receive chain: receiver input -> gain -> ADC
% Signals at the receiver input are in sqrt(W): a sine of amplitude A has power A^2/2.
sys.N0 = cfg.const.kB * cfg.const.T0 * 10^(cfg.rx.noiseFigure_dB/10);   % noise density [W/Hz]
sys.rx.sigmaTh   = sqrt(sys.N0 * fs/2);                       % thermal noise per real sample [sqrt(W)]
sys.rx.voltGain  = 10^(cfg.rx.gain_dB/20);                    % receiver path gain (amplitude)
sys.rx.leakAmp   = sqrt(2 * cfg.tx.peakPower_W * 10^(-cfg.rx.txIsolation_dB/10));  % leakage [sqrt(W)]
if cfg.adc.nBits > 16
    error('cfg.adc.nBits = %d: the ADC samples are stored as int16 (at most 16 bits)', cfg.adc.nBits);
end
sys.adc.fullScale = 2^(cfg.adc.nBits - 1);
fullScale_W      = 1e-3 * 10^(cfg.adc.fullScale_dBm/10);      % full-scale sine power at the ADC [W]
sys.adc.lsbPerUnit = sys.adc.fullScale / sqrt(2*fullScale_W); % ADC input [sqrt(W)] -> LSB
sinad      = 10^((6.02*cfg.adc.enob + 1.76)/10);              % full-scale sine / ADC noise
adcNoise2  = sys.adc.fullScale^2/2 / sinad;                   % total ADC noise power [LSB^2]
% the noise added after the 16-bit quantization; with the two roundings
% (1/12 LSB^2 each) the total is adcNoise2
sys.adc.sigmaAdc = sqrt(max(adcNoise2 - 2/12, 0));            % [LSB]
sys.adc.sigmaTh  = sys.rx.sigmaTh * sys.rx.voltGain * sys.adc.lsbPerUnit;   % thermal noise at the ADC [LSB]
sys.N0eff = sys.N0 * (sys.adc.sigmaTh^2 + adcNoise2) / sys.adc.sigmaTh^2;    % incl. ADC noise

%% PRI sequence and MTI weights
switch lower(cfg.prf.mode)
    case 'fixed'
        sys.pri = cfg.prf.pri_s;
        sys.mtiCoeffs = [1 -2 1];                             % three-pulse canceller
    case 'staggered'
        st = load_input(cfg, 'stagger', []);
        sys.pri = st.priVector_s(:);
        sys.mtiCoeffs = st.mtiCoeffs;                         % one row per stagger position
    otherwise
        error('Unknown PRF mode "%s"', cfg.prf.mode);
end
sys.staggered = strcmpi(cfg.prf.mode, 'staggered');
sys.nStagger  = numel(sys.pri);
if size(sys.mtiCoeffs, 1) == 1
    sys.mtiCoeffs = repmat(sys.mtiCoeffs, sys.nStagger, 1);
end

%% DDC filters (signed 16-bit integers with coefFracBits fractional bits) and NCO
D1 = ddcDecim(1);
D2 = ddcDecim(2);
sys.ddc.h1 = double(ddcF.h1(:));                              % integer coefficients
sys.ddc.h2 = double(ddcF.h2(:));
sys.ddc.shift = cfg.ddc.shift;
sys.ddc.decim = [D1 D2];
L1 = numel(sys.ddc.h1);
L2 = numel(sys.ddc.h2);
fs1 = fs / D1;
sys.fsOut = fs1 / D2;
% Real-valued filters (unit DC gain shape) and the extra gain of each stage
% (2^(fracBits - shift) = 1 with the default shifts)
q = 2^cfg.ddc.coefFracBits;
sys.ddc.gain = 2.^(cfg.ddc.coefFracBits - cfg.ddc.shift);
sys.ddc.h1real = sys.ddc.h1 / q;
sys.ddc.h2real = sys.ddc.h2 / q;
% Single filter at the ADC rate equivalent to the two stages (for the noise bandwidth)
h2Up = zeros(D1*(L2-1)+1, 1);
h2Up(1:D1:end) = sys.ddc.h2real;
sys.ddc.hEq = conv(sys.ddc.h1real, h2Up);
sys.ddc.noiseBw_Hz = fs * sum(sys.ddc.hEq.^2);                % complex noise bandwidth [Hz]

%% PRI timeline: short pulse, receive, long pulse, receive
% One ADC record per PRI runs from the PRI start to the end of the long-pulse
% receive window. While a pulse is transmitted the receiver sees the
% transmitter leakage (these samples are not used by the processing).
names = {'short', 'long'};
codes = [];                                                   % phase codes and mismatch filters
if any(strcmpi({cfg.wave.short.type, cfg.wave.long.type}, 'phasecode'))
    nChips = [round(cfg.wave.short.width_s / cfg.wave.short.chip_s), round(cfg.wave.long.width_s / cfg.wave.long.chip_s)];
    codes = load_input(cfg, 'phaseCode', struct('shortChips', nChips(1), 'longChips', nChips(2)));
end
for i = 1:2
    wv = cfg.wave.(names{i});
    wv.name     = names{i};
    wv.lambda_m = sys.lambda;                                 % both sub-pulses are on fc
    wv.compression = lower(wv.compression);
    if strcmpi(wv.type, 'phasecode')
        wv.code = codes.([names{i} 'Code'])(:);
        if numel(wv.code) ~= round(wv.width_s / wv.chip_s)
            error('The %s code of %s has %d chips, but width_s/chip_s = %d.', names{i}, cfg.files.phaseCode, ...
                  numel(wv.code), round(wv.width_s / wv.chip_s));
        end
        if strcmp(wv.compression, 'mismatch')
            wv.mismatch = codes.([names{i} 'Mismatch'])(:);
            wv.mismatchRate_Hz = 1 / wv.chip_s;                % one coefficient per chip (default)
            if isfield(codes, [names{i} 'MismatchRate_Hz'])
                wv.mismatchRate_Hz = codes.([names{i} 'MismatchRate_Hz']);
            end
        end
    end
    % Pulse compression filter (unit noise gain) and the samples by which its
    % peak comes earlier than that of the matched filter
    [wv.mfTaps, wv.mfShift] = compression_filter(wv, cfg, sys);
    sys.wave.(names{i}) = wv;
end
ws = sys.wave.short;
wl = sys.wave.long;
if wl.tTx_s < ws.tTx_s + ws.width_s
    error('The long pulse must start after the end of the short pulse.');
end
filterSpan    = (L1-1)/fs + (L2-1)/fs1;                       % input span of one DDC output sample
sys.rxStart_s = 0;                                            % the record starts with the PRI
% (the CFAR tests the range cells up to the end of the window: its range
% window is shifted inwards at the edges, so no extra cells are needed)
% a mismatch filter also needs the samples after the end of the echo
mfTail = max(0, numel(wl.mfTaps) - round(wl.width_s*sys.fsOut) - wl.mfShift) / sys.fsOut;
rxEnd = wl.tTx_s + wl.width_s + 2*cfg.rx.maxRange_m/c + mfTail + filterSpan + 2/sys.fsOut;
if rxEnd > min(sys.pri)
    error(['Receive window (%.1f us) does not fit in the shortest PRI (%.1f us): reduce ', ...
           'cfg.rx.maxRange_m or start the long pulse earlier.'], rxEnd*1e6, min(sys.pri)*1e6);
end
sys.nAdc = round((rxEnd - sys.rxStart_s) * fs);
sys.tAdc = sys.rxStart_s + (0:sys.nAdc-1)' / fs;             % ADC sample times in the PRI [s]

% NCO of the DDC: the real line is multiplied by cos(2*pi*fIF*t), the imaginary
% line by -sin(2*pi*fIF*t). At fIF = fs/4 these are exactly 1, 0, -1, 0, ...
sys.ddc.cos      = round(cos(2*pi*sys.fIF*sys.tAdc) * 2^15) / 2^15;   % 16-bit NCO values
sys.ddc.minusSin = round(-sin(2*pi*sys.fIF*sys.tAdc) * 2^15) / 2^15;

% DDC output samples: their time (a valid output is centred on its input window, linear-phase
% FIR) and the first/last ADC sample each of them depends on
t1 = sys.tAdc(1:D1:sys.nAdc-L1+1) + (L1-1)/2/fs;
sys.tDdc = t1(1:D2:numel(t1)-L2+1) + (L2-1)/2/fs1;
sys.nDdc = numel(sys.tDdc);
k = (0:sys.nDdc-1)';
firstAdc = k*D2*D1 + 1;
lastAdc  = (k*D2 + L2 - 1)*D1 + L1;

%% Receive interval, DDC segment and delay / range axes of every sub-pulse
rxInterval.short = [ws.tTx_s + ws.width_s, wl.tTx_s];         % [start end) of the receive interval
rxInterval.long  = [wl.tTx_s + wl.width_s, inf];
for i = 1:2
    wv = sys.wave.(names{i});
    % DDC samples whose whole input window lies inside the receive interval of this sub-pulse
    iv = rxInterval.(names{i});
    wv.ddcIdx = find(sys.tAdc(firstAdc) >= iv(1) & sys.tAdc(lastAdc) < iv(2)).';
    L = numel(wv.mfTaps);                                     % compression filter length [samples]
    if numel(wv.ddcIdx) < L + 2*(sys.cfar.guardRange + sys.cfar.trainRange) + 1
        error(['The %s-pulse receive interval is too short for its compression filter and the CFAR ', ...
               'window: move the long pulse later or increase cfg.rx.maxRange_m.'], wv.name);
    end
    % Output k of the valid convolution of the matched filter correlates the
    % reference with the segment starting at sample k. The reference holds
    % the pulse at the mid-sample times, so output k peaks for an echo that
    % starts half a sample before tDdc(k). The peak of a mismatch filter
    % comes mfShift samples earlier: output k belongs to the echo starting at
    % sample k + mfShift.
    wv.nRange  = numel(wv.ddcIdx) - L + 1;
    wv.delay_s = sys.tDdc(wv.ddcIdx((1:wv.nRange) + wv.mfShift)).' - 0.5/sys.fsOut - wv.tTx_s;   % delay axis [s]
    wv.range_m = c * wv.delay_s / 2;                              % range axis [m]
    sys.wave.(names{i}) = wv;
end
% Handover: the short pulse is tested up to where the long pulse starts (its
% first range cell, tested thanks to the shifted CFAR window), so a target is
% not reported twice. The short pulse keeps oR more cells as training cells
% (and at least one whole CFAR range window).
oR = sys.cfar.guardRange + sys.cfar.trainRange;
handover_m = sys.wave.long.range_m(1);
ws = sys.wave.short;
if ws.range_m(end) + abs(diff(ws.range_m(1:2))) < handover_m
    warning('Gap in the range coverage: short pulse up to %.1f km, long pulse from %.1f km.', ...
        ws.range_m(end)/1e3, handover_m/1e3);
end
nKeep = max(find(ws.range_m < handover_m, 1, 'last') + oR, 2*oR + 1);   % + cells only used as CFAR training
if nKeep < ws.nRange
    Lp = numel(ws.mfTaps);
    ws.nRange  = nKeep;
    ws.ddcIdx  = ws.ddcIdx(1:nKeep + Lp - 1);                 % the compression filter then gives nKeep cells
    ws.delay_s = ws.delay_s(1:nKeep);
    ws.range_m = ws.range_m(1:nKeep);
    sys.wave.short = ws;
end
sys.handover_m = handover_m;

% Range cells tested by the CFAR (both sub-pulses): cfg.cfar.<pulse>.min/maxRange_m
% and, for the short pulse, the handover. testIdx = [first last] (1-based).
for i = 1:2
    wv = sys.wave.(names{i});
    lim = cfg.cfar.(names{i});
    ok = wv.range_m >= lim.minRange_m & wv.range_m <= lim.maxRange_m;
    if i == 1, ok = ok & wv.range_m < handover_m; end
    first = find(ok, 1, 'first');
    last  = find(ok, 1, 'last');
    if isempty(first)
        warning('No %s-pulse range cell is tested: check cfg.cfar.%s.minRange_m / maxRange_m.', names{i}, names{i});
        first = 1;  last = 0;
    end
    wv.testIdx = [first last];
    sys.wave.(names{i}) = wv;
end

%% Digital beamforming (disabled with a single receive element)
sys.dbf.enabled = sys.nCh > 1;
if sys.dbf.enabled
    dbfF = load_input(cfg, 'dbf', struct('fc', cfg.fc, 'nElements', sys.nCh, 'tilt_deg', cfg.ant.tilt_deg, ...
                                         'spacing', cfg.ant.elementSpacing_lambda));
    if size(dbfF.W, 2) ~= sys.nCh
        error('The DBF matrix has %d columns but cfg.ant.nElements = %d.', size(dbfF.W, 2), sys.nCh);
    end
    sys.dbf.W = dbfF.W;
    sys.dbf.beamEl_deg = dbfF.beamEl_deg(:);
    sys.dbf.beamU = sind(sys.dbf.beamEl_deg - sys.ant.tilt_deg);
else
    sys.dbf.W = [];
    sys.dbf.beamEl_deg = NaN;                                 % no elevation measurement
    sys.dbf.beamU = NaN;
end
sys.nBeams = numel(sys.dbf.beamEl_deg);

% All clutter arrives from (about) the same elevation: one steering vector
sys.clutterSteer = steering_vector(sys, 0, cfg.clutter.elevation_deg);

%% CPI buffer and Doppler processing
% CPI = time a target stays inside the 3 dB azimuth beam; the buffer holds
% the pulses that arrive during one CPI.
sys.cpi.time_s = sys.ant.azBeamwidth_deg / sys.ant.rate_dps;
if isempty(cfg.proc.cpiLength)
    sys.cpi.length = round(sys.cpi.time_s / mean(sys.pri));
else
    sys.cpi.length = cfg.proc.cpiLength;
end
if sys.staggered
    % Staggered PRF: a whole number of stagger periods in the buffer
    sys.cpi.length = sys.nStagger * floor(sys.cpi.length / sys.nStagger);
end
nTaps = size(sys.mtiCoeffs, 2);
nOut  = sys.cpi.length - (nTaps - 1);                         % MTI output samples
if sys.staggered
    nRows  = sys.nStagger;                                    % one row per stagger position
    rowLen = floor(nOut / nRows);
    rowPri = sum(sys.pri);                                    % sampling period inside a row
else
    nRows  = 1;                                               % one FFT over the whole CPI
    rowLen = nOut;
    rowPri = sys.pri;
end
nfft = cfg.proc.nfft;
if isempty(nfft) && sys.staggered
    nfft = 2^nextpow2(rowLen);                                % 28 samples per row -> 32
elseif isempty(nfft)
    nfft = 2^nextpow2(sys.cpi.length);                        % 114 pulses -> 128
end
rowLen = min(rowLen, nfft);                                   % a shorter FFT keeps the newest samples
if rowLen < 2
    error('The CPI of %d pulses is too short for the MTI and the Doppler processing (cfg.proc.cpiLength).', ...
          sys.cpi.length);
end
sys.cpi.nRows  = nRows;
sys.cpi.rowLen = rowLen;
w = hamming_window(rowLen);
sys.dop.window = w / norm(w);
sys.dop.nfft   = nfft;
sys.dop.rowPri_s = rowPri;
sys.dop.f_Hz   = ((0:nfft-1) - floor(nfft/2)) / (nfft*rowPri); % after fftshift
sys.dop.v_mps  = -sys.dop.f_Hz * sys.lambda / 2;              % radial velocity (+ = away)

% Noise power of every Doppler bin for white noise at the MTI input. The MTI
% colours the noise (e.g. +4 dB near PRF/2 for [1 -2 1]); doppler_process
% divides by this gain so that the CFAR sees a flat noise floor.
n = sys.cpi.length;
sys.dop.noiseGain = 1;
impulses = reshape(eye(n), 1, n, n);                          % one impulse per "range cell"
st = mod(0:n-1, sys.nStagger) + 1;
[pN, Y] = doppler_process(mti_filter(impulses, sys.mtiCoeffs, st), sys);
sys.dop.noiseGain = sum(pN, 2);
% Effective number of independent looks of the non-coherent integration.
% De-interleaved rows (staggered PRF) are built from overlapping MTI outputs,
% so their noise is correlated. The tail of the averaged noise power decays
% with the largest eigenvalue of the row covariance, so the number of looks
% is trace(C) / max(eig(C)): 1 for a single FFT, 1.5 for the default stagger.
nLooks = zeros(1, nfft);
for k = 1:nfft
    Yk = reshape(Y(1, :, k, :), n, nRows);
    C  = Yk' * Yk;                                            % noise covariance of the rows
    nLooks(k) = real(trace(C)) / max(real(eig(C)));
end
% MTI notch: Doppler bins around zero (clutter residue) that the CFAR does not test.
% The residue is as wide as the Doppler window mainlobe: +-2*nfft/rowLen bins for Hamming.
sys.dop.notchBins = round(cfg.cfar.notchWidth * 2 * nfft / rowLen);
zeroBin = floor(nfft/2) + 1;
tested = abs((1:nfft) - zeroBin) > sys.dop.notchBins;
if ~any(tested)
    error('The MTI notch (cfg.cfar.notchWidth = %g: +-%d bins) covers every Doppler bin (FFT length %d).', ...
          cfg.cfar.notchWidth, sys.dop.notchBins, nfft);
end
sys.dop.nLooks = min(nLooks(tested));
% Median / mean of the noise power (used to estimate the noise floor of a beam)
L = sys.dop.nLooks;
sys.dop.medianToMean = fzero(@(x) gammainc(L*x, L) - 0.5, [0.1 2]);

%% CFAR threshold of the selected method
c = sys.cfar;
nOuter = (2*(c.guardRange + c.trainRange) + 1) * (2*(c.guardDoppler + c.trainDoppler) + 1);
nGuard = (2*c.guardRange + 1) * (2*c.guardDoppler + 1);
sys.cfar.nTrain = nOuter - nGuard;
if isempty(c.threshold_dB)
    sys.cfar.alpha = cfar_alpha(c.pfa, sys.cfar.nTrain, sys.dop.nLooks);
else
    sys.cfar.alpha = 10^(c.threshold_dB/10);
end
sys.cfar.pfa = cfar_pfa(sys.cfar.alpha, sys.cfar.nTrain, sys.dop.nLooks);   % of the threshold in use
if 2*(c.guardDoppler + c.trainDoppler) + 1 > nfft
    error('CFAR Doppler window (%d bins) is larger than the FFT length (%d).', ...
          2*(c.guardDoppler + c.trainDoppler) + 1, nfft);
end

%% Tracking
sys.track = track_setup(cfg, sys);
end

function [taps, shift] = compression_filter(wv, cfg, sys)
% Pulse compression filter of a sub-pulse at the DDC rate (row vector, unit
% noise gain) and the number of samples by which its peak comes earlier than
% that of the matched filter.
%   matched  : conj(flip(reference)), the reference sampled at the mid-sample
%              times (Taylor window for LFM); shift 0
%   mismatch : phase code only: every coefficient of wv.mismatch (rate
%              wv.mismatchRate_Hz, by default one per chip) is held for the
%              DDC samples of its period (with one per chip: a chip-matched
%              sum followed by the chip-rate filter). shift: from the peaks of
%              both filters for the pulse after the DDC filters (the filter
%              may have been designed on that pulse or on the ideal chips)
fsOut = sys.fsOut;
Lp  = round(wv.width_s * fsOut);                              % pulse length [samples]
ref = waveform_envelope(wv, ((0:Lp-1)' + 0.5) / fsOut);       % mid-sample times (not on chip edges)
win = ones(Lp, 1);
if strcmpi(wv.type, 'lfm')
    win = taylor_window(Lp, cfg.proc.taylorNbar, cfg.proc.taylorSll_dB);
end
hMatched = conj(flipud(ref .* win));
switch wv.compression
    case 'matched'
        h = hMatched;
        shift = 0;
    case 'mismatch'
        if ~strcmpi(wv.type, 'phasecode')
            error('cfg.wave.%s.compression = ''mismatch'' needs a phase code (type ''phasecode'').', wv.name);
        end
        spc = fsOut / wv.mismatchRate_Hz;                     % DDC samples per filter coefficient
        if abs(spc - round(spc)) > 1e-9 || spc < 0.5
            error(['The %s-pulse mismatch filter (%g MHz) does not fit the DDC output rate (%g MHz): ' ...
                   'the rate must be a whole multiple of the filter rate.'], wv.name, ...
                  wv.mismatchRate_Hz/1e6, fsOut/1e6);
        end
        h = kron(wv.mismatch(:), ones(round(spc), 1));
        % valid output index of the peak = echo start + peak of the full
        % convolution - (filter length - 1), for both filters
        y = ddc_pulse(sys, wv);                               % the pulse after the DDC filters
        qMatched = peak_position(conv(y, hMatched));
        q = peak_position(conv(y, h));
        shift = round((numel(h) - 1 - q) - (Lp - 1 - qMatched));
    otherwise
        error('cfg.wave.%s.compression must be ''matched'' or ''mismatch''.', wv.name);
end
taps = (h / norm(h)).';
end

function q = peak_position(y)
% Position (0-based, with parabolic interpolation) of the largest |y|.
a = abs(y(:));
[~, k] = max(a);
q = k - 1;
if k > 1 && k < numel(a)
    d = a(k-1) - 2*a(k) + a(k+1);
    if d < 0, q = q + 0.5*(a(k-1) - a(k+1)) / d; end
end
end

function w = cfar_window(cfg)
% Window and threshold settings of the selected CFAR method, as a range x Doppler window.
c = cfg.cfar;
w.method = lower(c.method);
w.notchWidth = c.notchWidth;
switch w.method
    case '2d'
        m = c.twoD;
        w.guardRange = m.guardRange;    w.trainRange = m.trainRange;
        w.guardDoppler = m.guardDoppler;  w.trainDoppler = m.trainDoppler;
        w.label = '2-D';
    case 'range'
        m = c.range;
        w.guardRange = m.guard;  w.trainRange = m.train;  w.guardDoppler = 0;  w.trainDoppler = 0;
        w.label = '1-D range';
    case 'doppler'
        m = c.doppler;
        w.guardRange = 0;  w.trainRange = 0;  w.guardDoppler = m.guard;  w.trainDoppler = m.train;
        w.label = '1-D Doppler';
    otherwise
        error('cfg.cfar.method must be ''2d'', ''range'' or ''doppler''.');
end
w.pfa = m.pfa;
w.threshold_dB = m.threshold_dB;
end

function bw = beamwidth_3dB(az_deg, gain_dB)
% 3 dB width of a pattern, interpolated between the table points.
g = gain_dB - max(gain_dB);
[~, iMax] = max(g);
iL = find(g(1:iMax) < -3.0103, 1, 'last');                    % last point below -3 dB on the left
iR = iMax - 1 + find(g(iMax:end) < -3.0103, 1, 'first');      % first point below -3 dB on the right
azL = interp1(g(iL:iL+1), az_deg(iL:iL+1), -3.0103);
azR = interp1(g(iR-1:iR), az_deg(iR-1:iR), -3.0103);
bw = azR - azL;
end

function data = load_input(cfg, which, madeFor)
% Load input file cfg.files.(which). A missing file, or a default file built by
% make_input_files for other configuration values (its 'madeFor'), is rebuilt.
file = cfg.files.(which);
if exist(file, 'file') ~= 2
    make_input_files(cfg, which);
end
data = load(file);
if isfield(data, 'madeFor') && ~isequal(data.madeFor, madeFor)
    make_input_files(cfg, which);
    data = load(file);
end
end

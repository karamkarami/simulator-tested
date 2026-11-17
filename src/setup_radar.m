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
sys.report.snrNoise = 'median';                               % noise floor of the SNR of the reports
sys.report.dopplerInterp = false;                             % peak interpolation (peak_interpolation.m)
sys.report.elevationInterp = true;
if isfield(cfg, 'report')
    if isfield(cfg.report, 'snrNoise'), sys.report.snrNoise = lower(cfg.report.snrNoise); end
    if isfield(cfg.report, 'dopplerInterp'), sys.report.dopplerInterp = logical(cfg.report.dopplerInterp); end
    if isfield(cfg.report, 'elevationInterp'), sys.report.elevationInterp = logical(cfg.report.elevationInterp); end
end
if ~any(strcmp(sys.report.snrNoise, {'median', 'cfar'}))
    error('cfg.report.snrNoise must be ''median'' or ''cfar''.');
end

%% Input files (defaults are built when missing or when made for other cfg values)
ant  = load_input(cfg, 'antenna', struct('fc', cfg.fc, 'spacing', cfg.ant.elementSpacing_lambda));
[ddcDecim, ddcFile] = ddc_decimation(cfg);                    % 2 or 4 MHz output: its filter file
ddcF = load_input(cfg, ddcFile, struct('fs_Hz', fs, 'decimation', ddcDecim, ...
                                       'coefFracBits', cfg.ddc.coefFracBits));

%% Antenna
sys.ant = ant;
sys.ant.height_m = cfg.ant.height_m;
sys.ant.tilt_deg = cfg.ant.tilt_deg;
[sys.ant.rpm, sys.ant.availableRpm] = rpm_choice(cfg.ant);    % rotation speed of the mode [rev/min]
sys.ant.rate_dps = sys.ant.rpm * 360 / 60;                    % rotation rate [deg/s]
sys.ant.azBeamwidth_deg = beamwidth_3dB(ant.az_deg, ant.azGain_dB);
sys.ant.azSpan_deg = 3 * sys.ant.azBeamwidth_deg;             % echoes outside +-azSpan are ignored
sys.nCh = cfg.ant.nElements;
% channels of the DDC output the processing uses: the first cfg.rx.nChannelsUsed ([] or missing: all)
sys.nChUsed = sys.nCh;
if isfield(cfg.rx, 'nChannelsUsed') && ~isempty(cfg.rx.nChannelsUsed), sys.nChUsed = cfg.rx.nChannelsUsed; end
if ~isscalar(sys.nChUsed) || ~isreal(sys.nChUsed) || sys.nChUsed ~= round(sys.nChUsed) || sys.nChUsed < 1 || ...
        sys.nChUsed > sys.nCh
    error('cfg.rx.nChannelsUsed must be a whole number 1 ... cfg.ant.nElements = %d (or [] = all).', sys.nCh);
end

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

%% Radar mode (cfg.modes, apply_radar_mode.m)
sys.mode = radar_mode_check(cfg);

%% PRI sequence and MTI weights
sys.scan2scan = false;
switch lower(cfg.prf.mode)
    case 'fixed'
        sys.pri = cfg.prf.pri_s;
        sys.mtiCoeffs = [1 -2 1];                             % three-pulse canceller
    case 'staggered'
        % several (PRI sequence T, MTI weights w) pairs in cfg.files.stagger: the mode may use the
        % pairs cfg.prf.availableStaggerSets and selects cfg.prf.staggerSet (stagger_set below)
        [setIdx, sets] = stagger_choice(cfg);
        st = load_input(cfg, 'stagger', struct('nSets', max(sets)));
        [sys.pri, sys.mtiCoeffs] = stagger_set(st, setIdx, sets, cfg.files.stagger);
        sys.availableStaggerSets = sets;
        sys.staggerSet = setIdx;                              % weights: one row per stagger position
    case 'scan2scan'
        % one PRF per antenna scan (constant during the scan), used in turn: scan s has
        % prf_Hz(mod(s - 1, numel(prf_Hz)) + 1). Inside a scan as the fixed PRF.
        if ~isfield(cfg.files, 'scan2scan')
            error('cfg.prf.mode = ''scan2scan'' needs the PRF file cfg.files.scan2scan.');
        end
        s2 = load_input(cfg, 'scan2scan', []);
        if ~isfield(s2, 'prf_Hz') || isempty(s2.prf_Hz) || ~isnumeric(s2.prf_Hz) || ~isreal(s2.prf_Hz) || ...
                any(~isfinite(s2.prf_Hz(:))) || any(s2.prf_Hz(:) <= 0)
            error('%s must hold the PRFs of the scans [Hz] in the variable ''prf_Hz'' (positive, finite).', ...
                  cfg.files.scan2scan);
        end
        sys.prfScan_Hz = double(s2.prf_Hz(:)).';
        sys.priScan_s = 1 ./ sys.prfScan_Hz;
        sys.pri = sys.priScan_s(1);                          % the slow-time settings below: of scan 1
        sys.mtiCoeffs = [1 -2 1];
        sys.scan2scan = true;
    otherwise
        error('Unknown PRF mode "%s" (''fixed'', ''staggered'' or ''scan2scan'')', cfg.prf.mode);
end
sys.staggered = strcmpi(cfg.prf.mode, 'staggered');
if ~sys.staggered, sys.staggerSet = 0;  sys.availableStaggerSets = []; end   % (no stagger pair in use)
sys.nStagger  = numel(sys.pri);
% Clutter filter: the canceller above, or an FIR filter from cfg.files.mtiFir
sys.mti.type = 'canceller';
if isfield(cfg, 'mti') && isfield(cfg.mti, 'type'), sys.mti.type = lower(cfg.mti.type); end
switch sys.mti.type
    case 'canceller'
    case 'fir'
        if ~isfield(cfg.files, 'mtiFir')
            error('cfg.mti.type = ''fir'' needs the coefficient file cfg.files.mtiFir.');
        end
        fir = load_input(cfg, 'mtiFir', []);
        if ~isfield(fir, 'coeffs') || isempty(fir.coeffs) || ~isnumeric(fir.coeffs)
            error('%s must hold the FIR coefficients in the variable ''coeffs''.', cfg.files.mtiFir);
        end
        h = double(fir.coeffs);
        if ~isreal(h) || any(~isfinite(h(:)))
            error('The FIR coefficients of %s must be real and finite.', cfg.files.mtiFir);
        end
        if isvector(h), h = h(:).'; end                       % a column is one filter too
        if size(h, 1) ~= 1 && size(h, 1) ~= sys.nStagger
            error('%s: ''coeffs'' has %d rows: 1 (one filter) or %d (one per stagger position).', ...
                  cfg.files.mtiFir, size(h, 1), sys.nStagger);
        end
        if ~any(h(:))
            error('The FIR coefficients of %s are all zero.', cfg.files.mtiFir);
        end
        sys.mtiCoeffs = h;
    otherwise
        error('cfg.mti.type must be ''canceller'' or ''fir''.');
end
if size(sys.mtiCoeffs, 1) == 1
    sys.mtiCoeffs = repmat(sys.mtiCoeffs, sys.nStagger, 1);
end
% Range limit of the clutter filter (cfg.mti.maxRange_m, Inf = every cell). The cells beyond it
% get the "bypass" weights: only the last tap (output m = pulse m + L - 1), so their CPI has the
% same nOut samples and the same time as the filtered cells.
sys.mti.maxRange_m = Inf;
if isfield(cfg, 'mti') && isfield(cfg.mti, 'maxRange_m') && ~isempty(cfg.mti.maxRange_m)
    sys.mti.maxRange_m = cfg.mti.maxRange_m;
end
if ~isscalar(sys.mti.maxRange_m) || ~isreal(sys.mti.maxRange_m) || isnan(sys.mti.maxRange_m) || sys.mti.maxRange_m < 0
    error('cfg.mti.maxRange_m must be a non-negative scalar range in metres (Inf is allowed).');
end
sys.mtiBypass = zeros(size(sys.mtiCoeffs));
sys.mtiBypass(:, end) = 1;

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
% Noise at the DDC output per I or Q component [LSB^2]: receiver and ADC
% noise through both filters, and the rounding of both stages to 16 bit (the
% stage-1 rounding filtered by stage 2). When the receiver noise is not well
% above 1 LSB the rounding costs SNR: sys.N0ddc is N0eff with it.
adcLevel = (sys.adc.sigmaTh^2 + adcNoise2) * sum(sys.ddc.hEq.^2) * prod(sys.ddc.gain)^2 / 2;
rounding = cfg.ddc.fixedPoint * (1 + sum(sys.ddc.h2real.^2) * sys.ddc.gain(2)^2) / 12;
sys.ddc.noiseLsb = sqrt(adcLevel + rounding);                 % rms of I or Q [LSB]
sys.N0ddc = sys.N0eff * (adcLevel + rounding) / adcLevel;      % noise density after the DDC

%% PRI timeline: short pulse, receive, long pulse, receive
% One ADC record per PRI runs from the PRI start to the end of the long-pulse
% receive window. While a pulse is transmitted the receiver sees the
% transmitter leakage (these samples are not used by the processing).
names = {'short', 'long'};
codes = [];                                                   % phase codes and mismatch filters
% Code agility: cfg.files.phaseCode holds several code / mismatch-filter pairs per sub-pulse
% (<pulse>Code_<k>, <pulse>Mismatch_<k>; a file with only <pulse>Code, <pulse>Mismatch holds pair 1);
% the mode may use the pairs cfg.wave.<pulse>.availableCodes and selects cfg.wave.<pulse>.codeIndex.
if isfield(cfg.wave, 'maxCodes')
    error(['cfg.wave.maxCodes no longer exists: cfg.wave.<pulse>.availableCodes lists the code / mismatch ' ...
           'pairs of cfg.files.phaseCode a sub-pulse may use (cfg.wave.<pulse>.codeIndex is one of them).']);
end
avCodes = {available_numbers(cfg.wave.short, 'availableCodes', 'wave.short.availableCodes', 'code / mismatch pairs'), ...
           available_numbers(cfg.wave.long, 'availableCodes', 'wave.long.availableCodes', 'code / mismatch pairs')};
maxCodes = max([avCodes{:}]);                                 % the default file holds pairs 1 ... maxCodes
if any(strcmpi({cfg.wave.short.type, cfg.wave.long.type}, 'phasecode'))
    nChips = [round(cfg.wave.short.width_s / cfg.wave.short.chip_s), round(cfg.wave.long.width_s / cfg.wave.long.chip_s)];
    codes = load_input(cfg, 'phaseCode', struct('shortChips', nChips(1), 'shortMismatchLength', ...
                       mismatch_length(cfg.wave.short), 'longChips', nChips(2), 'longMismatchLength', ...
                       mismatch_length(cfg.wave.long), 'nCodes', maxCodes));
end
for i = 1:2
    wv = cfg.wave.(names{i});
    wv.name     = names{i};
    % sub-pulse on or off (cfg.wave.<pulse>.enable) and its carrier (cfg.wave.<pulse>.fc_Hz, [] = cfg.fc).
    % Every sub-pulse reaches the ADC on the same IF cfg.adc.fIF_Hz: the receiver LO is tuned to
    % fc_Hz - fIF while the sub-pulse is transmitted and received (it switches to the LO of the long
    % pulse when the long pulse is transmitted), and the NCO of the DDC stays at fIF.
    if ~isfield(wv, 'enable') || isempty(wv.enable), wv.enable = true; end
    wv.enable = logical(wv.enable);
    if ~isfield(wv, 'fc_Hz') || isempty(wv.fc_Hz), wv.fc_Hz = cfg.fc; end
    if ~isscalar(wv.fc_Hz) || ~(wv.fc_Hz > sys.fIF) || ~isfinite(wv.fc_Hz)
        error('cfg.wave.%s.fc_Hz must be a carrier frequency [Hz] above cfg.adc.fIF_Hz (or [] = cfg.fc).', names{i});
    end
    wv.lambda_m = c / wv.fc_Hz;
    wv.fIF_Hz   = sys.fIF;                                   % IF at the ADC (the same for both)
    wv.fLO_Hz   = wv.fc_Hz - sys.fIF;                        % receiver LO of this sub-pulse
    wv.radarConst = sys.radarConst * (wv.lambda_m / sys.lambda)^2;   % received power ~ lambda^2
    wv.compression = lower(wv.compression);
    if strcmpi(wv.type, 'phasecode')
        [wv.code, wv.codeIndex, sfx] = code_pair(codes, names{i}, wv, avCodes{i}, cfg.files.phaseCode);
        if numel(wv.code) ~= round(wv.width_s / wv.chip_s)
            error('The %s code %d of %s has %d chips, but width_s/chip_s = %d.', names{i}, wv.codeIndex, ...
                  cfg.files.phaseCode, numel(wv.code), round(wv.width_s / wv.chip_s));
        end
        if wv.enable && wv.chip_s * sys.fsOut < 1 - 1e-9
            % every chip must hold at least one DDC sample: the echo is sampled per chip
            % (waveform_envelope: chip floor(t/chip_s) at every ADC sample) and the filters
            % hold every chip / mismatch coefficient for its DDC samples (2 per 1 us chip at 2 MHz)
            error(['cfg.wave.%s.chip_s = %.3g us is shorter than the DDC sample period %.3g us ' ...
                   '(%.3g MHz): the code bandwidth 1/chip_s does not fit the DDC output.'], ...
                  names{i}, wv.chip_s*1e6, 1e6/sys.fsOut, sys.fsOut/1e6);
        end
        if strcmp(wv.compression, 'mismatch')
            if ~isfield(codes, [names{i} 'Mismatch' sfx])
                error('%s has no mismatch filter ''%sMismatch%s'' for the %s code %d.', cfg.files.phaseCode, ...
                      names{i}, sfx, names{i}, wv.codeIndex);
            end
            wv.mismatch = codes.([names{i} 'Mismatch' sfx])(:);
            wv.mismatchRate_Hz = 1 / wv.chip_s;                % one coefficient per chip (default)
            if isfield(codes, [names{i} 'MismatchRate_Hz'])
                wv.mismatchRate_Hz = codes.([names{i} 'MismatchRate_Hz']);
            end
        end
    end
    wv.availableCodes = avCodes{i};
    if ~strcmpi(wv.type, 'phasecode'), wv.codeIndex = 0;  wv.availableCodes = []; end   % (0: no phase code)
    % Pulse compression filter (unit noise gain) and the samples by which its
    % peak comes earlier than that of the matched filter
    [wv.mfTaps, wv.mfShift] = compression_filter(wv, cfg, sys);
    sys.wave.(names{i}) = wv;
end
ws = sys.wave.short;
wl = sys.wave.long;
if ~ws.enable && ~wl.enable
    error('Both sub-pulses are switched off (cfg.wave.short.enable, cfg.wave.long.enable).');
end
if ws.enable && wl.enable && wl.tTx_s < ws.tTx_s + ws.width_s
    error('The long pulse must start after the end of the short pulse.');
end
sys.pulseOn = [ws.enable, wl.enable];                         % sub-pulses transmitted and processed
filterSpan    = (L1-1)/fs + (L2-1)/fs1;                       % input span of one DDC output sample
sys.rxStart_s = 0;                                            % the record starts with the PRI
% Receive window of every sub-pulse [start end) in the PRI: cfg.wave.<pulse>.rxStart_s / rxEnd_s
% (times from the PRI start, as tTx_s). A set start / end selects the DDC output samples by their
% time (the centre of their input); [] = the DDC samples whose whole input lies after the end of
% the pulse and before the transmission of the long pulse (short pulse) or the end of the record
% (no transmitter leakage). The record (ADC samples, DDC output) runs from the PRI start to the
% end of the last receive window: rxEnd_s of the last transmitted sub-pulse (plus half the input
% span of a DDC sample), or without it the end of its echoes from cfg.rx.maxRange_m. The
% processing takes only the receive windows out of the DDC output (the CFAR tests the range
% cells up to the end of the window: its range window is shifted inwards at the edges, so no
% extra cells are needed). The record, and so the DDC output, depends only on the timing, the
% DDC filters and the decimation: not on the pulse compression filters or the CFAR (a filter
% longer than the pulse sees zeros outside the receive window).
wLast = wl;                                                   % the last transmitted sub-pulse ends the record
if ~wl.enable, wLast = ws; end
[rxIv, rxEnd, rxSet] = receive_windows(ws, wl, wLast.tTx_s + wLast.width_s + 2*cfg.rx.maxRange_m/c + filterSpan + 2/sys.fsOut, ...
                                     filterSpan/2 + 1/sys.fsOut);
priMin = min(sys.pri);
if sys.scan2scan, priMin = min(sys.priScan_s); end
manualDdc = isfield(cfg.rx, 'manualDdc') && ~isempty(cfg.rx.manualDdc) && logical(cfg.rx.manualDdc);
if ~manualDdc && rxEnd > priMin + 1e-12
    error(['Receive window (%.1f us) does not fit in the shortest PRI (%.1f us): reduce ', ...
           'cfg.rx.maxRange_m (or cfg.wave.%s.rxEnd_s) or start the long pulse earlier.'], ...
           rxEnd*1e6, priMin*1e6, wLast.name);
end
sys.nAdc = round((rxEnd - sys.rxStart_s) * fs);
% Manual DDC record of the mode (cfg.rx.manualDdc): cfg.rx.nDdc DDC samples, the first one at
% cfg.rx.ddcT0_s ([] = the automatic time, half the input span of a DDC sample after the PRI
% start). The record then starts half that span before ddcT0_s (on the ADC sample grid: ddcT0_s
% is rounded to it; before the PRI start too, e.g. ddcT0_s = 0) and holds the ADC samples of
% exactly nDdc DDC samples. The receive windows are taken out of it as with the automatic record.
m0 = 0;                                                       % first ADC sample of the record
halfSpan = (L1-1)/2/fs + (L2-1)/2/fs1;                        % ddcT0_s of the automatic record
sys.ddcRecord.manual = manualDdc;
if sys.ddcRecord.manual
    nSet = [];  if isfield(cfg.rx, 'nDdc'), nSet = cfg.rx.nDdc; end
    if ~isscalar(nSet) || ~isreal(nSet) || nSet < 1 || nSet ~= round(nSet)
        error('cfg.rx.manualDdc = true needs cfg.rx.nDdc: the number of DDC samples per PRI (a whole number >= 1).');
    end
    t0 = halfSpan;
    if isfield(cfg.rx, 'ddcT0_s') && ~isempty(cfg.rx.ddcT0_s), t0 = cfg.rx.ddcT0_s; end
    if ~isscalar(t0) || ~isreal(t0) || ~isfinite(t0)
        error('cfg.rx.ddcT0_s must be the time of the first DDC sample in the PRI [s] (or [] = automatic).');
    end
    m0 = round((t0 - halfSpan) * fs);
    sys.nAdc = (nSet - 1)*D1*D2 + (L2 - 1)*D1 + L1;
    sys.rxStart_s = m0 / fs;
    recEnd = (m0 + sys.nAdc) / fs;
    if recEnd > priMin + 1e-12
        error(['The manual DDC record (cfg.rx.nDdc = %d samples from cfg.rx.ddcT0_s = %.2f us) ends at %.1f us, ', ...
               'after the shortest PRI (%.1f us).'], nSet, (m0/fs + halfSpan)*1e6, recEnd*1e6, priMin*1e6);
    end
    for i = 1:2
        nm = {'short', 'long'};
        iv = rxIv.(nm{i});
        if sys.pulseOn(i) && rxSet.(nm{i})(2) && iv(2) > recEnd - halfSpan + 1e-12
            error(['cfg.wave.%s.rxEnd_s = %.2f us is after the last DDC sample of the manual record ', ...
                   '(%.2f us: cfg.rx.nDdc, cfg.rx.ddcT0_s).'], nm{i}, iv(2)*1e6, (recEnd - halfSpan)*1e6);
        end
    end
end
sys.ddcRecord.m0 = m0;                                        % (ADC samples before the PRI start: -m0)
sys.tAdc = (m0 + (0:sys.nAdc-1)') / fs;                       % ADC sample times in the PRI [s]
shortEnd = inf;  if wl.enable, shortEnd = wl.tTx_s; end      % (without the long pulse the short pulse
sys.shortRxEnd_s = shortEnd;                                  %  is followed by no other one)
sys.rxWindow_s.short = rxIv.short;                            % receive windows [start end) [s]
sys.rxWindow_s.long  = rxIv.long;
sys.adc.rxSamples = find((ws.enable & sys.tAdc >= rxIv.short(1) & sys.tAdc < rxIv.short(2)) | ...
                         (wl.enable & sys.tAdc >= rxIv.long(1) & sys.tAdc < rxIv.long(2)));  % (adc_overload.m)
sys.rx.stcGain = single(stc_gain(sys, sys.tAdc));            % STC in front of the receiver (1 = off)

% Receiver LO at every ADC sample: that of the short pulse, from the transmission of the long
% pulse on that of the long pulse. An echo of sub-pulse p received while the LO is that of the
% other sub-pulse lands on fIF + (fc_p - fc_other) (the DDC filters reject it when the carriers are
% more than the bandwidth apart): sys.rx.loShift(:, p) = exp(-1j*2*pi*(fLO(t) - fLO_p)*(t - tTx_p))
% multiplies its complex IF signal ([] = one LO for the whole PRI).
fLO = repmat(ws.fLO_Hz, sys.nAdc, 1);
if ~ws.enable, fLO(:) = wl.fLO_Hz; end
if ws.enable && wl.enable, fLO(sys.tAdc >= wl.tTx_s) = wl.fLO_Hz; end
sys.rx.fLO_Hz = fLO;
sys.rx.loShift = [];
if ws.enable && wl.enable && ws.fc_Hz ~= wl.fc_Hz
    sys.rx.loShift = single([exp(-1j*2*pi*(fLO - ws.fLO_Hz).*(sys.tAdc - ws.tTx_s)), ...
                             exp(-1j*2*pi*(fLO - wl.fLO_Hz).*(sys.tAdc - wl.tTx_s))]);
end

% NCO of the DDC: the real line is multiplied by cos(2*pi*fIF*t), the imaginary
% line by -sin(2*pi*fIF*t). At fIF = fs/4 these are exactly 1, 0, -1, 0, ...
sys.ddc.fNco_Hz  = sys.fIF;
sys.ddc.cos      = round(cos(2*pi*sys.fIF*sys.tAdc) * 2^15) / 2^15;   % 16-bit NCO values
sys.ddc.minusSin = round(-sin(2*pi*sys.fIF*sys.tAdc) * 2^15) / 2^15;

% DDC output samples: their time (a valid output is centred on its input window, linear-phase
% FIR) and the first/last ADC sample each of them depends on
t1 = sys.tAdc(1:D1:sys.nAdc-L1+1) + (L1-1)/2/fs;
sys.tDdc = t1(1:D2:numel(t1)-L2+1) + (L2-1)/2/fs1;
sys.nDdc = numel(sys.tDdc);
sys.ddcT0_s = sys.tDdc(1);                                    % time of the first DDC sample
k = (0:sys.nDdc-1)';
firstAdc = k*D2*D1 + 1;
lastAdc  = (k*D2 + L2 - 1)*D1 + L1;

%% Receive interval, DDC segment and delay / range axes of every sub-pulse
% rxIdx : DDC samples whose whole input window lies inside the receive
%         interval of the sub-pulse (sent in the PRI header). They depend only
%         on the timing, the DDC filters and the decimation.
% ddcIdx: input of the pulse compression (valid convolution), nRange + L - 1
%         samples. Output k belongs to the echo starting at rxIdx(k), for
%         every echo that lies entirely inside the receive interval, whatever
%         the length of the compression filter: a filter longer than the
%         pulse (mismatch filter) needs samples before and after the interval,
%         which are taken as zeros (ddc_segment.m), so its range cells are the
%         same as those of the matched filter.
rxInterval = rxIv;                                            % [start end) of the receive interval
for i = 1:2
    wv = sys.wave.(names{i});
    iv = rxInterval.(names{i});
    if ~wv.enable                                             % switched off: no receive interval, no cells
        wv.rxIdx = zeros(1, 0);  wv.nRange = 0;  wv.ddcIdx = zeros(1, 0);
        wv.delay_s = zeros(1, 0);  wv.range_m = zeros(1, 0);
        sys.wave.(names{i}) = wv;
        continue
    end
    % a set window start / end (cfg.wave.<pulse>.rxStart_s / rxEnd_s): the DDC samples whose time
    % (the centre of their input) lies in it; [] = the DDC samples whose whole input lies between
    % the end of the pulse and the next transmission / the end of the record (no leakage)
    set = rxSet.(names{i});
    inStart = sys.tAdc(firstAdc) >= iv(1);  if set(1), inStart = sys.tDdc >= iv(1); end
    inEnd   = sys.tAdc(lastAdc) < iv(2);    if set(2), inEnd = sys.tDdc < iv(2); end
    wv.rxIdx = find(inStart & inEnd).';
    Lp = round(wv.width_s * sys.fsOut);                       % pulse length [samples]
    L  = numel(wv.mfTaps);                                    % compression filter length [samples]
    wv.nRange = numel(wv.rxIdx) - Lp + 1;
    if wv.nRange < 2*(sys.cfar.guardRange + sys.cfar.trainRange) + 1
        error(['The %s-pulse receive interval is too short for its pulse and the CFAR ', ...
               'window: move the long pulse later or increase cfg.rx.maxRange_m (or widen ', ...
               'cfg.wave.%s.rxStart_s / rxEnd_s).'], wv.name, wv.name);
    end
    % Output k of the valid convolution of the matched filter correlates the
    % reference with the segment starting at sample k. The reference holds
    % the pulse at the mid-sample times, so output k peaks for an echo that
    % starts half a sample before tDdc of that sample. The peak of a mismatch
    % filter comes mfShift samples earlier: output k belongs to the echo
    % starting at segment sample k + mfShift, so the segment starts mfShift
    % samples before the receive interval.
    wv.ddcIdx  = wv.rxIdx(1) - wv.mfShift + (0:wv.nRange + L - 2);
    wv.delay_s = sys.tDdc(wv.rxIdx(1:wv.nRange)).' - 0.5/sys.fsOut - wv.tTx_s;   % delay axis [s]
    wv.range_m = c * wv.delay_s / 2;                              % range axis [m]
    sys.wave.(names{i}) = wv;
end
% Handover: the short pulse is tested up to where the long pulse starts (its
% first range cell, tested thanks to the shifted CFAR window), so a target is
% not reported twice. The short pulse keeps oR more cells as training cells
% (and at least one whole CFAR range window). Only its processing is
% shortened (ddcIdx, range cells): its receive interval rxIdx, and so the PRI
% header, stays the same for every CFAR method.
oR = sys.cfar.guardRange + sys.cfar.trainRange;
handover_m = inf;                                             % (no long pulse: the short pulse to the end)
if sys.wave.long.enable, handover_m = sys.wave.long.range_m(1); end
ws = sys.wave.short;
if ws.enable && wl.enable && ws.range_m(end) + abs(diff(ws.range_m(1:2))) < handover_m
    warning('Gap in the range coverage: short pulse up to %.1f km, long pulse from %.1f km.', ...
        ws.range_m(end)/1e3, handover_m/1e3);
end
nKeep = max(find(ws.range_m < handover_m, 1, 'last') + oR, 2*oR + 1);   % + cells only used as CFAR training
if ws.enable && wl.enable && nKeep < ws.nRange
    ws.nRange  = nKeep;
    ws.ddcIdx  = ws.ddcIdx(1:nKeep + numel(ws.mfTaps) - 1);   % the compression filter then gives nKeep cells
    ws.delay_s = ws.delay_s(1:nKeep);
    ws.range_m = ws.range_m(1:nKeep);
    sys.wave.short = ws;
end
sys.handover_m = handover_m;

% Range cells tested by the CFAR (both sub-pulses): cfg.cfar.<pulse>.min/maxRange_m
% and, for the short pulse, the handover. testIdx = [first last] (1-based).
for i = 1:2
    wv = sys.wave.(names{i});
    if ~wv.enable
        wv.testIdx = [1 0];  wv.mtiCells = 0;
        sys.wave.(names{i}) = wv;
        continue
    end
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
    wv.mtiCells = sum(wv.range_m <= sys.mti.maxRange_m);      % cells 1 ... mtiCells get the clutter filter
    sys.wave.(names{i}) = wv;
end

%% Digital beamforming (disabled with a single receive element, or with one channel used)
% With cfg.rx.nChannelsUsed < cfg.ant.nElements the DBF uses only the first nChannelsUsed elements:
% the first nChannelsUsed columns of the DBF matrix, every beam scaled to unit noise gain again (a
% matrix of your own may also have nChannelsUsed columns: then it is used as it is).
sys.dbf.enabled = sys.nChUsed > 1;
if sys.dbf.enabled
    dbfF = load_input(cfg, 'dbf', struct('fc', cfg.fc, 'nElements', sys.nCh, 'tilt_deg', cfg.ant.tilt_deg, ...
                                         'spacing', cfg.ant.elementSpacing_lambda));
    if size(dbfF.W, 2) ~= sys.nCh && size(dbfF.W, 2) ~= sys.nChUsed
        error('The DBF matrix has %d columns but cfg.ant.nElements = %d.', size(dbfF.W, 2), sys.nCh);
    end
    sys.dbf.W = dbfF.W;
    if size(sys.dbf.W, 2) ~= sys.nChUsed                      % the sub-array of the first elements
        sys.dbf.W = sys.dbf.W(:, 1:sys.nChUsed);
        sys.dbf.W = sys.dbf.W ./ sqrt(sum(abs(sys.dbf.W).^2, 2));
    end
    sys.dbf.beamEl_deg = dbfF.beamEl_deg(:);
    sys.dbf.beamU = sind(sys.dbf.beamEl_deg - sys.ant.tilt_deg);
else
    sys.dbf.W = [];
    sys.dbf.beamEl_deg = NaN;                                 % no elevation measurement
    sys.dbf.beamU = NaN;
end
sys.nBeams = numel(sys.dbf.beamEl_deg);
% The DBF weights steer a beam with fixed element phases: on another carrier the beam
% points elsewhere, u = beamU * cfg.fc / fc (beam elevation of every sub-pulse, NaN without DBF)
for i = 1:2
    wv = sys.wave.(names{i});
    wv.beamU = sys.dbf.beamU * (cfg.fc / wv.fc_Hz);
    wv.beamU(abs(wv.beamU) > 1) = sign(wv.beamU(abs(wv.beamU) > 1));
    wv.beamEl_deg = sys.ant.tilt_deg + asind(wv.beamU);
    sys.wave.(names{i}) = wv;
end

%% Ground clutter: reflectivity from the CNR, elevation groups of the terrain
% cfg.clutter.cnr_dB is the CNR of land at cfg.clutter.cnrRange_m at the
% matched filter output of one element (clutter_cnr.m): it sets the
% reflectivity sigma0 of land, the other terrain types are relative to it.
% Without cnr_dB ([]) sigma0 is cfg.clutter.sigma0_dB. Every elevation group
% of the terrain (clutter_terrain.m: flat ground, each hill or mountain
% region) has its own steering vector (columns of sys.clutterSteer).
cc = cfg.clutter;
if ~isfield(cc, 'cnrRange_m') || isempty(cc.cnrRange_m), cc.cnrRange_m = 5000; end
if ~(cc.cnrRange_m > 0), error('cfg.clutter.cnrRange_m must be a positive range [m].'); end
cnrPerSigma0 = clutter_cnr(sys, 1, cc.cnrRange_m, cc.elevation_deg);   % CNR of sigma0 = 1 (0 dB)
if isfield(cc, 'cnr_dB') && ~isempty(cc.cnr_dB)
    sys.clutter.sigma0_dB = cc.cnr_dB - 10*log10(cnrPerSigma0);
else
    sys.clutter.sigma0_dB = cc.sigma0_dB;
end
sys.clutter.cnr_dB = sys.clutter.sigma0_dB + 10*log10(cnrPerSigma0);
sys.clutter.cnrRange_m = cc.cnrRange_m;
sys.clutter.groupEl_deg = cc.elevation_deg;
if cc.enable
    T = clutter_terrain(cfg, sys);
    sys.clutter.groupEl_deg = T.groupEl_deg;
    sys.clutter.typeNames = T.typeNames;
    sys.clutter.maxRange_m = T.range_m(end);
    clear T
end
sys.clutterSteer = zeros(sys.nCh, numel(sys.clutter.groupEl_deg));
for g = 1:numel(sys.clutter.groupEl_deg)
    sys.clutterSteer(:, g) = steering_vector(sys, 0, sys.clutter.groupEl_deg(g));
end
for i = 1:2                                                   % on the carrier of every sub-pulse
    wv = sys.wave.(names{i});
    wv.clutterSteer = sys.clutterSteer;
    if wv.fc_Hz ~= cfg.fc
        for g = 1:numel(sys.clutter.groupEl_deg)
            wv.clutterSteer(:, g) = steering_vector(sys, 0, sys.clutter.groupEl_deg(g), wv.lambda_m);
        end
    end
    sys.wave.(names{i}) = wv;
end

%% CPI buffer and Doppler processing (slow_time_setup below; with the scan-to-scan PRF one
% set per PRF of the scans: sys.prfProfiles, select_prf_profile.m)
sys = slow_time_setup(sys, cfg);
if sys.scan2scan
    prof = struct('pri', {}, 'cpi', {}, 'dop', {}, 'v_mps', {});
    for j = 1:numel(sys.priScan_s)
        sj = sys;
        sj.pri = sys.priScan_s(j);
        sj = slow_time_setup(sj, cfg);
        prof(j).pri = sj.pri;
        prof(j).cpi = sj.cpi;
        prof(j).dop = sj.dop;
        prof(j).v_mps = {sj.wave.short.v_mps, sj.wave.long.v_mps};
    end
    sys.prfProfiles = prof;
    sys.prfProfile = 1;
end

%% CFAR threshold of the selected detector and window
c = sys.cfar;
nOuter = (2*(c.guardRange + c.trainRange) + 1) * (2*(c.guardDoppler + c.trainDoppler) + 1);
nGuard = (2*c.guardRange + 1) * (2*c.guardDoppler + 1);
sys.cfar.nTrain = nOuter - nGuard;
% GO / SO: cells per half window (2-D and range: leading / lagging in range,
% without the training cells on the range of the cell; Doppler: below / above)
oR = c.guardRange + c.trainRange;  oD = c.guardDoppler + c.trainDoppler;
if strcmp(c.method, 'doppler')
    sys.cfar.nHalf = c.trainDoppler;
else
    sys.cfar.nHalf = oR * (2*oD + 1) - c.guardRange * (2*c.guardDoppler + 1);
end
sys.cfar.osK = min(sys.cfar.nTrain, max(1, round(c.osRank * sys.cfar.nTrain)));   % OS: rank of the statistic
% CE: the training cells kept (the round(ceCensor * nTrain) largest are left out of the mean)
sys.cfar.ceKeep = sys.cfar.nTrain - round(c.ceCensor * sys.cfar.nTrain);
if strcmp(c.type, 'ce') && sys.cfar.ceKeep < 1
    error('cfg.cfar.ceCensor = %g leaves no training cell of the %d (CE-CFAR).', c.ceCensor, sys.cfar.nTrain);
end
if any(strcmp(c.type, {'go', 'so'})) && sys.cfar.nHalf < 1
    error('cfg.cfar.type = ''%s'' needs training cells on both sides of the cell (cfg.cfar.%s).', c.type, c.method);
end
if isempty(c.threshold_dB)
    sys.cfar.alpha = cfar_alpha(c.pfa, sys.cfar, sys.dop.nLooks);
else
    sys.cfar.alpha = 10^(c.threshold_dB/10);
end
sys.cfar.pfa = cfar_pfa(sys.cfar.alpha, sys.cfar, sys.dop.nLooks);   % of the threshold in use
% lowest threshold: the factor of the largest allowed Pfa (cfg.cfar.maxPfa)
maxPfa = 1e-4;
if isfield(cfg.cfar, 'maxPfa') && ~isempty(cfg.cfar.maxPfa), maxPfa = cfg.cfar.maxPfa; end
if ~isscalar(maxPfa) || ~isreal(maxPfa) || ~(maxPfa > 0 && maxPfa < 1)
    error('cfg.cfar.maxPfa must be a scalar in (0, 1).');
end
sys.cfar.maxPfa = maxPfa;
sys.cfar.alphaMin = cfar_alpha(maxPfa, sys.cfar, sys.dop.nLooks);
if sys.cfar.alpha < sys.cfar.alphaMin * (1 - 1e-9)
    error(['CFAR threshold %.2f dB (Pfa %.2e) is below the lowest threshold %.2f dB of ' ...
           'cfg.cfar.maxPfa = %.1e.'], 10*log10(sys.cfar.alpha), sys.cfar.pfa, ...
          10*log10(sys.cfar.alphaMin), maxPfa);
end
% mean of the statistic for noise of mean 1: the noise estimate of a cell is statistic / noiseScale
sys.cfar.noiseScale = 1;
if ~strcmp(c.type, 'ca')
    [zq, wq] = cfar_quadrature(sys.dop.nLooks);
    sys.cfar.noiseScale = sum(wq .* (1 - cfar_statistic(zq, sys.cfar, sys.dop.nLooks)));
end
nfftMin = sys.dop.nfft;
if sys.scan2scan, nfftMin = min(arrayfun(@(q) q.dop.nfft, sys.prfProfiles)); end
if 2*(c.guardDoppler + c.trainDoppler) + 1 > nfftMin
    error('CFAR Doppler window (%d bins) is larger than the FFT length (%d).', ...
          2*(c.guardDoppler + c.trainDoppler) + 1, nfftMin);
end

%% VBDR stream and board settings (cfg.udp.format, cfg.udp.vbdr, cfg.board, cfg.fcMin_Hz / fcMax_Hz)
sys.vbdr = vbdr_setup(cfg, sys.nDdc, sys.nChUsed);

%% Tracking
sys.track = track_setup(cfg, sys);
end

function sys = slow_time_setup(sys, cfg)
% CPI length, MTI output, Doppler rows, FFT length, window, Doppler and velocity axes, noise
% gain and looks of every Doppler bin, MTI notch and median / mean of the noise: everything that
% depends on the PRI (sys.pri). setup_radar runs it once, and with the scan-to-scan PRF once per
% PRF of the scans (sys.prfProfiles).
names = {'short', 'long'};
% CPI = time a target stays inside the 3 dB azimuth beam; the buffer holds
% the pulses that arrive during one CPI.
sys.cpi.time_s = sys.ant.azBeamwidth_deg / sys.ant.rate_dps;
% cfg.proc.cpiMode: 'auto' (from the beam dwell) or 'manual' (cfg.proc.cpiLength); without
% the field (older configurations) cpiLength = [] means auto
cpiMode = 'manual';
if isfield(cfg.proc, 'cpiMode') && ~isempty(cfg.proc.cpiMode)
    cpiMode = lower(cfg.proc.cpiMode);
elseif ~isfield(cfg.proc, 'cpiLength') || isempty(cfg.proc.cpiLength)
    cpiMode = 'auto';
end
switch cpiMode
    case 'auto'
        sys.cpi.length = round(sys.cpi.time_s / mean(sys.pri));
    case 'manual'
        if ~isfield(cfg.proc, 'cpiLength') || isempty(cfg.proc.cpiLength) || ~isscalar(cfg.proc.cpiLength) || ...
                cfg.proc.cpiLength < 1 || cfg.proc.cpiLength ~= round(cfg.proc.cpiLength)
            error('cfg.proc.cpiMode = ''manual'' needs cfg.proc.cpiLength: a whole number of pulses.');
        end
        sys.cpi.length = cfg.proc.cpiLength;
    otherwise
        error('cfg.proc.cpiMode must be ''auto'' or ''manual''.');
end
sys.cpi.mode = cpiMode;
if sys.staggered
    % Staggered PRF: a whole number of stagger periods in the buffer
    sys.cpi.length = sys.nStagger * floor(sys.cpi.length / sys.nStagger);
end
nTaps = size(sys.mtiCoeffs, 2);
nOut  = sys.cpi.length - (nTaps - 1);                         % MTI output samples
% Doppler processing of the staggered PRF (cfg.proc.staggerDoppler): 'deinterleave' (rows of one
% stagger position, FFT, non-coherent sum of the rows) or 'nufft' (the whole CPI as one non-uniformly
% sampled sequence: non-uniform DFT on a velocity grid, nufft_setup below)
sys.dop.method = 'fft';
if sys.staggered && isfield(cfg.proc, 'staggerDoppler') && ~isempty(cfg.proc.staggerDoppler)
    switch lower(cfg.proc.staggerDoppler)
        case 'deinterleave'
        case 'nufft', sys.dop.method = 'nufft';
        otherwise, error('cfg.proc.staggerDoppler must be ''deinterleave'' or ''nufft''.');
    end
end
if strcmp(sys.dop.method, 'nufft')
    sys = nufft_setup(sys, cfg, nOut);
    return
end
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
[w, sys.dop.mainlobe] = doppler_window(cfg, rowLen);         % cfg.proc.dopplerWindow
sys.dop.window = w / norm(w);
sys.dop.nfft   = nfft;
sys.dop.rowPri_s = rowPri;
sys.dop.f_Hz   = ((0:nfft-1) - floor(nfft/2)) / (nfft*rowPri); % after fftshift
sys.dop.v_mps  = -sys.dop.f_Hz * sys.lambda / 2;              % radial velocity (+ = away) at cfg.fc
for i = 1:2                                                   % and at the carrier of every sub-pulse
    sys.wave.(names{i}).v_mps = -sys.dop.f_Hz * sys.wave.(names{i}).lambda_m / 2;
end

% Noise power of every Doppler bin for white noise at the MTI input. The MTI
% colours the noise (e.g. +4 dB near PRF/2 for [1 -2 1]); doppler_process
% divides by this gain so that the CFAR sees a flat noise floor.
n = sys.cpi.length;
sys.dop.noiseGain = 1;
impulses = reshape(eye(n), 1, n, n);                          % one impulse per "range cell"
st = mod(0:n-1, sys.nStagger) + 1;
[pN, Y] = doppler_process(mti_filter(impulses, sys.mtiCoeffs, st), sys);
sys.dop.noiseGain = sum(pN, 2);
% the same for the range cells beyond cfg.mti.maxRange_m (no clutter filter)
sys.dop.noiseGainBypass = sys.dop.noiseGain;
bypassUsed = sys.wave.short.mtiCells < sys.wave.short.nRange || sys.wave.long.mtiCells < sys.wave.long.nRange;
if bypassUsed
    sysB = sys;  sysB.dop.noiseGain = 1;                      % (doppler_process divides by it)
    [pB, YB] = doppler_process(mti_filter(impulses, sys.mtiBypass, st), sysB);
    sys.dop.noiseGainBypass = sum(pB, 2);
end
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
    if bypassUsed                                             % the unfiltered cells: the fewer looks count
        Yk = reshape(YB(1, :, k, :), n, nRows);
        C  = Yk' * Yk;
        nLooks(k) = min(nLooks(k), real(trace(C)) / max(real(eig(C))));
    end
end
% MTI notch: Doppler bins around zero (clutter residue) that the CFAR does not test.
% The residue is as wide as the Doppler window mainlobe: +-m*nfft/rowLen bins (m = 2 for Hamming
% and Hann, 3 for Blackman, 1 for the rectangular window).
sys.dop.notchBins = round(cfg.cfar.notchWidth * sys.dop.mainlobe * nfft / rowLen);
zeroBin = floor(nfft/2) + 1;
sys.dop.zeroBin = zeroBin;                                    % zero Doppler after fftshift (1-based)
if sys.cfar.notchMaxRange_m < cfg.rx.maxRange_m
    tested = true(1, nfft);                  % beyond the notch range, every Doppler bin is tested
else
    tested = abs((1:nfft) - zeroBin) > sys.dop.notchBins;
end
if ~any(tested)
    error('The MTI notch (cfg.cfar.notchWidth = %g: +-%d bins) covers every Doppler bin (FFT length %d).', ...
          cfg.cfar.notchWidth, sys.dop.notchBins, nfft);
end
sys.dop.nLooks = min(nLooks(tested));
% Median / mean of the noise power (used to estimate the noise floor of a beam)
L = sys.dop.nLooks;
sys.dop.medianToMean = fzero(@(x) gammainc(L*x, L) - 0.5, [0.1 2]);
end

function sys = nufft_setup(sys, cfg, nOut)
% Staggered PRF with cfg.proc.staggerDoppler = 'nufft': the MTI outputs of the CPI (nOut, oldest
% first; output m at the time of its oldest pulse) are one non-uniformly sampled sequence. The
% Doppler spectrum is its non-uniform DFT on cfg.proc.nufft.nPoints frequencies,
%   Y(f_j) = sum_m a_m y_m exp(-j 2 pi f_j t_m),   P = |Y|^2,
% with a_m the window cfg.proc.nufft.window evaluated at the sample times (t_m scaled to [0, 1]:
% a non-uniformly sampled window) times the density compensation (the mean spacing of the
% neighbouring samples, cfg.proc.nufft.densityComp) and the bins f_j = -2 v_j / lambda of the velocity
% grid from cfg.proc.nufft.vMax_mps down to vMin_mps (nPoints bins of (vMax - vMin)/nPoints, ordered by
% increasing Doppler frequency like the FFT after fftshift; vMin itself is the next bin). The times
% depend on the stagger position of the oldest output (the phase of the CPI): one matrix A(:, :, p) per
% position p, every row scaled so that white noise at the MTI input gives unit noise power in every
% bin (no further noise-gain equalisation; the range cells beyond cfg.mti.maxRange_m, without clutter
% filter, are divided by their own gain bypassGain(:, p)). Every bin is a linear combination of the
% samples: one look of the noise (nLooks = 1).
nu = struct();
if isfield(cfg.proc, 'nufft') && ~isempty(cfg.proc.nufft), nu = cfg.proc.nufft; end
K = 2^nextpow2(nOut);
if isfield(nu, 'nPoints') && ~isempty(nu.nPoints), K = nu.nPoints; end
if ~isscalar(K) || ~isreal(K) || K < 2 || K ~= round(K) || K > 65536
    error('cfg.proc.nufft.nPoints must be a whole number 2 ... 65536 (velocity bins of the NUFFT).');
end
if nOut < 2
    error('The CPI of %d pulses is too short for the MTI and the Doppler processing (cfg.proc.cpiLength).', ...
          sys.cpi.length);
end
lam = sys.lambda;
vMax = lam / (4 * mean(sys.pri));                             % default: +- half the span of the mean PRF
vMin = -vMax;
if isfield(nu, 'vMax_mps') && ~isempty(nu.vMax_mps), vMax = nu.vMax_mps; end
if isfield(nu, 'vMin_mps') && ~isempty(nu.vMin_mps), vMin = nu.vMin_mps; end
if ~isscalar(vMax) || ~isscalar(vMin) || ~isfinite(vMax) || ~isfinite(vMin) || ~(vMax > vMin)
    error('cfg.proc.nufft.vMin_mps must be smaller than vMax_mps (finite velocities [m/s]).');
end
dens = true;
if isfield(nu, 'densityComp') && ~isempty(nu.densityComp), dens = logical(nu.densityComp); end
[wname, wp] = window_spec(nu, 'window', 'hamming', 'cfg.proc.nufft.');
dv = (vMax - vMin) / K;
df = 2 * dv / lam;
f = -2 * vMax / lam + (0:K-1) * df;                           % increasing frequency, decreasing velocity
n = sys.cpi.length;
L = size(sys.mtiCoeffs, 2);
nS = sys.nStagger;
A = zeros(K, nOut, nS);
gB = zeros(K, nS);
for p = 1:nS                                                  % stagger position of the oldest pulse
    pos = mod(p - 1 + (0:n-1), nS) + 1;                       % positions of the pulses of the CPI
    tp = [0, cumsum(reshape(sys.pri(pos(1:n-1)), 1, []))];    % their times; pulse k is followed by pri(pos(k))
    t = tp(1:nOut);                                           % output m: its oldest pulse
    w = window_function(wname, nOut, wp, t / t(end));
    if dens
        d = diff(t);
        dd = [d(1), (d(1:end-1) + d(2:end)) / 2, d(end)];
        w = w .* (dd(:) / mean(dd));
    end
    a = exp(-2i*pi * f(:) * t) .* w(:).';                     % [K x nOut]
    H = zeros(nOut, n);                                       % the MTI of this phase
    for m = 1:nOut
        H(m, m:m+L-1) = sys.mtiCoeffs(pos(m), :);
    end
    g = sum(abs(a * H).^2, 2);                                % noise gain of every bin
    g = max(g, eps * max(g));
    a = a ./ sqrt(g);
    A(:, :, p) = a;
    gB(:, p) = sum(abs(a).^2, 2);                             % white noise without the clutter filter
end
[w1, mainlobe] = window_function(wname, nOut, wp);
sys.cpi.nRows = 1;
sys.cpi.rowLen = nOut;
sys.dop.nfft = K;
sys.dop.window = w1 / norm(w1);
sys.dop.mainlobe = mainlobe;
sys.dop.rowPri_s = 1 / (K * df);                              % velocity span of the grid: lambda/(2*rowPri)
sys.dop.f_Hz = f;
sys.dop.v_mps = -f * lam / 2;
names = {'short', 'long'};
for i = 1:2
    sys.wave.(names{i}).v_mps = -f * sys.wave.(names{i}).lambda_m / 2;
end
sys.dop.noiseGain = ones(K, 1);                               % (in the matrices)
sys.dop.noiseGainBypass = ones(K, 1);
sys.dop.nufft = struct('A', A, 'bypassGain', gB, 'vMin_mps', vMin, 'vMax_mps', vMax, 'window', wname, ...
                       'densityComp', dens);
% MTI notch: the window mainlobe of the CPI duration around zero velocity
sys.dop.notchBins = round(cfg.cfar.notchWidth * mainlobe / (nOut * mean(sys.pri)) / df);
j0 = round(2 * vMax / lam / df) + 1;                          % the bin of zero velocity (1-based)
sys.dop.zeroBin = 0;                                          % 0: zero velocity is not on the grid (no notch)
if j0 >= 1 && j0 <= K && abs(f(j0)) <= df/2 * (1 + 1e-9), sys.dop.zeroBin = j0; end
if sys.dop.zeroBin > 0 && ~(sys.cfar.notchMaxRange_m < cfg.rx.maxRange_m) && 2*sys.dop.notchBins + 1 >= K
    error('The MTI notch (cfg.cfar.notchWidth = %g: +-%d bins) covers every Doppler bin (NUFFT %d points).', ...
          cfg.cfar.notchWidth, sys.dop.notchBins, K);
end
sys.dop.nLooks = 1;
sys.dop.medianToMean = log(2);                                % median / mean of an exponential power
end

function [rpm, list] = rpm_choice(ant)
% cfg.ant.availableRpm: the rotation speeds the mode may use (1 ... 64 different positive numbers,
% missing or empty: cfg.ant.rpm alone); cfg.ant.rpm (missing or empty: the first of them): the one in use.
rpm = [];
if isfield(ant, 'rpm'), rpm = ant.rpm; end
list = rpm;
if isfield(ant, 'availableRpm') && ~isempty(ant.availableRpm), list = ant.availableRpm; end
if isempty(list)
    error('cfg.ant.availableRpm (rotation speeds of the mode) and cfg.ant.rpm are not set (radar mode file).');
end
if ~isnumeric(list) || ~isreal(list) || numel(list) > 64 || any(~isfinite(list(:))) || any(list(:) <= 0) || ...
        numel(unique(list(:))) ~= numel(list)
    error('cfg.ant.availableRpm must hold 1 to 64 different positive rotation speeds [rev/min] (is %s).', mat2str(list));
end
list = double(list(:)).';
if isempty(rpm), rpm = list(1); end
if ~isnumeric(rpm) || ~isscalar(rpm) || ~any(rpm == list)
    error('cfg.ant.rpm = %s is not one of cfg.ant.availableRpm = %s.', mat2str(rpm), mat2str(list));
end
rpm = double(rpm);
end

function [k, sets] = stagger_choice(cfg)
% cfg.prf.availableStaggerSets (default 1): the (T, w) pairs of the stagger file the mode may use;
% cfg.prf.staggerSet (default: the first of them): the pair in use.
if isfield(cfg.prf, 'maxStaggerSets')
    error(['cfg.prf.maxStaggerSets no longer exists: cfg.prf.availableStaggerSets lists the (T, w) pairs of ' ...
           'cfg.files.stagger the mode may use (cfg.prf.staggerSet is one of them).']);
end
sets = available_numbers(cfg.prf, 'availableStaggerSets', 'prf.availableStaggerSets', '(T, w) pairs');
k = chosen_number(cfg.prf, 'staggerSet', sets, 'prf', 'availableStaggerSets');
end

function v = available_numbers(s, f, path, what)
% s.(f): 1 ... 64 different whole numbers >= 1 (missing or empty: 1), as a row.
v = 1;
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
if ~isnumeric(v) || ~isreal(v) || numel(v) > 64 || any(~isfinite(v(:))) || any(v(:) < 1) || ...
        any(v(:) ~= round(v(:))) || numel(unique(v(:))) ~= numel(v)
    error('cfg.%s must hold 1 to 64 different whole numbers >= 1 (the %s the mode may use).', path, what);
end
v = double(v(:)).';
end

function k = chosen_number(s, f, list, prefix, listName)
% s.(f) (missing or empty: the first of list) must be one of list (cfg.<prefix>.<listName>).
k = list(1);
if isfield(s, f) && ~isempty(s.(f)), k = s.(f); end
if ~isnumeric(k) || ~isscalar(k) || ~any(k == list)
    error('cfg.%s.%s = %s is not one of cfg.%s.%s = %s.', prefix, f, mat2str(k), prefix, listName, mat2str(list));
end
k = double(k);
end

function [pri, w] = stagger_set(st, k, sets, file)
% Pair k of a stagger file: priVector_s_<k> and mtiCoeffs_<k>, or priVector_s and mtiCoeffs (a file
% with one pair: pair 1). Every pair the mode may use (sets) must be in the file, its MTI weights
% with one row per stagger position: the time-varying filter has length(T) weight sequences.
for j = sets
    [T, wj, sfx] = stagger_pair(st, j, file);
    if isempty(T) || any(~isfinite(T)) || any(T <= 0)
        error('%s: the PRIs ''priVector_s%s'' must be positive [s].', file, sfx);
    end
    if size(wj, 1) ~= numel(T) || isempty(wj) || any(~isfinite(wj(:)))
        error(['%s: ''mtiCoeffs%s'' needs %d rows, one weight sequence per stagger position of ' ...
               '''priVector_s%s'' (the MTI is time varying).'], file, sfx, numel(T), sfx);
    end
    if j == k, pri = T;  w = wj; end
end
end

function [T, w, sfx] = stagger_pair(st, j, file)
% Variables of pair j: priVector_s_<j>, mtiCoeffs_<j> (or priVector_s, mtiCoeffs for pair 1).
sfx = sprintf('_%d', j);
if ~isfield(st, ['priVector_s' sfx]) && j == 1 && isfield(st, 'priVector_s'), sfx = ''; end
if ~isfield(st, ['priVector_s' sfx])
    error('cfg.prf.availableStaggerSets holds pair %d, but %s has no ''priVector_s_%d''.', j, file, j);
end
if ~isfield(st, ['mtiCoeffs' sfx])
    error('%s: the MTI weights ''mtiCoeffs%s'' of pair %d are missing.', file, sfx, j);
end
T = double(st.(['priVector_s' sfx])(:));
w = double(st.(['mtiCoeffs' sfx]));
end

function [code, k, sfx] = code_pair(codes, name, wv, list, file)
% Code cfg.wave.<name>.codeIndex (default: the first of cfg.wave.<name>.availableCodes) of the
% phase-code file: <name>Code_<k>, or <name>Code (a file with one code: pair 1); sfx is the suffix
% of its variables. Every pair of availableCodes must be in the file with width_s/chip_s chips
% (and its mismatch filter, for mismatch compression).
k = chosen_number(wv, 'codeIndex', list, ['wave.' name], 'availableCodes');
for j = list
    s = sprintf('_%d', j);
    if ~isfield(codes, [name 'Code' s]) && j == 1 && isfield(codes, [name 'Code']), s = ''; end
    if ~isfield(codes, [name 'Code' s])
        error('cfg.wave.%s.availableCodes holds pair %d, but %s has no ''%sCode_%d''.', name, j, file, name, j);
    end
    if numel(codes.([name 'Code' s])) ~= round(wv.width_s / wv.chip_s)
        error('The %s code %d of %s has %d chips, but width_s/chip_s = %d.', name, j, file, ...
              numel(codes.([name 'Code' s])), round(wv.width_s / wv.chip_s));
    end
    if strcmpi(wv.compression, 'mismatch') && ~isfield(codes, [name 'Mismatch' s])
        error('%s has no mismatch filter ''%sMismatch%s'' for the %s code %d.', file, name, s, name, j);
    end
    if j == k, sfx = s; end
end
code = codes.([name 'Code' sfx])(:);
end

function [w, mainlobe] = doppler_window(cfg, N)
% Window of the Doppler FFT (cfg.proc.dopplerWindow, default 'hamming'; N x 1, not normalised) and
% the half-width of its mainlobe in units of nfft/N bins (window_function.m).
[name, p] = window_spec(cfg.proc, 'dopplerWindow', 'hamming');
[w, mainlobe] = window_function(name, N, p);
end

function [name, p] = window_spec(s, key, default, path)
% Window s.(key) (default: default) and its parameters s.(key)Sll_dB, s.(key)Nbar, s.(key)Beta
% (window_function.m). The range window also takes the older cfg.proc.taylorNbar / taylorSll_dB.
% path: the prefix of the fields in the messages (default 'cfg.proc.').
if nargin < 4, path = 'cfg.proc.'; end
name = default;
if isfield(s, key) && ~isempty(s.(key)), name = s.(key); end
if ~ischar(name), error('%s%s must be the name of a window.', path, key); end
name = lower(name);
p = struct('sll_dB', -35, 'nbar', 4, 'beta', 6);
if strcmp(key, 'rangeWindow')
    if isfield(s, 'taylorSll_dB') && ~isempty(s.taylorSll_dB), p.sll_dB = s.taylorSll_dB; end
    if isfield(s, 'taylorNbar') && ~isempty(s.taylorNbar), p.nbar = s.taylorNbar; end
end
f = {'Sll_dB', 'sll_dB'; 'Nbar', 'nbar'; 'Beta', 'beta'};
for i = 1:size(f, 1)
    if isfield(s, [key f{i, 1}]) && ~isempty(s.([key f{i, 1}])), p.(f{i, 2}) = s.([key f{i, 1}]); end
end
known = {'rect', 'rectangular', 'none', 'hamming', 'hann', 'blackman', 'blackmanharris', 'taylor', 'chebyshev', 'kaiser'};
if ~any(strcmp(name, known))
    error(['%s%s must be ''rect'', ''hamming'', ''hann'', ''blackman'', ''blackmanharris'', ''taylor'', ' ...
           '''chebyshev'' or ''kaiser'' (is ''%s'').'], path, key, name);
end
end

function [iv, recEnd, isSet] = receive_windows(ws, wl, defaultEnd, halfSpan)
% Receive windows [start end) of the sub-pulses in the PRI (cfg.wave.<pulse>.rxStart_s / rxEnd_s,
% [] = from the end of the pulse to the next transmission or to the end of the record), whether
% each start / end is set (isSet.<pulse> = [start end]) and the end of the record: with rxEnd_s
% of the last transmitted sub-pulse set, halfSpan (half the input span of a DDC sample and one
% sample) after it, so that every DDC sample up to rxEnd_s is in the record; else defaultEnd.
% An end of inf means the end of the record. A switched-off sub-pulse has the window [0 0].
w = {ws, wl};
names = {'short', 'long'};
last = 2;  if ~wl.enable, last = 1; end
iv.short = [0 0];  iv.long = [0 0];
isSet.short = [false false];  isSet.long = [false false];
recEnd = defaultEnd;
for i = 1:2
    v = w{i};
    if ~v.enable, continue; end
    if ~isscalar(v.tTx_s) || ~isreal(v.tTx_s) || ~(v.tTx_s >= 0)
        error('cfg.wave.%s.tTx_s must be a time in the PRI (>= 0) [s].', names{i});
    end
    pulseEnd = v.tTx_s + v.width_s;
    a = pulseEnd;  b = inf;
    if i == 1 && wl.enable, b = wl.tTx_s; end
    if isfield(v, 'rxStart_s') && ~isempty(v.rxStart_s), a = v.rxStart_s;  isSet.(names{i})(1) = true; end
    if isfield(v, 'rxEnd_s') && ~isempty(v.rxEnd_s)
        b = v.rxEnd_s;  isSet.(names{i})(2) = true;
        if i == last, recEnd = b + halfSpan; end
    end
    if ~isscalar(a) || ~isreal(a) || ~isfinite(a) || a < pulseEnd - 1e-12
        error(['cfg.wave.%s.rxStart_s = %.2f us: the receive window starts at the end of the pulse ', ...
               '(tTx_s + width_s = %.2f us) or later.'], names{i}, a*1e6, pulseEnd*1e6);
    end
    if ~isscalar(b) || ~isreal(b) || isnan(b) || b <= a
        error('cfg.wave.%s.rxEnd_s must be later than the start of its receive window (%.2f us).', names{i}, a*1e6);
    end
    if i == 1 && wl.enable && b > wl.tTx_s + 1e-12
        error(['cfg.wave.short.rxEnd_s = %.2f us: the short-pulse receive window ends at the ', ...
               'transmission of the long pulse (%.2f us) or earlier.'], b*1e6, wl.tTx_s*1e6);
    end
    iv.(names{i}) = [a b];
end
end

function m = radar_mode_check(cfg)
% Checks cfg.modes (radar_config.m, apply_radar_mode.m): at most 8 modes, numbers 0 ... 7, the
% active one among them. Without cfg.modes (older configurations): mode 0.
m = struct('available', 0, 'active', 0, 'name', '');
if ~isfield(cfg, 'modes'), return; end
md = cfg.modes;
if isfield(md, 'available'), m.available = md.available; end
if isfield(md, 'active'), m.active = md.active; end
if isfield(md, 'name'), m.name = md.name; end
if isfield(md, 'subMode') || (isfield(cfg, 'wave') && isfield(cfg.wave, 'singlePulse'))
    error(['cfg.modes.subMode and cfg.wave.singlePulse no longer exist: a mode with one sub-pulse ' ...
           'per PRI switches the other one off (cfg.wave.<pulse>.enable = false in its mode file).']);
end
a = m.available;
if isempty(a) || ~isnumeric(a) || numel(a) > 8 || any(a(:) ~= round(a(:))) || any(a(:) < 0 | a(:) > 7) || ...
        numel(unique(a)) ~= numel(a)
    error('cfg.modes.available must hold 1 to 8 different mode numbers 0 ... 7.');
end
if ~isscalar(m.active) || ~any(a == m.active)
    error('cfg.modes.active = %s is not one of cfg.modes.available = [%s].', mat2str(m.active), num2str(a(:).'));
end
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
if strcmpi(wv.type, 'lfm')                                    % cfg.proc.rangeWindow (default Taylor)
    [wname, wp] = window_spec(cfg.proc, 'rangeWindow', 'taylor');
    win = window_function(wname, Lp, wp);
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

function B = waveform_bandwidth(wv)
% Bandwidth of a sub-pulse [Hz]: the LFM sweep, or 1/chip of a phase code.
if strcmpi(wv.type, 'phasecode'), B = 1 / wv.chip_s; else, B = wv.bandwidth_Hz; end
end

function K = mismatch_length(wv)
% cfg.wave.<pulse>.mismatchLength ([] = default length of make_input_files; also for a cfg without the field).
K = [];
if isfield(wv, 'mismatchLength'), K = wv.mismatchLength; end
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
% Detector, window and threshold settings of the selected CFAR method, as a range x Doppler window.
c = cfg.cfar;
w.method = lower(c.method);
w.notchWidth = c.notchWidth;
w.notchMaxRange_m = Inf;
if isfield(c, 'notchMaxRange_m') && ~isempty(c.notchMaxRange_m)
    w.notchMaxRange_m = c.notchMaxRange_m;
end
if ~isscalar(w.notchMaxRange_m) || ~isreal(w.notchMaxRange_m) || isnan(w.notchMaxRange_m) || w.notchMaxRange_m < 0
    error('cfg.cfar.notchMaxRange_m must be a non-negative scalar range in metres (Inf is allowed).');
end
w.type = 'ca';                                     % detector: cfg.cfar.type
if isfield(c, 'type') && ~isempty(c.type), w.type = lower(c.type); end
if ~any(strcmp(w.type, {'ca', 'go', 'so', 'os', 'ce'}))
    error('cfg.cfar.type must be ''ca'', ''go'', ''so'', ''os'' or ''ce''.');
end
w.osRank = 0.9;
if isfield(c, 'osRank') && ~isempty(c.osRank), w.osRank = c.osRank; end
if ~(w.osRank > 0 && w.osRank <= 1)
    error('cfg.cfar.osRank must lie in (0, 1].');
end
w.ceCensor = 0.1;                                  % CE: fraction of the training cells left out (the largest)
if isfield(c, 'ceCensor') && ~isempty(c.ceCensor), w.ceCensor = c.ceCensor; end
if ~isscalar(w.ceCensor) || ~isreal(w.ceCensor) || ~(w.ceCensor >= 0 && w.ceCensor < 1)
    error('cfg.cfar.ceCensor must lie in [0, 1).');
end
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
w.windowLabel = w.label;
w.label = [upper(w.type) ' ' w.label];             % e.g. 'CA 2-D', 'OS 1-D range'
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

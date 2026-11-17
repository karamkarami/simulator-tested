function ok = test_processing_options()
%TEST_PROCESSING_OPTIONS  Windows of range / Doppler / NUFFT, the NUFFT of the staggered PRF and the channels used.
%
%   1) window_function.m: the Chebyshev and Kaiser windows agree with the definitions of
%      MATLAB's chebwin / kaiser (reference values), every window has its peak sidelobe
%      level where expected (Chebyshev and Taylor at their sll_dB, Blackman-Harris below
%      -90 dB), the continuous evaluation (NUFFT) at the uniform positions equals the
%      uniform window, and an unknown name stops with a message (range, Doppler, NUFFT)
%   2) cfg.proc.rangeWindow (LFM): the range sidelobes of the compressed pulse follow the
%      window (Taylor -35 dB, Chebyshev -45 dB lower than with the rectangular window)
%   3) cfg.proc.staggerDoppler = 'nufft': for every stagger phase white noise at the MTI
%      input gives a flat noise floor of 1 (the scaling of the matrices), the bins equal the
%      direct non-uniform DFT of the windowed samples (doppler_process), a target at
%      +250 m/s (beyond the +-112 m/s span of the de-interleaved rows) peaks in its own
%      velocity bin, with the de-interleaving it is folded, and zero velocity is the notch
%      bin (none when the grid does not hold it)
%   4) cfg.rx.nChannelsUsed: 4 of 48 channels: a DBF matrix of 4 columns, unit noise gain
%      per beam; 1: DBF off (one beam); 0 or 49 stop with a message
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'), fullfile(root, 'tests'));
ok = true;

%% 1) windows
w = window_function('chebyshev', 15, struct('sll_dB', -50));
ref = [0.050612 0.131406 0.265439 0.442893];                 % chebwin(15, 50) (normalised to 1)
ok = check(max(abs(w(1:4).' - ref)) < 1e-5 && abs(max(w) - 1) < 1e-12, 'Chebyshev window = chebwin') && ok;
w = window_function('kaiser', 16, struct('beta', 6));
ok = check(abs(w(1) - 0.014873337104763) < 1e-12 && abs(w(8) - 0.987896888460194) < 1e-12, 'Kaiser window = kaiser') && ok;
names = {'rect', 'hamming', 'hann', 'blackman', 'blackmanharris', 'taylor', 'chebyshev', 'kaiser'};
expPsl = [-13.3 -42.7 -31.5 -58.1 -92 -35 -50 -44];          % peak sidelobe levels [dB] (N = 64)
p = struct('sll_dB', -50, 'nbar', 5, 'beta', 6);
for i = 1:numel(names)
    q = p;  if strcmp(names{i}, 'taylor'), q.sll_dB = -35; end
    [w, ml] = window_function(names{i}, 64, q);
    W = abs(fft(w, 64*64));  W = 20*log10(W / max(W) + 1e-300);
    k = find(diff(W(1:end/2)) > 0, 1);                       % first null: end of the mainlobe
    psl = max(W(k:end/2));
    x = window_function(names{i}, 64, q, ((0:63).' + 0.5*strcmp(names{i}, 'taylor')) / (64 - ~strcmp(names{i}, 'taylor')));
    good = abs(psl - expPsl(i)) < 1.5 + 2*strcmp(names{i}, 'kaiser') && max(abs(x - w)) < 1e-9 && ...
           abs(k / 64 - ml) < 0.6;
    ok = check(good, sprintf('window %-15s PSL %6.1f dB (expected %6.1f), mainlobe %.2f bins (first null %.2f)', ...
                             names{i}, psl, expPsl(i), ml, k / 64)) && ok;
end
c = radar_config();  c.proc.dopplerWindow = 'gauss';
ok = check(contains_text(error_of(@() setup_radar(c)), 'cfg.proc.dopplerWindow'), 'unknown Doppler window refused') && ok;
c = lfm_mode_config();  c.proc.rangeWindow = 'gauss';
ok = check(contains_text(error_of(@() setup_radar(c)), 'cfg.proc.rangeWindow'), 'unknown range window refused') && ok;
c = radar_config();  c.prf.mode = 'staggered';  c.proc.staggerDoppler = 'nufft';  c.proc.nufft.window = 'gauss';
ok = check(contains_text(error_of(@() setup_radar(c)), 'cfg.proc.nufft.window'), 'unknown NUFFT window refused') && ok;

%% 2) range window of the LFM matched filter
c = lfm_mode_config();
c.ant.nElements = 1;
psl = zeros(1, 3);
wins = {'rect', 'taylor', 'chebyshev'};
for i = 1:3
    c.proc.rangeWindow = wins{i};  c.proc.rangeWindowSll_dB = -35 - 10*(i == 3);
    s = setup_radar(c);
    wv = s.wave.long;
    y = abs(conv(waveform_envelope(wv, ((0:round(wv.width_s*s.fsOut)-1).' + 0.5) / s.fsOut), wv.mfTaps(:)));
    [m, k] = max(y);
    y(max(1, k-3):min(end, k+3)) = 0;                       % mainlobe (time-bandwidth 125: a few samples)
    psl(i) = 20*log10(max(y) / m);
end
ok = check(psl(2) < psl(1) - 10 && psl(3) < psl(2) - 3, sprintf(['range window (LFM, 2 MHz samples): range sidelobes ' ...
           'rect %.1f, Taylor %.1f, Chebyshev %.1f dB'], psl)) && ok;

%% 3) NUFFT of the staggered PRF
c = radar_config();
c.prf.mode = 'staggered';
c.proc.staggerDoppler = 'nufft';
c.ant.nElements = 1;
s = setup_radar(c);
K = s.dop.nfft;  n = s.cpi.length;  nS = s.nStagger;
flat = true;  direct = true;
rng(5);
for ph = 1:nS
    st = mod(ph - 1 + (0:n-1), nS) + 1;
    g = sum(doppler_process(mti_filter(reshape(eye(n), 1, n, n), s.mtiCoeffs, st), s, [], st), 2);
    flat = flat && max(abs(g(:) - 1)) < 1e-9;
    x = complex(randn(1, 3, s.cpi.rowLen), randn(1, 3, s.cpi.rowLen));
    [P, Y] = doppler_process(x, s, [], st);
    pri = reshape(s.pri(st), 1, []);  t = [0, cumsum(pri(1:end-1))];  t = t(1:s.cpi.rowLen);
    A = s.dop.nufft.A(:, :, ph);
    yd = A * reshape(x(1, 2, :), [], 1);
    direct = direct && max(abs(yd - reshape(Y(1, 2, :), [], 1))) < 1e-9 * max(abs(yd)) && ...
             max(abs(abs(yd).^2 - reshape(P(1, 2, :), [], 1))) < 1e-9 * max(abs(yd).^2);
    % the phase factors of A are those of the sample times
    direct = direct && abs(angle(A(1, 5) / A(1, 4)) - angle(exp(-2i*pi*s.dop.f_Hz(1)*(t(5) - t(4))))) < 1e-9;
end
ok = check(flat, sprintf('NUFFT: flat noise floor (1) in all %d bins for every stagger phase (%d)', K, nS)) && ok;
ok = check(direct, 'NUFFT: doppler_process = the direct non-uniform DFT of the windowed samples') && ok;
v0 = 250;                                                    % beyond the +-112 m/s of the rows
pri = reshape(s.pri(mod(0:n-1, nS) + 1), 1, []);  t = [0, cumsum(pri(1:end-1))];
x = reshape(exp(-4i*pi*v0*t / s.lambda), 1, 1, n);
y = mti_filter(x, s.mtiCoeffs, mod(0:n-1, nS) + 1);
P = squeeze(doppler_process(y, s, [], mod(0:n-1, nS) + 1));
[~, k] = max(P);
sD = c;  sD.proc.staggerDoppler = 'deinterleave';  sD = setup_radar(sD);
PD = squeeze(doppler_process(y, sD));
[~, kD] = max(PD);
dv = abs(s.dop.v_mps(2) - s.dop.v_mps(1));
ok = check(abs(s.dop.v_mps(k) - v0) <= dv/2 + 1e-9 && abs(sD.dop.v_mps(kD) - v0) > 50, ...
           sprintf(['NUFFT: a target at %g m/s peaks at %.1f m/s (bin %.1f m/s); de-interleaved rows: %.1f m/s ' ...
                    '(span +-%.0f m/s)'], v0, s.dop.v_mps(k), dv, sD.dop.v_mps(kD), max(abs(sD.dop.v_mps)))) && ok;
ok = check(s.dop.zeroBin == K/2 + 1 && abs(s.dop.v_mps(s.dop.zeroBin)) < 1e-9 && s.dop.nLooks == 1, ...
           sprintf('NUFFT: zero velocity in bin %d of %d, %d notch bins, one look', s.dop.zeroBin, K, s.dop.notchBins)) && ok;
c2 = c;  c2.proc.nufft.vMin_mps = 20;  c2.proc.nufft.vMax_mps = 400;  c2.proc.nufft.nPoints = 95;
s2 = setup_radar(c2);
ok = check(s2.dop.zeroBin == 0 && s2.dop.nfft == 95 && abs(max(s2.dop.v_mps) - 400) < 1e-9, ...
           'NUFFT: a grid without zero velocity has no notch (95 points, 20 ... 400 m/s)') && ok;

%% 4) channels used
c = radar_config();
c.rx.nChannelsUsed = 4;
s = setup_radar(c);
ok = check(s.dbf.enabled && isequal(size(s.dbf.W), [s.nBeams 4]) && max(abs(sum(abs(s.dbf.W).^2, 2) - 1)) < 1e-12 && ...
           s.nCh == 48, 'nChannelsUsed = 4: DBF of the first 4 elements, unit noise gain') && ok;
c.rx.nChannelsUsed = 1;
s = setup_radar(c);
ok = check(~s.dbf.enabled && s.nBeams == 1, 'nChannelsUsed = 1: no DBF, one beam') && ok;
for bad = [0 49]
    c.rx.nChannelsUsed = bad;
    ok = check(contains_text(error_of(@() setup_radar(c)), 'nChannelsUsed'), sprintf('nChannelsUsed = %d refused', bad)) && ok;
end
fprintf('Processing options (windows, NUFFT, channels used): %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function msg = error_of(f)
msg = '';
try
    f();
catch err
    msg = err.message;
end
end

function t = contains_text(s, p)
t = ~isempty(strfind(s, p));
end

function ok = check(cond, msg)
ok = logical(cond);
fprintf('  %-110s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

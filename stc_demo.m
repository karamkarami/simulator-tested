function out = stc_demo(cfg)
%STC_DEMO  The STC (sensitivity time control) of the receiver: why it is needed and what it does.
%
%   out = STC_DEMO() uses radar_config(); STC_DEMO(cfg) another configuration
%   (its receiver gain and STC curves). Five figures, also saved as PNG files
%   in <cfg.io.outFolder>/stc_demo:
%     1. stc_curve.png: the STC attenuation of the short- and of the
%        long-pulse window against the equivalent range c*dt/2 (dt = time
%        since the end of the transmitted pulse) and the attenuation the
%        clutter needs: its peak level at the ADC over 24 azimuths (no STC)
%        above full scale, plus 6 dB margin. The long-pulse window starts
%        with the eclipsed echoes of the targets nearer than its blind range
%        (22.5 km), so its curve goes on to there.
%     2. adc_pri.png: the ADC samples of one element over one PRI (clutter,
%        six targets, noise) without and with the STC.
%     3. clipping.png: the fraction of the receive samples at full scale for
%        24 azimuths: cfg.rx.gain_dB without STC, with STC, and 30 dB without
%        STC (the receiver gain that avoids clipping without STC).
%     4. snr_range.png: the expected SNR after the Doppler processing
%        (snr_budget.m) of a 1 m^2 target against range, with and without
%        the STC and with 30 dB without STC: the STC lowers the sensitivity
%        at the start of each receive window (short pulse up to ~12 km, long
%        pulse up to ~35 km), where the targets are strong anyway.
%     5. eclipse.png: one CPI on a 10 m^2 target at 5 km (one element): the
%        long-pulse power in its Doppler bin against range, without STC (at
%        30 dB) and with it, and the CFAR detections. The end of the long
%        pulse echoed by this target arrives at the start of the long-pulse
%        window; without STC its compressed response reaches the CFAR
%        threshold (false plots up to 27.5 km), the long-window STC
%        attenuates it below the noise.
%   out holds the numbers of the figures (equivalent range, needed and
%   applied attenuation, clipped fractions, SNR) and out.suggested, an STC
%   curve [range_m atten_dB] derived from the clutter of cfg (for another
%   clutter, gain or carrier).

root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'));
if nargin < 1, cfg = radar_config(); end
cfg.console.enable = false;
cfg.display.enable = false;
folder = fullfile(cfg.io.outFolder, 'stc_demo');
if ~exist(folder, 'dir'), mkdir(folder); end
visible = 'on';
if ~usejava('awt') && ~exist('OCTAVE_VERSION', 'builtin'), visible = 'off'; end
if exist('OCTAVE_VERSION', 'builtin') && isempty(getenv('DISPLAY')), visible = 'off'; end
c = cfg.const.c;
if isempty(cfg.rx.stc.curve), error('cfg.rx.stc.curve is empty.'); end

%% 1) attenuation the clutter needs (peak level without STC, 1 element: the same clutter per element)
c1 = cfg;  c1.ant.nElements = 1;  c1.rx.stc.enable = false;
c1.rx.addThermalNoise = false;  c1.rx.addLeakage = false;  c1.adc.addAdcNoise = false;  c1.adc.quantize = false;
sys1 = setup_radar(c1);
rng(1);
clut = clutter_init(sys1);
pk = zeros(sys1.nAdc, 1);
for az = 0:15:345
    x = adc_sample(generate_rx_signal(sys1, [], 0, az, clutter_signal(clut, sys1, az, 0)), sys1);
    pk = max(pk, abs(double(x(:, 1))));
end
ws = sys1.wave.short;  wl = sys1.wave.long;
Re = (0:100:25e3).';                                          % equivalent range grid [m]
need = zeros(numel(Re), 1);
for w = {ws, wl}
    t0 = w{1}.tTx_s + w{1}.width_s;                           % start of the receive window
    for k = 1:numel(Re)
        n = round((t0 + 2*Re(k)/c) * sys1.fs) + 1 + (-60:60); % +-0.5 us
        n = n(n >= 1 & n <= sys1.nAdc);
        need(k) = max(need(k), 20*log10(max(pk(n)) / sys1.adc.fullScale) + 6);
    end
end
need = max(need, 0);
sysS = setup_radar(cfg);                                      % the configured STC
applied = -20*log10(stc_gain(sysS, ws.tTx_s + ws.width_s + 2*Re/c));
appliedLong = -20*log10(stc_gain(sysS, wl.tTx_s + wl.width_s + 2*Re/c));
out.range_m = Re;  out.needed_dB = need;  out.stc_dB = applied;  out.stcLong_dB = appliedLong;
% a suggested curve: the needed attenuation on a 1 km grid (upper envelope), 0 dB beyond the clutter
rGrid = (0:1000:ceil(Re(find(need > 0, 1, 'last'))/1000)*1000).';
sug = arrayfun(@(r) max([need(Re >= r - 1000 & Re <= r + 1000); 0]), rGrid);
out.suggested = [rGrid, ceil(sug); rGrid(end) + 500, 0];

f1 = figure('Visible', visible, 'Name', 'STC curve');
plot(Re/1e3, need, 'r-', Re/1e3, applied, 'b-', Re/1e3, appliedLong, 'g--', 'LineWidth', 1.5);
grid on;  xlabel('equivalent range c\Deltat/2 after the end of the pulse [km]');  ylabel('attenuation [dB]');
title(sprintf('STC: needed by the clutter at %g dB gain (peak + 6 dB) and applied', cfg.rx.gain_dB));
legend('needed by the clutter (peak above full scale + 6 dB)', 'short-pulse window (cfg.rx.stc.curve)', ...
       'long-pulse window (cfg.rx.stc.curveLong)', 'Location', 'northeast');
print(f1, fullfile(folder, 'stc_curve.png'), '-dpng');

%% 2) ADC samples of one PRI, without and with the STC
c2 = cfg;  c2.ant.nElements = 1;
tg = draw_target_rcs(cfg.targets, 2, 0);
x = cell(1, 2);
for k = 1:2
    c2.rx.stc.enable = k == 2;
    s2 = setup_radar(c2);
    rng(2);
    cl = clutter_init(s2);
    x{k} = double(adc_sample(generate_rx_signal(s2, tg, 0, 30, clutter_signal(cl, s2, 30, 0)), s2));
end
f2 = figure('Visible', visible, 'Name', 'ADC samples');
labels = {'without STC', 'with STC'};
for k = 1:2
    subplot(2, 1, k);
    plot(s2.tAdc*1e6, x{k}, 'k-');  hold on;
    plot(s2.tAdc([1 end])*1e6, s2.adc.fullScale*[1 1], 'r--', s2.tAdc([1 end])*1e6, -s2.adc.fullScale*[1 1], 'r--');
    grid on;  ylim(1.1*s2.adc.fullScale*[-1 1]);  ylabel('ADC [LSB]');
    title(sprintf('one element, azimuth 30 deg, %g dB gain, %s: %.2f%% of the receive samples at full scale', ...
                  cfg.rx.gain_dB, labels{k}, 100*adc_overload(int16(x{k}), s2)));
end
xlabel('time in the PRI [\mus]');
print(f2, fullfile(folder, 'adc_pri.png'), '-dpng');

%% 3) clipped receive samples for 24 azimuths
cases = {cfg.rx.gain_dB, false; cfg.rx.gain_dB, true; 30, false};
names = cell(1, 3);
az = 0:15:345;
clip = zeros(3, numel(az));
for k = 1:3
    c3 = cfg;  c3.ant.nElements = 1;  c3.rx.gain_dB = cases{k, 1};  c3.rx.stc.enable = cases{k, 2};
    s3 = setup_radar(c3);
    rng(1);
    cl = clutter_init(s3);
    for j = 1:numel(az)
        clip(k, j) = adc_overload(adc_sample(generate_rx_signal(s3, tg, 0, az(j), ...
                                  clutter_signal(cl, s3, az(j), 0)), s3), s3);
    end
    names{k} = sprintf('%g dB, STC %s: mean %.3f%%', cases{k, 1}, ternary(cases{k, 2}, 'on', 'off'), 100*mean(clip(k, :)));
end
out.clipped = clip;  out.clipCases = cases;
f3 = figure('Visible', visible, 'Name', 'Clipping');
semilogy(az, max(100*clip, 1e-4).', 'o-', 'LineWidth', 1.2);
grid on;  xlabel('antenna azimuth [deg]');  ylabel('receive samples at full scale [%] (1e-4 = none)');
title('ADC clipping in the receive windows');  legend(names, 'Location', 'east');
print(f3, fullfile(folder, 'clipping.png'), '-dpng');

%% 4) expected SNR of a 1 m^2 target against range
R = (4.5:0.5:80) * 1e3;
snr = zeros(3, numel(R));
for k = 1:3
    c4 = cfg;  c4.rx.gain_dB = cases{k, 1};  c4.rx.stc.enable = cases{k, 2};
    s4 = setup_radar(c4);
    for i = 1:numel(R)
        b = snr_budget(s4, make_target('S', 1, R(i), 0, 1500, 200, 0), 0);
        snr(k, i) = b.doppler;
    end
end
out.snrRange_m = R;  out.snr_dB = snr;
f4 = figure('Visible', visible, 'Name', 'SNR');
plot(R/1e3, snr.', 'LineWidth', 1.5);
grid on;  xlabel('range [km]');  ylabel('SNR after Doppler processing [dB]');
title('1 m^2 target at 1500 m height, 200 m/s');
legend(cellfun(@(n) n(1:strfind(n, ':') - 1), names, 'UniformOutput', false), 'Location', 'northeast');
print(f4, fullfile(folder, 'snr_range.png'), '-dpng');

%% 5) the eclipsed long-pulse echo of a near target, without and with the STC
Rn = 5e3;                                                       % range of the near target [m]
tg5 = draw_target_rcs(make_target('N', 10, Rn, 40, 1000, -150, 0), 2, 0);
f5 = figure('Visible', visible, 'Name', 'Eclipsed echo');
styles = {'r', 'b'};
hl = [];  leg = {};
for k = 1:2
    c5 = cfg;  c5.ant.nElements = 1;  c5.rx.stc.enable = k == 2;
    if k == 1, c5.rx.gain_dB = 30; end                      % no STC: the gain that does not clip
    s5 = setup_radar(c5);
    [prof, det, rng_m, nGhost] = eclipse_profile(c5, s5, tg5, 40, Rn);
    hl(end+1) = plot(rng_m/1e3, prof, [styles{k} '-'], 'LineWidth', 1.2);  hold on; %#ok<AGROW>
    x = rng_m(det)/1e3;  y = prof(det);
    if isempty(x), x = NaN;  y = NaN; end                   % keeps the legend entry
    hl(end+1) = plot(x, y, [styles{k} 'o'], 'MarkerSize', 7, 'LineWidth', 1.5); %#ok<AGROW>
    leg = [leg, {sprintf('%s (%g dB gain)', ternary(k == 2, 'STC', 'no STC'), c5.rx.gain_dB), ...
                 sprintf('CFAR detections (%d)', nGhost)}]; %#ok<AGROW>
    out.eclipseGhosts(k) = nGhost;
end
grid on;  xlabel('long-pulse range [km]');  ylabel('power / noise floor [dB]');  xlim([22 (Rn + 30e3)/1e3]);
title(sprintf('10 m^2 at %g km: its eclipsed long-pulse echo (one element)', Rn/1e3));
legend(hl, leg, 'Location', 'northeast');
print(f5, fullfile(folder, 'eclipse.png'), '-dpng');

fprintf(['STC demo (%g dB gain): clipped receive samples %s; ghosts of the eclipsed echo of 10 m^2 at %g km: ', ...
         '%d without, %d with the STC; figures in %s\n'], cfg.rx.gain_dB, strjoin(names, ' | '), Rn/1e3, ...
        out.eclipseGhosts, folder);
fprintf('suggested STC curve for this clutter [range_m atten_dB]:\n');
fprintf('   %6.0f %4.0f\n', out.suggested.');
end

function [prof, det, rng_m, nGhost] = eclipse_profile(cfg, sys, tg, az0, Rn)
% One CPI centred on azimuth az0, long pulse: power in the Doppler bin of the
% target (at range Rn) over range (dB above the median of the map), the CFAR
% detections in that bin (+-1) and their number up to Rn + the blind range.
rng(11);
clut = clutter_init(sys);
n = sys.cpi.length;
t0 = (az0 - cfg.ant.az0_deg) / sys.ant.rate_dps;
wv = sys.wave.long;
buf = cpi_buffer_init(sys.nBeams, wv.nRange, n);
for k = 1:n
    t = t0 + (k - 1 - n/2) * sys.pri(1);
    az = mod(cfg.ant.az0_deg + sys.ant.rate_dps * t, 360);
    [I, Q] = ddc(adc_sample(generate_rx_signal(sys, tg, t, az, clutter_signal(clut, sys, az, t)), sys), sys);
    seg = ddc_segment(complex(single(I), single(Q)), wv);
    buf = cpi_buffer_push(buf, matched_filter(seg.', wv), 1, az, t);
end
[cpi, info] = cpi_buffer_read(buf);
rd = doppler_process(mti_filter(cpi, sys.mtiCoeffs, info.staggerIdx), sys);
isDet = cfar_detect(rd, sys, wv);
[~, ~, ~, vr] = target_geometry(tg, t0, cfg.ant.height_m);
[~, d] = min(abs(sys.dop.v_mps - vr));                    % Doppler bin of the target
floorP = median(rd(:));
prof = 10*log10(double(squeeze(rd(1, :, d))) / double(floorP));
det = squeeze(any(isDet(1, :, max(1, d-1):min(end, d+1)), 3)).';
rng_m = wv.range_m;
nGhost = nnz(det(:) & rng_m(:) <= Rn + cfg.const.c * wv.width_s / 2 + 1e3);
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

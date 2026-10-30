function results = run_radar_sim(cfg)
%RUN_RADAR_SIM  Simulation of a rotating S-band search radar with elevation DBF.
%
%   results = RUN_RADAR_SIM()    runs the simulation with radar_config().
%   results = RUN_RADAR_SIM(cfg) runs it with a modified configuration.
%
%   Processing chain of every PRI (a short and then a long sub-pulse):
%     1. Signal generation : radar equation per target (Swerling 1), clutter, Tx
%                            leakage and thermal noise on the 30 MHz IF, receiver
%                            gain, 48 elements, fs = 120 MHz, 16-bit ADC (12.7 ENOB)
%     2. Recording         : raw ADC data of every PRI -> <rawFolder>/pulse_<n>.mat
%                            (optional, cfg.io.saveRaw)
%     3. DDC               : I/Q mixer (cos / -sin 30 MHz), two fixed-point filter
%                            stages with decimation 10 and 6 (valid convolutions),
%                            signed 16-bit I/Q output (optional recording to
%                            <ddcFolder>/ddc_<n>.mat, cfg.io.saveDdc, and UDP streaming)
%     4. DBF               : 48 elements -> 21 elevation beams (off for 1 element)
%     5. Matched filter    : per beam (valid convolution, Taylor window for LFM)
%     6. CPI buffers       : circular buffers holding the pulses of one CPI (the
%                            time a target stays in the 3 dB azimuth beam)
%     7. MTI + Doppler     : 3-pulse canceller + FFT over the CPI (fixed PRF), or
%                            stagger MTI + de-interleaving + FFT + NCI (staggered PRF)
%     8. CFAR              : on every range-Doppler map (2-D, 1-D range or 1-D
%                            Doppler window, cfg.cfar.method)
%     9. Pruning           : local peaks only, then azimuth peak picking (plots)
%    10. Tracking          : EKF on the plots, 2-D or 3-D, constant velocity or
%                            constant acceleration, GNN association (optional, cfg.track)
%    11. PPI display       : updated while the radar is running: plots and/or the
%                            CFAR detections before the extractor, tracks, mean
%                            processing time of every block (also of the display
%                            and of the console output); the CFAR threshold can
%                            be changed in the PPI panel
%
%   results.plots holds the extracted plots (time, azimuth, range, elevation,
%   ambiguous radial velocity, SNR, pulse type). results.tracker holds the
%   tracks at the end and the log of all track updates (track_update.m).
%   results.timing_us is the mean time per PRI of DBF, matched filter, MTI,
%   Doppler, detection (CFAR + pruning), plot extractor, tracking, console
%   output (progress lines) and display (PPI updates) [us].

if nargin < 1, cfg = radar_config(); end
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'src'));
rng(cfg.sim.seed);
tStart = tic;

%% Setup
sys   = setup_radar(cfg);
sched = pulse_schedule(sys);
targets = draw_target_rcs(cfg.targets, max(sched.scan) + 1, cfg.sim.swerling);   % RCS of every antenna pass
names = {'short', 'long'};
if cfg.sim.verbose
    print_summary(sys, sched);
    results.budget = print_snr_budget(sys, cfg.targets);
end

clut = [];
if cfg.clutter.enable
    clut = clutter_init(sys);
end
if cfg.io.saveRaw && ~exist(cfg.io.rawFolder, 'dir')
    mkdir(cfg.io.rawFolder);
end
if cfg.io.saveDdc && ~exist(cfg.io.ddcFolder, 'dir')
    mkdir(cfg.io.ddcFolder);
end
if ~exist(cfg.io.outFolder, 'dir')
    mkdir(cfg.io.outFolder);
end
udp = [];
if cfg.udp.enable
    udp = udp_sender_open(cfg.udp);
    closeUdp = onCleanup(@() udp_sender_close(udp));   %#ok<NASGU> closes it at the end, also after an error
end
ppi = [];
detLog = [];
showPpi = cfg.display.enable;
if showPpi && exist('OCTAVE_VERSION', 'builtin') && isempty(available_graphics_toolkits())
    warning('run_radar_sim: no graphics toolkit (e.g. octave-cli): the PPI display is switched off.');
    showPpi = false;
end
if showPpi
    ppi = ppi_init(sys);
    detLog = detection_log_init(cfg.display.maxDetections);   % CFAR detections for the PPI
end

for i = 1:2
    p = names{i};
    buf.(p) = cpi_buffer_init(sys.nBeams, sys.wave.(p).nRange, sys.cpi.length);
    ext.(p) = plot_extractor_init(sys.dop.nfft);
end
plots = [];
tracker = track_init(sys.track);
blockUs = zeros(1, 9);                     % accumulated time of the blocks, console, display [us]
nProcessed = 0;
wasActive = false;
% Plot extraction options: minimum azimuth extent -> minimum number of CPI detections
extractOpt = cfg.extract;
extractOpt.minHits = max(1, ceil(cfg.extract.minExtent_deg / ...
                     (sys.ant.rate_dps * mean(sys.pri) * cfg.proc.updateEvery)));
extractOpt.flush = false;
flushOpt = extractOpt;
flushOpt.flush = true;

%% Pulse loop
for k = 1:sched.nPulses
    t  = sched.time_s(k);
    az = sched.az_deg(k);

    % Skipped sector: restart the CPI buffers and close the open plots
    if ~sched.active(k)
        if wasActive
            gapPlots = [];
            for i = 1:2
                p = names{i};
                buf.(p) = cpi_buffer_init(sys.nBeams, sys.wave.(p).nRange, sys.cpi.length);
                [ext.(p), newPlots] = plot_extractor(ext.(p), [], flushOpt);
                plots = [plots, newPlots]; %#ok<AGROW>
                gapPlots = [gapPlots, newPlots]; %#ok<AGROW>
            end
            if sys.track.enable                     % at the last PRI before the gap (as the C receiver)
                tracker = track_update(tracker, gapPlots, sched.time_s(k-1), true);
            end
        end
        wasActive = false;
        continue
    end
    wasActive = true;

    % 1) Received signal of all elements and ADC
    clutterIf = [];
    if ~isempty(clut)
        clutterIf = clutter_signal(clut, sys, az, t);
    end
    adcData = adc_sample(generate_rx_signal(sys, targets, t, az, clutterIf), sys);

    % 2) Raw data recording
    hdr = pulse_header(sys, sched, k);
    if cfg.io.saveRaw
        save_pulse_file(cfg.io.rawFolder, cfg.io.matFormat, hdr, adcData);
    end

    % 3) DDC of the whole record (both sub-pulses), recording and streaming of its output
    [I, Q] = ddc(adcData, sys);                         % signed 16 bit [nDdc x nElements]
    if cfg.io.saveDdc
        save_ddc_file(cfg.io.ddcFolder, cfg.io.matFormat, hdr, I, Q);
    end
    if ~isempty(udp)
        udp_send_pulse(udp, hdr, I, Q);
    end
    bb = complex(single(I), single(Q));

    % 4-9) Processing chain, separately for the short and the long sub-pulse
    priPlots = [];
    for i = 1:2
        p  = names{i};
        wv = sys.wave.(p);
        seg = bb(wv.ddcIdx, :);                         % DDC samples of this sub-pulse
        tb = tic;
        if sys.dbf.enabled
            bm = digital_beamforming(seg, sys.dbf.W);   % [nBeams x nSamples]
        else
            bm = seg.';                                 % single element: no beamforming
        end
        blockUs(1) = blockUs(1) + 1e6*toc(tb);
        tb = tic;
        rp = matched_filter(bm, wv);                    % [nBeams x nRange]
        blockUs(2) = blockUs(2) + 1e6*toc(tb);
        if k == cfg.display.rangeProfilePulse
            plot_range_profile(rp, wv, cfg.display.rangeProfileAxis);
        end

        tb = tic;                                       % CPI buffer and clutter filter
        buf.(p) = cpi_buffer_push(buf.(p), rp, sched.staggerIdx(k), az, t);
        if buf.(p).count == sys.cpi.length && mod(buf.(p).nPushed, cfg.proc.updateEvery) == 0
            [cpi, info] = cpi_buffer_read(buf.(p));
            y  = mti_filter(cpi, sys.mtiCoeffs, info.staggerIdx);
            blockUs(3) = blockUs(3) + 1e6*toc(tb);
            tb = tic;
            rd = doppler_process(y, sys);               % [nBeams x nRange x nfft]
            blockUs(4) = blockUs(4) + 1e6*toc(tb);
            tb = tic;
            [isDet, noise] = cfar_detect(rd, sys, wv);
            dets = prune_peaks(isDet, rd, noise, sys, wv);
            reports = detection_reports(dets, rd, sys, wv, info);
            blockUs(5) = blockUs(5) + 1e6*toc(tb);
            if ~isempty(detLog)
                detLog = detection_log_add(detLog, reports);
            end
            tb = tic;
            [ext.(p), newPlots] = plot_extractor(ext.(p), reports, extractOpt);
            blockUs(6) = blockUs(6) + 1e6*toc(tb);
            plots = [plots, newPlots]; %#ok<AGROW>
            priPlots = [priPlots, newPlots]; %#ok<AGROW>
        else
            blockUs(3) = blockUs(3) + 1e6*toc(tb);
        end
    end
    if sys.track.enable                                 % 10) tracking (also deletes old tracks)
        tb = tic;
        tracker = track_update(tracker, priPlots, t);
        blockUs(7) = blockUs(7) + 1e6*toc(tb);
    end
    nProcessed = nProcessed + 1;

    % 11) Display
    if ~isempty(ppi) && mod(nProcessed, cfg.display.updateEvery) == 0
        tb = tic;
        sys = ppi_update(ppi, sys, az, t, sched.scan(k), plots, detLog.M(1:detLog.n, :), ...
                         track_report(tracker), blockUs / nProcessed);   % may set a new CFAR threshold
        blockUs(9) = blockUs(9) + 1e6*toc(tb);
    end
    if cfg.sim.verbose && mod(nProcessed, 500) == 0
        tb = tic;
        fprintf('  pulse %6d (az %6.1f deg): %d plots so far, %.0f s elapsed\n', k, az, numel(plots), toc(tStart));
        blockUs(8) = blockUs(8) + 1e6*toc(tb);
    end
end

%% Close the remaining plots, associate the last GNN batch
endPlots = [];
for i = 1:2
    p = names{i};
    [ext.(p), newPlots] = plot_extractor(ext.(p), [], flushOpt);
    plots = [plots, newPlots]; %#ok<AGROW>
    endPlots = [endPlots, newPlots]; %#ok<AGROW>
end
lastActive = find(sched.active, 1, 'last');
if sys.track.enable && ~isempty(lastActive)
    tracker = track_update(tracker, endPlots, sched.time_s(lastActive), true);
end
if ~isempty(ppi) && ishandle(ppi.fig)               % (the user may have closed the window)
    ppi_update(ppi, sys, sched.az_deg(end), sched.time_s(end), sched.scan(end), plots, detLog.M(1:detLog.n, :), ...
               track_report(tracker), blockUs / max(nProcessed, 1));
    print(ppi.fig, fullfile(cfg.io.outFolder, 'ppi.png'), '-dpng');
end

results.plots      = plots;
results.tracker    = tracker;
results.timing_us  = blockUs / max(nProcessed, 1);
results.sys        = sys;
results.targets    = targets;                       % incl. the RCS of every scan (rcsScan)
results.nProcessed = nProcessed;
results.runTime_s  = toc(tStart);
if cfg.sim.verbose
    print_plots(plots);
    if sys.track.enable
        print_tracks(tracker);
    end
    fprintf(['Mean processing time per PRI [us]: DBF %.0f, matched filter %.0f, MTI %.0f, Doppler %.0f, ' ...
             'CFAR %s + pruning %.0f, extractor %.0f, tracking %.0f, console %.0f, display %.0f\n'], ...
            results.timing_us(1:4), sys.cfar.label, results.timing_us(5:9));
    fprintf('%d PRIs simulated and processed in %.1f s (%.1f ms per PRI)\n', nProcessed, ...
        results.runTime_s, 1e3*results.runTime_s/max(nProcessed, 1));
end
end

function print_summary(sys, sched)
% Short description of the simulated radar.
cfg = sys.cfg;
fprintf('\n=== Rotating search radar simulation ===\n');
if sys.dbf.enabled, dbfText = sprintf('%d beams', sys.nBeams); else, dbfText = 'no beamforming'; end
fprintf('fc = %.2f GHz, Pt = %.0f kW, NF = %.1f dB, %d element(s) -> %s\n', cfg.fc/1e9, ...
    cfg.tx.peakPower_W/1e3, cfg.rx.noiseFigure_dB, sys.nCh, dbfText);
fprintf('PRF mode: %s, mean PRI = %.0f us, %d rpm, azimuth beamwidth %.2f deg\n', cfg.prf.mode, ...
    mean(sys.pri)*1e6, cfg.ant.rpm, sys.ant.azBeamwidth_deg);
fprintf('ADC record: %.1f ... %.1f us of the PRI, IF %.0f MHz (%d samples, %.1f MB per PRI)\n', ...
    sys.rxStart_s*1e6, sys.tAdc(end)*1e6, sys.fIF/1e6, sys.nAdc, sys.nAdc*sys.nCh*2/2^20);
noiseDdc = sys.adc.sigmaTh^2 * sum(sys.ddc.hEq.^2) * prod(sys.ddc.gain)^2;   % per complex sample
fprintf(['Receiver gain %.0f dB: thermal noise at the ADC %.0f LSB rms (%.1f dBFS), ADC noise %.1f LSB, ', ...
         'Tx leakage %.0f dBm at the ADC, DDC output noise %.1f LSB rms (I or Q)\n'], cfg.rx.gain_dB, ...
    sys.adc.sigmaTh, 20*log10(sys.adc.sigmaTh*sqrt(2)/sys.adc.fullScale), sys.adc.sigmaAdc, ...
    10*log10(sys.rx.leakAmp^2/2 * sys.rx.voltGain^2 * 1e3), sqrt(noiseDdc/2));
names = {'short', 'long'};
for i = 1:2
    wv = sys.wave.(names{i});
    fprintf('%-5s pulse: %s %3.0f us at t = %3.0f us, range %.2f ... %.2f km (%d cells)\n', wv.name, wv.type, ...
        wv.width_s*1e6, wv.tTx_s*1e6, wv.range_m([1 end])/1e3, wv.nRange);
end
fprintf('CPI: %.1f ms in the beam -> %d pulses = %d row(s) x %d samples, FFT %d, Doppler bin %.2f m/s, unambiguous +-%.1f m/s\n', ...
    sys.cpi.time_s*1e3, sys.cpi.length, sys.cpi.nRows, sys.cpi.rowLen, sys.dop.nfft, abs(diff(sys.dop.v_mps(1:2))), ...
    sys.lambda/(4*sys.dop.rowPri_s));
fprintf('CFAR %s: Pfa = %.0e, %d training cells, %.1f looks, alpha = %.2f (%.1f dB)\n', sys.cfar.label, ...
    sys.cfar.pfa, sys.cfar.nTrain, sys.dop.nLooks, sys.cfar.alpha, 10*log10(sys.cfar.alpha));
fprintf('Simulated pulses: %d of %d (%.1f s of radar time)\n', sum(sched.active), sched.nPulses, sched.time_s(end));
end

function print_tracks(tracker)
% List the tracks at the end of the run.
tk = tracker.tracks;
fprintf('\nTracks at the end (%d-D EKF, %s model, GNN): %d (%d confirmed), %d track updates\n', tracker.p.dim, ...
        upper(tracker.p.model), numel(tk), ...
        sum([tk.confirmed]), size(tracker.log, 1));
fprintf('%5s %4s %9s %8s %8s %8s %8s %5s\n', 'id', 'conf', 'R[km]', 'az[deg]', 'h[m]', 'v[m/s]', 'hdg[deg]', 'hits');
for i = 1:numel(tk)
    x = tk(i).x;  d = tracker.p.dim;
    z = 0;  if d == 3, z = x(3); end
    v = x(d+1:d+2);
    fprintf('%5d %4d %9.3f %8.2f %8.0f %8.1f %8.1f %5d\n', tk(i).id, tk(i).confirmed, norm(x(1:d))/1e3, ...
            mod(atan2(x(1), x(2))*180/pi, 360), z + tracker.p.height_m, norm(v), mod(atan2(v(1), v(2))*180/pi, 360), tk(i).hits);
end
end

function print_plots(plots)
% List the extracted plots.
fprintf('\nExtracted plots: %d\n', numel(plots));
fprintf('%8s %8s %9s %7s %8s %7s %6s %5s\n', 'time[s]', 'az[deg]', 'R[km]', 'el[deg]', 'vr[m/s]', ...
    'SNR[dB]', 'pulse', 'hits');
for i = 1:numel(plots)
    p = plots(i);
    fprintf('%8.3f %8.2f %9.3f %7.2f %8.2f %7.1f %6s %5d\n', p.time_s, p.az_deg, p.range_m/1e3, ...
        p.el_deg, p.vr_mps, p.snr_dB, p.pulse, p.hits);
end
end

function ok = test_slc(quick)
%TEST_SLC  Sidelobe canceller (cfg.slc), DBF / SLC elements and noise jammers (cfg.jammer).
%
%   ok = TEST_SLC() runs every part; TEST_SLC(true) skips the long simulation of part 6.
%
%   1) sidelobe_canceller.m on synthetic data (4 auxiliary elements, 3 beams): a noise jammer
%      30 dB above the noise is cancelled by QRD-SMI and QRD-RLS to within 1 dB of the noise
%      (NLMS: 2 dB, its misadjustment with the step 0.1), a target outside the training samples keeps its amplitude (0.5 dB), without the
%      jammer the SLC adds less than 0.5 dB of noise; QRD-SMI (qr) equals the QRD-RLS Givens
%      rotations with lambda = 1 and the same loading; the weights of NLMS and QRD-RLS go on
%      from PRI to PRI (a second PRI cancels better than or as well as the first); two jammers
%      from different directions are cancelled with 4 auxiliary elements
%   2) the jammer model of generate_rx_signal.m (no clutter, leakage, STC): its power at the
%      DDC output over the noise of an element equals sys.jammer.jnr_dB with the antenna toward
%      it (0.5 dB), it follows the azimuth pattern of the element (10 deg off: 1 dB), and its
%      phase from element to element is that of its direction (0.02 rad)
%   3) the C generator radar_generate_ddc (when built): the same jammer power (0.5 dB) and phase
%      step (0.02 rad) as MATLAB; radar_live_sim (fidelity 'full' and 'fast', over UDP into the
%      C receiver, which logs the DDC output): the same power (1 dB) and phase step
%   4) the C receiver (when built) on the DDC files of run_radar_sim with the jammer: its logged
%      DBF and SLC outputs equal those of MATLAB (DBF 1e-5, SLC 1e-3 relative) for QRD-SMI,
%      QRD-RLS and NLMS, DBF of the elements 1 ... 44, SLC with 45 ... 48; its PPI modes: the
%      profile of the SLC output exists only with the SLC
%   5) DBF of the elements [3 10 17 40] and SLC with [1 48] in run_radar_sim (the selection by
%      element number), the single-element main channel with an SLC (no DBF), the 'slc' display
%      of the PPI only with the SLC
%   6) (not quick) a weak target behind the elevation sidelobe of a jammer: detected without the
%      jammer, missed with it, detected again with the SLC (each algorithm) with an SNR within
%      6 dB of the case without the jammer
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'), fullfile(root, 'tests'));
if nargin < 1, quick = false; end
ok = true;
build = fullfile(receiver_dir(), 'build');
exeSfx = '';  if ispc, exeSfx = '.exe'; end
out = fullfile(root, 'output', 'test_slc');
if ~exist(out, 'dir'), mkdir(out); end

%% 1) sidelobe_canceller.m on synthetic data
randn('state', 7);  rand('state', 7); %#ok<RAND>
nA = 4;  nB = 3;  n = 3000;  K = 256;
vA = exp(1j * 0.9 * (0:nA-1)).';                              % jammer on the auxiliary elements
gB = [0.3; 0.05 * exp(1j); 0.6 * exp(-0.4j)] * 2;              % ... and in the beams (sidelobes)
noise = @(r, c) (randn(r, c) + 1j*randn(r, c)) / sqrt(2);
jam = 10^(30/20) * noise(1, n);                                % 30 dB above the noise
tgt = zeros(1, n);  tgt(500:520) = 20;                         % a target (only in the beams)
trainIdx = n - K + 1:n;
names = {'qrd-smi', 'qrd-rls', 'nlms'};
for m = 0:2
    slc = slc_settings(m);
    st = [];
    for pri = 1:2
        aux = vA * jam + noise(nA, n);
        main = gB * jam + noise(nB, n) + [tgt; tgt; tgt];
        [y, st] = sidelobe_canceller(st, main, aux, main(:, trainIdx), aux(:, trainIdx), slc);
        res(pri) = mean(mean(abs(y(:, 1000:n-K)).^2, 2)); %#ok<AGROW>
    end
    tg = mean(abs(y(:, 505)).^2 ./ abs(main(:, 505) - gB * jam(505)).^2);
    lim = 1 + (m == 2);                                        % NLMS: its misadjustment (step 0.1)
    ok = check(10*log10(res(2)) < lim && res(2) <= res(1) * 1.2 && abs(10*log10(tg)) < 0.5, ...
               sprintf('%s: jammer 30 dB -> residual %.2f dB above the noise (first PRI %.2f dB), target %.2f dB', ...
                       upper(names{m+1}), 10*log10(res(2)), 10*log10(res(1)), 10*log10(tg))) && ok;
    st = [];  aux = noise(nA, n);  main = noise(nB, n);
    [y, st] = sidelobe_canceller(st, main, aux, main(:, trainIdx), aux(:, trainIdx), slc);
    [y, ~] = sidelobe_canceller(st, main, aux, main(:, trainIdx), aux(:, trainIdx), slc);
    inc = mean(abs(y(:)).^2) / mean(abs(main(:)).^2);
    ok = check(10*log10(inc) < 0.5, sprintf('%s without jammer: %.2f dB more noise', upper(names{m+1}), 10*log10(inc))) && ok;
end
aux = vA * jam + noise(nA, n);  main = gB * jam + noise(nB, n);
s0 = slc_settings(0);  s1 = slc_settings(1);  s1.rlsLambda = 1;  s1.rlsDelta_dB = s0.smiLoading_dB;
[y0, a] = sidelobe_canceller([], main, aux, main(:, trainIdx), aux(:, trainIdx), s0);
[y1, b] = sidelobe_canceller([], main, aux, main(:, trainIdx), aux(:, trainIdx), s1);
ok = check(max(abs(a.W(:) - b.W(:))) / max(abs(a.W(:))) < 1e-9 && max(abs(y0(:) - y1(:))) < 1e-6 * max(abs(y0(:))), ...
           'QRD-SMI (qr) = QRD-RLS Givens rotations with lambda 1 and the same loading') && ok;
vA2 = exp(-1j * 2.1 * (0:nA-1)).';  jam2 = 10^(25/20) * noise(1, n);
aux = vA * jam + vA2 * jam2 + noise(nA, n);  main = gB * jam + [0.2; 0.4; 0.1] * jam2 + noise(nB, n);
y = sidelobe_canceller([], main, aux, main(:, trainIdx), aux(:, trainIdx), s0);
r2 = mean(abs(y(:, 1:n-K)).^2, 2);
ok = check(all(10*log10(r2) < 1.5), sprintf('two jammers (30, 25 dB), QRD-SMI: residual %.2f dB above the noise', ...
                                            10*log10(max(r2)))) && ok;

%% 2) jammer model of generate_rx_signal.m
cfg = radar_config();
cfg.clutter.enable = false;  cfg.rx.addLeakage = false;  cfg.rx.stc.enable = false;
cfg.jammer.enable = true;  cfg.jammer.az_deg = 60;  cfg.jammer.el_deg = 2;
sys = setup_radar(cfg);
c0 = cfg;  c0.jammer.enable = false;  s0sys = setup_radar(c0);
idx = sys.wave.long.rxIdx;
[I, Q] = ddc(adc_sample(generate_rx_signal(s0sys, [], 0, 60, []), s0sys), s0sys);
x0 = complex(double(I(idx, :)), double(Q(idx, :)));
pN = mean(abs(x0(:)).^2);
mat = struct();
for azA = [60 70]
    [I, Q] = ddc(adc_sample(generate_rx_signal(sys, [], 0, azA, []), sys), sys);
    x = complex(double(I(idx, :)), double(Q(idx, :)));
    [jnr, ph] = jam_stats(x, pN);
    azRel = wrap180(60 - azA);
    [~, u] = steering_vector(sys, azRel, 2, sys.jammer.lambda_m);
    expJ = sys.jammer.jnr_dB + 10*log10(pattern_gain(sys.ant.az_deg, sys.ant.azGain_dB, azRel));
    expPh = angle(exp(1j*2*pi*sys.ant.elementSpacing_m/sys.jammer.lambda_m*u));
    lim = 0.5 + 0.5*(azA ~= 60);
    ok = check(abs(jnr - expJ) < lim && abs(angle(exp(1j*(ph - expPh)))) < 0.02, ...
               sprintf('MATLAB jammer, antenna at %g deg: %.2f dB above the noise (expected %.2f), phase step %.3f (%.3f)', ...
                       azA, jnr, expJ, ph, expPh)) && ok;
    mat.(sprintf('az%d', azA)) = [jnr ph];
end

%% 3) C generator and real-time simulator
gen = fullfile(build, ['radar_generate_ddc' exeSfx]);
if exist(gen, 'file')
    cfgM = fullfile(out, 'jam_cfg.m');  tgtCsv = fullfile(out, 'jam_targets.csv');
    c = cfg;  c.sim.sectors_deg = [59.8 60.2];  c.targets = c.targets(1);   % T1 is far from 60 deg
    write_config_m(c, cfgM);  export_targets_csv(c.targets, tgtCsv);
    for jamOn = [false true]
        folder = fullfile(out, sprintf('gen_%d', jamOn));
        if exist(folder, 'dir'), rmdir(folder, 's'); end
        [st, msg] = system(sprintf('"%s" --config-m "%s" --targets "%s" --out "%s" --no-raw --quiet --set "cfg.jammer.enable = %s;"', ...
                                   gen, cfgM, tgtCsv, folder, ternary(jamOn, 'true', 'false')));
        if st ~= 0, fprintf('%s', msg); end
        f = dir(fullfile(folder, 'ddc_*.mat'));
        x = [];
        for k = 1:numel(f)
            d = load(fullfile(folder, f(k).name));
            x = [x; complex(double(d.I(idx, :)), double(d.Q(idx, :)))]; %#ok<AGROW>
        end
        if ~jamOn, pNc = mean(abs(x(:)).^2);
        else, [jnrC, phC] = jam_stats(x, pNc); end
    end
    ok = check(abs(jnrC - mat.az60(1)) < 0.5 && abs(angle(exp(1j*(phC - mat.az60(2))))) < 0.02, ...
               sprintf('radar_generate_ddc: jammer %.2f dB above the noise (MATLAB %.2f), phase step %.3f (%.3f)', ...
                       jnrC, mat.az60(1), phC, mat.az60(2))) && ok;
    % radar_live_sim over UDP into the C receiver, which logs the DDC output of every PRI
    live = fullfile(build, ['radar_live_sim' exeSfx]);
    recv = fullfile(build, ['radar_receiver' exeSfx]);
    if exist(live, 'file') && exist(recv, 'file')
        port = 50371;
        cl = cfg;  cl.ant.az0_deg = 59.9;  cl.sim.sectors_deg = [];
        binL = fullfile(out, 'live_receiver.bin');
        evalc('export_receiver_config(cl, binL)');
        logDir = fullfile(out, 'logs');
        for fid = {'full', 'fast'}
            pw = zeros(1, 2);  phL = 0;
            for jamOn = [false true]
                logName = sprintf('live_%s_%d', fid{1}, jamOn);
                if exist(fullfile(logDir, logName), 'dir'), rmdir(fullfile(logDir, logName), 's'); end
                system(sprintf(['("%s" --config "%s" --port %d --log-dir "%s" --log-name %s --log-blocks ddc ' ...
                                '--idle-timeout 3 --quiet > "%s" 2>&1) &'], recv, binL, port, logDir, logName, ...
                               fullfile(out, 'receiver.log')));
                pause(2);
                % (the radar_config.m of the project: the real-time simulator sets the modes up itself)
                sets = ['--set "cfg.clutter.enable = false;" --set "cfg.rx.addLeakage = false;" ' ...
                        '--set "cfg.rx.stc.enable = false;" --set "cfg.jammer.az_deg = 60;" --set "cfg.jammer.el_deg = 2;" ' ...
                        '--set "cfg.ant.az0_deg = 59.9;" --set "cfg.sim.sectors_deg = [];" ' ...
                        sprintf('--set "cfg.jammer.enable = %s;"', ternary(jamOn, 'true', 'false'))];
                system(sprintf('"%s" --targets "%s" --no-board --quiet --port %d --fidelity %s --max-pulses 6 --speed 0.05 %s > "%s" 2>&1', ...
                               live, tgtCsv, port, fid{1}, sets, fullfile(out, 'live.log')));
                pause(5);                                       % (the receiver stops 3 s after the data)
                f = dir(fullfile(logDir, logName, 'ddc', 'ddc_*.mat'));
                x = [];
                for k = 1:numel(f)
                    d = load(fullfile(logDir, logName, 'ddc', f(k).name));
                    x = [x; complex(double(d.I(idx, :)), double(d.Q(idx, :)))]; %#ok<AGROW>
                end
                if isempty(x), break; end
                if ~jamOn, pw(1) = mean(abs(x(:)).^2);
                else, [pw(2), phL] = jam_stats(x, pw(1)); end
            end
            if pw(2) ~= 0
                ok = check(abs(pw(2) - mat.az60(1)) < 1 && abs(angle(exp(1j*(phL - mat.az60(2))))) < 0.02, ...
                           sprintf('radar_live_sim (%s): jammer %.2f dB above the noise (MATLAB %.2f), phase step %.3f', ...
                                   fid{1}, pw(2), mat.az60(1), phL)) && ok;
            else
                ok = check(false, sprintf('radar_live_sim (%s): no PRIs logged by the receiver', fid{1})) && ok;
            end
        end
    end
else
    fprintf('  radar_generate_ddc not built (%s): C generator checks skipped\n', gen);
end

%% 4) C receiver against MATLAB: DBF and SLC outputs
recv = fullfile(build, ['radar_receiver' exeSfx]);
send = fullfile(build, ['radar_stream_send' exeSfx]);
c = radar_config();
c.display.enable = false;  c.console.enable = false;  c.io.saveRaw = false;
c.io.saveDdc = true;  c.io.outFolder = out;  c.io.ddcFolder = fullfile(out, 'ddc_pulses');
c.sim.sectors_deg = [59.85 60.15];  c.targets = c.targets(1);
c.rx.nChannelsUsedDbf = 1:44;  c.rx.nChannelsUsedSlc = 45:48;  c.slc.enable = true;  c.jammer.enable = true;
if exist(c.io.ddcFolder, 'dir'), rmdir(c.io.ddcFolder, 's'); end
run_radar_sim(c);
f = dir(fullfile(c.io.ddcFolder, 'ddc_*.mat'));
if exist(recv, 'file') && exist(send, 'file')
    for m = 0:2
        cm = c;  cm.slc.algorithm = names{m+1};
        sysM = setup_radar(cm);
        cfgFile = fullfile(out, 'slc_receiver.bin');
        evalc('export_receiver_config(cm, cfgFile)');
        logDir = fullfile(out, 'logs');  logName = names{m+1};
        if exist(fullfile(logDir, logName), 'dir'), rmdir(fullfile(logDir, logName), 's'); end
        port = 50381 + m;
        system(sprintf(['("%s" --config "%s" --port %d --log-dir "%s" --log-name %s --log-blocks dbf,slc --threads 2 ' ...
                        '--idle-timeout 3 --quiet > "%s" 2>&1) &'], recv, cfgFile, port, logDir, logName, ...
                       fullfile(out, 'receiver.log')));
        pause(2);
        system(sprintf('"%s" "%s" --port %d --speed 0.2 --quiet', send, c.io.ddcFolder, port));
        t0 = tic;
        while toc(t0) < 30 && numel(dir(fullfile(logDir, logName, 'slc', 'slc_*.mat'))) < numel(f), pause(1); end
        pause(5);                                               % (the receiver stops after 3 s without data)
        st = struct('short', [], 'long', []);
        errD = 0;  errS = 0;  nCmp = 0;
        for k = 1:numel(f)
            d = load(fullfile(c.io.ddcFolder, f(k).name));
            lf = fullfile(logDir, logName, 'slc', sprintf('slc_%07d.mat', d.hdr.pulseNumber));
            df = fullfile(logDir, logName, 'dbf', sprintf('dbf_%07d.mat', d.hdr.pulseNumber));
            bb = complex(single(d.I(:, sysM.chDbf)), single(d.Q(:, sysM.chDbf)));
            bbAux = complex(single(d.I(:, sysM.slc.ch)), single(d.Q(:, sysM.slc.ch)));
            for p = {'short', 'long'}
                wv = sysM.wave.(p{1});
                bm = digital_beamforming(ddc_segment(bb, wv), sysM.dbf.W);
                mainTrain = digital_beamforming(bb(wv.slcTrainIdx, :), sysM.dbf.W);
                [y, st.(p{1})] = sidelobe_canceller(st.(p{1}), bm, ddc_segment(bbAux, wv).', mainTrain, ...
                                                    bbAux(wv.slcTrainIdx, :).', sysM.slc);
                if exist(lf, 'file') && exist(df, 'file')
                    L = load(lf);  D = load(df);
                    errD = max(errD, norm(double(D.(p{1})(:)) - double(bm(:))) / norm(double(bm(:))));
                    errS = max(errS, norm(double(L.(p{1})(:)) - double(y(:))) / norm(double(y(:))));
                    nCmp = nCmp + 1;
                end
            end
        end
        ok = check(nCmp == 2*numel(f) && errD < 1e-5 && errS < 1e-3, ...
                   sprintf('C receiver = MATLAB, %s: %d of %d sub-pulse outputs, DBF %.1e, SLC %.1e relative', ...
                           upper(names{m+1}), nCmp, 2*numel(f), errD, errS)) && ok;
    end
else
    fprintf('  radar_receiver / radar_stream_send not built: C receiver comparison skipped\n');
end

%% 5) element selection, single main element, PPI modes
c5 = radar_config();
c5.rx.nChannelsUsedDbf = [3 10 17 40];  c5.rx.nChannelsUsedSlc = [1 48];  c5.slc.enable = true;
s5 = setup_radar(c5);
modes = ppi_display_modes(s5);
ok = check(isequal(s5.chDbf, [3 10 17 40]) && isequal(s5.slc.ch, [1 48]) && any(strcmp(modes, 'slc')) && ...
           any(strcmp(modes, 'dbf')), 'DBF [3 10 17 40], SLC [1 48]: the display modes have DBF and SLC out') && ok;
c5.rx.nChannelsUsedDbf = 5;  c5.rx.nChannelsUsedSlc = [6 7];
s5 = setup_radar(c5);
modes = ppi_display_modes(s5);
ok = check(~s5.dbf.enabled && s5.nBeams == 1 && any(strcmp(modes, 'slc')) && ~any(strcmp(modes, 'dbf')), ...
           'one main element (no DBF) with an SLC: SLC out, no DBF out') && ok;
s5 = setup_radar(radar_config());
modes = ppi_display_modes(s5);
ok = check(~any(strcmp(modes, 'slc')) && any(strcmp(modes, 'dbf')), 'without the SLC: no SLC out') && ok;
% the single element with the SLC through the processing (a jammer: the SLC output below the element)
c5 = radar_config();
c5.display.enable = false;  c5.console.enable = false;  c5.io.saveRaw = false;  c5.clutter.enable = false;
c5.sim.sectors_deg = [59.9 60.1];  c5.targets = c5.targets(1);
c5.rx.nChannelsUsedDbf = 5;  c5.rx.nChannelsUsedSlc = [6 7];  c5.slc.enable = true;  c5.jammer.enable = true;
try
    r5 = run_radar_sim(c5);
    ok = check(r5.nProcessed > 0, sprintf('run_radar_sim with one main element and an SLC: %d PRIs', r5.nProcessed)) && ok;
catch err
    ok = check(false, ['run_radar_sim with one main element and an SLC: ' err.message]) && ok;
end

%% 6) a weak target behind the elevation sidelobe of a jammer
if ~quick
    c6 = radar_config();
    c6.display.enable = false;  c6.console.enable = false;  c6.io.saveRaw = false;  c6.console.plots = false;
    c6.targets = c6.targets(3);                                % T3: 34 km, el 10 deg, at az 160 deg
    c6.targets.rcs_m2 = 0.05;  c6.sim.swerling = 0;
    c6.sim.sectors_deg = [153 167];
    c6.rx.nChannelsUsedDbf = 1:44;  c6.rx.nChannelsUsedSlc = 45:48;
    c6.jammer.az_deg = 160;  c6.jammer.el_deg = 2;  c6.jammer.erp_W = 2;
    runs = {'no jammer', false, false, ''; 'jammer, no SLC', true, false, ''; 'jammer, QRD-SMI', true, true, 'qrd-smi'; ...
            'jammer, QRD-RLS', true, true, 'qrd-rls'; 'jammer, NLMS', true, true, 'nlms'};
    snr = nan(1, size(runs, 1));
    for k = 1:size(runs, 1)
        cc = c6;  cc.jammer.enable = runs{k, 2};  cc.slc.enable = runs{k, 3};
        if runs{k, 3}, cc.slc.algorithm = runs{k, 4}; end
        r = run_radar_sim(cc);
        p = r.plots;
        if ~isempty(p)
            near = abs([p.range_m] - 34.2e3) < 600;
            if any(near), snr(k) = max([p(near).snr_dB]); end
        end
        fprintf('  %-18s: %d plots, target SNR %.1f dB\n', runs{k, 1}, numel(p), snr(k));
    end
    ok = check(isfinite(snr(1)) && isnan(snr(2)), sprintf(['weak target: detected without the jammer (%.1f dB), ' ...
               'missed with it'], snr(1))) && ok;
    for k = 3:5
        ok = check(isfinite(snr(k)) && snr(k) > snr(1) - 6, sprintf('%s: the target again, SNR %.1f dB (no jammer %.1f dB)', ...
                   runs{k, 1}, snr(k), snr(1))) && ok;
    end
end
fprintf('Sidelobe canceller, element lists and jammers: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function slc = slc_settings(method)
names = {'qrd-smi', 'qrd-rls', 'nlms'};
slc = struct('method', method, 'algorithm', names{method+1}, 'smiLoading_dB', -30, 'rlsLambda', 0.999, ...
             'rlsDelta_dB', -20, 'nlmsMu', 0.1, 'nlmsEps', 1);
end

function [jnr, ph] = jam_stats(x, pN)
% Jammer power over the noise pN [dB] and the mean phase step from element to element.
jnr = 10*log10(mean(abs(x(:)).^2) / pN - 1);
ph = angle(sum(sum(x(:, 2:end) .* conj(x(:, 1:end-1)))));
end

function ok = check(cond, msg)
ok = logical(cond);
fprintf('  %-110s %s\n', msg, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

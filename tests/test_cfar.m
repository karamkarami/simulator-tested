function ok = test_cfar()
%TEST_CFAR  CFAR detectors (cfar_detect.m): CA, GO, SO, OS with the 2-D, 1-D range and 1-D Doppler windows.
%   For every detector (cfg.cfar.type) and window (cfg.cfar.method), on maps of
%   exponential noise (one look) with the threshold for Pfa = 1e-3
%   (cfar_alpha.m for that detector):
%   - the false alarm rate of the interior cells is close to 1e-3 (the OS
%     detector tests every cell here: sys.cfar.osAllCells), and the one of
%     the range-edge cells (shifted window) is reported (CA: close to 1e-3;
%     GO / SO: one half is smaller or empty there),
%   - a target in the first and in the last range cell is detected,
%   - only the cells of the tested range interval wv.testIdx are tested,
%   - the noise estimate equals the statistic computed directly from the
%     training cells of the (shifted) window (mean, mean of the larger /
%     smaller half, k-th smallest) / sys.cfar.noiseScale, also when the whole
%     tested interval lies within the first or last guard + training cells,
%   - OS with only the local maxima tested (the default) detects exactly the
%     local maxima that every-cell testing detects.
%   cfar_pfa.m also agrees with the closed forms of the single-look Pfa of
%   CA, GO, SO (Gandhi and Kassam) and OS (Rohling).

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = closed_forms();
ok = range_limited_notch() && ok;
rng(11);
for type = {'ca', 'go', 'so', 'os'}
    for method = {'2d', 'range', 'doppler'}
        cfg = radar_config();
        cfg.console.enable = false;
        cfg.cfar.type = type{1};
        cfg.cfar.method = method{1};
        sys = setup_radar(cfg);
        sys.dop.notchBins = 0;                             % test every Doppler bin
        sys.cfar.alpha = cfar_alpha(1e-3, sys.cfar, 1);
        sys.cfar.noiseScale = noise_scale(sys.cfar, 1);
        sys.cfar.osAllCells = true;
        oR = sys.cfar.guardRange + sys.cfar.trainRange;
        nR = 240;  nD = 64;  nB = 4;
        nTrials = 12;
        if strcmp(type{1}, 'os') && strcmp(method{1}, '2d'), nTrials = 4; end   % 296 cells per test
        wv.testIdx = [1 nR];
        wv.range_m = (0:nR-1) * 100;
        fa = 0;  faEdge = 0;  nEdge = 0;  nInner = 0;
        edge = false(1, nR);
        edge([1:oR, nR-oR+1:nR]) = true;
        for k = 1:nTrials
            P = single(-log(rand(nB, nR, nD)));            % exponential noise, mean 1
            isDet = cfar_detect(P, sys, wv);
            fa = fa + nnz(isDet(:, ~edge, :));
            faEdge = faEdge + nnz(isDet(:, edge, :));
            nInner = nInner + nB*nD*sum(~edge);
            nEdge = nEdge + nB*nD*sum(edge);
        end
        pfaIn = fa / nInner;
        pfaEdge = faEdge / max(nEdge, 1);
        % targets in the first and last range cells (20 dB above the noise)
        P = single(-log(rand(nB, nR, nD)));
        P(2, 1, 10) = 100;  P(2, nR, 20) = 100;
        isDet = cfar_detect(P, sys, wv);
        edgeOk = isDet(2, 1, 10) && isDet(2, nR, 20);
        % tested interval: cells 50 ... 120 only
        wv2.testIdx = [50 120];  wv2.range_m = wv.range_m;
        P(:, 30, 5) = 1e3;  P(:, 80, 5) = 1e3;  P(:, 150, 5) = 1e3;
        isDet = cfar_detect(P, sys, wv2);
        hit = squeeze(any(any(isDet, 1), 3));
        limitOk = ~any(hit([1:49, 121:end])) && hit(80);
        % noise estimate = the statistic of the training cells of the (shifted) window
        noiseOk = true;
        for iv = [1 3; 1 oR+1; 4 9; nR-2 nR; nR-oR nR; nR-12 nR-6; 100 140; 1 nR].'
            wv3.testIdx = iv.';  wv3.range_m = wv.range_m;
            [isD, noise] = cfar_detect(P, sys, wv3);
            ref = direct_noise(P, sys, iv(1):iv(2));
            got = double(noise(:, iv(1):iv(2), :));
            nD = size(P, 3);  z = floor(nD/2) + 1;        % the MTI notch bins are not tested (noise 0)
            tested = repmat(reshape(abs((1:nD) - z) > sys.dop.notchBins, 1, 1, nD), size(got, 1), size(got, 2), 1);
            if strcmp(type{1}, 'os')                       % OS: the estimate only at the detections
                m = isD(:, iv(1):iv(2), :) & tested;
            else
                m = tested;
            end
            got = got(m);  ref = ref(m);
            err = abs(got - ref) ./ ref;
            noiseOk = noiseOk && (isempty(err) || max(err(:)) < 1e-4);
        end
        % OS: testing only the local maxima detects the same local maxima
        peaksOk = true;
        if strcmp(type{1}, 'os')
            P = single(-log(rand(nB, nR, nD)));
            P(3, 100, 30) = 40;  P(3, 101, 30) = 35;  P(1, 7, 3) = 60;
            all = cfar_detect(P, sys, wv);
            sys2 = sys;  sys2.cfar.osAllCells = false;
            pk = cfar_detect(P, sys2, wv);
            isMax = P >= local_max(P);
            peaksOk = isequal(pk, all & isMax) && pk(3, 100, 30) && pk(1, 7, 3);
        end
        tol = 0.25;  if nTrials < 12, tol = 0.4; end
        good = abs(pfaIn/1e-3 - 1) < tol && edgeOk && limitOk && noiseOk && peaksOk;
        if oR > 0 && strcmp(type{1}, 'ca'), good = good && abs(pfaEdge/1e-3 - 1) < 0.5; end
        edgeText = 'no range window';
        if oR > 0, edgeText = sprintf('%.2e', pfaEdge); end
        fprintf(['  %-16s %3d training cells, alpha %5.2f dB: Pfa interior %.2e, range edges %s (design 1e-3); ' ...
                 'edge targets %s; interval %s; noise = direct %s%s: %s\n'], sys.cfar.label, sys.cfar.nTrain, ...
                10*log10(sys.cfar.alpha), pfaIn, edgeText, ternary(edgeOk, 'detected', 'MISSED'), ...
                ternary(limitOk, 'ok', 'FAIL'), ternary(noiseOk, 'ok', 'FAIL'), ...
                ternary(strcmp(type{1}, 'os'), ternary(peaksOk, '; local maxima only: same', '; local maxima only: DIFFERENT'), ''), ...
                ternary(good, 'ok', 'FAIL'));
        ok = ok && good;
    end
end
fprintf('CFAR detectors (CA, GO, SO, OS), windows, edges and range interval: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ok = range_limited_notch()
% The MTI notch suppresses its Doppler bins only through the configured range.
cfg = radar_config();
cfg.console.enable = false;
cfg.cfar.method = '2d';
cfg.cfar.notchMaxRange_km = 0.5;
sys = setup_radar(cfg);
sys.dop.notchBins = 2;
sys.cfar.alpha = 10;
nR = 100;  nD = 64;  zeroBin = floor(nD/2) + 1;
wv.testIdx = [1 nR];
wv.range_m = (0:nR-1) * 100;
P = ones(1, nR, nD, 'single');
P(1, 3, zeroBin) = 100;                         % 200 m: notch still applies
P(1, 80, zeroBin) = 100;                        % 7.9 km: full CFAR applies
isDet = cfar_detect(P, sys, wv);
ok = ~isDet(1, 3, zeroBin) && isDet(1, 80, zeroBin);
fprintf('  Range-limited MTI notch (500 m): %s\n', ternary(ok, 'ok', 'FAIL'));
end

function ok = closed_forms()
% cfar_pfa against the single-look closed forms.
N = 32;  n = 16;  k = 24;
c = struct('type', 'ca', 'nTrain', N, 'nHalf', n, 'osK', k);
T = 10;
s = 0;
for j = 0:n-1, s = s + nchoosek(n+j-1, j) * (2 + T/n)^(-(n+j)); end
ref = [(1 + T/N)^(-N), 2*(1 + T/n)^(-n) - 2*s, 2*s, prod((N - (0:k-1)) ./ (N - (0:k-1) + T))];
got = zeros(1, 4);
types = {'ca', 'go', 'so', 'os'};
for i = 1:4
    c.type = types{i};
    got(i) = cfar_pfa(T, c, 1);
end
err = abs(got ./ ref - 1);
ok = all(err < 1e-6);
fprintf('  cfar_pfa against the closed forms (one look, alpha = 10): CA %.3e, GO %.3e, SO %.3e, OS %.3e, largest relative error %.1e: %s\n', ...
        got, max(err), ternary(ok, 'ok', 'FAIL'));
end

function s = noise_scale(c, L)
% Mean of the statistic for unit noise (as setup_radar.m).
s = 1;
if ~strcmp(c.type, 'ca')
    [z, w] = cfar_quadrature(L);
    s = sum(w .* (1 - cfar_statistic(z, c, L)));
end
end

function ref = direct_noise(P, sys, r)
% The statistic of the training cells of every cell of the range cells r (all
% beams and Doppler bins) / noiseScale: the outer window shifted inwards at the
% range edges, the guard cells inside the data excluded, Doppler circular.
c = sys.cfar;
gR = c.guardRange;  oR = gR + c.trainRange;
gD = c.guardDoppler;  oD = gD + c.trainDoppler;
[nB, nR, nD] = size(P);
ref = zeros(nB, numel(r), nD);
for i = 1:numel(r)
    rc = min(max(r(i), oR + 1), nR - oR);
    [RR, DD] = ndgrid(rc-oR:rc+oR, -oD:oD);
    train = ~(abs(RR - r(i)) <= gR & abs(DD) <= gD);
    for d = 1:nD
        di = mod(d - 1 + DD, nD) + 1;
        for b = 1:nB
            v = double(P(b + nB*(RR - 1) + nB*nR*(di - 1)));
            switch c.type
                case 'ca'
                    z = mean(v(train));
                case {'go', 'so'}
                    if strcmp(c.method, 'doppler')
                        A = v(train & DD < 0);  B = v(train & DD > 0);
                    else
                        A = v(train & RR < r(i));  B = v(train & RR > r(i));
                    end
                    m = [mean(A), mean(B)];
                    m = m(~isnan(m));
                    if strcmp(c.type, 'go'), z = max(m); else, z = min(m); end
                case 'os'
                    t = sort(v(train));
                    kr = max(1, min(numel(t), round(c.osK * numel(t) / c.nTrain)));
                    z = t(kr);
            end
            ref(b, i, d) = z / c.noiseScale;
        end
    end
end
end

function M = local_max(P)
[nB, nR, nD] = size(P);
M = P;
M = max(M, cat(1, M(2:end, :, :), -inf(1, nR, nD, class(P))));
M = max(M, cat(1, -inf(1, nR, nD, class(P)), M(1:end-1, :, :)));
M2 = max(M, cat(2, M(:, 2:end, :), -inf(nB, 1, nD, class(P))));
M = max(M2, cat(2, -inf(nB, 1, nD, class(P)), M(:, 1:end-1, :)));
M = max(M, max(circshift(M, 1, 3), circshift(M, -1, 3)));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

function ok = test_cfar()
%TEST_CFAR  CA-CFAR (cfar_detect.m): the three methods, the edges and the tested range interval.
%   For every method (2-D, 1-D range, 1-D Doppler) on maps of exponential
%   noise (one look) with the threshold for Pfa = 1e-3:
%   - the false alarm rate of the interior cells and of the range-edge cells
%     (shifted window) is close to 1e-3,
%   - a target in the first and in the last range cell is detected,
%   - only the cells of the tested range interval wv.testIdx are tested,
%   - the noise estimate equals the mean of the training cells of the
%     (shifted) window, also when the whole tested interval lies within the
%     first or last guard + training cells of the range axis.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'src'));
ok = true;
rng(11);
for method = {'2d', 'range', 'doppler'}
    cfg = radar_config();
    cfg.sim.verbose = false;
    cfg.cfar.method = method{1};
    sys = setup_radar(cfg);
    sys.dop.notchBins = 0;                             % test every Doppler bin
    sys.cfar.alpha = cfar_alpha(1e-3, sys.cfar.nTrain, 1);
    oR = sys.cfar.guardRange + sys.cfar.trainRange;
    nR = 240;  nD = 64;  nB = 4;  nTrials = 12;
    wv.testIdx = [1 nR];
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
    wv2.testIdx = [50 120];
    P(:, 30, 5) = 1e3;  P(:, 80, 5) = 1e3;  P(:, 150, 5) = 1e3;
    isDet = cfar_detect(P, sys, wv2);
    hit = squeeze(any(any(isDet, 1), 3));
    limitOk = ~any(hit([1:49, 121:end])) && hit(80);
    % noise estimate = mean of the training cells of the (shifted) window, also
    % when the whole tested interval lies within oR cells of an edge
    noiseOk = true;
    for iv = [1 3; 1 oR+1; 4 9; nR-2 nR; nR-oR nR; nR-12 nR-6; 100 140; 1 nR].'
        wv3.testIdx = iv.';
        [~, noise] = cfar_detect(P, sys, wv3);
        ref = direct_noise(P, sys, iv(1):iv(2));
        err = abs(double(noise(:, iv(1):iv(2), :)) - ref) ./ ref;
        noiseOk = noiseOk && max(err(:)) < 1e-4;
    end
    good = abs(pfaIn/1e-3 - 1) < 0.25 && (oR == 0 || abs(pfaEdge/1e-3 - 1) < 0.5) && edgeOk && limitOk && noiseOk;
    edgeText = 'no range window';
    if oR > 0, edgeText = sprintf('%.2e', pfaEdge); end
    fprintf(['  %-12s %3d training cells: Pfa interior %.2e, range edges %s (design 1e-3); edge targets %s; ' ...
             'interval %s; noise = direct sums %s: %s\n'], sys.cfar.label, sys.cfar.nTrain, pfaIn, edgeText, ...
            ternary(edgeOk, 'detected', 'MISSED'), ternary(limitOk, 'ok', 'FAIL'), ternary(noiseOk, 'ok', 'FAIL'), ...
            ternary(good, 'ok', 'FAIL'));
    ok = ok && good;
end
fprintf('CFAR methods, edges and range interval: %s\n', ternary(ok, 'PASS', 'FAIL'));
end

function ref = direct_noise(P, sys, r)
% Mean of the training cells of every cell of the range cells r (all beams and
% Doppler bins): the outer window shifted inwards at the range edges, the guard
% cells inside the data excluded, Doppler circular.
c = sys.cfar;
gR = c.guardRange;  oR = gR + c.trainRange;
gD = c.guardDoppler;  oD = gD + c.trainDoppler;
[nB, nR, nD] = size(P);
ref = zeros(nB, numel(r), nD);
for i = 1:numel(r)
    rc = min(max(r(i), oR + 1), nR - oR);
    for d = 1:nD
        dOut = mod(d - 1 + (-oD:oD), nD) + 1;
        dGuard = mod(d - 1 + (-gD:gD), nD) + 1;
        rGuard = max(1, r(i) - gR):min(nR, r(i) + gR);
        sOut = sum(sum(double(P(:, rc-oR:rc+oR, dOut)), 2), 3);
        sGuard = sum(sum(double(P(:, rGuard, dGuard)), 2), 3);
        ref(:, i, d) = (sOut - sGuard) / ((2*oR + 1)*numel(dOut) - numel(rGuard)*numel(dGuard));
    end
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end

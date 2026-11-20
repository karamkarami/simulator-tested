function [y, st] = sidelobe_canceller(st, main, aux, mainTrain, auxTrain, slc)
%SIDELOBE_CANCELLER  Adaptive sidelobe canceller after the DBF (one sub-pulse of one PRI).
%   [y, st] = SIDELOBE_CANCELLER(st, main, aux, mainTrain, auxTrain, slc) subtracts the
%   auxiliary elements from every main channel (the beams of the DBF, or the single element):
%       y(b, :) = main(b, :) - w_b' * aux
%   main [nBeams x n], aux [nAux x n]: the pulse compression input of the main channels and of
%   the auxiliary elements; mainTrain [nBeams x K], auxTrain [nAux x K]: the training samples
%   (setup_radar.m: wv.slcTrainIdx), in the order of the samples. The weights W = [w_1 ... w_nBeams]
%   [nAux x nBeams] minimise the residual |main - W' * aux|^2 over the training samples:
%     slc.method 0 'qrd-smi'  least squares of this PRI by the QR decomposition of the training
%                             samples, diagonal loading 10^(smiLoading_dB/10) * P * K (P = mean
%                             power of an auxiliary element over the training samples)
%     slc.method 1 'qrd-rls'  the same recursively with Givens rotations, carried over the PRIs:
%                             before every sample R and z are scaled by sqrt(rlsLambda); the first
%                             PRI starts from R = sqrt(10^(rlsDelta_dB/10) * P * K) * I, z = 0
%     slc.method 2 'nlms'     W = W + nlmsMu * a * e' / (nlmsEps + a' * a), e = m - W' * a for
%                             every training sample, carried over the PRIs (W = 0 at the start)
%   st is the state of the sub-pulse ([] or struct() at the start): st.R, st.z (QRD-RLS), st.W.
%   The C receiver (processor.c) does the same in double precision.
nAux = size(aux, 1);
nB = size(main, 1);
K = size(auxTrain, 2);
if isempty(st) || ~isfield(st, 'W'), st = struct('W', zeros(nAux, nB), 'R', [], 'z', []); end
A = double(auxTrain);
M = double(mainTrain);
P = sum(abs(A(:)).^2) / max(1, nAux * K);                     % mean power of an auxiliary element
switch slc.method
    case 0                                                    % QRD-SMI: R, z of this PRI
        % (the QR decomposition of [loading; A'] at once: the same R, z as the Givens rotations)
        X = [sqrt(10^(slc.smiLoading_dB/10) * P * K) * eye(nAux); A'];
        [Q, R] = qr(X, 0);
        z = Q' * [zeros(nAux, nB); M'];
        st.W = tri_solve(R, z);
    case 1                                                    % QRD-RLS: carried over the PRIs
        if isempty(st.R)
            st.R = sqrt(10^(slc.rlsDelta_dB/10) * P * K) * eye(nAux);
            st.z = zeros(nAux, nB);
        end
        [st.R, st.z] = qrd_update(st.R, st.z, A, M, slc.rlsLambda);
        st.W = tri_solve(st.R, st.z);
    otherwise                                                 % NLMS
        W = st.W;
        for k = 1:K
            a = A(:, k);
            e = M(:, k) - W' * a;                             % [nB x 1] errors of all beams
            W = W + (slc.nlmsMu / (slc.nlmsEps + real(a' * a))) * a * e';
        end
        st.W = W;
end
y = main - cast(st.W', 'like', main) * aux;
end

function [R, z] = qrd_update(R, z, A, M, lambda)
% Rows [a_k', m_k'] of the training samples rotated into the upper triangular [R z] by Givens
% rotations (R(i,i) real >= 0); before every row [R z] is scaled by sqrt(lambda).
nAux = size(R, 1);
sl = sqrt(lambda);
for k = 1:size(A, 2)
    if lambda ~= 1, R = sl * R;  z = sl * z; end
    x = A(:, k)';                                             % row a'
    d = M(:, k)';                                             % row m' (one value per beam)
    for i = 1:nAux
        rho = real(R(i, i));  xi = x(i);
        nu = sqrt(rho^2 + abs(xi)^2);
        if nu == 0, continue; end
        cc = rho / nu;  s = xi / nu;
        Ri = R(i, i:end);  zi = z(i, :);
        R(i, i:end) = cc * Ri + conj(s) * x(i:end);
        x(i:end) = -s * Ri + cc * x(i:end);
        z(i, :) = cc * zi + conj(s) * d;
        d = -s * zi + cc * d;
        R(i, i) = nu;  x(i) = 0;
    end
end
end

function W = tri_solve(R, z)
% W = R \ z by back substitution (R upper triangular; a zero diagonal gives a zero weight row).
n = size(R, 1);
W = zeros(size(z));
for i = n:-1:1
    if R(i, i) == 0, continue; end
    W(i, :) = (z(i, :) - R(i, i+1:n) * W(i+1:n, :)) / R(i, i);
end
end

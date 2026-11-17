function T = cfar_ce_table(N, m, L)
%CFAR_CE_TABLE  Distribution function of the censored-mean statistic of the CE-CFAR (table).
%   T = CFAR_CE_TABLE(N, m, L) returns T = [log(z); log(G(z))] (2 x nz) on a logarithmic grid of z,
%   with G(z) = P(Z <= z) for the statistic Z of the CE (censored) CFAR: the mean of the m smallest of
%   N training cells, every cell Gamma(L) noise of mean 1 (L looks, not necessarily whole). It is
%   computed by conditioning on the m-th smallest cell t (its distribution: u = F(t) ~ Beta(m, N-m+1);
%   the other m-1 cells kept are then independent cells conditioned to lie below t) and the
%   distribution of their sum by the FFT of the (m-1)-fold convolution of their probabilities on a
%   lattice (cell masses at the centres of K1 cells of [0 t]):
%     G(z) = int_0^1 P(t + sum of the m-1 cells <= m z | t) dv,  v = Beta distribution function of u,
%   integrated by 8-point Gauss-Legendre panels in log(v) from -70 to 0 (v near 0 gives the small z of
%   the lower tail). The C tools compute the same table (radar_export_config: export_config.c); the C
%   receiver interpolates the exported table for the false alarm probability (cfar.c).
%   Below the z where G reaches 1e-12 (where the FFT rounding would dominate) log G continues on the
%   straight line (in log-log) through the points where G is 1e-12 and 1e-9 (the true tail falls even
%   faster: the Pfa from this part is negligible and conservative).
%   The table is kept for the last (N, m, L) (persistent), see cfar_statistic.m.
persistent key tab
if ~isempty(key) && isequal(key, [N m L])
    T = tab;
    return
end
x8 = [-0.9602898564975363 -0.7966664774136267 -0.5255324099163290 -0.1834346424956498 ...
       0.1834346424956498  0.5255324099163290  0.7966664774136267  0.9602898564975363];
g8 = [ 0.1012285362903763  0.2223810344533745  0.3137066458778873  0.3626837833783620 ...
       0.3626837833783620  0.3137066458778873  0.2223810344533745  0.1012285362903763];
nz = 1500;
lz = linspace(log(1e-4), log(30), nz);
z  = exp(lz);
G  = zeros(1, nz);
nPan = 35;  e = linspace(-70, 0, nPan + 1);
n = m - 1;                                       % cells kept besides the m-th smallest
Kmax = 16384;
for k = 1:nPan
    h = e(k+1) - e(k);
    for i = 1:8
        lv = (e(k) + e(k+1))/2 + x8(i) * h/2;
        v  = exp(lv);
        wv = g8(i) * h/2 * v;                    % dv = v d(log v)
        u  = betaincinv(v, m, N - m + 1);
        t  = gammaincinv(u, L) / L;
        if n == 0
            G = G + wv * (m*z >= t);
            continue
        end
        K1 = max(4, min(512, floor((Kmax - 1) / n)));
        hc = t / K1;
        c  = gammainc(L * hc * (0:K1), L);
        p  = diff(c) / c(end);                   % masses of the K1 cells of [0 t]
        Mv = n * (K1 - 1) + 1;                   % lattice points of the sum
        M  = 2^nextpow2(Mv);
        q  = real(ifft(fft(p, M).^n));
        C  = cumsum(q(1:Mv));
        pos = ((0:Mv-1) + n/2 + 0.5) * hc;       % CDF just after lattice point j (the sum (j + n/2) hc)
        y  = m*z - t;
        Gy = interp1(pos, C, y, 'linear');
        Gy(y < pos(1)) = 0;
        Gy(y > pos(end)) = 1;
        G  = G + wv * Gy;
    end
end
G = min(max(G, 0), 1);
lG = log(max(G, realmin));
i12 = find(G >= 1e-12, 1);
i9 = find(G >= 1e-9, 1);
if ~isempty(i12) && i12 > 1 && ~isempty(i9) && i9 > i12
    s = (lG(i9) - lG(i12)) / (lz(i9) - lz(i12));
    lG(1:i12-1) = lG(i12) + s * (lz(1:i12-1) - lz(i12));
end
T = [lz; lG];
key = [N m L];
tab = T;
end

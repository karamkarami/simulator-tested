function pfa = cfar_pfa(alpha, c, nLooks)
%CFAR_PFA  False alarm probability of a CFAR detector (CA, GO, SO or OS) with threshold factor alpha.
%   pfa = CFAR_PFA(alpha, c, nLooks) for the detector c (sys.cfar: c.type =
%   'ca', 'go', 'so' or 'os', c.nTrain training cells, c.nHalf cells per half
%   window (GO, SO), c.osK rank of the order statistic (OS)); c may also be
%   the number of training cells (CA). The noise power of a cell is the mean
%   of nLooks exponential (square-law) samples, a Gamma(nLooks) variable with
%   mean 1 (nLooks need not be an integer). A cell S is detected when
%   S > alpha * Z with the statistic Z of the training cells:
%     CA  the mean of the nTrain cells: Pfa = I(1/(1 + alpha/nTrain);
%         nTrain*nLooks, nLooks) (regularized incomplete beta function)
%     GO  the greater of the means of the two halves (nHalf cells each)
%     SO  the smaller of them
%     OS  the osK-th smallest of the nTrain cells
%   For GO, SO and OS Pfa = P(Z < S/alpha) = int_0^inf G(s/alpha) f_S(s) ds
%   with the distribution function G of Z (cfar_statistic.m) and the
%   density f_S of S, integrated by Gauss-Legendre panels (the C receiver
%   uses the same rule: cfar.c).
%   alpha may be an array.
if isnumeric(c), c = struct('type', 'ca', 'nTrain', c); end
L = nLooks;
if strcmp(c.type, 'ca')
    N = c.nTrain;
    pfa = betainc(1 ./ (1 + alpha/N), N*L, L);
    return
end
[s, w] = cfar_quadrature(L);
lf = L*log(L) + (L - 1)*log(s) - L*s - gammaln(L);   % log density of S (mean 1)
pfa = zeros(size(alpha));
for i = 1:numel(alpha)
    pfa(i) = sum(w .* cfar_statistic(s / alpha(i), c, L) .* exp(lf));
end
end

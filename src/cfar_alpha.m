function alpha = cfar_alpha(pfa, nTrain, nLooks)
%CFAR_ALPHA  Threshold factor of a cell-averaging CFAR.
%   alpha = CFAR_ALPHA(pfa, nTrain, nLooks) returns the factor that multiplies
%   the mean of nTrain training cells. The noise power of every cell is the
%   average of nLooks independent exponential (square-law) samples, i.e. a
%   Gamma(nLooks) variable (nLooks = 1: no integration; non-integer values
%   approximate correlated rows). With S the cell under test and Q the sum of
%   the training cells, Pfa = P(S > (alpha/nTrain) Q) and S/(S+Q) is Beta
%   distributed. For nLooks = 1: alpha = nTrain*(pfa^(-1/nTrain) - 1).

N = nTrain;
L = nLooks;
pfaOf  = @(a) betainc(1 ./ (1 + a/N), N*L, L);
alpha0 = N * (pfa^(-1/N) - 1);                           % single-look value
alpha  = fzero(@(a) log(pfaOf(a)) - log(pfa), [1e-3 10] * alpha0);
end

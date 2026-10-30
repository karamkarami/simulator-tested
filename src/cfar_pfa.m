function pfa = cfar_pfa(alpha, nTrain, nLooks)
%CFAR_PFA  False alarm probability of a cell-averaging CFAR with factor alpha.
%   pfa = CFAR_PFA(alpha, nTrain, nLooks) is the inverse of cfar_alpha.m:
%   Pfa = P(S > (alpha/nTrain) Q) = I(1/(1 + alpha/nTrain); nTrain*nLooks, nLooks)
%   with I the regularized incomplete beta function (see cfar_alpha.m).
pfa = betainc(1 ./ (1 + alpha/nTrain), nTrain*nLooks, nLooks);
end

function [s, w] = cfar_quadrature(L)
%CFAR_QUADRATURE  Nodes and weights for the integrals over the cell power S (Gamma(L), mean 1).
%   [s, w] = CFAR_QUADRATURE(L): 8-point Gauss-Legendre rule on 2000 equal
%   panels of [0, sHi], sHi = (L + 10*sqrt(L) + 75) / L, beyond which the
%   density of S is below 1e-30 (cfar_pfa.m; the C receiver: cfar.c).
x = [-0.9602898564975363 -0.7966664774136267 -0.5255324099163290 -0.1834346424956498 ...
      0.1834346424956498  0.5255324099163290  0.7966664774136267  0.9602898564975363];
g = [ 0.1012285362903763  0.2223810344533745  0.3137066458778873  0.3626837833783620 ...
      0.3626837833783620  0.3137066458778873  0.2223810344533745  0.1012285362903763];
sHi = (L + 10*sqrt(L) + 75) / L;
nP = 2000;
h = sHi / nP;
mid = ((0:nP-1).' + 0.5) * h;
s = reshape((mid + x * h/2).', [], 1);
w = reshape(repmat(g * h/2, nP, 1).', [], 1);
end

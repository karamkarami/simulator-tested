function alpha = cfar_alpha(pfa, c, nLooks)
%CFAR_ALPHA  Threshold factor of a CFAR detector (CA, GO, SO or OS) for a false alarm probability.
%   alpha = CFAR_ALPHA(pfa, c, nLooks) inverts cfar_pfa.m: the factor that
%   multiplies the noise statistic of the training cells (c = sys.cfar, or
%   the number of training cells for CA) so that a cell of pure noise
%   (Gamma(nLooks), nLooks looks) exceeds it with probability pfa.
%   CA, one look: alpha = nTrain*(pfa^(-1/nTrain) - 1).
if isnumeric(c), c = struct('type', 'ca', 'nTrain', c); end
N = c.nTrain;
alpha0 = N * (pfa^(-1/N) - 1);                           % CA, single look
f = @(a) log(cfar_pfa(a, c, nLooks)) - log(pfa);
lo = 1e-3 * alpha0;  hi = 10 * alpha0;
while f(hi) > 0, hi = 4 * hi; end                        % SO, OS with a high rank: larger factors
while f(lo) < 0, lo = lo / 4; end
alpha = fzero(f, [lo hi]);
end

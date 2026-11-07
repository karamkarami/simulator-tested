function G = cfar_statistic(z, c, L)
%CFAR_STATISTIC  Distribution function of the noise statistic of a CFAR detector.
%   G = CFAR_STATISTIC(z, c, L) = P(Z <= z) for the statistic Z of the
%   training cells of detector c (sys.cfar; see cfar_pfa.m) when every cell
%   holds Gamma(L) noise of mean 1 (L looks):
%     CA  mean of N cells:  gammainc(N*L*z, N*L)
%     GO  max of the means of two halves of n cells: F_n(z)^2
%     SO  min of them: 1 - (1 - F_n(z))^2,  F_n(z) = gammainc(n*L*z, n*L)
%     OS  k-th smallest of N cells: betainc(F_1(z), k, N - k + 1)
%         (at least k of the N cells below z)
z = max(z, 0);
switch c.type
    case 'ca'
        N = c.nTrain;
        G = gammainc(N*L*z, N*L);
    case 'go'
        n = c.nHalf;
        G = gammainc(n*L*z, n*L).^2;
    case 'so'
        n = c.nHalf;
        F = gammainc(n*L*z, n*L);
        G = F .* (2 - F);                            % = 1 - (1 - F)^2, also for a tiny F
    case 'os'
        N = c.nTrain;  k = c.osK;
        G = betainc(gammainc(L*z, L), k, N - k + 1);
    otherwise
        error('Unknown CFAR type %s.', c.type);
end
end

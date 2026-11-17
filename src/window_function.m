function [w, mainlobe] = window_function(name, N, p, x)
%WINDOW_FUNCTION  Window of the range (LFM matched filter), Doppler FFT and NUFFT processing.
%   [w, mainlobe] = WINDOW_FUNCTION(name, N, p) returns the N x 1 window name on uniformly spaced
%   samples (not normalised) and the half-width of its mainlobe (first null) in DFT bins of N
%   samples. The windows (no toolbox needed; the C receiver tools use the same definitions):
%     'rect'           rectangular (no taper)                         mainlobe 1
%     'hamming'        0.54 - 0.46 cos                                 2
%     'hann'           0.5 - 0.5 cos                                   2
%     'blackman'       0.42 - 0.5 cos + 0.08 cos(2.)                   3
%     'blackmanharris' 4-term Blackman-Harris (-92 dB sidelobes)       4
%     'taylor'         Taylor, p.nbar nearly constant sidelobes of p.sll_dB (taylor_window.m)
%     'chebyshev'      Dolph-Chebyshev, all sidelobes at p.sll_dB (as MATLAB's chebwin)
%     'kaiser'         Kaiser with shape p.beta (I0 Bessel function)
%   The cosine windows and the Kaiser window are symmetric over k/(N-1), k = 0 ... N-1; the Taylor
%   window is sampled at (k + 0.5)/N (taylor_window.m, taylorwin).
%   p (struct, optional): sll_dB (default -35), nbar (4), beta (6).
%
%   w = WINDOW_FUNCTION(name, N, p, x) evaluates the window at the positions x (in [0, 1], any
%   number, e.g. non-uniform sample times scaled to [0, 1]: the window of the NUFFT): the cosine,
%   Taylor and Kaiser windows from their continuous definition, the Chebyshev window (defined only
%   on N uniform samples) by linear interpolation of its N-sample design.
if nargin < 3 || isempty(p), p = struct(); end
sll = -35;  nbar = 4;  beta = 6;
if isfield(p, 'sll_dB') && ~isempty(p.sll_dB), sll = p.sll_dB; end
if isfield(p, 'nbar') && ~isempty(p.nbar), nbar = p.nbar; end
if isfield(p, 'beta') && ~isempty(p.beta), beta = p.beta; end
name = lower(name);
if any(strcmp(name, {'rectangular', 'none'})), name = 'rect'; end
if any(strcmp(name, {'taylor', 'chebyshev'})) && ~(isscalar(sll) && isreal(sll) && sll < 0)
    error('The %s window needs a negative sidelobe level [dB] (is %s).', name, mat2str(sll));
end
if strcmp(name, 'taylor') && ~(isscalar(nbar) && nbar >= 1 && nbar == round(nbar))
    error('The Taylor window needs nbar: a whole number >= 1 (is %s).', mat2str(nbar));
end
if strcmp(name, 'kaiser') && ~(isscalar(beta) && isreal(beta) && beta >= 0)
    error('The Kaiser window needs beta >= 0 (is %s).', mat2str(beta));
end
uniform = nargin < 4;
if uniform
    if strcmp(name, 'taylor')
        x = ((0:N-1)' + 0.5) / N;
    else
        x = (0:N-1)' / max(N - 1, 1);
    end
end
x = x(:);
c = 2*pi*x;
switch name
    case 'rect',     w = ones(size(x));  mainlobe = 1;
    case 'hamming',  w = 0.54 - 0.46*cos(c);  mainlobe = 2;
    case 'hann',     w = 0.5 - 0.5*cos(c);  mainlobe = 2;
    case 'blackman', w = 0.42 - 0.5*cos(c) + 0.08*cos(2*c);  mainlobe = 3;
    case 'blackmanharris'
        w = 0.35875 - 0.48829*cos(c) + 0.14128*cos(2*c) - 0.01168*cos(3*c);  mainlobe = 4;
    case 'taylor'
        A = acosh(10^(-sll/20)) / pi;
        sp2 = nbar^2 / (A^2 + (nbar - 0.5)^2);
        w = ones(size(x));
        n = (1:nbar-1)';
        for m = 1:nbar-1
            q = n(n ~= m);
            num = prod(1 - (m^2/sp2) ./ (A^2 + (n - 0.5).^2));
            den = prod(1 - m^2 ./ q.^2);
            w = w + 2 * (-1)^(m+1) * num / (2*den) * cos(2*pi*m*(x - 0.5));
        end
        mainlobe = sqrt(sp2) * sqrt(A^2 + 0.25);           % first null: sigma * sqrt(A^2 + (1 - 0.5)^2)
    case 'kaiser'
        r = 1 - (2*x - 1).^2;
        w = besseli(0, beta * sqrt(max(r, 0))) / besseli(0, beta);
        mainlobe = sqrt(1 + (beta/pi)^2);
    case 'chebyshev'
        wc = chebyshev_window(N, sll);
        if uniform
            w = wc;
        else
            w = interp1((0:N-1)' / max(N - 1, 1), wc, min(max(x, 0), 1), 'linear');
            if N == 1, w = ones(size(x)); end
        end
        r = 10^(-sll/20);
        if N > 1
            x0 = cosh(acosh(r) / (N - 1));
            mainlobe = N / pi * acos(cos(pi / (2*(N - 1))) / x0);
        else
            mainlobe = 1;
        end
    otherwise
        error(['Unknown window ''%s'': ''rect'', ''hamming'', ''hann'', ''blackman'', ''blackmanharris'', ' ...
               '''taylor'', ''chebyshev'' or ''kaiser''.'], name);
end
if uniform && N == 1, w = 1; end
if ~any(w), w(:) = 1; end                                   % (Hann/Blackman of 2 samples are zero)
end

function w = chebyshev_window(N, sll)
% Dolph-Chebyshev window of N samples with all sidelobes at sll [dB] (the algorithm of chebwin:
% the DFT of the Chebyshev polynomial T_{N-1}(x0 cos(pi k / N))), peak 1.
if N == 1, w = 1; return; end
r  = 10^(-sll/20);
x0 = cosh(acosh(r) / (N - 1));
k  = (0:N-1)';
xc = x0 * cos(pi * k / N);
P  = zeros(N, 1);
hi = xc > 1;  lo = xc < -1;  mid = ~hi & ~lo;
P(hi)  = cosh((N - 1) * acosh(xc(hi)));
P(lo)  = (1 - 2*rem(N - 1, 2)) * cosh((N - 1) * acosh(-xc(lo)));
P(mid) = cos((N - 1) * acos(xc(mid)));
E = exp(-2i*pi * k * k.' / N);                              % DFT matrix (N is small)
if rem(N, 2)
    v = real(E * P);
    M = (N + 1) / 2;
    v = v(1:M) / v(1);
    w = [v(M:-1:2); v];
else
    v = real(E * (P .* exp(1i*pi/N * k)));
    M = N/2 + 1;
    v = v(2:M) / v(2);
    w = [v(M-1:-1:1); v];
end
w = w / max(w);
end

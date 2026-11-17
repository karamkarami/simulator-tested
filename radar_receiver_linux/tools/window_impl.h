/*
 * window_impl.h - the windows of window_function.m (range window of the LFM matched filter,
 * Doppler FFT window, NUFFT window) for the C tools radar_export_config and
 * radar_generate_ddc. Included after their helpers fail, xcalloc, has, empty, num, str and
 * same_text (and PI): static functions of the including file.
 */
#ifndef WINDOW_IMPL_H
#define WINDOW_IMPL_H

/* ------------------------------------------------------------------------ */
/* Windows (window_function.m)                                               */
/* ------------------------------------------------------------------------ */

enum { WIN_RECT, WIN_HAMMING, WIN_HANN, WIN_BLACKMAN, WIN_BLACKMANHARRIS, WIN_TAYLOR, WIN_CHEBYSHEV, WIN_KAISER };
typedef struct { int type; double sll_dB, nbar, beta; } WinSpec;

/* Modified Bessel function I0 (series). */
static double bessel_i0(double x)
{
    double sum = 1.0, term = 1.0;
    const double q = x * x / 4.0;
    for (int k = 1; k < 500; k++) {
        term *= q / ((double)k * k);
        sum += term;
        if (term < sum * 1e-17) break;
    }
    return sum;
}

/* Dolph-Chebyshev window of N samples, all sidelobes at sll [dB], peak 1 (chebyshev_window of
 * window_function.m: the DFT of T_{N-1}(x0 cos(pi k / N))). */
static void chebyshev_window(int N, double sll, double *w)
{
    if (N == 1) { w[0] = 1.0; return; }
    const double r = pow(10.0, -sll / 20.0), x0 = cosh(acosh(r) / (N - 1));
    double *P = xcalloc((size_t)N, sizeof(double)), *v = xcalloc((size_t)N, sizeof(double));
    for (int k = 0; k < N; k++) {
        const double xc = x0 * cos(PI * k / N);
        if (xc > 1) P[k] = cosh((N - 1) * acosh(xc));
        else if (xc < -1) P[k] = (1 - 2 * ((N - 1) % 2)) * cosh((N - 1) * acosh(-xc));
        else P[k] = cos((N - 1) * acos(xc));
    }
    for (int n = 0; n < N; n++) {                         /* real part of the DFT */
        double s = 0.0;
        for (int k = 0; k < N; k++) {
            const double ph = -2.0 * PI * (double)n * k / N;
            if (N % 2) s += P[k] * cos(ph);
            else {                                        /* P .* exp(1i pi k / N) */
                const double a = PI * k / N;
                s += P[k] * (cos(a) * cos(ph) - sin(a) * sin(ph));
            }
        }
        v[n] = s;
    }
    if (N % 2) {
        const int M = (N + 1) / 2;
        for (int i = 0; i < M; i++) w[M - 1 + i] = v[i] / v[0];          /* [v(M:-1:2); v(1:M)] */
        for (int i = 1; i < M; i++) w[M - 1 - i] = v[i] / v[0];
    } else {
        const int M = N / 2 + 1;                          /* v = v(2:M) / v(2); w = [v(M-1:-1:1); v] */
        for (int i = 0; i < M - 1; i++) w[N / 2 + i] = v[1 + i] / v[1];
        for (int i = 0; i < M - 1; i++) w[N / 2 - 1 - i] = v[1 + i] / v[1];
    }
    double mx = 0;
    for (int k = 0; k < N; k++) mx = w[k] > mx ? w[k] : mx;
    for (int k = 0; k < N; k++) w[k] /= mx;
    free(P);
    free(v);
}

/* window_function(name, N, p) (x == NULL: N uniform samples) or window_function(name, N, p, x)
 * (n positions x in [0 1]); returns the mainlobe half-width in bins of N samples. */
static double window_eval(const WinSpec *s, int N, const double *x, int n, double *w)
{
    const int uniform = x == NULL;
    if (uniform) n = N;
    double mainlobe = 1.0;
    double A = 0, sp2 = 0;
    if (s->type == WIN_TAYLOR) {
        A = acosh(pow(10.0, -s->sll_dB / 20.0)) / PI;
        sp2 = s->nbar * s->nbar / (A * A + (s->nbar - 0.5) * (s->nbar - 0.5));
    }
    double *cheb = NULL;
    if (s->type == WIN_CHEBYSHEV) {
        cheb = xcalloc((size_t)N, sizeof(double));
        chebyshev_window(N, s->sll_dB, cheb);
    }
    for (int k = 0; k < n; k++) {
        double xk;
        if (!uniform) xk = x[k];
        else if (s->type == WIN_TAYLOR) xk = (k + 0.5) / N;
        else xk = (double)k / (N > 1 ? N - 1 : 1);
        const double c = 2.0 * PI * xk;
        switch (s->type) {
        case WIN_RECT:     w[k] = 1.0; break;
        case WIN_HAMMING:  w[k] = 0.54 - 0.46 * cos(c); break;
        case WIN_HANN:     w[k] = 0.5 - 0.5 * cos(c); break;
        case WIN_BLACKMAN: w[k] = 0.42 - 0.5 * cos(c) + 0.08 * cos(2 * c); break;
        case WIN_BLACKMANHARRIS:
            w[k] = 0.35875 - 0.48829 * cos(c) + 0.14128 * cos(2 * c) - 0.01168 * cos(3 * c); break;
        case WIN_TAYLOR: {
            double v = 1.0;
            const int nbar = (int)s->nbar;
            for (int m = 1; m <= nbar - 1; m++) {
                double num = 1.0, den = 1.0;
                for (int q = 1; q <= nbar - 1; q++) {
                    num *= 1.0 - (m * (double)m / sp2) / (A * A + (q - 0.5) * (q - 0.5));
                    if (q != m) den *= 1.0 - (double)m * m / ((double)q * q);
                }
                v += 2.0 * ((m + 1) % 2 ? -1.0 : 1.0) * num / (2.0 * den) * cos(2.0 * PI * m * (xk - 0.5));
            }
            w[k] = v;
            break;
        }
        case WIN_KAISER: {
            const double r = 1.0 - (2.0 * xk - 1.0) * (2.0 * xk - 1.0);
            w[k] = bessel_i0(s->beta * sqrt(r > 0 ? r : 0)) / bessel_i0(s->beta);
            break;
        }
        default: {                                        /* Chebyshev: linear interpolation of the design */
            if (uniform || N == 1) { w[k] = uniform ? cheb[k] : 1.0; break; }
            double t = (xk < 0 ? 0 : xk > 1 ? 1 : xk) * (N - 1);
            int i0 = (int)floor(t);
            if (i0 >= N - 1) i0 = N - 2;
            w[k] = cheb[i0] + (t - i0) * (cheb[i0 + 1] - cheb[i0]);
            break;
        }
        }
    }
    switch (s->type) {
    case WIN_RECT: mainlobe = 1; break;
    case WIN_HAMMING: case WIN_HANN: mainlobe = 2; break;
    case WIN_BLACKMAN: mainlobe = 3; break;
    case WIN_BLACKMANHARRIS: mainlobe = 4; break;
    case WIN_TAYLOR: mainlobe = sqrt(sp2) * sqrt(A * A + 0.25); break;
    case WIN_KAISER: mainlobe = sqrt(1.0 + (s->beta / PI) * (s->beta / PI)); break;
    default:
        if (N > 1) {
            const double r = pow(10.0, -s->sll_dB / 20.0), x0 = cosh(acosh(r) / (N - 1));
            mainlobe = N / PI * acos(cos(PI / (2.0 * (N - 1))) / x0);
        }
        break;
    }
    if (uniform && N == 1) w[0] = 1.0;
    int any = 0;
    for (int k = 0; k < n; k++) any |= w[k] != 0;
    if (!any) for (int k = 0; k < n; k++) w[k] = 1.0;
    free(cheb);
    return mainlobe;
}

/* The window cfg.<prefix><key> and its parameters <key>Sll_dB, <key>Nbar, <key>Beta (window_spec of
 * setup_radar.m; the range window also takes cfg.proc.taylorNbar / taylorSll_dB). */
static WinSpec window_spec(const char *prefix, const char *key, const char *def)
{
    char p[128];
    WinSpec s = { WIN_HAMMING, -35, 4, 6 };
    snprintf(p, sizeof p, "%s%s", prefix, key);
    const char *name = has(p) && !empty(p) ? str(p) : def;
    static const char *const names[] = { "rect", "hamming", "hann", "blackman", "blackmanharris", "taylor", "chebyshev", "kaiser" };
    s.type = -1;
    for (int i = 0; i < 8; i++) if (same_text(name, names[i])) s.type = i;
    if (same_text(name, "rectangular") || same_text(name, "none")) s.type = WIN_RECT;
    if (s.type < 0)
        fail("cfg.%s must be 'rect', 'hamming', 'hann', 'blackman', 'blackmanharris', 'taylor', 'chebyshev' or 'kaiser' "
             "(is '%s').", p, name);
    if (!strcmp(key, "rangeWindow")) {
        if (has("proc.taylorSll_dB") && !empty("proc.taylorSll_dB")) s.sll_dB = num("proc.taylorSll_dB");
        if (has("proc.taylorNbar") && !empty("proc.taylorNbar")) s.nbar = num("proc.taylorNbar");
    }
    snprintf(p, sizeof p, "%s%sSll_dB", prefix, key);  if (has(p) && !empty(p)) s.sll_dB = num(p);
    snprintf(p, sizeof p, "%s%sNbar", prefix, key);    if (has(p) && !empty(p)) s.nbar = num(p);
    snprintf(p, sizeof p, "%s%sBeta", prefix, key);    if (has(p) && !empty(p)) s.beta = num(p);
    if ((s.type == WIN_TAYLOR || s.type == WIN_CHEBYSHEV) && !(s.sll_dB < 0))
        fail("The %s window needs a negative sidelobe level [dB] (is %g).", names[s.type], s.sll_dB);
    if (s.type == WIN_TAYLOR && !(s.nbar >= 1 && s.nbar == floor(s.nbar)))
        fail("The Taylor window needs nbar: a whole number >= 1 (is %g).", s.nbar);
    if (s.type == WIN_KAISER && !(s.beta >= 0)) fail("The Kaiser window needs beta >= 0 (is %g).", s.beta);
    return s;
}

#endif

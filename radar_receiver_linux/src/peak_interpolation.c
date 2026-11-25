/* Sub-bin interpolation of radial velocity and beam elevation. */
#include "processor_internal.h"

#include <math.h>

/* Vertex offset of the parabola through the powers [dB] p0, p1 (peak), p2;
 * 0 when it has no maximum (or a power is zero: -Inf dB). */
static int parabola_vertex(double p0, double p1, double p2, double *delta)
{
    const double l0 = 10 * log10(p0), l1 = 10 * log10(p1), l2 = 10 * log10(p2);
    const double den = l0 - 2 * l1 + l2;
    *delta = 0.5 * (l0 - l2) / den;
    return den < 0 && isfinite(*delta);
}

/* Peak interpolation of the reports of the latest CPI update (as
 * peak_interpolation.m): the radial velocity from the parabola through the
 * Doppler bins d-1, d, d+1 (circular), the elevation from the parabola
 * through the beams b-1, b, b+1 (fractional beam -> u by linear interpolation
 * of beamU); otherwise the centres of the bin and of the beam stay. */
void peak_interpolation(const Processor *p, int pulse, PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const int nfft = cfg->nfft, nB = cfg->nBeams;
    const int doDop = cfg->dopplerInterp && nfft > 2;
    const int doEl = cfg->elevationInterp && cfg->dbfEnabled && nB > 2;
    if (!doDop && !doEl) return;
    const double df = cfg->dopFreq_Hz[1] - cfg->dopFreq_Hz[0];      /* Doppler bin width [Hz] */
    const double lambda = cfg->pulse[pulse].lambda_m;
    const double *beamU = cfg->pulse[pulse].beamU ? cfg->pulse[pulse].beamU : cfg->beamU;
    for (int i = 0; i < s->nReports; i++) {
        Detection *d = &s->reports[i];
        const int b = d->beam, r = d->rIdx, k = d->dIdx;
        double delta;
        if (doDop && parabola_vertex(P_at(s, (k + nfft - 1) % nfft, b, r), P_at(s, k, b, r),
                                     P_at(s, (k + 1) % nfft, b, r), &delta))
            d->vr_mps = -(cfg->dopFreq_Hz[k] + delta * df) * lambda / 2;
        if (doEl && b > 0 && b < nB - 1 &&
            parabola_vertex(P_at(s, k, b - 1, r), P_at(s, k, b, r), P_at(s, k, b + 1, r), &delta)) {
            const double pos = b + delta;                              /* fractional beam index */
            const int i0 = (int)floor(pos);
            const double f = pos - i0;
            const double u = (i0 + 1 < nB) ? beamU[i0] * (1 - f) + beamU[i0 + 1] * f : beamU[i0];
            d->el_deg = cfg->tilt_deg + asin(u) * 180 / PROCESSOR_PI;
        }
    }
}

#include "processor_internal.h"

#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------------ */
/* Sidelobe canceller (sidelobe_canceller.m)                                */
/* ------------------------------------------------------------------------ */

/* Givens rotations of the row [x d] (x = a^H [nA], d = m^H [nB], interleaved re/im) into the upper
 * triangular [R z] (R(i,i) real >= 0): QR decomposition of the training samples one row at a time. */
static void qrd_row(double *Rr, double *Ri, double *Zr, double *Zi, double *x, double *d, int nA, int nB)
{
    for (int i = 0; i < nA; i++) {
        const double rho = Rr[i * nA + i], xr = x[2 * i], xi = x[2 * i + 1];
        const double nu = sqrt(rho * rho + xr * xr + xi * xi);
        if (nu == 0.0) continue;
        const double c = rho / nu, sr = xr / nu, si = xi / nu;        /* s = x_i / nu */
        for (int j = i + 1; j < nA; j++) {
            const double rr = Rr[i * nA + j], ri = Ri[i * nA + j], ar = x[2 * j], ai = x[2 * j + 1];
            Rr[i * nA + j] = c * rr + (sr * ar + si * ai);            /* c R + conj(s) x */
            Ri[i * nA + j] = c * ri + (sr * ai - si * ar);
            x[2 * j]     = -(sr * rr - si * ri) + c * ar;              /* -s R + c x */
            x[2 * j + 1] = -(sr * ri + si * rr) + c * ai;
        }
        for (int b = 0; b < nB; b++) {
            const double zr = Zr[i * nB + b], zi = Zi[i * nB + b], dr = d[2 * b], di = d[2 * b + 1];
            Zr[i * nB + b] = c * zr + (sr * dr + si * di);
            Zi[i * nB + b] = c * zi + (sr * di - si * dr);
            d[2 * b]     = -(sr * zr - si * zi) + c * dr;
            d[2 * b + 1] = -(sr * zi + si * zr) + c * di;
        }
        Rr[i * nA + i] = nu;
        Ri[i * nA + i] = 0.0;
        x[2 * i] = x[2 * i + 1] = 0.0;
    }
}

/* W = R \ z by back substitution (a zero diagonal gives a zero weight row). */
static void tri_solve(const double *Rr, const double *Ri, const double *Zr, const double *Zi, double *Wr, double *Wi,
                      int nA, int nB)
{
    for (int i = nA - 1; i >= 0; i--) {
        const double rii = Rr[i * nA + i];
        for (int b = 0; b < nB; b++) {
            double ar = Zr[i * nB + b], ai = Zi[i * nB + b];
            for (int j = i + 1; j < nA; j++) {
                const double rr = Rr[i * nA + j], ri = Ri[i * nA + j], wr = Wr[j * nB + b], wi = Wi[j * nB + b];
                ar -= rr * wr - ri * wi;
                ai -= rr * wi + ri * wr;
            }
            Wr[i * nB + b] = rii != 0.0 ? ar / rii : 0.0;
            Wi[i * nB + b] = rii != 0.0 ? ai / rii : 0.0;
        }
    }
}

/* Smallest eigenvalue of the covariance A A^H / K of the auxiliary training samples (A [nA][K],
 * complex): the cyclic Jacobi method on the real symmetric matrix [Re -Im; Im Re] (2nA x 2nA, every
 * eigenvalue of the Hermitian matrix twice). work: 4 nA^2 doubles. */
static double min_eigenvalue(const double *ar, const double *ai, int nA, int K, double *work)
{
    const int n = 2 * nA;
    double *M = work;
    for (int i = 0; i < nA; i++)
        for (int j = 0; j < nA; j++) {
            double cr = 0, ci = 0;                          /* C_ij = sum_k a_ik conj(a_jk) / K */
            for (int k = 0; k < K; k++) {
                const double xr = ar[(size_t)i * K + k], xi = ai[(size_t)i * K + k];
                const double yr = ar[(size_t)j * K + k], yi = ai[(size_t)j * K + k];
                cr += xr * yr + xi * yi;
                ci += xi * yr - xr * yi;
            }
            cr /= K;
            ci /= K;
            M[i * n + j] = cr;               M[i * n + j + nA] = -ci;
            M[(i + nA) * n + j] = ci;        M[(i + nA) * n + j + nA] = cr;
        }
    for (int sweep = 0; sweep < 50; sweep++) {
        double off = 0, diag = 0;
        for (int p = 0; p < n; p++) {
            diag += M[p * n + p] * M[p * n + p];
            for (int q = p + 1; q < n; q++) off += M[p * n + q] * M[p * n + q];
        }
        if (off <= 1e-24 * diag || off == 0) break;
        for (int p = 0; p < n - 1; p++)
            for (int q = p + 1; q < n; q++) {
                const double apq = M[p * n + q];
                if (apq == 0) continue;
                const double theta = (M[q * n + q] - M[p * n + p]) / (2 * apq);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1));
                const double c = 1 / sqrt(t * t + 1), sn = t * c;
                for (int k = 0; k < n; k++) {                 /* columns p, q */
                    const double mkp = M[k * n + p], mkq = M[k * n + q];
                    M[k * n + p] = c * mkp - sn * mkq;
                    M[k * n + q] = sn * mkp + c * mkq;
                }
                for (int k = 0; k < n; k++) {                 /* rows p, q */
                    const double mpk = M[p * n + k], mqk = M[q * n + k];
                    M[p * n + k] = c * mpk - sn * mqk;
                    M[q * n + k] = sn * mpk + c * mqk;
                }
            }
    }
    double lam = M[0];
    for (int p = 1; p < n; p++) if (M[p * n + p] < lam) lam = M[p * n + p];
    return lam > 0 ? lam : 0.0;
}

/* The weights of the sub-pulse from the training samples of this PRI (QRD-SMI; QRD-RLS and NLMS
 * go on from the state of the earlier PRIs), then every beam minus the weighted auxiliary
 * elements, in place in beamRe / beamIm. */
void sidelobe_canceller(const RadarConfig *cfg, const PulseConfig *pc, const PriData *pri, PulseState *s,
                        int nThreads)
{
    const int nA = cfg->nAux, nB = cfg->nBeams, K = pc->slcTrainCount, n = pc->ddcCount;
    const int nD = cfg->nChUsed > 0 ? cfg->nChUsed : cfg->nCh;
    (void)nThreads;
    /* training samples: the auxiliary elements and the DBF output of the main channels */
    double P = 0.0;
    for (int a = 0; a < nA; a++) {
        const int16_t *x = pri->iq + 2 * ((size_t)cfg->chSlc[a] * pri->nSamples + pc->slcTrainFirst);
        for (int k = 0; k < K; k++) {
            const double re = x[2 * k], im = x[2 * k + 1];
            s->trARe[(size_t)a * K + k] = re;
            s->trAIm[(size_t)a * K + k] = im;
            P += re * re + im * im;
        }
    }
    P /= (double)nA * K;                                    /* mean power of an auxiliary element */
    /* noise power of an auxiliary element: the smallest eigenvalue of the covariance of the training
     * samples (one auxiliary element: 30 dB below its power) */
    const double Pn = nA > 1 ? min_eigenvalue(s->trARe, s->trAIm, nA, K, s->slcEig) : P * 1e-3;
    const int j0 = pc->slcTrainFirst - pc->ddcFirst;      /* the training samples in the segment */
    if (j0 >= 0 && j0 + K <= n) {                           /* ... the DBF output has them already */
        for (int b = 0; b < nB; b++)
            for (int k = 0; k < K; k++) {
                s->trMRe[(size_t)b * K + k] = s->beamRe[(size_t)b * n + j0 + k];
                s->trMIm[(size_t)b * K + k] = s->beamIm[(size_t)b * n + j0 + k];
            }
    } else {                                                /* (short pulse: segment cut at the handover) */
        memset(s->trMRe, 0, sizeof(double) * (size_t)nB * K);
        memset(s->trMIm, 0, sizeof(double) * (size_t)nB * K);
        for (int c = 0; c < nD; c++) {
            const int16_t *x = pri->iq + 2 * ((size_t)cfg->chDbf[c] * pri->nSamples + pc->slcTrainFirst);
            for (int b = 0; b < nB; b++) {
                const double wr = cfg->dbfEnabled ? cfg->wRe[b * nD + c] : 1.0, wi = cfg->dbfEnabled ? cfg->wIm[b * nD + c] : 0.0;
                double *mr = s->trMRe + (size_t)b * K, *mi = s->trMIm + (size_t)b * K;
                for (int k = 0; k < K; k++) {
                    const double re = x[2 * k], im = x[2 * k + 1];
                    mr[k] += wr * re - wi * im;
                    mi[k] += wr * im + wi * re;
                }
            }
        }
    }
    double *x = s->slcX, *d = s->slcD;
    if (cfg->slcMethod == SLC_NLMS) {
        if (!s->slcStarted) {
            memset(s->slcWRe, 0, sizeof(double) * (size_t)nA * nB);
            memset(s->slcWIm, 0, sizeof(double) * (size_t)nA * nB);
        }
        const double leak = isfinite(cfg->slcNlmsLoading_dB) ? pow(10.0, cfg->slcNlmsLoading_dB / 10.0) * Pn : 0.0;
        for (int k = 0; k < K; k++) {
            double norm = cfg->slcNlmsEps;
            for (int a = 0; a < nA; a++) norm += s->trARe[(size_t)a * K + k] * s->trARe[(size_t)a * K + k] +
                                                 s->trAIm[(size_t)a * K + k] * s->trAIm[(size_t)a * K + k];
            const double g = norm > 0 ? cfg->slcNlmsMu / norm : 0.0;
            for (int b = 0; b < nB; b++) {                  /* e_b = m_b - sum_a conj(W_ab) a_a */
                double er = s->trMRe[(size_t)b * K + k], ei = s->trMIm[(size_t)b * K + k];
                for (int a = 0; a < nA; a++) {
                    const double wr = s->slcWRe[a * nB + b], wi = s->slcWIm[a * nB + b];
                    const double ar = s->trARe[(size_t)a * K + k], ai = s->trAIm[(size_t)a * K + k];
                    er -= wr * ar + wi * ai;
                    ei -= wr * ai - wi * ar;
                }
                d[2 * b] = er;
                d[2 * b + 1] = ei;
            }
            for (int a = 0; a < nA; a++) {                  /* W_ab += g (a_a conj(e_b) - leak W_ab) */
                const double ar = s->trARe[(size_t)a * K + k], ai = s->trAIm[(size_t)a * K + k];
                for (int b = 0; b < nB; b++) {
                    const double wr = s->slcWRe[a * nB + b], wi = s->slcWIm[a * nB + b];
                    s->slcWRe[a * nB + b] += g * (ar * d[2 * b] + ai * d[2 * b + 1] - leak * wr);
                    s->slcWIm[a * nB + b] += g * (ai * d[2 * b] - ar * d[2 * b + 1] - leak * wi);
                }
            }
        }
    } else {
        const int rls = cfg->slcMethod == SLC_QRD_RLS;
        if (!rls || !s->slcStarted) {                       /* R = sqrt(loading P K) I, z = 0 */
            const double r0 = sqrt(pow(10.0, (rls ? cfg->slcRlsDelta_dB : cfg->slcSmiLoading_dB) / 10.0) * Pn * K);
            memset(s->slcRRe, 0, sizeof(double) * (size_t)nA * nA);
            memset(s->slcRIm, 0, sizeof(double) * (size_t)nA * nA);
            memset(s->slcZRe, 0, sizeof(double) * (size_t)nA * nB);
            memset(s->slcZIm, 0, sizeof(double) * (size_t)nA * nB);
            for (int a = 0; a < nA; a++) s->slcRRe[a * nA + a] = r0;
        }
        const double sl = rls ? sqrt(cfg->slcRlsLambda) : 1.0;
        for (int k = 0; k < K; k++) {
            if (sl != 1.0) {
                for (int i = 0; i < nA * nA; i++) { s->slcRRe[i] *= sl; s->slcRIm[i] *= sl; }
                for (int i = 0; i < nA * nB; i++) { s->slcZRe[i] *= sl; s->slcZIm[i] *= sl; }
            }
            for (int a = 0; a < nA; a++) {                  /* the row a^H */
                x[2 * a] = s->trARe[(size_t)a * K + k];
                x[2 * a + 1] = -s->trAIm[(size_t)a * K + k];
            }
            for (int b = 0; b < nB; b++) {                  /* m^H */
                d[2 * b] = s->trMRe[(size_t)b * K + k];
                d[2 * b + 1] = -s->trMIm[(size_t)b * K + k];
            }
            qrd_row(s->slcRRe, s->slcRIm, s->slcZRe, s->slcZIm, x, d, nA, nB);
        }
        tri_solve(s->slcRRe, s->slcRIm, s->slcZRe, s->slcZIm, s->slcWRe, s->slcWIm, nA, nB);
    }
    /* the weights applied: a second-order mean over about slcWeightMemory PRIs (two exponential means
     * of slcWeightMemory/2 PRIs one after the other; the first PRI: W), smooth enough for the MTI */
    const double bm = 1.0 - 2.0 / cfg->slcWeightMemory > 0 ? 1.0 - 2.0 / cfg->slcWeightMemory : 0.0;
    if (!s->slcStarted || bm == 0.0) {
        memcpy(s->slcBRe, s->slcWRe, sizeof(double) * (size_t)nA * nB);
        memcpy(s->slcBIm, s->slcWIm, sizeof(double) * (size_t)nA * nB);
        memcpy(s->slcARe, s->slcWRe, sizeof(double) * (size_t)nA * nB);
        memcpy(s->slcAIm, s->slcWIm, sizeof(double) * (size_t)nA * nB);
    } else {
        for (int i = 0; i < nA * nB; i++) {
            s->slcBRe[i] = bm * s->slcBRe[i] + (1.0 - bm) * s->slcWRe[i];
            s->slcBIm[i] = bm * s->slcBIm[i] + (1.0 - bm) * s->slcWIm[i];
            s->slcARe[i] = bm * s->slcARe[i] + (1.0 - bm) * s->slcBRe[i];
            s->slcAIm[i] = bm * s->slcAIm[i] + (1.0 - bm) * s->slcBIm[i];
        }
    }
    s->slcStarted = 1;
    /* y_b = x_b - sum_a conj(A_ab) aux_a over the whole segment */
    load_channels(pc, pri, cfg->chSlc, nA, s->auxRe, s->auxIm);
    #pragma omp parallel for schedule(static) num_threads(team_size(nThreads, nB)) if(nB > 1)
    for (int b = 0; b < nB; b++) {
        float *restrict yr = s->beamRe + (size_t)b * n, *restrict yi = s->beamIm + (size_t)b * n;
        for (int a = 0; a < nA; a++) {
            const float wr = (float)s->slcARe[a * nB + b], wi = (float)-s->slcAIm[a * nB + b];   /* conj(A) */
            const float *restrict ar = s->auxRe + (size_t)a * n, *restrict ai = s->auxIm + (size_t)a * n;
            IVDEP
            for (int k = 0; k < n; k++) {
                yr[k] -= wr * ar[k] - wi * ai[k];
                yi[k] -= wr * ai[k] + wi * ar[k];
            }
        }
    }
}

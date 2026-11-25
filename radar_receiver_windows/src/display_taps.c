/* Range profiles and range-Doppler maps for the receiver display. */
#include "processor_internal.h"

#include <math.h>
#include <string.h>

static void rdmap_update_velocity(Processor *p)
{
    const RadarConfig *cfg = p->cfg;
    const int vp = cfg->pulse[PULSE_LONG].enable ? PULSE_LONG : PULSE_SHORT;   /* velocity on its carrier */
    for (int d = 0; d < cfg->nfft; d++) p->rdVel[d] = -cfg->dopFreq_Hz[d] * cfg->pulse[vp].lambda_m / 2;
}

static void rdmap_clear(Processor *p)
{
    memset(p->rd, 0, sizeof(float) * (size_t)p->nBins * p->cfg->nfft);
    p->nRdDets = 0;
    p->rdLast = -1e300;
}

void rdmap_prf_changed(Processor *p)
{
    if (p->rd) {
        rdmap_update_velocity(p);
        rdmap_clear(p);
    }
}

static void taps_alloc(Processor *p)
{
    if (p->profile) return;
    const RadarConfig *cfg = p->cfg;
    p->dR = 299792458.0 / (2.0 * cfg->fsOut_Hz);
    p->nBins = (int)ceil(cfg->maxRange_m / p->dR);
    if (p->nBins < 1) p->nBins = 1;
    p->profile = proc_calloc((size_t)p->nBins, sizeof(float));
    p->profileCnt = proc_calloc((size_t)p->nBins, sizeof(int));
    const int nfftMax = cfg->maxNfft > cfg->nfft ? cfg->maxNfft : cfg->nfft;   /* (scan-to-scan PRF: the largest) */
    p->rd = proc_calloc((size_t)p->nBins * nfftMax, sizeof(float));
    p->rdVel = proc_calloc((size_t)nfftMax, sizeof(double));
    p->rdWork = proc_calloc((size_t)p->nBins * nfftMax + 1, sizeof(float));
    rdmap_update_velocity(p);
}

void processor_set_profile(Processor *p, int kind, int beam)
{
    if (beam < 0 || beam >= p->cfg->nBeams) beam = 0;
    p->profileKind = kind == PROFILE_DBF || kind == PROFILE_MF || (kind == PROFILE_SLC && p->cfg->slcEnabled) ? kind
                   : PROFILE_NONE;
    p->profileBeam = beam;
    if (p->profileKind) taps_alloc(p);
}

void processor_set_rdmap(Processor *p, int enabled, int beam, double minPeriod_s)
{
    if (beam < 0 || beam >= p->cfg->nBeams) beam = 0;
    if (enabled && (!p->rdEnabled || beam != p->rdBeam)) {   /* a new map: no old cells */
        taps_alloc(p);
        rdmap_clear(p);
    }
    p->rdEnabled = enabled != 0;
    p->rdBeam = beam;
    p->rdPeriod_s = minPeriod_s;
}

/* Range bin of a range (-1 outside the taps). */
static int tap_bin(const Processor *p, double r)
{
    const int b = (int)floor(r / p->dR);
    return b >= 0 && b < p->nBins ? b : -1;
}

/* Adds |x|^2 of the selected beam after DBF or after the SLC (pulse compression input inside the
 * receive interval) or after the pulse compression to the profile. */
void profile_add(Processor *p, int pulse, const PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const PulseConfig *pc = &cfg->pulse[pulse];
    const double handover = cfg->handover_m;
    const int b = p->profileBeam;
    if (p->profileKind == PROFILE_DBF || p->profileKind == PROFILE_SLC) {   /* (the SLC: the same samples) */
        const int n = pc->ddcCount;
        const float *xr = s->beamRe + (size_t)b * n, *xi = s->beamIm + (size_t)b * n;
        for (int k = 0; k < pc->rxCount; k++) {
            const int j = pc->rxFirst + k - pc->ddcFirst;        /* segment index of rx sample k */
            const double r = pc->range_m[0] + k * p->dR;
            if (pulse == PULSE_SHORT && !(r < handover)) break;
            const int bin = tap_bin(p, r);
            if (bin < 0 || j < 0 || j >= n) continue;
            p->profile[bin] += xr[j] * xr[j] + xi[j] * xi[j];
            p->profileCnt[bin]++;
        }
    } else {
        const float *yr = s->cpiRe + (size_t)s->mfLast * s->nCells + (size_t)b * s->nR;
        const float *yi = s->cpiIm + (size_t)s->mfLast * s->nCells + (size_t)b * s->nR;
        for (int r = 0; r < s->nR; r++) {
            if (pulse == PULSE_SHORT && !(pc->range_m[r] < handover)) break;
            const int bin = tap_bin(p, pc->range_m[r]);
            if (bin < 0) continue;
            p->profile[bin] += yr[r] * yr[r] + yi[r] * yi[r];
            p->profileCnt[bin]++;
        }
    }
}

void profile_post(Processor *p, int profileKind, const PriData *pri, const ProcessorCallbacks *cb)
{
    for (int b = 0; b < p->nBins; b++) {                   /* mean of the samples in every bin */
        p->profile[b] = p->profileCnt[b] ? p->profile[b] / p->profileCnt[b] : 0.0f;
    }
    cb->profile(cb->user, profileKind, p->profile, p->nBins, p->dR, pri->az_deg, pri->time_s);
}

/* The cells of the selected beam of a sub-pulse's map into the joined range-Doppler map. */
void rdmap_add(Processor *p, int pulse, const PulseState *s)
{
    const RadarConfig *cfg = p->cfg;
    const PulseConfig *pc = &cfg->pulse[pulse];
    const double handover = cfg->handover_m;
    const int b = p->rdBeam, nfft = cfg->nfft;
    for (int r = 0; r < s->nR; r++) {
        if (pulse == PULSE_SHORT && !(pc->range_m[r] < handover)) break;
        const int bin = tap_bin(p, pc->range_m[r]);
        if (bin < 0) continue;
        for (int d = 0; d < nfft; d++) p->rd[(size_t)d * p->nBins + bin] = P_at(s, d, b, r);
    }
    for (int i = 0; i < s->nReports; i++) {                /* the reports of this beam */
        if (s->reports[i].beam != b) continue;
        if (p->nRdDets == p->capRdDets) {
            p->capRdDets = p->capRdDets ? 2 * p->capRdDets : 64;
            p->rdDets = proc_realloc(p->rdDets, sizeof(Detection) * p->capRdDets);
        }
        p->rdDets[p->nRdDets++] = s->reports[i];
    }
    p->rdUpdated = 1;
}

void rdmap_post(Processor *p, const PriData *pri, const ProcessorCallbacks *cb)
{
    const RadarConfig *cfg = p->cfg;
    const size_t n = (size_t)p->nBins * cfg->nfft;
    size_t m = 0;
    for (size_t i = 0; i < n; i++) if (p->rd[i] > 0) p->rdWork[m++] = p->rd[i];
    RdMap map;
    memset(&map, 0, sizeof map);
    map.beam = p->rdBeam;
    map.nBins = p->nBins;
    map.nfft = cfg->nfft;
    map.dR_m = p->dR;
    map.power = p->rd;
    map.vel_mps = p->rdVel;
    map.notchBins = cfg->notchBins;
    map.zeroBin = cfg->zeroBin;
    map.noise = m ? processor_select_kth(p->rdWork, (int)m, (int)(m / 2)) : 0.0;
    map.time_s = pri->time_s;
    map.az_deg = pri->az_deg;
    map.dets = p->rdDets;
    map.nDets = p->nRdDets;
    cb->rdmap(cb->user, &map);
}

/* Receiver lifecycle, PRI validation and scheduling of the processing stages.
 * Each stage implementation is listed in the receiver README. */
#include "processor_internal.h"

#include <math.h>
#include <string.h>

Processor *processor_create(const RadarConfig *cfg, int nThreads)
{
#ifdef _OPENMP
    /* default: OMP_NUM_THREADS if set, else all cores but one (left for the
     * UDP input and the display) */
    if (nThreads <= 0 && !getenv("OMP_NUM_THREADS"))
        nThreads = omp_get_num_procs() > 1 ? omp_get_num_procs() - 1 : 1;
    if (nThreads > 0) omp_set_num_threads(nThreads);
    const int nThr = omp_get_max_threads();
#else
    (void)nThreads;
    const int nThr = 1;
#endif
    Processor *p = proc_calloc(1, sizeof(Processor));
    p->cfg = cfg;
    p->cfgPrf = (RadarConfig *)cfg;
    /* buffers for the largest slow-time settings (scan-to-scan PRF: of all PRFs) */
    const int maxCpi = cfg->maxCpiLength > cfg->cpiLength ? cfg->maxCpiLength : cfg->cpiLength;
    const int maxNfft = cfg->maxNfft > cfg->nfft ? cfg->maxNfft : cfg->nfft;
    const int maxRow = cfg->maxRowLen > cfg->rowLen ? cfg->maxRowLen : cfg->rowLen;
    const int nChUsed = cfg->nChUsed > 0 ? cfg->nChUsed : cfg->nCh;
    /* (the NUFFT has no FFT: its length need not be a power of 2) */
    if (cfg->dopplerMethod == DOPPLER_FFT && fft_plan_init(&p->fft, cfg->nfft) != 0) { free(p); return NULL; }
    int maxNR = 1;
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        PulseState *s = &p->ps[i];
        if (!pc->enable) {                         /* sub-pulse switched off: no buffers */
            extractor_init(&s->ext, cfg->nfft, cfg->minHits, cfg->maxMiss);
            s->ext.nBeams = cfg->nBeams;
            cpi_buffer_reset(s);
            continue;
        }
        s->nR = pc->nRange;
        if (s->nR > maxNR) maxNR = s->nR;
        s->nCells = cfg->nBeams * s->nR;
        s->nSlots = cfg->cpiLength;
        s->nMtiOut = cfg->cpiLength - cfg->nMtiTaps + 1;
        size_t seg = (size_t)nChUsed * pc->ddcCount, beams = (size_t)cfg->nBeams * pc->ddcCount;
        size_t cpi = (size_t)cfg->nMtiTaps * s->nCells;
        size_t mti = (size_t)(s->nCells + CELL_BLOCK - 1) / CELL_BLOCK * (maxCpi - cfg->nMtiTaps + 1) * CELL_BLOCK;
        size_t map = (size_t)maxNfft * s->nCells;
        s->segRe = proc_calloc(seg, sizeof(float));   s->segIm = proc_calloc(seg, sizeof(float));
        s->beamRe = proc_calloc(beams, sizeof(float)); s->beamIm = proc_calloc(beams, sizeof(float));
        if (cfg->slcEnabled) {
            const size_t nA = (size_t)cfg->nAux, nB = (size_t)cfg->nBeams, K = (size_t)pc->slcTrainCount;
            s->auxRe = proc_calloc(nA * pc->ddcCount, sizeof(float)); s->auxIm = proc_calloc(nA * pc->ddcCount, sizeof(float));
            s->trMRe = proc_calloc(nB * K, sizeof(double)); s->trMIm = proc_calloc(nB * K, sizeof(double));
            s->trARe = proc_calloc(nA * K, sizeof(double)); s->trAIm = proc_calloc(nA * K, sizeof(double));
            s->slcRRe = proc_calloc(nA * nA, sizeof(double)); s->slcRIm = proc_calloc(nA * nA, sizeof(double));
            s->slcZRe = proc_calloc(nA * nB, sizeof(double)); s->slcZIm = proc_calloc(nA * nB, sizeof(double));
            s->slcWRe = proc_calloc(nA * nB, sizeof(double)); s->slcWIm = proc_calloc(nA * nB, sizeof(double));
            s->slcARe = proc_calloc(nA * nB, sizeof(double)); s->slcAIm = proc_calloc(nA * nB, sizeof(double));
            s->slcBRe = proc_calloc(nA * nB, sizeof(double)); s->slcBIm = proc_calloc(nA * nB, sizeof(double));
            s->slcX = proc_calloc(2 * nA, sizeof(double)); s->slcD = proc_calloc(2 * nB, sizeof(double));
            s->slcEig = proc_calloc(4 * nA * nA, sizeof(double));
        }
        s->cpiRe = proc_calloc(cpi, sizeof(float));   s->cpiIm = proc_calloc(cpi, sizeof(float));
        s->mtiRe = proc_calloc(mti, sizeof(float));   s->mtiIm = proc_calloc(mti, sizeof(float));
        s->staggerIdx = proc_calloc(maxCpi, sizeof(uint16_t));
        s->mtiStag = proc_calloc(maxCpi, sizeof(uint16_t));
        s->az = proc_calloc(maxCpi, sizeof(double));
        s->time = proc_calloc(maxCpi, sizeof(double));
        s->P = proc_calloc(map, sizeof(float));
        s->det = proc_calloc(map, 1);
        s->cntInv = proc_calloc(s->nR, sizeof(double));
        s->thrScale = proc_calloc(s->nR, sizeof(double));
        if (cfg->snrNoise == SNR_NOISE_CFAR) s->noise = proc_calloc(map, sizeof(float));
        s->lastMedian = proc_calloc(cfg->nBeams, sizeof(float));
        cfar_init_counts(cfg, s);
        if (mf_fft_init(cfg, pc, s, cfg->mfMethod, nThr) != 0) { processor_destroy(p); return NULL; }
        extractor_init(&s->ext, cfg->nfft, cfg->minHits, cfg->maxMiss);
        s->ext.nBeams = cfg->nBeams;
        cpi_buffer_reset(s);
    }
    p->nThreads = nThr;
    p->slowTime = 1;
    p->alpha = cfg->alpha;
    if (cfg->trackEnable) tracker_init(&p->tracker, cfg);
    init_pack_table();
    p->scratch = proc_calloc(nThr, sizeof(Scratch));
    for (int t = 0; t < nThr; t++) {
        Scratch *w = &p->scratch[t];
        const int nx = maxNfft > maxRow ? maxNfft : maxRow;      /* (NUFFT: the rowLen inputs of a block) */
        w->xr = proc_calloc((size_t)nx * CELL_BLOCK, sizeof(float));
        w->xi = proc_calloc((size_t)nx * CELL_BLOCK, sizeof(float));
        cfar_scratch_init(&w->cfar, cfg, maxNR);
        w->work = proc_calloc((size_t)maxNfft * maxNR + 8, sizeof(float));
    }
    return p;
}

void processor_destroy(Processor *p)
{
    if (!p) return;
    free(p->profile); free(p->profileCnt); free(p->rd); free(p->rdVel); free(p->rdDets); free(p->rdWork);
    for (int i = 0; i < N_PULSES; i++) {
        PulseState *s = &p->ps[i];
        free(s->segRe); free(s->segIm); free(s->beamRe); free(s->beamIm);
        free(s->auxRe); free(s->auxIm); free(s->trMRe); free(s->trMIm); free(s->trARe); free(s->trAIm);
        free(s->slcRRe); free(s->slcRIm); free(s->slcZRe); free(s->slcZIm); free(s->slcWRe); free(s->slcWIm);
        free(s->slcARe); free(s->slcAIm); free(s->slcBRe); free(s->slcBIm);
        free(s->slcX); free(s->slcD); free(s->slcEig);
        free(s->cpiRe); free(s->cpiIm); free(s->mtiRe); free(s->mtiIm);
        free(s->staggerIdx); free(s->mtiStag); free(s->az); free(s->time);
        free(s->P); free(s->det); free(s->cntInv); free(s->thrScale); free(s->noise); free(s->lastMedian);
        free(s->pkM); free(s->pkA); free(s->pkE); free(s->pkD); free(s->pkT);
        free(s->reports);
        if (s->mfFft) fft_plan_free(&s->mfPlan);
        if (s->mfFour) { fft_plan_free(&s->mfPlanQ); fft_plan_free(&s->mfPlan8); }
        free(s->mfHre); free(s->mfHim); free(s->mfTwRe); free(s->mfTwIm); free(s->mfWork);
        extractor_free(&s->ext);
    }
    for (int t = 0; p->scratch && t < p->nThreads; t++) {
        Scratch *w = &p->scratch[t];
        free(w->xr); free(w->xi); free(w->work);
        cfar_scratch_free(&w->cfar);
    }
    free(p->scratch);
    fft_plan_free(&p->fft);
    tracker_free(&p->tracker);
    free(p->newPlots);
    free(p->trackPlots);
    free(p);
}

void processor_set_slow_time(Processor *p, int enabled)
{
    p->slowTime = enabled;
}

void processor_set_log(Processor *p, int mask)
{
    p->logMask = mask & (LOG_DBF | LOG_SLC | LOG_MF | LOG_DOPPLER);
}

void processor_set_cfar_alpha(Processor *p, double alpha)
{
    if (alpha > 0.0 && isfinite(alpha)) p->alpha = alpha;
}

double processor_cfar_alpha(const Processor *p)
{
    return p->alpha;
}

void processor_timing(const Processor *p, ProcessorTiming *t)
{
    *t = p->timing;
}

int processor_threads(const Processor *p)
{
    return p->nThreads;
}

/* ------------------------------------------------------------------------ */
/* One PRI                                                                  */
/* ------------------------------------------------------------------------ */

static double elapsed_ms(double *t0)
{
    double t = now_seconds(), dt = (t - *t0) * 1e3;
    *t0 = t;
    return dt;
}

/* The extractor's plots are collected first, so that the extractor and the
 * tracker are timed separately (and the output callbacks are not timed). */
static void plot_callback(void *user, const Detection *plot)
{
    Processor *p = user;
    if (p->nNewPlots == p->capNewPlots) {
        p->capNewPlots = p->capNewPlots ? 2 * p->capNewPlots : 16;
        p->newPlots = proc_realloc(p->newPlots, sizeof(Detection) * p->capNewPlots);
    }
    p->newPlots[p->nNewPlots++] = *plot;
}

/* Track updates and deletions are collected first and reported after the
 * tracker has run (so that the output callbacks are not timed). */
typedef struct {
    Track *list;
    int n, cap;
} TrackReports;

static void track_callback(void *user, const Track *track, int plotIndex)
{
    TrackReports *r = user;
    (void)plotIndex;
    if (r->n == r->cap) {
        r->cap = r->cap ? 2 * r->cap : 16;
        r->list = proc_realloc(r->list, sizeof(Track) * r->cap);
    }
    r->list[r->n++] = *track;
}

/* Hands the collected plots to the output; the tracker gets them at the end of the PRI. */
static void plots_out(Processor *p, const ProcessorCallbacks *cb)
{
    for (int k = 0; k < p->nNewPlots; k++) {
        const Detection *plot = &p->newPlots[k];
        if (cb && cb->plot) cb->plot(cb->user, plot);
        if (!p->cfg->trackEnable) continue;
        if (p->nTrackPlots == p->capTrackPlots) {
            p->capTrackPlots = p->capTrackPlots ? 2 * p->capTrackPlots : 16;
            p->trackPlots = proc_realloc(p->trackPlots, sizeof(Detection) * p->capTrackPlots);
        }
        p->trackPlots[p->nTrackPlots++] = *plot;
    }
    p->nNewPlots = 0;
}

/* The plots of this PRI to the tracker (GNN batches, deletions) at time t_s, as
 * track_update.m once per PRI; flush = 1 associates the open batch. */
static void run_tracker(Processor *p, double t_s, int flush, const ProcessorCallbacks *cb)
{
    if (!p->cfg->trackEnable) return;
    TrackReports r = { NULL, 0, 0 };
    double t0 = now_seconds();
    tracker_update(&p->tracker, p->trackPlots, p->nTrackPlots, t_s, flush, track_callback, &r);
    p->timing.track_ms += elapsed_ms(&t0);
    p->nTrackPlots = 0;
    for (int i = 0; i < r.n; i++)
        if (cb && cb->track) cb->track(cb->user, &r.list[i]);
    free(r.list);
}

/* The PRI has the layout of the configuration: channels, DDC samples, the
 * receive intervals of both sub-pulses (header, 1-based) and the DDC rate.
 * (None of them depends on the CFAR or on the pulse compression filters.)
 * The record may be longer than the configured one (e.g. the DDC output of
 * the whole PRI): only the receive windows of the sub-pulses are taken out
 * of it, at their samples from the PRI start. It may also have more channels
 * than those used (cfg.rx.nChannelsUsedDbf, nChannelsUsedSlc): it needs the
 * channels up to the largest element number used. */
static int pri_matches_config(const RadarConfig *cfg, const PriData *pri)
{
    const int need = cfg->chNeeded > 0 ? cfg->chNeeded : cfg->nCh;
    if (pri->nCh < need || pri->nSamples < cfg->nDdc || !pri->iq) return 0;
    const uint32_t *idx[N_PULSES] = { pri->shortIdx, pri->longIdx };
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        if (!pc->enable) continue;                          /* (its header interval is [0 0]) */
        if (idx[i][0] != (uint32_t)pc->rxFirst + 1 || idx[i][1] != (uint32_t)(pc->rxFirst + pc->rxCount))
            return 0;
    }
    return fabs(pri->ddcFs_Hz - cfg->fsOut_Hz) <= 1e-3 * cfg->fsOut_Hz;
}

void processor_flush(Processor *p, const ProcessorCallbacks *cb)
{
    for (int i = 0; i < N_PULSES; i++) {
        extractor_update(&p->ps[i].ext, NULL, 0, 1, plot_callback, p);
        plots_out(p, cb);
    }
    if (p->havePulse) run_tracker(p, p->lastTime, 1, cb);   /* at the newest PRI, as run_radar_sim.m */
}

/* Scan-to-scan PRF: the PRI belongs to a scan with PRF k (0-based, the stagger index of its
 * header). The open plots are closed and go to the tracker (with the Doppler axis of their
 * scan), then the slow-time settings of PRF k are put into the configuration and the CPI
 * buffers, extractors and the FFT start anew (as run_radar_sim.m). */
static void prf_switch(Processor *p, int k, const ProcessorCallbacks *cb)
{
    if (p->havePulse) {
        for (int i = 0; i < N_PULSES; i++) {
            extractor_update(&p->ps[i].ext, NULL, 0, 1, plot_callback, p);
            plots_out(p, cb);
        }
        if (p->nTrackPlots > 0) run_tracker(p, p->lastTime, 0, cb);   /* (as run_radar_sim.m) */
    }
    if (!radar_config_select_prf(p->cfgPrf, k)) return;
    const RadarConfig *cfg = p->cfg;
    fft_plan_free(&p->fft);
    if (cfg->dopplerMethod == DOPPLER_FFT && fft_plan_init(&p->fft, cfg->nfft) != 0) abort();
    for (int i = 0; i < N_PULSES; i++) {
        PulseState *s = &p->ps[i];
        extractor_free(&s->ext);
        extractor_init(&s->ext, cfg->nfft, cfg->minHits, cfg->maxMiss);
        s->ext.nBeams = cfg->nBeams;
        if (!cfg->pulse[i].enable) continue;
        s->nSlots = cfg->cpiLength;
        s->nMtiOut = cfg->cpiLength - cfg->nMtiTaps + 1;
        cpi_buffer_reset(s);
    }
    if (cfg->trackEnable) tracker_set_vr_span(&p->tracker, cfg);
    rdmap_prf_changed(p);
    p->timing.prfSwitches++;
}

void processor_process(Processor *p, const PriData *pri, const ProcessorCallbacks *cb)
{
    const RadarConfig *cfg = p->cfg;
    if (!pri_matches_config(cfg, pri)) {                    /* made for another configuration */
        p->timing.rejected++;
        return;
    }

    /* a gap in the pulse numbers (skipped sector): restart the CPIs, close the plots */
    if (p->havePulse && pri->pulseNumber != p->lastPulse + 1) {
        processor_flush(p, cb);
        for (int i = 0; i < N_PULSES; i++) cpi_buffer_reset(&p->ps[i]);
    }
    if (cfg->scan2scan) {                                   /* scan-to-scan PRF: the PRF of this scan */
        const int k = pri->staggerIndex >= 1 && pri->staggerIndex <= cfg->nPrf ? pri->staggerIndex - 1 : 0;
        if (k != cfg->prfIndex) prf_switch(p, k, cb);
    }
    p->lastPulse = pri->pulseNumber;
    p->lastTime = pri->time_s;
    p->havePulse = 1;
    if (cb && cb->antenna) cb->antenna(cb->user, pri->az_deg, pri->time_s, pri->scanNumber);

    const double tStart = now_seconds();
    const int profileKind = cb && cb->profile ? p->profileKind : PROFILE_NONE;
    const int rdOn = cb && cb->rdmap && p->rdEnabled && pri->time_s - p->rdLast >= p->rdPeriod_s;
    if (profileKind) {
        memset(p->profile, 0, sizeof(float) * (size_t)p->nBins);
        memset(p->profileCnt, 0, sizeof(int) * (size_t)p->nBins);
    }
    if (p->rdEnabled && pri->time_s < p->rdLast) p->rdLast = -1e300;   /* a new run */
    if (rdOn) { p->nRdDets = 0; p->rdUpdated = 0; }
    const int firstOn = cfg->pulse[PULSE_SHORT].enable ? PULSE_SHORT : PULSE_LONG;   /* counts the updates */
    for (int i = 0; i < N_PULSES; i++) {
        const PulseConfig *pc = &cfg->pulse[i];
        PulseState *s = &p->ps[i];
        if (!pc->enable) continue;                          /* sub-pulse switched off */
        double t0 = now_seconds();

        /* fast time into the next CPI slot */
        load_segment(cfg, pc, pri, s);
        beamform(cfg, pc, s, p->nThreads);
        p->timing.dbf_ms += elapsed_ms(&t0);
        if ((p->logMask & LOG_DBF) && cb && cb->block) {      /* (logging: the application's time) */
            const BlockOutput o = { LOG_DBF, i, cfg->nBeams, pc->nRange, cfg->nfft, pc->ddcCount, pc->ddcFirst,
                                    s->beamRe, s->beamIm, pri->pulseNumber };
            cb->block(cb->user, &o);
            t0 = now_seconds();
        }
        if (profileKind == PROFILE_DBF) profile_add(p, i, s);   /* (the display's time, not the block's) */
        if (cfg->slcEnabled) {                                  /* sidelobe canceller (in place) */
            t0 = now_seconds();
            sidelobe_canceller(cfg, pc, pri, s, p->nThreads);
            p->timing.slc_ms += elapsed_ms(&t0);
            if ((p->logMask & LOG_SLC) && cb && cb->block) {
                const BlockOutput o = { LOG_SLC, i, cfg->nBeams, pc->nRange, cfg->nfft, pc->ddcCount, pc->ddcFirst,
                                        s->beamRe, s->beamIm, pri->pulseNumber };
                cb->block(cb->user, &o);
                t0 = now_seconds();
            }
            if (profileKind == PROFILE_SLC) profile_add(p, i, s);
        }
        const int slot = cpi_buffer_begin(s, cfg->nMtiTaps);
        matched_filter(cfg, pc, s, s->mfLast, p->nThreads);
        if (profileKind == PROFILE_MF) profile_add(p, i, s);
        if ((p->logMask & LOG_MF) && cb && cb->block) {
            const double tl = now_seconds();
            const BlockOutput o = { LOG_MF, i, cfg->nBeams, pc->nRange, cfg->nfft, 0, 0,
                                    s->cpiRe + (size_t)s->mfLast * s->nCells, s->cpiIm + (size_t)s->mfLast * s->nCells,
                                    pri->pulseNumber };
            cb->block(cb->user, &o);
            t0 += now_seconds() - tl;                       /* (not the time of the block) */
        }
        cpi_buffer_push(s, cfg, pri, slot);
        p->timing.mf_ms += elapsed_ms(&t0);

        /* slow time */
        if (s->count >= cfg->nMtiTaps) mti_newest(cfg, pc, s);
        if (s->count == s->nSlots) s->sinceUpdate++;
        p->timing.mti_ms += elapsed_ms(&t0);
        if (s->count < s->nSlots || s->nPushed != 0) continue;
        if (!p->slowTime) {                                 /* overload: this update is skipped */
            if (i == firstOn) p->timing.skipped++;
            continue;
        }
        /* a plot needs detections over the same azimuth extent also when
         * updates were skipped: an update counts for the pulses it covers */
        s->ext.weight = s->sinceUpdate / cfg->updateEvery > 1 ? s->sinceUpdate / cfg->updateEvery : 1;
        s->sinceUpdate = 0;
        doppler(p, pc, s);
        p->timing.doppler_ms += elapsed_ms(&t0);
        if ((p->logMask & LOG_DOPPLER) && cb && cb->block) {
            const BlockOutput o = { LOG_DOPPLER, i, cfg->nBeams, pc->nRange, cfg->nfft, 0, 0, s->P, NULL, pri->pulseNumber };
            cb->block(cb->user, &o);
            t0 = now_seconds();
        }
        processor_cfar(p, s, pc);
        prune_and_report(p, i, s);
        p->timing.cfar_ms += elapsed_ms(&t0);
        peak_interpolation(p, i, s);
        p->timing.interp_ms += elapsed_ms(&t0);
        if (rdOn) rdmap_add(p, i, s);
        if (cb && cb->detections) cb->detections(cb->user, i, s->reports, s->nReports);
        t0 = now_seconds();
        extractor_update(&s->ext, s->reports, s->nReports, 0, plot_callback, p);
        p->timing.extract_ms += elapsed_ms(&t0);
        plots_out(p, cb);
        if (i == firstOn) p->timing.updates++;
    }
    if (profileKind) profile_post(p, profileKind, pri, cb);
    if (rdOn && p->rdUpdated) {
        rdmap_post(p, pri, cb);
        p->rdLast = pri->time_s;
    }
    run_tracker(p, pri->time_s, 0, cb);                    /* plots of this PRI, missed scans, deletion */
    p->timing.total_ms += (now_seconds() - tStart) * 1e3;
    p->timing.pris++;
}

/*
 * params_mat.c - radar_parameters.mat of a log (logger.c, written when the log starts): every
 * setting of the processing chain in use at that moment, as one MATLAB struct 'params':
 *
 *   params.info        start time, configuration file, log blocks, last PRI before the start
 *   params.antenna     channels (all / DBF elements chDbf), beams, DBF matrix W [nBeams x nChUsed]
 *                      (complex), beam directions, tilt, height, rotation, ranges
 *   params.slc         sidelobe canceller: on/off, auxiliary elements, algorithm and its settings,
 *                      training samples of each sub-pulse
 *   params.ddc         DDC record: samples, rate, first sample
 *   params.short, params.long   each sub-pulse: carrier, waveform, receive interval, pulse
 *                      compression filter (complex taps), range axis, CFAR cells tested, ...
 *   params.slowTime    CPI, stagger, clutter filter, Doppler processing (FFT / NUFFT with its
 *                      matrices [nfft x rowLen x nStagger]), window, noise gains, axes, notch
 *   params.cfar        detector, window, threshold factor in use (the PPI can change it) and
 *                      its Pfa, the design value of the configuration, CE table, ...
 *   params.report      peak interpolation, plot extraction
 *   params.tracker     tracker settings (EKF / PMBM / TOMHT)
 *   params.runtime     options of the receiver: matched filter method, threads, update rate,
 *                      Doppler/CFAR skipped (behind the stream), range-Doppler rate of the PPI
 *
 * Indices are 1-based (MATLAB), arrays are columns as in the exported configuration.
 */
#if defined(_WIN32) && !defined(_CRT_SECURE_NO_WARNINGS)
#  define _CRT_SECURE_NO_WARNINGS
#endif

#include "radar_rx.h"
#include "mat_write.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SC(st, name, v) mat_struct_scalar(st, name, (double)(v))

static void col_double(MatStruct *s, const char *name, const double *v, int n)
{
    if (v && n > 0) mat_struct_double(s, name, v, (size_t)n, 1);
    else mat_struct_double(s, name, NULL, 0, 0);
}

static void col_single(MatStruct *s, const char *name, const float *re, const float *im, int n)
{
    const size_t d[2] = { (size_t)(re && n > 0 ? n : 0), re && n > 0 ? 1u : 0u };
    mat_struct_single(s, name, re && n > 0 ? re : NULL, re && n > 0 ? im : NULL, d, 2);
}

/* row-major [rows][cols] -> MATLAB [rows x cols] */
static void matrix_double(MatStruct *s, const char *name, const double *v, int rows, int cols)
{
    if (!v || rows <= 0 || cols <= 0) { mat_struct_double(s, name, NULL, 0, 0); return; }
    double *t = malloc(sizeof(double) * (size_t)rows * cols);
    if (!t) abort();
    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++) t[r + (size_t)rows * c] = v[(size_t)r * cols + c];
    mat_struct_double(s, name, t, (size_t)rows, (size_t)cols);
    free(t);
}

/* row-major [d3][rows][cols] (complex) -> MATLAB [rows x cols x d3] */
static void array_single(MatStruct *s, const char *name, const float *re, const float *im, int d3, int rows, int cols)
{
    const size_t n = (size_t)d3 * rows * cols;
    if (!re || !n) { const size_t z[2] = { 0, 0 }; mat_struct_single(s, name, NULL, NULL, z, 2); return; }
    float *tr = malloc(sizeof(float) * n), *ti = im ? malloc(sizeof(float) * n) : NULL;
    if (!tr || (im && !ti)) abort();
    for (int k = 0; k < d3; k++)
        for (int r = 0; r < rows; r++)
            for (int c = 0; c < cols; c++) {
                const size_t src = ((size_t)k * rows + r) * cols + c, dst = r + (size_t)rows * (c + (size_t)cols * k);
                tr[dst] = re[src];
                if (ti) ti[dst] = im[src];
            }
    const size_t dims[3] = { (size_t)rows, (size_t)cols, (size_t)d3 };
    mat_struct_single(s, name, tr, ti, dims, 3);
    free(tr);
    free(ti);
}

static MatStruct *pulse_struct(const RadarConfig *cfg, int p)
{
    const PulseConfig *q = &cfg->pulse[p];
    MatStruct *s = mat_struct_new();
    SC(s, "enable", q->enable);
    SC(s, "fc_Hz", q->fc_Hz);
    SC(s, "lambda_m", q->lambda_m);
    mat_struct_char(s, "waveform", q->waveform == WAVE_PHASE_CODE ? "phase code" : "lfm");
    mat_struct_char(s, "compression", q->compression == COMPRESSION_MISMATCH ? "mismatch" : "matched");
    SC(s, "nChips", q->nChips);
    SC(s, "codeIndex", cfg->codeIndex[p]);
    SC(s, "rxFirst", q->rxFirst + 1);                    /* receive interval: DDC samples (1-based) */
    SC(s, "rxCount", q->rxCount);
    SC(s, "ddcFirst", q->ddcFirst + 1);                  /* pulse compression input */
    SC(s, "ddcCount", q->ddcCount);
    SC(s, "nTaps", q->nTaps);
    col_single(s, "taps", q->tapsRe, q->tapsIm, q->nTaps);
    SC(s, "mfShift", q->mfShift);
    SC(s, "nRange", q->nRange);
    col_double(s, "range_m", q->range_m, q->nRange);
    col_double(s, "delay_s", q->delay_s, q->nRange);
    SC(s, "testFirst", q->testFirst + 1);                /* range cells tested by the CFAR */
    SC(s, "testCount", q->testCount);
    SC(s, "mtiCells", q->mtiCells);
    SC(s, "vrSpan_mps", q->vrSpan_mps);
    col_double(s, "beamU", q->beamU, q->beamU ? cfg->nBeams : 0);
    if (cfg->nfft > 0 && cfg->dopFreq_Hz && q->lambda_m > 0) {   /* velocity of the Doppler bins */
        double *v = malloc(sizeof(double) * (size_t)cfg->nfft);
        if (!v) abort();
        for (int k = 0; k < cfg->nfft; k++) v[k] = -cfg->dopFreq_Hz[k] * q->lambda_m / 2.0;   /* + = away */
        col_double(s, "velocity_mps", v, cfg->nfft);
        free(v);
    }
    return s;
}

int radar_parameters_write(const char *path, const RadarConfig *cfg, const RuntimeParams *rt)
{
    MatWriter *w = mat_write_open(path, 0);
    if (!w) return -1;
    MatStruct *P = mat_struct_new();

    MatStruct *s = mat_struct_new();                     /* info */
    mat_struct_char(s, "started", rt->started ? rt->started : "");
    mat_struct_char(s, "configFile", rt->configPath ? rt->configPath : "");
    mat_struct_char(s, "logBlocks", rt->logBlocks ? rt->logBlocks : "");
    SC(s, "screenRecord", rt->screenRecord);
    SC(s, "lastPulseNumber", rt->havePri ? (double)rt->pulseNumber : NAN);
    SC(s, "lastPriTime_s", rt->havePri ? rt->time_s : NAN);
    SC(s, "lastAntennaAz_deg", rt->havePri ? rt->az_deg : NAN);
    SC(s, "radarMode", cfg->radarMode);
    mat_struct_struct(P, "info", s);

    s = mat_struct_new();                                /* antenna and DBF */
    SC(s, "nCh", cfg->nCh);
    SC(s, "nChUsed", cfg->nChUsed);
    {                                                    /* element numbers (1-based) */
        double *v = malloc(sizeof(double) * (size_t)(cfg->nChUsed > 0 ? cfg->nChUsed : 1));
        if (!v) abort();
        for (int k = 0; k < cfg->nChUsed; k++) v[k] = (cfg->chDbf ? cfg->chDbf[k] : k) + 1;
        mat_struct_double(s, "chDbf", v, 1, (size_t)cfg->nChUsed);
        free(v);
    }
    SC(s, "nBeams", cfg->nBeams);
    SC(s, "dbfEnabled", cfg->dbfEnabled);
    if (cfg->wRe) {                                      /* [nBeams x nChUsed] */
        const size_t d[2] = { (size_t)cfg->nBeams, (size_t)cfg->nChUsed };
        float *re = malloc(sizeof(float) * d[0] * d[1]), *im = malloc(sizeof(float) * d[0] * d[1]);
        if (!re || !im) abort();
        for (size_t b = 0; b < d[0]; b++)
            for (size_t c = 0; c < d[1]; c++) {
                re[b + d[0] * c] = cfg->wRe[b * d[1] + c];
                im[b + d[0] * c] = cfg->wIm ? cfg->wIm[b * d[1] + c] : 0.0f;
            }
        mat_struct_single(s, "W", re, im, d, 2);
        free(re);
        free(im);
    }
    col_double(s, "beamU", cfg->beamU, cfg->beamU ? cfg->nBeams : 0);
    if (cfg->beamU) {
        double *el = malloc(sizeof(double) * (size_t)cfg->nBeams);
        if (!el) abort();
        for (int b = 0; b < cfg->nBeams; b++) el[b] = asin(cfg->beamU[b]) * 57.29577951308232 + cfg->tilt_deg;
        col_double(s, "beamEl_deg", el, cfg->nBeams);
        free(el);
    }
    SC(s, "tilt_deg", cfg->tilt_deg);
    SC(s, "height_m", cfg->height_m);
    SC(s, "rate_dps", cfg->rate_dps);
    SC(s, "azBeamwidth_deg", cfg->azBeamwidth_deg);
    SC(s, "maxRange_m", cfg->maxRange_m);
    SC(s, "handover_m", cfg->handover_m);
    mat_struct_struct(P, "antenna", s);

    s = mat_struct_new();                                /* sidelobe canceller */
    SC(s, "enabled", cfg->slcEnabled);
    if (cfg->slcEnabled) {
        static const char *algo[3] = { "qrd-smi", "qrd-rls", "nlms" };
        double *v = malloc(sizeof(double) * (size_t)(cfg->nAux > 0 ? cfg->nAux : 1));
        if (!v) abort();
        for (int k = 0; k < cfg->nAux; k++) v[k] = cfg->chSlc[k] + 1;
        mat_struct_double(s, "chSlc", v, 1, (size_t)cfg->nAux);
        free(v);
        mat_struct_char(s, "algorithm", algo[cfg->slcMethod]);
        SC(s, "smiLoading_dB", cfg->slcSmiLoading_dB);
        SC(s, "rlsLambda", cfg->slcRlsLambda);
        SC(s, "rlsDelta_dB", cfg->slcRlsDelta_dB);
        SC(s, "nlmsMu", cfg->slcNlmsMu);
        SC(s, "nlmsEps", cfg->slcNlmsEps);
        for (int p = 0; p < N_PULSES; p++) {             /* training samples: 1-based DDC samples */
            char nm[40];
            snprintf(nm, sizeof nm, "%sTrainFirstSample", p ? "long" : "short");
            SC(s, nm, cfg->pulse[p].enable ? cfg->pulse[p].slcTrainFirst + 1 : NAN);
            snprintf(nm, sizeof nm, "%sTrainCount", p ? "long" : "short");
            SC(s, nm, cfg->pulse[p].enable ? cfg->pulse[p].slcTrainCount : 0);
        }
    }
    mat_struct_struct(P, "slc", s);

    s = mat_struct_new();                                /* DDC */
    SC(s, "nDdc", cfg->nDdc);
    SC(s, "fsOut_Hz", cfg->fsOut_Hz);
    SC(s, "manual", cfg->ddcManual);
    SC(s, "t0_s", cfg->ddcT0_s);
    mat_struct_struct(P, "ddc", s);

    mat_struct_struct(P, "short", pulse_struct(cfg, PULSE_SHORT));
    mat_struct_struct(P, "long", pulse_struct(cfg, PULSE_LONG));

    s = mat_struct_new();                                /* slow time: CPI, MTI, Doppler */
    SC(s, "cpiLength", cfg->cpiLength);
    SC(s, "staggered", cfg->staggered);
    SC(s, "nStagger", cfg->nStagger);
    SC(s, "staggerSet", cfg->staggerSet);
    SC(s, "meanPri_s", cfg->meanPri_s);
    SC(s, "scan2scan", cfg->scan2scan);
    SC(s, "nPrf", cfg->nPrf);
    SC(s, "prfIndex", cfg->prfIndex + 1);
    if (cfg->scan2scan && cfg->prf) {
        double *pri = malloc(sizeof(double) * (size_t)cfg->nPrf);
        if (!pri) abort();
        for (int k = 0; k < cfg->nPrf; k++) pri[k] = cfg->prf[k].pri_s;
        col_double(s, "prfPri_s", pri, cfg->nPrf);
        free(pri);
    }
    mat_struct_char(s, "mtiType", cfg->mtiType == MTI_FIR ? "fir" : "canceller");
    matrix_double(s, "mtiCoeffs", cfg->mtiCoeffs, cfg->nStagger, cfg->nMtiTaps);
    SC(s, "mtiMaxRange_m", cfg->mtiMaxRange_m);
    mat_struct_char(s, "dopplerMethod", cfg->dopplerMethod == DOPPLER_NUFFT ? "nufft" : "fft");
    SC(s, "nRows", cfg->nRows);
    SC(s, "rowLen", cfg->rowLen);
    SC(s, "nfft", cfg->nfft);
    col_single(s, "dopWindow", cfg->dopWindow, NULL, cfg->dopplerMethod == DOPPLER_NUFFT ? 0 : cfg->rowLen);
    col_single(s, "noiseGain", cfg->noiseGain, NULL, cfg->nfft);
    col_single(s, "noiseGainBypass", cfg->noiseGainBypass, NULL, cfg->nfft);
    col_double(s, "dopFreq_Hz", cfg->dopFreq_Hz, cfg->nfft);
    if (cfg->dopplerMethod == DOPPLER_NUFFT) {
        array_single(s, "nufftA", cfg->nufftARe, cfg->nufftAIm, cfg->nStagger, cfg->nfft, cfg->rowLen);
        if (cfg->nufftBypassGain) {                      /* [nStagger][nfft] -> [nfft x nStagger] */
            const size_t d[2] = { (size_t)cfg->nfft, (size_t)cfg->nStagger };
            float *g = malloc(sizeof(float) * d[0] * d[1]);
            if (!g) abort();
            for (size_t p = 0; p < d[1]; p++)
                for (size_t k = 0; k < d[0]; k++) g[k + d[0] * p] = cfg->nufftBypassGain[p * d[0] + k];
            mat_struct_single(s, "nufftBypassGain", g, NULL, d, 2);
            free(g);
        }
    }
    SC(s, "notchBins", cfg->notchBins);
    SC(s, "zeroBin", cfg->zeroBin >= 0 ? cfg->zeroBin + 1 : 0);     /* 0: no zero-Doppler bin */
    SC(s, "notchMaxRange_m", cfg->notchMaxRange_m);
    SC(s, "medianToMean", cfg->medianToMean);
    SC(s, "updateEvery", cfg->updateEvery);
    mat_struct_struct(P, "slowTime", s);

    s = mat_struct_new();                                /* CFAR */
    const double alpha = rt->alpha > 0 ? rt->alpha : cfg->alpha;
    mat_struct_char(s, "type", cfar_type_name(cfg->cfarType));
    mat_struct_char(s, "method", cfar_method_name(cfg->cfarMethod));
    SC(s, "guardRange", cfg->guardRange);
    SC(s, "trainRange", cfg->trainRange);
    SC(s, "guardDoppler", cfg->guardDoppler);
    SC(s, "trainDoppler", cfg->trainDoppler);
    SC(s, "nTrain", cfg->nTrain);
    SC(s, "nHalf", cfg->cfarNHalf);
    SC(s, "osK", cfg->cfarOsK);
    SC(s, "ceKeep", cfg->cfarCeKeep);
    SC(s, "alpha", alpha);                               /* in use (the PPI can change it) */
    SC(s, "threshold_dB", 10.0 * log10(alpha));
    SC(s, "pfa", cfar_pfa(cfg, alpha));
    SC(s, "alphaDesign", cfg->alpha);                    /* from cfg.cfar.pfa */
    SC(s, "thresholdDesign_dB", 10.0 * log10(cfg->alpha));
    SC(s, "pfaDesign", cfar_pfa(cfg, cfg->alpha));
    SC(s, "alphaMin", cfg->alphaMin);
    SC(s, "nLooks", cfg->nLooks);
    SC(s, "noiseScale", cfg->cfarNoiseScale);
    SC(s, "sidelobeMargin_dB", cfg->sidelobeMargin_dB);
    if (cfg->nCeTable > 0) {
        col_double(s, "ceLogZ", cfg->ceLogZ, cfg->nCeTable);
        col_double(s, "ceLogG", cfg->ceLogG, cfg->nCeTable);
    }
    mat_struct_struct(P, "cfar", s);

    s = mat_struct_new();                                /* reports and plots */
    mat_struct_char(s, "snrNoise", cfg->snrNoise == SNR_NOISE_CFAR ? "cfar" : "median");
    SC(s, "dopplerInterp", cfg->dopplerInterp);
    SC(s, "elevationInterp", cfg->elevationInterp);
    SC(s, "minHits", cfg->minHits);
    SC(s, "maxMiss", cfg->maxMiss);
    SC(s, "minExtent_deg", cfg->minExtent_deg);
    mat_struct_struct(P, "report", s);

    s = mat_struct_new();                                /* tracker */
    char name[96];
    tracker_name(cfg, name, sizeof name);
    SC(s, "enable", cfg->trackEnable);
    mat_struct_char(s, "name", cfg->trackEnable ? name : "off");
    if (cfg->trackEnable) {
        mat_struct_char(s, "method", cfg->trackMethod == TRACK_METHOD_PMBM ? "pmbm" :
                                     cfg->trackMethod == TRACK_METHOD_TOMHT ? "tomht" : "ekf");
        SC(s, "dim", cfg->trackDim);
        SC(s, "az0_deg", cfg->trackAz0_deg);
        SC(s, "scanTime_s", cfg->scanTime_s);
        mat_struct_char(s, "model", cfg->trackModel == TRACK_CA ? "ca" : "cv");
        mat_struct_double(s, "sigma", cfg->trackSigma, 1, 3);
        mat_struct_double(s, "q", cfg->trackQ, 1, 3);
        SC(s, "sigmaV0", cfg->trackSigmaV0);
        SC(s, "sigmaA0", cfg->trackSigmaA0);
        SC(s, "batchGap_s", cfg->trackBatchGap_s);
        SC(s, "batchMax_s", cfg->trackBatchMax_s);
        SC(s, "gate", cfg->trackGate);
        SC(s, "useVr", cfg->trackUseVr);
        SC(s, "vrUnambiguous", cfg->trackVrUnambiguous);
        SC(s, "sigmaVr", cfg->trackSigmaVr);
        SC(s, "vrSpan", cfg->trackVrSpan);
        SC(s, "gateVr", cfg->trackGateVr);
        SC(s, "confirmHits", cfg->trackConfirmHits);
        SC(s, "maxMissTentative", cfg->trackMaxMissTentative);
        SC(s, "maxMiss", cfg->trackMaxMiss);
        SC(s, "minSpeed", cfg->trackMinSpeed);
        if (cfg->trackMethod != TRACK_METHOD_EKF) {
            MatStruct *m = mat_struct_new();             /* IMM-UKF (PMBM, TOMHT) */
            SC(m, "pd", cfg->pmbmPd);
            SC(m, "ps", cfg->pmbmPs);
            SC(m, "logLambdaC", cfg->pmbmLogLambdaC);
            SC(m, "logLambdaU", cfg->pmbmLogLambdaU);
            SC(m, "logVrDensity", cfg->pmbmLogVrDensity);
            SC(m, "maxHyp", cfg->pmbmMaxHyp);
            SC(m, "logHypPrune", cfg->pmbmLogHypPrune);
            SC(m, "rPrune", cfg->pmbmRPrune);
            SC(m, "rConfirm", cfg->pmbmRConfirm);
            SC(m, "rDelete", cfg->pmbmRDelete);
            double models[3];
            for (int k = 0; k < cfg->pmbmNModels && k < 3; k++) models[k] = cfg->pmbmModels[k];
            mat_struct_double(m, "models", models, 1, (size_t)(cfg->pmbmNModels < 3 ? cfg->pmbmNModels : 3));
            mat_struct_double(m, "mu0", cfg->pmbmMu0, 1, (size_t)(cfg->pmbmNModels < 3 ? cfg->pmbmNModels : 3));
            SC(m, "switchLambda", cfg->pmbmSwitchLambda);
            mat_struct_double(m, "qCv", cfg->pmbmQCv, 1, 3);
            mat_struct_double(m, "qCa", cfg->pmbmQCa, 1, 3);
            mat_struct_double(m, "qCt", cfg->pmbmQCt, 1, 3);
            SC(m, "qTurn", cfg->pmbmQTurn);
            SC(m, "sigmaTurn0", cfg->pmbmSigmaTurn0);
            SC(m, "sigmaClimb0", cfg->pmbmSigmaClimb0);
            SC(m, "varUnusedA", cfg->pmbmVarUnusedA);
            SC(m, "varUnusedW", cfg->pmbmVarUnusedW);
            SC(m, "ukfScale", cfg->pmbmUkfScale);
            SC(m, "latency_s", cfg->pmbmLatency_s);
            SC(m, "maintStep_s", cfg->pmbmMaintStep_s);
            mat_struct_struct(s, "imm", m);
        }
        if (cfg->trackMethod == TRACK_METHOD_TOMHT) {
            MatStruct *m = mat_struct_new();
            SC(m, "pd", cfg->tomhtPd);
            SC(m, "logHypPrune", cfg->tomhtLogHypPrune);
            SC(m, "confirmScore", cfg->tomhtConfirmScore);
            SC(m, "deleteScore", cfg->tomhtDeleteScore);
            SC(m, "pDelete", cfg->tomhtPDelete);
            SC(m, "maxHyp", cfg->tomhtMaxHyp);
            SC(m, "kBest", cfg->tomhtKBest);
            SC(m, "nScan", cfg->tomhtNScan);
            mat_struct_struct(s, "tomht", m);
        }
    }
    mat_struct_struct(P, "tracker", s);

    s = mat_struct_new();                                /* receiver options */
    mat_struct_char(s, "mfMethod", cfg->mfMethod == MF_DIRECT ? "direct" : cfg->mfMethod == MF_FFT ? "fft" : "auto");
    SC(s, "threads", rt->threads);
    SC(s, "slowTimeSkipping", rt->slowTimeSkipping);
    SC(s, "rdRate_Hz", cfg->rdRate_Hz);
    SC(s, "consolePlots", cfg->consolePlots);
    mat_struct_struct(P, "runtime", s);

    mat_write_tree(w, "params", P);
    return mat_write_close(w);
}

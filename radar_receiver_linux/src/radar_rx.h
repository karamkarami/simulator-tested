/*
 * radar_rx.h - public types and functions of the radar receiver.
 *
 * The receiver gets the signed 16-bit DDC output of every PRI over UDP (sent by
 * the MATLAB simulator, see ../README.md for the packet format), and runs the
 * processing chain of the MATLAB simulator: DBF, matched filter, CPI buffers,
 * MTI, Doppler, 2-D CFAR, pruning and plot extraction. Plots can be shown on a
 * Qt PPI and the CFAR detections can be sent over a serial port.
 *
 * All quantities that MATLAB derives (DBF matrix, matched filters, axes, CPI
 * length, CFAR factor, ...) are read from a configuration file written by
 * export_receiver_config.m, so both implementations use the same numbers.
 */
#ifndef RADAR_RX_H
#define RADAR_RX_H

#include <signal.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------ */
/* Configuration (config.c)                                                 */
/* ------------------------------------------------------------------------ */

enum { PULSE_SHORT = 0, PULSE_LONG = 1, N_PULSES = 2 };
enum { MF_AUTO = 0, MF_DIRECT = 1, MF_FFT = 2 };
enum { CFAR_2D = 0, CFAR_RANGE = 1, CFAR_DOPPLER = 2 };        /* window (cfg.cfar.method)  */
enum { CFAR_CA = 0, CFAR_GO = 1, CFAR_SO = 2, CFAR_OS = 3, CFAR_CE = 4 };   /* detector (cfg.cfar.type) */
enum { DOPPLER_FFT = 0, DOPPLER_NUFFT = 1 };   /* Doppler processing (cfg.proc.staggerDoppler)    */
enum { TRACK_CV = 0, TRACK_CA = 1 };                          /* EKF motion model          */
enum { TRACK_METHOD_EKF = 0, TRACK_METHOD_PMBM = 1, TRACK_METHOD_TOMHT = 2 };   /* tracker (cfg.track.method) */
enum { PMBM_CV = 0, PMBM_CA = 1, PMBM_CT = 2 };               /* IMM models of the PMBM     */
enum { WAVE_LFM = 0, WAVE_PHASE_CODE = 1 };
enum { COMPRESSION_MATCHED = 0, COMPRESSION_MISMATCH = 1 };
enum { SNR_NOISE_MEDIAN = 0, SNR_NOISE_CFAR = 1 };   /* noise floor of the SNR of the reports */
enum { MTI_CANCELLER = 0, MTI_FIR = 1 };
enum { UDP_FORMAT_RDDC = 0, UDP_FORMAT_VBDR = 1 };   /* packet format of the DDC stream (cfg.udp.format) */             /* clutter filter (cfg.mti.type; information) */
const char *cfar_method_name(int method);   /* "2-D", "1-D range", "1-D Doppler" */
const char *cfar_type_name(int type);       /* "CA", "GO", "SO", "OS", "CE" */

typedef struct {
    int     enable;        /* the sub-pulse is transmitted and processed (cfg.wave.<pulse>.enable;
                              optional in the file: 1). A sub-pulse that is off has no cells
                              (nRange = rxCount = 0) and is skipped by the processing      */
    double  fc_Hz;         /* its carrier (cfg.wave.<pulse>.fc_Hz; optional: 0 = not given) */
    double *beamU;         /* [nBeams] beam directions on its carrier (optional: cfg->beamU) */
    double  vrSpan_mps;    /* span of its Doppler axis lambda_m / (2 row PRI) (tracker)   */
    int     rxFirst;       /* receive interval of this sub-pulse: DDC samples rxFirst ... */
    int     rxCount;       /* rxFirst+rxCount-1 (0-based; the PRI header holds it 1-based) */
    int     ddcFirst;      /* first DDC sample of the pulse compression input (0-based); the
                              input may start before or end after the receive interval:
                              those samples are zero (a filter longer than the pulse)  */
    int     ddcCount;      /* samples of the pulse compression input (= nRange + nTaps - 1) */
    int     nTaps;         /* pulse compression filter length                         */
    int     nRange;        /* range cells                                             */
    float  *tapsRe;        /* compression filter taps [nTaps]: y[r] = sum_j h[j]*x[r+nTaps-1-j]
                              (matched filter, or the mismatch filter of a phase code;
                              range_m/delay_s already hold the earlier peak of the latter) */
    float  *tapsIm;
    double *range_m;       /* range axis [nRange]                                     */
    double *delay_s;       /* delay axis [nRange]                                     */
    double  lambda_m;      /* wavelength used for the radial velocity                 */
    int     testFirst;     /* range cells tested by the CFAR: testFirst ... testFirst+testCount-1 */
    int     testCount;     /* (the others are only training cells)                     */
    int     waveform;      /* WAVE_LFM or WAVE_PHASE_CODE (information)                */
    int     compression;   /* COMPRESSION_MATCHED or COMPRESSION_MISMATCH (information)  */
    int     nChips;        /* phase code: chips (0 for LFM)                           */
    int     mfShift;       /* samples by which the mismatch filter peaks earlier (0; < 0: later) */
    int     mtiCells;      /* range cells 0 ... mtiCells-1 get the clutter filter (range_m <=
                              mtiMaxRange_m); the others pass it unfiltered (newest pulse) */
    int     nPulse;        /* the transmitted pulse at the DDC rate (for the tests)    */
    float  *pulseRe, *pulseIm;
} PulseConfig;

/* Slow-time settings of one PRF of the scan-to-scan PRF (cfg.prf.mode = 'scan2scan'):
 * records s2s<j>_... of the configuration file */
typedef struct {
    double  pri_s;
    int     cpiLength, rowLen, nfft, notchBins, minHits;
    float  *dopWindow;       /* [rowLen] */
    float  *noiseGain;       /* [nfft]   */
    float  *noiseGainBypass; /* [nfft]   */
    double *dopFreq_Hz;      /* [nfft]   */
    double  trackVrSpan, trackVrSpanPulse[N_PULSES];
} PrfProfile;

typedef struct {
    /* antenna and beamforming */
    int     nCh;           /* receive elements: channels of the DDC output (cfg.ant.nElements) */
    int     nChUsed;       /* channels processed: the first nChUsed of every PRI (cfg.rx.nChannelsUsed;
                              optional in the file: nCh). A PRI needs at least nChUsed channels */
    int     nBeams;        /* elevation beams (1 when DBF is disabled)                */
    int     dbfEnabled;    /* 0 with a single element                                 */
    float  *wRe, *wIm;     /* DBF matrix, row-major [nBeams x nChUsed]: beam b = sum_ch w[b*nChUsed+ch]*x[ch] */
    double *beamU;         /* [nBeams] sine of the beam angle from broadside (NaN without DBF) */
    double  tilt_deg;
    double  height_m;
    double  rate_dps;      /* antenna rotation rate [deg/s]                           */
    double  azBeamwidth_deg;
    double  maxRange_m;    /* instrumented range (for the display)                    */
    double  handover_m;    /* the short pulse is used nearer than this range (the first long-pulse
                              cell; Inf without the long pulse)                       */

    /* DDC output */
    int     nDdc;          /* DDC samples per PRI and channel (a longer record is accepted) */
    double  fsOut_Hz;      /* DDC output rate                                         */
    PulseConfig pulse[N_PULSES];

    /* slow time */
    int     cpiLength;     /* pulses in a CPI buffer                                  */
    int     nStagger;      /* PRIs of the stagger (1 = fixed PRF)                     */
    int     staggered;
    int     nMtiTaps;
    double *mtiCoeffs;     /* row-major [nStagger x nMtiTaps]: the canceller or the FIR
                              filter of cfg.files.mtiFir (cfg.mti.type)               */
    int     mtiType;       /* MTI_CANCELLER or MTI_FIR (optional in the file: canceller) */
    int     nRows;         /* rows of the Doppler processing (de-interleaved for staggered PRF) */
    int     rowLen;        /* samples per row                                         */
    int     nfft;          /* Doppler FFT length (power of 2)                         */
    float  *dopWindow;     /* [rowLen] Hamming window, unit energy                    */
    float  *noiseGain;     /* [nfft] noise gain of every Doppler bin (after fftshift) */
    double  mtiMaxRange_m; /* clutter filter only up to this range (optional in the file: Inf) */
    float  *noiseGainBypass; /* [nfft] noise gain of the cells beyond it (no clutter filter) */
    double *dopFreq_Hz;    /* [nfft] Doppler axis (after fftshift)                    */
    int     notchBins;     /* bins on each side of zero Doppler that are not tested   */
    int     zeroBin;       /* bin of zero Doppler (FFT: nfft/2; NUFFT: from the file, -1 = not on
                              the grid: no notch)                                      */
    /* Doppler processing (optional in the file: FFT): DOPPLER_FFT (the FFT of the CPI or of the
     * de-interleaved rows) or DOPPLER_NUFFT (staggered PRF, cfg.proc.staggerDoppler = 'nufft': the
     * non-uniform DFT of all rowLen MTI outputs on the nfft velocity bins, P = |A y|^2 with the matrix
     * of the stagger position of the oldest output, already scaled for a flat noise floor) */
    int     dopplerMethod;
    float  *nufftARe, *nufftAIm;   /* [nStagger][nfft][rowLen] */
    float  *nufftBypassGain;       /* [nStagger][nfft]: noise gain of the cells without clutter filter */
    double  notchMaxRange_m; /* apply the MTI notch only at range_m <= this limit      */
    double  medianToMean;  /* median / mean of the noise power (noise floor estimate) */

    /* CFAR, pruning, plot extraction */
    int     cfarType;      /* detector CFAR_CA, CFAR_GO, CFAR_SO or CFAR_OS (cfar.c; optional
                              in the file: CA)                                            */
    int     cfarMethod;    /* CFAR_2D, CFAR_RANGE or CFAR_DOPPLER: the window below   */
    int     guardRange, trainRange, guardDoppler, trainDoppler, nTrain;
    int     cfarNHalf;     /* GO, SO: training cells per half window                  */
    int     cfarOsK;       /* OS: rank of the statistic (of nTrain)                    */
    int     cfarCeKeep;    /* CE: training cells averaged (the nTrain - cfarCeKeep largest are
                              censored; at the range edges in proportion)                  */
    int     nCeTable;      /* CE: distribution of the statistic for the Pfa: log z, log G(z)  */
    double *ceLogZ, *ceLogG;
    double  cfarNoiseScale;  /* mean of the statistic for unit noise (noise estimate =
                                statistic / this; 1 for CA)                            */
    double  alpha;         /* CFAR threshold factor (design value from cfg.cfar.pfa)  */
    double  alphaMin;      /* lowest factor that can be set (cfg.cfar.maxPfa; optional
                              in the file: 0)                                          */
    double  nLooks;        /* looks of the noise power in a cell (for cfar_pfa)       */
    double  sidelobeMargin_dB;
    int     snrNoise;      /* SNR = power / noise - 1 of the reports with the noise of
                              SNR_NOISE_MEDIAN: the median of the beam's map / medianToMean, or
                              SNR_NOISE_CFAR: the CFAR estimate of the cell (no median computed);
                              cfg.report.snrNoise, optional in the file (median)            */
    int     dopplerInterp;   /* peak interpolation (cfg.report.dopplerInterp / elevationInterp;  */
    int     elevationInterp; /* optional in the file: 0 / 1): parabola through the powers [dB] of
                                the neighbouring Doppler bins / beams of every reported peak  */
    int     minHits, maxMiss, updateEvery;
    int     consolePlots;  /* print every plot and confirmed track (cfg.console.plots; --quiet, --print-plots) */
    double  rdRate_Hz;     /* range-Doppler maps per second handed to the PPI tab (cfg.display.rdRate_Hz;
                              optional in the file: 10; --ppi-rd-rate)                          */
    double  minExtent_deg;  /* azimuth extent of a plot (minHits = this / azimuth step per update) */
    double  meanPri_s;      /* mean PRI                                               */

    /* radar mode (information; optional in the file): cfg.modes.active (0 ... 7) */
    int     radarMode;
    /* (information; optional) DDC record set by hand (cfg.rx.manualDdc) and the time of its first
     * sample, code / mismatch pair of each sub-pulse (0: no phase code), (T, w) pair of the
     * staggered PRF (0: not staggered) */
    int     ddcManual;
    double  ddcT0_s;
    int     codeIndex[N_PULSES], staggerSet;
    /* the DDC stream (optional in the file: UDP_FORMAT_RDDC) and the PRIs of the stagger (optional: none;
     * the VBDR stream finds the stagger position of a PRI from its PRT with them) */
    int     udpFormat;
    int     nPriSequence;
    double *priSequence_s;
    /* scan-to-scan PRF (optional in the file: off): one set of slow-time settings per PRF of the
     * scans. The fields above (cpiLength, rowLen, nfft, dopWindow, noiseGain, noiseGainBypass,
     * dopFreq_Hz, notchBins, minHits, meanPri_s, trackVrSpan, pulse[].vrSpan_mps) hold those of
     * the PRF in use: radar_config_select_prf() (the processor calls it when the stagger index of
     * the PRI header, the number of the PRF of the scan, changes). The arrays have room for the
     * largest set (maxCpiLength, maxNfft, maxRowLen). */
    int     scan2scan;
    int     nPrf;
    PrfProfile *prf;        /* [nPrf]                                                 */
    int     prfIndex;       /* the PRF in use (0-based)                               */
    int     maxCpiLength, maxNfft, maxRowLen;

    /* tracker (optional in the file: off if missing); the common settings, the
     * EKF ones (trackModel ... trackSigmaA0) and the PMBM ones (pmbm...) */
    int     trackEnable, trackDim;   /* dimension 2 (range, azimuth) or 3 (+ elevation)  */
    int     trackMethod;    /* TRACK_METHOD_EKF or TRACK_METHOD_PMBM (optional: EKF)        */
    double  trackAz0_deg;   /* antenna azimuth at t = 0 (PMBM: look times)                */
    int     trackModel;     /* TRACK_CV (constant velocity) or TRACK_CA (constant acceleration) */
    double  trackSigma[3];  /* measurement errors: range [m], azimuth, elevation [rad]   */
    double  trackQ[3];      /* process noise per axis: white acceleration [m^2/s^3] (CV)
                               or white jerk [m^2/s^5] (CA) spectral density              */
    double  trackSigmaV0;   /* velocity std of a new track [m/s]                          */
    double  trackSigmaA0;   /* acceleration std of a new track (CA) [m/s^2]               */
    double  trackBatchGap_s, trackBatchMax_s;   /* GNN batch: time without a new plot,
                                                   longest collection time [s]            */
    double  trackGate;      /* chi-square gate of d^2                                     */
    int     trackUseVr;     /* the radial velocity of the plots is a measurement too      */
    int     trackVrUnambiguous;   /* no target faster than half the span: the plot velocity is true */
    double  trackSigmaVr;   /* its error [m/s]                                            */
    double  trackVrSpan;    /* the plot velocity is known modulo this (lambda / (2 row PRI)) [m/s] */
    double  trackGateVr;    /* chi-square gate of d^2 with the radial velocity            */
    int     trackConfirmHits, trackMaxMissTentative, trackMaxMiss;
    double  trackMinSpeed;  /* delete the tracks with 2+ updates slower than this (|v| of the state)
                               [m/s]; 0 = off (cfg.track.minSpeedEnable; optional in the file: off) */
    double  scanTime_s;     /* antenna revolution                                         */
    /* PMBM tracker with IMM-UKF (track_setup.m, tracker_pmbm.c) */
    double  pmbmPd, pmbmPs; /* detection probability, survival probability per scan       */
    double  pmbmLogLambdaC, pmbmLogLambdaU, pmbmLogVrDensity;   /* false plots, undetected
                               targets per unit measurement volume, 1/vrSpan (logs)         */
    int     pmbmMaxHyp;     /* global hypotheses kept                                     */
    double  pmbmLogHypPrune, pmbmRPrune, pmbmRConfirm, pmbmRDelete;
    int     pmbmNModels, pmbmModels[3];    /* PMBM_CV, PMBM_CA, PMBM_CT                   */
    double  pmbmMu0[3], pmbmSwitchLambda;  /* model probabilities of a new track; second
                                              eigenvalue of the switching per scan         */
    double  pmbmQCv[3], pmbmQCa[3], pmbmQCt[3];   /* process noise per axis of the models  */
    double  pmbmQTurn, pmbmSigmaTurn0, pmbmSigmaClimb0;   /* CT turn rate noise, new track:
                                                         turn rate and climb rate std   */
    double  pmbmVarUnusedA, pmbmVarUnusedW; /* acceleration / turn rate of the models without */
    double  pmbmUkfScale, pmbmWm0, pmbmWc0, pmbmWi;  /* unscented transform              */
    double  pmbmLatency_s, pmbmMaintStep_s; /* plots in the tracker this long after the look;
                                               missed looks checked every               */
    /* TOMHT tracker (track_tomht.m, tracker_tomht.c): the IMM-UKF, densities and timing
     * above (pmbm...), and its hypothesis management */
    double  tomhtPd, tomhtLogHypPrune, tomhtConfirmScore, tomhtDeleteScore, tomhtPDelete;
    int     tomhtMaxHyp, tomhtKBest, tomhtNScan;

    /* run-time option (not in the file): matched filter method,
     * MF_AUTO (the faster one for the thread count), MF_DIRECT or MF_FFT */
    int     mfMethod;
} RadarConfig;

/* Loads the file written by export_receiver_config.m. Returns 0 on success. */
int  radar_config_load(const char *path, RadarConfig *cfg, char *err, size_t errLen);
void radar_config_free(RadarConfig *cfg);
/* Scan-to-scan PRF: puts the slow-time settings of PRF k (0-based) into cfg. Returns 1 when
 * the PRF changed, 0 when it was in use already or cfg has no scan-to-scan PRF. */
int  radar_config_select_prf(RadarConfig *cfg, int k);

/* ------------------------------------------------------------------------ */
/* DDC output of one PRI (udp_input.c)                                       */
/* ------------------------------------------------------------------------ */

typedef struct {
    uint32_t pulseNumber;
    uint32_t scanNumber;
    float    az_deg;         /* antenna azimuth                                       */
    uint16_t staggerIndex;   /* 1-based stagger position (1 with a fixed PRF)         */
    double   time_s;         /* PRI start time                                        */
    float    ddcFs_Hz;
    float    ddcT0_s;        /* time of the first DDC sample in the PRI              */
    uint32_t shortIdx[2];    /* first/last DDC sample of the short pulse (1-based)    */
    uint32_t longIdx[2];
    int      nSamples;       /* DDC samples per channel                               */
    int      nCh;            /* channels                                              */
    int16_t *iq;             /* I,Q pairs, channel after channel:
                                I = iq[2*(ch*nSamples + n)], Q = iq[2*(ch*nSamples + n) + 1] */
    size_t   capacity;       /* allocated int16 values in iq (managed by udp_input)   */
    /* VBDR stream (format UDP_FORMAT_VBDR): pulseNumber = stream number, az_deg, nSamples, nCh = 4
     * and the fields below come from the stream; vbdr_complete_pri() derives the others (time,
     * scan, stagger index, DDC rate and layout) from the stream and the configuration */
    int      format;         /* UDP_FORMAT_RDDC or UDP_FORMAT_VBDR                    */
    uint16_t prfNumber;      /* PRF number of the system parameters                   */
    double   prevPrt_s;      /* previous PRT (the PRI before this one)                */
    double   curPrt_s;       /* current PRT (this PRI)                                */
    uint32_t loss;           /* "Loss=" of the label: packets the board lost so far   */
} PriData;

/* State of vbdr_complete_pri() (zero-initialised before the first PRI). */
typedef struct {
    int      have;           /* a PRI was completed before                            */
    uint32_t lastStream;
    double   lastTime_s, lastCurPrt_s;
    float    lastAz_deg;
    int      lastStagger;    /* 1-based stagger position of the last PRI              */
    uint32_t scan;
    uint64_t missing;        /* streams that never arrived (gaps of the stream numbers) */
} VbdrTiming;
/* Completes a PRI of the VBDR stream: PRI start time from the sum of the PRTs (the first PRI at
 * 0; streams that never arrived, found from the stream numbers, are estimated with the stagger
 * sequence or the last PRT), scan number from the azimuth wrap (first scan 1), stagger index from
 * the current PRT (staggered PRF: the position in cfg->priSequence_s; scan-to-scan PRF: the PRF of
 * the scan), DDC rate, time of the first sample and the receive intervals of the configuration. */
void vbdr_complete_pri(VbdrTiming *st, const RadarConfig *cfg, PriData *pri);

typedef struct UdpInput UdpInput;

typedef struct {
    uint64_t datagrams;      /* datagrams received                                    */
    uint64_t foreign;        /* datagrams with a wrong magic number or size           */
    uint64_t prisComplete;   /* PRIs delivered                                        */
    uint64_t prisIncomplete; /* PRIs dropped because packets were missing             */
    uint32_t boardLoss;      /* VBDR: "Loss=" of the newest label (packets the board lost) */
} UdpStats;

/* Opens the UDP input on 'port' (all interfaces). recordPath (may be NULL)
 * stores every received datagram in a file; replayPath (may be NULL) reads
 * datagrams from such a file instead of the network. Returns NULL on error
 * (message in err). */
UdpInput *udp_input_open(int port, const char *recordPath, const char *replayPath, char *err, size_t errLen);
/* The packet format: UDP_FORMAT_RDDC (default) or UDP_FORMAT_VBDR (VBDR board: a stream of N packets
 * with a 28-byte label per PRI, CH1 & CH2 in the first half, CH3 & CH4 in the second; the PRIs have
 * 4 channels, complete them with vbdr_complete_pri()). Call before the first udp_input_next(). */
void udp_input_set_format(UdpInput *in, int format);
/* Waits for the next complete PRI. Returns 1 when 'pri' was filled, 0 on
 * timeout, -1 at the end of a replay file or on an error. 'pri' must be
 * zero-initialised before the first call; its buffer is reused. */
int  udp_input_next(UdpInput *in, PriData *pri, int timeout_ms);
void udp_input_stats(const UdpInput *in, UdpStats *stats);
void udp_input_close(UdpInput *in);
void pri_data_free(PriData *pri);

/* Reception in a thread of its own (input_thread.c): the PRIs from
 * udp_input_next() go into a queue of 'capacity' PRIs, so the socket is
 * emptied while the processing runs. dropWhenFull = 1 (live stream): a PRI
 * that finds the queue full is dropped and counted; 0 (replay): the reader
 * waits. input_thread_next() returns 1 with *pri valid until the next call,
 * 0 on timeout, -1 when the input has ended and the queue is empty. */
typedef struct InputThread InputThread;
InputThread *input_thread_start(UdpInput *in, int capacity, int dropWhenFull);
int      input_thread_next(InputThread *t, PriData **pri, int timeout_ms);
int      input_thread_backlog(InputThread *t);     /* PRIs waiting in the queue */
uint64_t input_thread_dropped(InputThread *t);
void     input_thread_stop(InputThread *t);        /* stops and joins the thread */

/* ------------------------------------------------------------------------ */
/* Detections and plots (processor.c, extractor.c)                           */
/* ------------------------------------------------------------------------ */

typedef struct {
    double   time_s;       /* time of the CPI centre                                  */
    double   az_deg;       /* azimuth of the CPI centre                               */
    double   range_m;
    double   el_deg;       /* NaN without beamforming                                 */
    double   vr_mps;       /* ambiguous radial velocity (+ = away)                    */
    double   snr_dB;       /* against the noise floor of the beam                     */
    double   power;        /* range-Doppler power of the cell                         */
    int      beam, rIdx, dIdx;   /* 0-based cell indices                              */
    int      pulse;        /* PULSE_SHORT or PULSE_LONG                               */
    int      hits;         /* plots: CPIs with a detection (0 for CFAR detections)    */
    uint32_t pulseNumber;  /* newest pulse of the CPI                                 */
} Detection;

/* A track of the tracker (tracker.c: EKF or PMBM, as track_update.m of the
 * simulator). Position and velocity relative to the antenna: x East, y North,
 * z up (z = 0 in 2-D). Reported after every update, and once more with
 * deleted = 1. */
typedef struct {
    uint32_t id;
    int      confirmed;    /* 0: tentative                                           */
    int      deleted;      /* 1: the track has been deleted (last report)            */
    int      dim;          /* 2 or 3                                                  */
    double   time_s;       /* time of the last update                                 */
    double   pos_m[3], vel_mps[3];
    double   acc_mps2[3];  /* acceleration (constant-acceleration model, else 0)     */
    double   alt_m;        /* height above ground (3-D), NaN in 2-D                  */
    int      hits;         /* updates                                                 */
    int      misses;       /* scans without a plot since the last update             */
    int      method;       /* TRACK_METHOD_EKF, TRACK_METHOD_PMBM or TRACK_METHOD_TOMHT */
    double   existence;    /* PMBM: probability that a target made its newest plot (NaN: EKF) */
    int      nModes;       /* PMBM: IMM models (0: EKF)                              */
    int      modes[3];     /* PMBM: PMBM_CV, PMBM_CA, PMBM_CT                        */
    double   modeProb[3];  /* PMBM: their probabilities                              */
} Track;

/* Name of the configured tracker, e.g. "EKF 3-D CV", "PMBM 2-D IMM CV/CA/CT" (tracker.c). */
void tracker_name(const RadarConfig *cfg, char *buf, size_t len);

/* Range-Doppler map of one beam (processor_set_rdmap, callback rdmap): the
 * power maps of both sub-pulses after the CFAR's noise-gain equalisation,
 * joined in range like the video (short-pulse cells nearer than the first
 * long-pulse cell, then the long-pulse cells), in range bins of dR_m from 0
 * (power[d * nBins + b], 0 = no cell in the bin) and the nfft Doppler bins
 * (after fftshift). dets: the reports of that beam in the same update. */
typedef struct {
    int           beam, nBins, nfft;
    double        dR_m;
    const float  *power;       /* [nfft][nBins]                                    */
    const double *vel_mps;     /* [nfft] radial velocity of every bin (+ = away)  */
    int           notchBins;   /* bins on each side of zero Doppler not tested     */
    int           zeroBin;     /* bin of zero Doppler (-1: none, NUFFT grid)         */
    double        noise;       /* median power of the cells (noise floor estimate) */
    double        time_s, az_deg;
    const Detection *dets;
    int           nDets;
} RdMap;

/* Profile of the selected beam (processor_set_profile): PROFILE_DBF the
 * beamformer output |x|^2 of both receive intervals, PROFILE_MF the pulse
 * compression output |y|^2 of both sub-pulses, joined like the video (short
 * pulse up to the first long-pulse range cell, then the long pulse). */
enum { PROFILE_NONE = 0, PROFILE_DBF = 1, PROFILE_MF = 2 };

/* Block outputs for logging (processor_set_log): the bits of the mask. LOG_DDC and LOG_DETECTION
 * are written by the application (the PRI as received, the detections callback). */
enum { LOG_DDC = 1, LOG_DBF = 2, LOG_MF = 4, LOG_DOPPLER = 8, LOG_DETECTION = 16 };
/* LOG_SCREEN: the PPI records its window (with the microphone) to <log>/screen_record.mp4 (ppi_qt.cpp);
 * not a block of the processing (a log with only this bit holds radar_parameters.mat and the video) */
enum { LOG_SCREEN = 32, LOG_BLOCKS = LOG_DDC | LOG_DBF | LOG_MF | LOG_DOPPLER | LOG_DETECTION };

typedef struct {
    int          kind;         /* LOG_DBF, LOG_MF or LOG_DOPPLER                              */
    int          pulse;        /* PULSE_SHORT or PULSE_LONG                                   */
    int          nBeams, nRange, nfft;
    /* LOG_DBF:     beamformer output [nBeams][nSamples] (complex): the pulse compression input,
     *              DDC samples first ... first+nSamples-1 (0-based; zero outside the receive window)
     * LOG_MF:      pulse compression output [nBeams][nRange] (complex)
     * LOG_DOPPLER: range-Doppler power after the noise-gain equalisation [nfft][nBeams][nRange]
     *              (im = NULL), Doppler axis dopFreq_Hz */
    int          nSamples, first;
    const float *re, *im;
    uint32_t     pulseNumber;
} BlockOutput;

typedef struct {
    /* called for the pruned CFAR detections of every CPI update (n may be 0) */
    void (*detections)(void *user, int pulse, const Detection *dets, int n);
    /* called for every new plot */
    void (*plot)(void *user, const Detection *plot);
    /* called once per PRI with the antenna azimuth (display) */
    void (*antenna)(void *user, double az_deg, double time_s, uint32_t scan);
    void *user;
    /* called for every track update and deletion (tracking enabled in the configuration) */
    void (*track)(void *user, const Track *track);
    /* called once per PRI with the profile of processor_set_profile(): row[nBins], range
     * bins of dR_m from 0 (0 = no sample in the bin) */
    void (*profile)(void *user, int kind, const float *row, int nBins, double dR_m, double az_deg, double time_s);
    /* called with the range-Doppler map of processor_set_rdmap() */
    void (*rdmap)(void *user, const RdMap *map);
    /* called with the block outputs of processor_set_log() (logging) */
    void (*block)(void *user, const BlockOutput *out);
} ProcessorCallbacks;

typedef struct {
    double dbf_ms, mf_ms, mti_ms, doppler_ms, cfar_ms, extract_ms, track_ms, total_ms;
    double interp_ms;    /* peak interpolation of the reports (Doppler, elevation) */
    /* time of the application in its output callbacks (main.c adds it; the
     * processor leaves it 0): console output (printf of plots, tracks and
     * messages) and the data handed to the PPI display */
    double console_ms, display_ms;
    double log_ms;       /* logging (main.c adds it): the copies of the block outputs handed
                            to the log writer thread (the files are written by that thread) */
    uint64_t pris, updates;
    uint64_t rejected;   /* PRIs whose layout (channels, DDC samples, sub-pulse samples,
                            DDC rate) differs from the configuration: not processed */
    uint64_t skipped;    /* CPI updates skipped (processor_set_slow_time(p, 0))      */
    uint64_t prfSwitches; /* scan-to-scan PRF: changes of the PRF (CPIs restarted)    */
} ProcessorTiming;     /* accumulated times */

typedef struct Processor Processor;
Processor *processor_create(const RadarConfig *cfg, int nThreads);
void processor_process(Processor *p, const PriData *pri, const ProcessorCallbacks *cb);
void processor_flush(Processor *p, const ProcessorCallbacks *cb);   /* close all open plots */
/* enabled = 0: the next PRIs get only the fast-time processing and the MTI
 * (so the CPI buffers stay complete), their Doppler/CFAR/extraction updates
 * are skipped. Used by the input loop when it falls behind the stream. */
void processor_set_slow_time(Processor *p, int enabled);
/* Profile of one beam (0-based) for the display, every PRI: kind PROFILE_DBF,
 * PROFILE_MF or PROFILE_NONE (default: nothing is computed). */
void processor_set_profile(Processor *p, int kind, int beam);
/* Range-Doppler map of one beam (0-based) for the display: enabled = 0
 * (default) computes nothing; else it is handed to the callback after the
 * CPI updates of a PRI, at most every minPeriod_s seconds. */
void processor_set_rdmap(Processor *p, int enabled, int beam, double minPeriod_s);
/* Block outputs handed to the block callback (logging): mask of LOG_DBF, LOG_MF, LOG_DOPPLER
 * (0 = none, the default). */
void processor_set_log(Processor *p, int mask);
/* CFAR threshold factor (of the configured detector) used from the next
 * CPI update on (starts at cfg->alpha). threshold [dB] = 10*log10(alpha).
 * Call from the processing thread. */
void   processor_set_cfar_alpha(Processor *p, double alpha);
double processor_cfar_alpha(const Processor *p);
/* False alarm probability of the configured CFAR detector (cfg->cfarType,
 * window, nLooks) with factor alpha (cfar.c, as cfar_pfa.m; the inverse of
 * cfar_alpha.m): CA I(1/(1 + alpha/nTrain); nTrain*nLooks, nLooks), GO, SO
 * and OS by numerical integration. NaN if nLooks <= 0. */
double cfar_pfa(const RadarConfig *cfg, double alpha);
void processor_timing(const Processor *p, ProcessorTiming *t);
int  processor_threads(const Processor *p);    /* processing threads (OpenMP team size) */
void processor_destroy(Processor *p);

/* ------------------------------------------------------------------------ */
/* Logging of the block outputs (logger.c)                                   */
/* ------------------------------------------------------------------------ */

/* logger_start: the log <folder>/<name> (name NULL or "": log_<date>_<time>) with a folder per
 * block of mask (LOG_DDC ddc/, LOG_DBF dbf/, LOG_MF mf/, LOG_DOPPLER doppler/, LOG_DETECTION
 * detection/), one MAT-file per PRI, written by a thread of its own. Per PRI the processing
 * thread calls logger_begin_pri (copies the DDC output), hands the block outputs (processor
 * callback block), the detections and plots, then logger_end_pri. logger_stop writes what is
 * queued and ends the log. */
typedef struct Logger Logger;
typedef struct {
    uint64_t pris, files, dropped, bytes, errors;
    size_t   queuedBytes;          /* waiting for the writer                          */
    double   seconds;              /* since the start                                 */
    double   writeSeconds;         /* time of the writer thread in the files           */
    char     lastError[256];
} LoggerStats;
/* The state of the receiver when the log starts, for radar_parameters.mat (params_mat.c) */
typedef struct {
    const char *started;          /* date and time                                     */
    const char *configPath;       /* configuration file                                */
    const char *logBlocks;        /* filled by logger_start                            */
    int         screenRecord;     /* LOG_SCREEN requested                              */
    double      alpha;            /* CFAR threshold factor in use (0: the configuration's) */
    int         threads, slowTimeSkipping;
    int         havePri;          /* the last PRI processed before the start:          */
    uint32_t    pulseNumber;
    double      time_s, az_deg;
} RuntimeParams;
/* Writes every setting of the processing chain (cfg) and the receiver state (rt) as struct 'params'
 * to a MAT-file (-v6). Returns 0 on success. */
int radar_parameters_write(const char *path, const RadarConfig *cfg, const RuntimeParams *rt);

/* logger_start also writes <log>/radar_parameters.mat (rt may be NULL: the configuration's values)
 * and log_info.txt; logger_note appends a line (e.g. a parameter changed during the log) to it. */
Logger *logger_start(const char *folder, const char *name, int mask, const RadarConfig *cfg, const RuntimeParams *rt,
                     char *err, size_t errLen);
void logger_note(Logger *L, const char *text);
void logger_begin_pri(Logger *L, const PriData *pri);
void logger_block(Logger *L, const BlockOutput *out);
void logger_detections(Logger *L, const Detection *dets, int n);
void logger_plot(Logger *L, const Detection *plot);
void logger_end_pri(Logger *L);
void logger_stats(Logger *L, LoggerStats *s);
const char *logger_dir(const Logger *L);
int  logger_mask(const Logger *L);
void logger_stop(Logger *L, LoggerStats *final);   /* final: the statistics at the end (may be NULL) */

/* ------------------------------------------------------------------------ */
/* Raw video of the PPI (video.c)                                            */
/* ------------------------------------------------------------------------ */

/* The DDC output of a PRI averaged non-coherently over the elements used (mean of
 * I^2 + Q^2), the receive intervals of the short pulse (up to the first range
 * cell of the long pulse) and of the long pulse joined, in range bins of
 * c / (2 fsOut) from 0 to maxRange_m (as ppi_video_add.m). */
typedef struct {
    int     nBins;        /* range bins: bin b is b*dR_m ... (b+1)*dR_m           */
    double  dR_m;         /* range per DDC sample                                  */
    int     nSamples;     /* DDC samples of the joined video                       */
    int    *sample;       /* [nSamples] DDC sample index (0-based)                 */
    int    *bin;          /* [nSamples] its range bin                              */
    float  *binInv;       /* [nBins] 1 / (samples x channels) of the bin, 0 = none */
    float  *acc;          /* work buffer [nSamples]                                */
    int     nCh;          /* channels averaged: the first nChUsed                  */
} VideoLayout;
int  video_layout_init(VideoLayout *v, const RadarConfig *cfg, double maxRange_m);   /* 0 on success */
void video_row(const VideoLayout *v, const PriData *pri, float *row);              /* row[nBins] */
void video_layout_free(VideoLayout *v);

/* ------------------------------------------------------------------------ */
/* Serial output of the CFAR detections (serial_out.c)                       */
/* ------------------------------------------------------------------------ */

typedef struct SerialOut SerialOut;
/* device: "/dev/ttyUSB0" (Linux) or "COM3" (Windows); baud e.g. 921600 */
SerialOut *serial_open(const char *device, int baud, char *err, size_t errLen);
/* Sends one frame with the detections of one CPI update (see README for the
 * frame format). Never blocks: when the port cannot take the frame it is
 * dropped and counted. Returns 1 if sent, 0 if dropped. */
int  serial_send_detections(SerialOut *s, int pulse, uint32_t pulseNumber, double az_deg,
                            const Detection *dets, int n);
/* Sends one plot (frame type 2). */
int  serial_send_plot(SerialOut *s, const Detection *plot);
void serial_stats(const SerialOut *s, uint64_t *framesSent, uint64_t *framesDropped, uint64_t *bytesSent);
void serial_close(SerialOut *s);
/* Frame encoder (also used by the tests): writes the frame to buf and returns
 * its length, or 0 if buf is too small. */
size_t serial_encode_detections(uint8_t *buf, size_t bufLen, int pulse, uint32_t pulseNumber,
                                double az_deg, const Detection *dets, int n);
size_t serial_encode_plot(uint8_t *buf, size_t bufLen, const Detection *plot);
uint16_t crc16_ccitt(const uint8_t *data, size_t len);

/* ------------------------------------------------------------------------ */
/* PPI display (ppi_qt.cpp, only with Qt)                                    */
/* ------------------------------------------------------------------------ */

/* Runs worker(arg) in a background thread and the Qt PPI window in the
 * calling (main) thread until the worker returns and the window is closed,
 * or the window is closed by the user (then *stopFlag is set to 1 so the
 * worker can stop). snapshotPath (may be NULL): the PPI is saved as a PNG
 * when the worker has finished (then the window closes by itself).
 * Returns the exit code of the application. */
int  ppi_run(int argc, char **argv, double maxRange_m, int (*worker)(void *), void *arg,
             volatile sig_atomic_t *stopFlag, const char *snapshotPath);
/* Thread-safe calls from the worker thread. */
void ppi_post_antenna(double az_deg, double time_s, uint32_t scan);
void ppi_post_plot(const Detection *plot);
/* The pruned CFAR detections of a CPI update (before the plot extractor).
 * The PPI draws each one as an arc of arcWidth_deg (the azimuth step between
 * CPI updates) when its display mode shows detections. */
void ppi_post_detections(const Detection *dets, int n, double arcWidth_deg);
/* CFAR threshold field of the PPI: the worker reports the threshold it uses
 * (and its Pfa, NaN if unknown; default_dB = the design value for the Default
 * button); ppi_take_cfar_request() returns 1 and the new threshold when the
 * user pressed Apply (or Default) since the last call. */
void ppi_post_cfar_threshold(double threshold_dB, double pfa, double default_dB, double min_dB);
int  ppi_take_cfar_request(double *threshold_dB);
/* Display mode at the start of ppi_run (call before it): plots after the plot
 * extractor (default), CFAR detections before it, both, or the raw video. */
enum { PPI_SHOW_PLOTS = 0, PPI_SHOW_DETECTIONS = 1, PPI_SHOW_BOTH = 2, PPI_SHOW_VIDEO = 3, PPI_SHOW_DBF = 4,
       PPI_SHOW_MF = 5 };
void ppi_set_display_mode(int mode);
/* Beams after the DBF for the beam selector of the panel (call before ppi_run):
 * their number and elevations (NaN without DBF); the beam (0-based) and the
 * tab (0 PPI, 1 range-Doppler map with a small PPI) and the plot below the
 * range-Doppler map (0 Doppler spectrum of a range bin, 1 range profile of a
 * Doppler bin) at the start. */
void ppi_set_beams(int nBeams, const double *beamEl_deg);
void ppi_set_beam(int beam);
void ppi_set_tab(int tab);
void ppi_set_rd_plot(int plot);
/* Display modes PPI_SHOW_DBF / PPI_SHOW_MF: ppi_profile_wanted() returns
 * PROFILE_DBF / PROFILE_MF (else PROFILE_NONE) and the selected beam; the worker
 * then hands the profile of every PRI to ppi_post_profile(). The range-Doppler
 * tab: ppi_rdmap_wanted() is always 1 (with the selected beam), whichever tab is shown, the
 * worker hands the maps to ppi_post_rdmap() (which copies them). */
int  ppi_profile_wanted(int *beam);
void ppi_post_profile(int kind, const float *row, int nBins, double dR_m, double az_deg, double time_s);
int  ppi_rdmap_wanted(int *beam);
void ppi_post_rdmap(const RdMap *map);
/* Raw video (display mode PPI_SHOW_VIDEO): ppi_video_wanted() is 1 while the
 * PPI shows the video (the worker then computes it, video.c); ppi_post_video()
 * hands the video of one PRI (row[nBins], range bins of dR_m from 0, the mean
 * power of nAvg complex samples: the receive elements) to the PPI, which keeps
 * the newest data of every azimuth bin (0.25 deg). nAvg sets the noise floor:
 * the 10th percentile of a mean of nAvg exponential powers is below its mean
 * (9.77 dB for one element, 0.86 dB for 48). */
int  ppi_video_wanted(void);
void ppi_post_video(const float *row, int nBins, int nAvg, double dR_m, double az_deg, double time_s);
/* Track updates and deletions (drawn as hollow triangles, Tracks checkbox). */
void ppi_post_track(const Track *track);
/* Range-Doppler tab before the first map: its empty axes show this Doppler axis [m/s] (call
 * before ppi_run; the range axis is 0 ... maxRange_m of ppi_run). */
void ppi_set_rd_axes(double vMin_mps, double vMax_mps);
/* Logging (tab "Controls"): the folder, name and blocks (LOG_* mask) the controls start with
 * (call before ppi_run); ppi_take_log_request() returns 1 when the user pressed Start (*start = 1,
 * with the folder, name and mask) or Stop (*start = 0); the worker reports the state of the log
 * (shown in the tab, the REC mark of the PPI) with ppi_post_log_status (error: NULL or a message). */
void ppi_set_log_defaults(const char *folder, const char *name, int mask);
int  ppi_take_log_request(int *start, char *folder, size_t folderLen, char *name, size_t nameLen, int *mask);
void ppi_post_log_status(int active, const char *dir, int mask, const LoggerStats *stats, const char *error);
/* Screen recording of a log with LOG_SCREEN (call before ppi_run; NULL / 0: the defaults): the
 * ffmpeg program (default: next to the program, else on PATH), the sound input ("auto": the
 * default microphone, "none", or "FORMAT:DEVICE" as ffmpeg's -f FORMAT -i DEVICE) and the frame
 * rate (default 10). The window (with the mouse pointer drawn) and the sound go to
 * <log>/screen_record.mp4 (H.264 + AAC) until the log stops. */
void ppi_set_screen_record(const char *ffmpeg, const char *audio, double fps);
/* Accumulated processing times: the panel shows the mean per PRI of every
 * block [us]. detectionName / trackingName label the CFAR and tracking rows
 * (e.g. "CFAR 2-D", "Track 3D CV"; NULL keeps the previous; the panel shortens
 * a label to 15 characters, "1-D" -> "1D"). */
void ppi_post_timing(const ProcessorTiming *timing, const char *detectionName, const char *trackingName);

/* ------------------------------------------------------------------------ */
/* Platform helpers (platform.c)                                             */
/* ------------------------------------------------------------------------ */

double now_seconds(void);   /* monotonic clock */

/* Real-time settings of the calling thread (Linux; elsewhere -1 when
 * something is requested): SCHED_FIFO with priority 1 ... 99 (0 = keep the
 * scheduling) and CPU affinity "2-5,7" (NULL or "" = keep; *nCpus gets the
 * number of CPUs, nCpus may be NULL). Threads created later by this thread
 * (the OpenMP workers, the input thread) inherit both. Returns 0, -2 for an
 * invalid CPU list (nothing changed), or -1 with the reason in err (the
 * settings that worked stay). */
int rt_thread_setup(int priority, const char *cpus, int *nCpus, char *err, size_t errLen);
/* Locks all present and future pages of the process in memory (mlockall)
 * and keeps freed memory in the heap, so the processing gets no page faults
 * (Linux). Returns 0, or -1 with the reason in err. */
int rt_lock_memory(char *err, size_t errLen);

#ifdef __cplusplus
}
#endif

#endif /* RADAR_RX_H */

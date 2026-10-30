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
enum { CFAR_2D = 0, CFAR_RANGE = 1, CFAR_DOPPLER = 2 };
enum { TRACK_CV = 0, TRACK_CA = 1 };
enum { WAVE_LFM = 0, WAVE_PHASE_CODE = 1 };
enum { COMPRESSION_MATCHED = 0, COMPRESSION_MISMATCH = 1 };
const char *cfar_method_name(int method);   /* "2-D", "1-D range", "1-D Doppler" */

typedef struct {
    int     ddcFirst;      /* first DDC sample of this sub-pulse (0-based)            */
    int     ddcCount;      /* DDC samples of this sub-pulse (= nRange + nTaps - 1)    */
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
    int     testCount;     /* (the others are only training cells; default: all)      */
    int     waveform;      /* WAVE_LFM or WAVE_PHASE_CODE (information; default LFM)   */
    int     compression;   /* COMPRESSION_MATCHED or COMPRESSION_MISMATCH (information)  */
    int     nChips;        /* phase code: chips (0 for LFM)                           */
    int     mfShift;       /* samples by which the mismatch filter peaks earlier (0)   */
    int     nPulse;        /* the transmitted pulse at the DDC rate (for tests; 0 if   */
    float  *pulseRe, *pulseIm;   /* the file does not have it)                       */
} PulseConfig;

typedef struct {
    /* antenna and beamforming */
    int     nCh;           /* receive elements                                        */
    int     nBeams;        /* elevation beams (1 when DBF is disabled)                */
    int     dbfEnabled;    /* 0 with a single element                                 */
    float  *wRe, *wIm;     /* DBF matrix, row-major [nBeams x nCh]: beam b = sum_ch w[b*nCh+ch]*x[ch] */
    double *beamU;         /* [nBeams] sine of the beam angle from broadside (NaN without DBF) */
    double  tilt_deg;
    double  height_m;
    double  rate_dps;      /* antenna rotation rate [deg/s]                           */
    double  azBeamwidth_deg;
    double  maxRange_m;    /* instrumented range (for the display)                    */

    /* DDC output */
    int     nDdc;          /* DDC samples per PRI and channel                         */
    double  fsOut_Hz;      /* DDC output rate                                         */
    PulseConfig pulse[N_PULSES];

    /* slow time */
    int     cpiLength;     /* pulses in a CPI buffer                                  */
    int     nStagger;      /* PRIs of the stagger (1 = fixed PRF)                     */
    int     staggered;
    int     nMtiTaps;
    double *mtiCoeffs;     /* row-major [nStagger x nMtiTaps]                         */
    int     nRows;         /* rows of the Doppler processing (de-interleaved for staggered PRF) */
    int     rowLen;        /* samples per row                                         */
    int     nfft;          /* Doppler FFT length (power of 2)                         */
    float  *dopWindow;     /* [rowLen] Hamming window, unit energy                    */
    float  *noiseGain;     /* [nfft] noise gain of every Doppler bin (after fftshift) */
    double *dopFreq_Hz;    /* [nfft] Doppler axis (after fftshift)                    */
    int     notchBins;     /* bins on each side of zero Doppler that are not tested   */
    double  medianToMean;  /* median / mean of the noise power (noise floor estimate) */

    /* CFAR, pruning, plot extraction */
    int     cfarMethod;    /* CFAR_2D, CFAR_RANGE or CFAR_DOPPLER: the window below   */
    int     guardRange, trainRange, guardDoppler, trainDoppler, nTrain;
    double  alpha;         /* CA-CFAR threshold factor (design value from cfg.cfar.pfa) */
    double  nLooks;        /* looks of the noise power in a cell (for cfar_pfa; 0 if the
                              file is older and does not have it)                     */
    double  sidelobeMargin_dB;
    int     minHits, maxMiss, updateEvery;
    double  minExtent_deg;  /* azimuth extent of a plot (minHits = this / azimuth step per update) */
    double  meanPri_s;      /* mean PRI                                               */

    /* EKF tracker (optional in the file: off if missing) */
    int     trackEnable, trackDim;   /* dimension 2 (range, azimuth) or 3 (+ elevation)  */
    int     trackModel;     /* TRACK_CV (constant velocity) or TRACK_CA (constant acceleration) */
    double  trackSigma[3];  /* measurement errors: range [m], azimuth, elevation [rad]   */
    double  trackQ[3];      /* process noise per axis: white acceleration [m^2/s^3] (CV)
                               or white jerk [m^2/s^5] (CA) spectral density              */
    double  trackSigmaV0;   /* velocity std of a new track [m/s]                          */
    double  trackSigmaA0;   /* acceleration std of a new track (CA) [m/s^2]               */
    double  trackBatchGap_s, trackBatchMax_s;   /* GNN batch: time without a new plot,
                                                   longest collection time [s]            */
    double  trackGate;      /* chi-square gate of d^2                                     */
    int     trackConfirmHits, trackMaxMissTentative, trackMaxMiss;
    double  scanTime_s;     /* antenna revolution                                         */

    /* run-time option (not in the file): matched filter method,
     * MF_AUTO (the faster one for the thread count), MF_DIRECT or MF_FFT */
    int     mfMethod;
} RadarConfig;

/* Loads the file written by export_receiver_config.m. Returns 0 on success. */
int  radar_config_load(const char *path, RadarConfig *cfg, char *err, size_t errLen);
void radar_config_free(RadarConfig *cfg);

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
} PriData;

typedef struct UdpInput UdpInput;

typedef struct {
    uint64_t datagrams;      /* datagrams received                                    */
    uint64_t foreign;        /* datagrams with a wrong magic number or size           */
    uint64_t prisComplete;   /* PRIs delivered                                        */
    uint64_t prisIncomplete; /* PRIs dropped because packets were missing             */
} UdpStats;

/* Opens the UDP input on 'port' (all interfaces). recordPath (may be NULL)
 * stores every received datagram in a file; replayPath (may be NULL) reads
 * datagrams from such a file instead of the network. Returns NULL on error
 * (message in err). */
UdpInput *udp_input_open(int port, const char *recordPath, const char *replayPath, char *err, size_t errLen);
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

/* A track of the EKF tracker (tracker.c, as track_update.m of the simulator).
 * Position and velocity relative to the antenna: x East, y North, z up (z = 0
 * in 2-D). Reported after every update, and once more with deleted = 1. */
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
} Track;

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
} ProcessorCallbacks;

typedef struct {
    double dbf_ms, mf_ms, mti_ms, doppler_ms, cfar_ms, extract_ms, track_ms, total_ms;
    /* time of the application in its output callbacks (main.c adds it; the
     * processor leaves it 0): console output (printf of plots, tracks and
     * messages) and the data handed to the PPI display */
    double console_ms, display_ms;
    uint64_t pris, updates;
    uint64_t rejected;   /* PRIs whose layout (channels, DDC samples, sub-pulse samples,
                            DDC rate) differs from the configuration: not processed */
    uint64_t skipped;    /* CPI updates skipped (processor_set_slow_time(p, 0))      */
} ProcessorTiming;     /* accumulated times */

typedef struct Processor Processor;
Processor *processor_create(const RadarConfig *cfg, int nThreads);
void processor_process(Processor *p, const PriData *pri, const ProcessorCallbacks *cb);
void processor_flush(Processor *p, const ProcessorCallbacks *cb);   /* close all open plots */
/* enabled = 0: the next PRIs get only the fast-time processing and the MTI
 * (so the CPI buffers stay complete), their Doppler/CFAR/extraction updates
 * are skipped. Used by the input loop when it falls behind the stream. */
void processor_set_slow_time(Processor *p, int enabled);
/* CA-CFAR threshold factor used from the next CPI update on (starts at
 * cfg->alpha). threshold [dB] = 10*log10(alpha). Call from the processing thread. */
void   processor_set_cfar_alpha(Processor *p, double alpha);
double processor_cfar_alpha(const Processor *p);
/* False alarm probability of the CA-CFAR with factor alpha (inverse of
 * cfar_alpha.m): I(1/(1 + alpha/nTrain); nTrain*nLooks, nLooks). NaN if
 * nLooks <= 0. */
double cfar_pfa(double alpha, int nTrain, double nLooks);
void processor_timing(const Processor *p, ProcessorTiming *t);
void processor_destroy(Processor *p);

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
void ppi_post_cfar_threshold(double threshold_dB, double pfa, double default_dB);
int  ppi_take_cfar_request(double *threshold_dB);
/* Display mode at the start of ppi_run (call before it): plots after the plot
 * extractor (default), CFAR detections before it, or both. */
enum { PPI_SHOW_PLOTS = 0, PPI_SHOW_DETECTIONS = 1, PPI_SHOW_BOTH = 2 };
void ppi_set_display_mode(int mode);
/* Track updates and deletions (drawn as hollow triangles, Tracks checkbox). */
void ppi_post_track(const Track *track);
/* Accumulated processing times: the panel shows the mean per PRI of every
 * block [us]. detectionName / trackingName label the CFAR and tracking rows
 * (e.g. "CFAR 2-D", "Tracking 3-D CV"; NULL keeps the previous). */
void ppi_post_timing(const ProcessorTiming *timing, const char *detectionName, const char *trackingName);

/* ------------------------------------------------------------------------ */
/* Platform helpers (platform.c)                                             */
/* ------------------------------------------------------------------------ */

double now_seconds(void);   /* monotonic clock */

#ifdef __cplusplus
}
#endif

#endif /* RADAR_RX_H */

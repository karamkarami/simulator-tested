#ifndef PROCESSOR_INTERNAL_H
#define PROCESSOR_INTERNAL_H

/* Shared storage and interfaces for receiver stages. The public processor API is in radar_rx.h. */
#include "radar_rx.h"
#include "extractor.h"
#include "tracker.h"
#include "cfar.h"
#include "fft.h"

#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define PROCESSOR_PI 3.14159265358979323846
#define CELL_BLOCK 32           /* cells per Doppler FFT block (32 KB: fits the L1 cache) */

#if defined(__clang__)
#  define IVDEP _Pragma("clang loop vectorize(assume_safety)")
#elif defined(__GNUC__)
#  define IVDEP _Pragma("GCC ivdep")
#elif defined(_MSC_VER)
#  define IVDEP __pragma(loop(ivdep))
#else
#  define IVDEP
#endif

/* Array layouts (float, real and imaginary parts in separate arrays):
 *   cell             = beam * nRange + rangeCell
 *   CPI buffer       = [slot][cell]
 *   MTI buffer       = [cellBlock][slot][CELL_BLOCK]   (cell = cellBlock * CELL_BLOCK + c)
 *   power map P      = [dopplerBin][cell]   (Doppler after fftshift)
 */
typedef struct {
    int nR, nCells, nSlots, nMtiOut;
    float *segRe, *segIm;       /* DDC samples of this sub-pulse [ch][n]          */
    float *beamRe, *beamIm;     /* beamformer output [beam][n] (after the SLC: its output) */
    /* sidelobe canceller (cfg->slcEnabled) */
    float  *auxRe, *auxIm;      /* auxiliary elements [aux][n] (zero outside the receive interval) */
    double *trMRe, *trMIm;      /* training samples of the main channels [beam][K] (DBF output)   */
    double *trARe, *trAIm;      /* training samples of the auxiliary elements [aux][K]           */
    double *slcRRe, *slcRIm;    /* QRD: upper triangular R [aux][aux] (QRD-RLS: kept over the PRIs) */
    double *slcZRe, *slcZIm;    /* QRD: z [aux][beam]                                           */
    double *slcWRe, *slcWIm;    /* weights W [aux][beam] of the algorithm                       */
    double *slcBRe, *slcBIm;    /* first exponential mean of W over the PRIs                     */
    double *slcARe, *slcAIm;    /* weights applied: the second mean (of B), y_b = x_b - sum_a
                                   conj(A[a][b]) aux_a                                        */
    double *slcX, *slcD;        /* Givens work: a row [2*aux], the beams [2*beam]                 */
    double *slcEig;             /* smallest eigenvalue: [2aux x 2aux]                            */
    int     slcStarted;         /* QRD-RLS / NLMS: state from an earlier PRI                    */
    float *cpiRe, *cpiIm;       /* matched filter outputs of the newest nMtiTaps
                                   pulses (all the MTI needs) [mfSlot][cell]      */
    float *mtiRe, *mtiIm;       /* MTI outputs [cell block][slot][CELL_BLOCK]     */
    int    last, count, nPushed;            /* CPI buffer state                   */
    int    mfLast;                          /* newest slot of cpiRe/cpiIm         */
    int    sinceUpdate;                     /* pulses since the last performed update */
    int    mtiLast, mtiCount;               /* MTI ring state                     */
    uint16_t *mtiStag;          /* stagger position of every MTI output (of its oldest pulse) */
    uint16_t *staggerIdx;       /* stagger position of every buffered pulse        */
    double   *az, *time;
    uint32_t newestPulse;
    float *P;                   /* power map [bin][cell]                          */
    unsigned char *det;         /* CFAR detections that are local peaks [bin][cell] */
    double *cntInv;             /* CFAR: 1 / training cells of every range cell (more at the edges) */
    double *thrScale;           /* CFAR: alpha * cntInv of every range cell       */
    float  *noise;              /* CFAR noise estimate of the detected cells [bin][cell]
                                   (cfg->snrNoise == SNR_NOISE_CFAR, else NULL)   */
    float  *lastMedian;         /* noise floor: median of every beam's map at its last
                                   computation (bracket of the next one; 0 = none) */
    float *pkM, *pkA;           /* dense sidelobe blanking: peak map, its maxima [bin][cell] */
    float *pkE, *pkD, *pkT;     /* maxima over the beams [bin][r], over Doppler [beam][r], temp */
    Detection *reports;         /* pruned detections of the latest CPI           */
    int nReports, capReports;
    Extractor ext;
    /* matched filter by fast convolution (used when cheaper than direct) */
    int    mfFft, mfN, mfChunks;
    int    mfFour;              /* one beam: FFT of length mfN as 8 FFTs of mfN/8 (four-step) */
    FftPlan mfPlan, mfPlanQ, mfPlan8;   /* lengths mfN; mfN/8 and 8 (four-step)  */
    float *mfHre, *mfHim;       /* FFT of the zero-padded taps, divided by mfN    */
    float *mfTwRe, *mfTwIm;     /* four-step twiddles exp(-2 pi j n1 k2 / mfN) [8][mfN/8] */
    float *mfWork;              /* per chunk: Xre, Xim, Yre, Yim [mfN][MF_LANES]
                                   (four-step: Are, Aim, Bre, Bim [mfN])           */
} PulseState;

typedef struct {                /* work buffers of one thread                      */
    float  *xr, *xi;            /* Doppler: one cell block [nfft][CELL_BLOCK]     */
    CfarScratch cfar;           /* CFAR (cfar.c)                                  */
    float  *work;               /* noise floor: median candidates [nfft*maxNR + 8] */
} Scratch;

struct Processor {
    const RadarConfig *cfg;
    RadarConfig *cfgPrf;        /* the same configuration: its slow-time settings follow the PRF of the
                                   scan (scan-to-scan PRF, radar_config_select_prf)                 */
    FftPlan fft;
    PulseState ps[N_PULSES];
    int nThreads;
    Scratch *scratch;           /* [nThreads] */
    int slowTime;               /* 0: skip the Doppler/CFAR updates (overload)    */
    int logMask;                /* block outputs for the block callback (LOG_DBF, LOG_SLC, LOG_MF, LOG_DOPPLER) */
    double alpha;               /* CFAR threshold factor in use                   */
    Tracker tracker;            /* EKF tracker (cfg->trackEnable)                 */
    Detection *newPlots;        /* plots of the latest extractor update           */
    int nNewPlots, capNewPlots;
    Detection *trackPlots;      /* plots of this PRI for the tracker              */
    int nTrackPlots, capTrackPlots;
    uint32_t lastPulse;
    double lastTime;            /* time of the newest PRI                          */
    int havePulse;
    ProcessorTiming timing;
    /* display taps (processor_set_profile, processor_set_rdmap) */
    int profileKind, profileBeam;
    int rdEnabled, rdBeam;
    double rdPeriod_s, rdLast;
    int nBins;                  /* range bins of the taps (dR from 0 to maxRange)  */
    double dR;
    float *profile;             /* [nBins]                                         */
    int *profileCnt;            /* [nBins]                                         */
    float *rd;                  /* [nfft][nBins]                                   */
    double *rdVel;              /* [nfft]                                          */
    Detection *rdDets;
    int nRdDets, capRdDets, rdUpdated;
    float *rdWork;
};

static inline void *proc_calloc(size_t n, size_t size)
{
    void *p = calloc(n ? n : 1, size);
    if (!p) abort();
    return p;
}

static inline void *proc_realloc(void *p, size_t size)
{
    p = realloc(p, size ? size : 1);
    if (!p) abort();
    return p;
}

static inline int thread_num(void)
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/* Threads for a parallel loop of nWork items: no more than there are items
 * (and no parallel region for one item: if(nWork > 1)).
 * Every thread of a team has to be woken and reach the closing barrier even
 * without work; when the other cores are busy (display, UDP input) that wait
 * costs more than the work. With one beam (no DBF) the loops over beams run
 * in the calling thread. */
static inline int team_size(int nThreads, int nWork)
{
    return nWork < nThreads ? (nWork > 1 ? nWork : 1) : nThreads;
}

/* Row of the MTI ring: the cells c0 ... c0+CELL_BLOCK-1 of cell block blk in
 * slot slot ([block][slot][CELL_BLOCK]: the Doppler processing of a block
 * reads contiguous memory). */
static inline float *mti_row(const PulseState *s, float *ring, int blk, int slot)
{
    return ring + ((size_t)blk * s->nMtiOut + slot) * CELL_BLOCK;
}

static inline float P_at(const PulseState *s, int d, int b, int r)
{
    return s->P[(size_t)d * s->nCells + (size_t)b * s->nR + r];
}

void load_channels(const PulseConfig *pc, const PriData *pri, const int *ch, int nCh, float *outRe, float *outIm);
void load_segment(const RadarConfig *cfg, const PulseConfig *pc, const PriData *pri, PulseState *s);
void beamform(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int nThreads);
void sidelobe_canceller(const RadarConfig *cfg, const PulseConfig *pc, const PriData *pri, PulseState *s,
                        int nThreads);
int mf_fft_init(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int method, int nThreads);
void matched_filter(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s, int slot, int nThreads);

void cpi_buffer_reset(PulseState *s);
int cpi_buffer_begin(PulseState *s, int nMtiTaps);
void cpi_buffer_push(PulseState *s, const RadarConfig *cfg, const PriData *pri, int slot);
void mti_newest(const RadarConfig *cfg, const PulseConfig *pc, PulseState *s);
void doppler(const Processor *p, const PulseConfig *pc, PulseState *s);

void cfar_init_counts(const RadarConfig *cfg, PulseState *s);
void processor_cfar(const Processor *p, PulseState *s, const PulseConfig *pc);
void prune_and_report(const Processor *p, int pulse, PulseState *s);
void peak_interpolation(const Processor *p, int pulse, PulseState *s);
double map_median(const float *base, int nRows, int rowLen, size_t stride, float *work, float *guess);
void init_pack_table(void);
float processor_select_kth(float *v, int n, int k);

void profile_add(Processor *p, int pulse, const PulseState *s);
void profile_post(Processor *p, int profileKind, const PriData *pri, const ProcessorCallbacks *cb);
void rdmap_prf_changed(Processor *p);
void rdmap_add(Processor *p, int pulse, const PulseState *s);
void rdmap_post(Processor *p, const PriData *pri, const ProcessorCallbacks *cb);

#endif

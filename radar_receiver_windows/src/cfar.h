/*
 * cfar.h - CFAR detectors of the receiver (internal): CA, GO, SO and OS with
 * the 2-D, 1-D range and 1-D Doppler window (as cfar_detect.m of the MATLAB
 * simulator), and their false alarm probability.
 */
#ifndef CFAR_H
#define CFAR_H

#include "radar_rx.h"

/* One power map of a sub-pulse and the outputs of the CFAR. */
typedef struct {
    const float *P;         /* power map [dopplerBin][cell], cell = beam * nR + rangeCell     */
    unsigned char *det;     /* out: detections that are also local peaks [bin][cell]          */
    float *noise;           /* out (NULL: not needed): noise estimate of every detection       */
    const double *cntInv;   /* CA: 1 / training cells of every range cell (more at the edges)  */
    const double *thrScale; /* CA: alpha * cntInv of every range cell                         */
    int nR, nCells, nBeams, nfft;
} CfarMap;

/* Work buffers of one thread (cfar_scratch_init: for range axes up to maxNR). */
typedef struct {
    float  *ringO, *ringG;  /* CA: range box sums of the Doppler rows d-oD ... d+1+oD [ringRows][maxNR] */
    float  *partial;        /* CA: partial range sums, 6 x (maxNR + 2*oR + 32)                */
    double *accO, *accG;    /* CA: Doppler running sums [maxNR]                                */
    float  *colMax;         /* local peaks: maxima over the neighbouring beams and bins [maxNR + 2] */
    double *cO, *cG;        /* GO/SO: Doppler sums of every range cell over +-oD, +-gD [maxNR] */
    double *sO, *sG;        /* GO/SO: their prefix sums along range [maxNR + 1]               */
    float  *win;            /* OS: the training cells of one cell [window size]                */
    int ringRows, maxNR;
} CfarScratch;

void cfar_scratch_init(CfarScratch *w, const RadarConfig *cfg, int maxNR);
void cfar_scratch_free(CfarScratch *w);

/* The tested range cells cf ... cl of beam b (all Doppler bins outside the
 * MTI notch) with the detector cfg->cfarType, the window of cfg and the
 * threshold factor alpha: m->det marks the cells that exceed alpha times the
 * statistic of their training cells and are local peaks (3 x 3 x 3 beam x
 * range x Doppler); m->noise gets the noise estimate of every detection. The
 * other cells of m->det are left as they are (cleared by the caller). */
void cfar_chunk(const RadarConfig *cfg, const CfarMap *m, double alpha, int b, int cf, int cl, CfarScratch *w);

#endif

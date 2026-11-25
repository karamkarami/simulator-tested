#ifndef DETECTION_REPORTS_H
#define DETECTION_REPORTS_H

#include "processor_internal.h"

typedef struct {
    int b, r, d;
    float pw;
} DetectionCell;

/* Convert retained cells to reports in their original scan order. */
void build_detection_reports(const Processor *p, int pulse, PulseState *s,
                             const DetectionCell *cells, const unsigned char *keep, int n);

#endif

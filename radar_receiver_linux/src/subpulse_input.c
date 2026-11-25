#include "processor_internal.h"

#include <string.h>

/* DDC samples of the sub-pulse of the channels ch[0 ... nCh-1], as float [k][n]; the samples outside
 * its receive interval (transmission, other sub-pulse, end of the record) are zero. */
void load_channels(const PulseConfig *pc, const PriData *pri, const int *ch, int nCh, float *outRe, float *outIm)
{
    const int n = pc->ddcCount;
    const int k0 = pc->rxFirst > pc->ddcFirst ? pc->rxFirst - pc->ddcFirst : 0;      /* first inside   */
    const int rxEnd = pc->rxFirst + pc->rxCount - pc->ddcFirst;
    const int k1 = rxEnd < n ? rxEnd : n;                                            /* first after it */
    for (int c = 0; c < nCh; c++) {
        const int16_t *src = pri->iq + 2 * ((size_t)ch[c] * pri->nSamples + pc->ddcFirst + k0);   /* sample k0 */
        float *re = outRe + (size_t)c * n, *im = outIm + (size_t)c * n;
        memset(re, 0, sizeof(float) * k0);
        memset(im, 0, sizeof(float) * k0);
        for (int k = k0; k < k1; k++) {
            re[k] = src[2 * (k - k0)];
            im[k] = src[2 * (k - k0) + 1];
        }
        memset(re + k1, 0, sizeof(float) * (n - k1));
        memset(im + k1, 0, sizeof(float) * (n - k1));
    }
}

void load_segment(const RadarConfig *cfg, const PulseConfig *pc, const PriData *pri, PulseState *s)
{
    load_channels(pc, pri, cfg->chDbf, cfg->nChUsed > 0 ? cfg->nChUsed : cfg->nCh, s->segRe, s->segIm);
}

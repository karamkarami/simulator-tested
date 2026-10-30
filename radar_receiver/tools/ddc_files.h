/*
 * ddc_files.h - the DDC files of the MATLAB simulator (ddc_<n>.mat, written
 * with cfg.io.saveDdc = true) and their UDP datagrams.
 *
 * A file holds I and Q (int16 [nDdc x nChannels]) and the header hdr of the
 * PRI. ddc_make_datagrams() splits a PRI into exactly the datagrams that the
 * simulator sends live (udp_make_packets.m): 64-byte header + int16 (I, Q)
 * pairs, channel after channel.
 */
#ifndef DDC_FILES_H
#define DDC_FILES_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint32_t pulseNumber, scanNumber;
    float    az_deg;
    uint16_t staggerIndex;
    double   time_s;            /* PRI start time                               */
    double   pri_s;             /* PRI (time to the next PRI)                   */
    float    ddcFs_Hz, ddcT0_s;
    uint32_t shortIdx[2], longIdx[2];
    int      nSamples, nCh;
    int16_t *iq;                /* I, Q pairs, sample after sample, channel after channel */
    size_t   capacity;          /* int16 values allocated in iq                 */
} DdcPri;

/* The ddc_*.mat files of a folder, sorted by name (= pulse number). Returns
 * the number of files (*paths is malloc'ed, free with ddc_free_list), or -1. */
int  ddc_list_files(const char *folder, char ***paths, char *err, size_t errLen);
void ddc_free_list(char **paths, int n);

/* Reads one DDC file into p (its buffer is reused). Returns 0 on success. */
int  ddc_load(const char *path, DdcPri *p, char *err, size_t errLen);
void ddc_free(DdcPri *p);

/* Number of datagrams and total bytes of a PRI with payloadBytes of samples
 * per datagram. */
int  ddc_datagram_count(const DdcPri *p, int payloadBytes, size_t *totalBytes);
/* Writes the datagrams back to back into out (totalBytes from
 * ddc_datagram_count) and their lengths into len[]. Returns their number. */
int  ddc_make_datagrams(const DdcPri *p, int payloadBytes, unsigned char *out, uint32_t *len);

#endif

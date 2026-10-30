/*
 * platform.c - small platform helpers (monotonic clock).
 */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200809L   /* clock_gettime with -std=c11 */
#endif

#include "radar_rx.h"

#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>

double now_seconds(void)
{
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    /* whole seconds and remainder separately to keep full precision */
    return (double)(count.QuadPart / freq.QuadPart)
         + (double)(count.QuadPart % freq.QuadPart) / (double)freq.QuadPart;
}

#else
#  include <time.h>

double now_seconds(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#endif

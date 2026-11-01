/*
 * platform.c - small platform helpers: monotonic clock; real-time settings
 * of threads and memory (Linux: SCHED_FIFO, CPU affinity, mlockall).
 */
#if defined(__linux__) && !defined(_GNU_SOURCE)
#  define _GNU_SOURCE               /* pthread_setaffinity_np, CPU_SET, mallopt */
#endif
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200809L   /* clock_gettime with -std=c11 */
#endif

#include "radar_rx.h"

#include <stdio.h>

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

/* ------------------------------------------------------------------------ */
/* Real-time settings                                                       */
/* ------------------------------------------------------------------------ */

#if defined(__linux__)
#  include <errno.h>
#  include <malloc.h>
#  include <pthread.h>
#  include <sched.h>
#  include <stdlib.h>
#  include <string.h>
#  include <sys/mman.h>

/* "2-5,7" -> CPU set; 0 on success. */
static int parse_cpus(const char *list, cpu_set_t *set)
{
    CPU_ZERO(set);
    const char *p = list;
    while (*p) {
        char *end;
        long a = strtol(p, &end, 10), b = a;
        if (end == p || a < 0 || a >= CPU_SETSIZE) return -1;
        p = end;
        if (*p == '-') {
            b = strtol(p + 1, &end, 10);
            if (end == p + 1 || b < a || b >= CPU_SETSIZE) return -1;
            p = end;
        }
        for (long c = a; c <= b; c++) CPU_SET((int)c, set);
        if (*p == ',') p++;
        else if (*p) return -1;
    }
    return CPU_COUNT(set) > 0 ? 0 : -1;
}

int rt_thread_setup(int priority, const char *cpus, int *nCpus, char *err, size_t errLen)
{
    int rc = 0;
    if (err && errLen) err[0] = '\0';
    if (cpus && cpus[0]) {
        cpu_set_t set;
        if (parse_cpus(cpus, &set) != 0) {
            if (err) snprintf(err, errLen, "invalid CPU list '%s' (e.g. 2-5,7)", cpus);
            return -2;
        }
        const int e = pthread_setaffinity_np(pthread_self(), sizeof set, &set);
        if (e != 0) {
            if (err) snprintf(err, errLen, "CPUs %s: %s", cpus, strerror(e));
            rc = -1;
        } else if (nCpus) {
            *nCpus = CPU_COUNT(&set);
        }
    }
    if (priority > 0) {
        struct sched_param sp;
        memset(&sp, 0, sizeof sp);
        sp.sched_priority = priority;
        const int e = pthread_setschedparam(pthread_self(), SCHED_FIFO, &sp);
        if (e != 0) {
            const size_t len = err ? strlen(err) : 0;
            if (err) snprintf(err + len, errLen - len, "%sSCHED_FIFO priority %d: %s%s", len ? "; " : "", priority,
                              strerror(e), e == EPERM ? " (allow it: rtprio in /etc/security/limits.conf, or "
                                                        "setcap cap_sys_nice,cap_ipc_lock+ep radar_receiver)" : "");
            rc = -1;
        }
    }
    return rc;
}

int rt_lock_memory(char *err, size_t errLen)
{
    /* freed memory stays in the heap and large blocks come from it as well,
     * so later allocations need no system call and no new pages */
    mallopt(M_TRIM_THRESHOLD, -1);
    mallopt(M_MMAP_MAX, 0);
    if (mlockall(MCL_CURRENT | MCL_FUTURE) != 0) {
        if (err) snprintf(err, errLen, "mlockall: %s%s", strerror(errno),
                          errno == ENOMEM || errno == EPERM ? " (allow it: memlock unlimited in "
                          "/etc/security/limits.conf, or setcap cap_ipc_lock+ep radar_receiver)" : "");
        return -1;
    }
    return 0;
}

#else

int rt_thread_setup(int priority, const char *cpus, int *nCpus, char *err, size_t errLen)
{
    (void)nCpus;
    if (priority <= 0 && (!cpus || !cpus[0])) return 0;
    if (err) snprintf(err, errLen, "real-time priority and CPU affinity are only implemented for Linux");
    return -1;
}

int rt_lock_memory(char *err, size_t errLen)
{
    if (err) snprintf(err, errLen, "locking the memory is only implemented for Linux");
    return -1;
}

#endif

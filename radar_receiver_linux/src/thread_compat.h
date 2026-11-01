/*
 * thread_compat.h - mutex, condition variable and thread helpers for Windows
 * and POSIX (used by input_thread.c and the stream sender).
 */
#ifndef THREAD_COMPAT_H
#define THREAD_COMPAT_H

#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200809L
#endif

#ifdef _WIN32
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
typedef CRITICAL_SECTION   Mutex;
typedef CONDITION_VARIABLE Cond;
#  define mutex_init(m)       InitializeCriticalSection(m)
#  define mutex_destroy(m)    DeleteCriticalSection(m)
#  define mutex_lock(m)       EnterCriticalSection(m)
#  define mutex_unlock(m)     LeaveCriticalSection(m)
#  define cond_init(c)        InitializeConditionVariable(c)
#  define cond_destroy(c)     ((void)0)
#  define cond_signal(c)      WakeConditionVariable(c)
#  define cond_wait_ms(c, m, ms) SleepConditionVariableCS(c, m, (DWORD)(ms))
#else
#  include <pthread.h>
#  include <time.h>
typedef pthread_mutex_t Mutex;
typedef pthread_cond_t  Cond;
#  define mutex_init(m)       pthread_mutex_init(m, NULL)
#  define mutex_destroy(m)    pthread_mutex_destroy(m)
#  define mutex_lock(m)       pthread_mutex_lock(m)
#  define mutex_unlock(m)     pthread_mutex_unlock(m)
#  define cond_init(c)        pthread_cond_init(c, NULL)
#  define cond_destroy(c)     pthread_cond_destroy(c)
#  define cond_signal(c)      pthread_cond_signal(c)
static inline void cond_wait_ms(Cond *c, Mutex *m, int ms)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_sec += ms / 1000;
    ts.tv_nsec += (long)(ms % 1000) * 1000000L;
    if (ts.tv_nsec >= 1000000000L) { ts.tv_sec++; ts.tv_nsec -= 1000000000L; }
    pthread_cond_timedwait(c, m, &ts);
}
#endif

#endif

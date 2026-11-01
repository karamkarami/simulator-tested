/*
 * input_thread.c - reception of the PRIs in a thread of its own.
 *
 * The thread calls udp_input_next() and puts every complete PRI into a
 * queue (a ring of PriData whose buffers are reused). The processing thread
 * takes them from the queue. The socket is thus emptied while a long
 * Doppler/CFAR update runs, and the receive buffer of the operating system
 * does not overflow. When the queue is full, a live stream drops the newest
 * PRI (counted); a replay waits instead, so nothing is lost.
 */
#if !defined(_WIN32) && !defined(_POSIX_C_SOURCE)
#  define _POSIX_C_SOURCE 200809L
#endif

#include "radar_rx.h"

#include <stdlib.h>
#include <string.h>

#include "thread_compat.h"

struct InputThread {
    UdpInput *in;
    PriData  *slots;            /* ring of capacity PRIs                          */
    PriData   spare;            /* receives the PRIs that are dropped             */
    int       capacity, head, count;   /* head: oldest queued PRI                 */
    int       held;             /* 1 while the consumer holds the oldest slot     */
    int       dropWhenFull, stop, ended;
    uint64_t  dropped;
    Mutex     lock;
    Cond      notEmpty, notFull;
#ifdef _WIN32
    HANDLE    thread;
#else
    pthread_t thread;
#endif
};

static void reader(InputThread *t)
{
    for (;;) {
        mutex_lock(&t->lock);
        while (!t->stop && t->count + t->held >= t->capacity && !t->dropWhenFull)
            cond_wait_ms(&t->notFull, &t->lock, 200);
        if (t->stop) { mutex_unlock(&t->lock); break; }
        const int full = t->count + t->held >= t->capacity;
        PriData *slot = full ? &t->spare : &t->slots[(t->head + t->held + t->count) % t->capacity];
        mutex_unlock(&t->lock);

        const int r = udp_input_next(t->in, slot, 200);       /* without the lock */

        mutex_lock(&t->lock);
        if (r < 0) t->ended = 1;
        if (r == 1) {
            if (full) t->dropped++;
            else t->count++;
        }
        cond_signal(&t->notEmpty);
        const int done = t->ended;
        mutex_unlock(&t->lock);
        if (done) break;
    }
}

#ifdef _WIN32
static DWORD WINAPI reader_main(LPVOID arg) { reader(arg); return 0; }
#else
static void *reader_main(void *arg) { reader(arg); return NULL; }
#endif

InputThread *input_thread_start(UdpInput *in, int capacity, int dropWhenFull)
{
    InputThread *t = calloc(1, sizeof *t);
    if (!t) return NULL;
    t->in = in;
    t->capacity = capacity >= 2 ? capacity : 2;
    t->dropWhenFull = dropWhenFull;
    t->slots = calloc((size_t)t->capacity, sizeof(PriData));
    if (!t->slots) { free(t); return NULL; }
    mutex_init(&t->lock);
    cond_init(&t->notEmpty);
    cond_init(&t->notFull);
#ifdef _WIN32
    t->thread = CreateThread(NULL, 0, reader_main, t, 0, NULL);
    const int ok = t->thread != NULL;
#else
    const int ok = pthread_create(&t->thread, NULL, reader_main, t) == 0;
#endif
    if (!ok) {
        mutex_destroy(&t->lock);
        cond_destroy(&t->notEmpty);
        cond_destroy(&t->notFull);
        free(t->slots);
        free(t);
        return NULL;
    }
    return t;
}

int input_thread_next(InputThread *t, PriData **pri, int timeout_ms)
{
    mutex_lock(&t->lock);
    if (t->held) {                                       /* release the previous PRI */
        t->held = 0;
        t->head = (t->head + 1) % t->capacity;
        cond_signal(&t->notFull);
    }
    if (t->count == 0 && !t->ended && timeout_ms != 0)
        cond_wait_ms(&t->notEmpty, &t->lock, timeout_ms < 0 ? 1000 : timeout_ms);
    int r = 0;
    if (t->count > 0) {
        *pri = &t->slots[t->head];
        t->count--;
        t->held = 1;
        r = 1;
    } else if (t->ended) {
        r = -1;
    }
    mutex_unlock(&t->lock);
    return r;
}

int input_thread_backlog(InputThread *t)
{
    mutex_lock(&t->lock);
    const int n = t->count;
    mutex_unlock(&t->lock);
    return n;
}

uint64_t input_thread_dropped(InputThread *t)
{
    mutex_lock(&t->lock);
    const uint64_t n = t->dropped;
    mutex_unlock(&t->lock);
    return n;
}

void input_thread_stop(InputThread *t)
{
    if (!t) return;
    mutex_lock(&t->lock);
    t->stop = 1;
    cond_signal(&t->notFull);
    mutex_unlock(&t->lock);
#ifdef _WIN32
    WaitForSingleObject(t->thread, INFINITE);
    CloseHandle(t->thread);
#else
    pthread_join(t->thread, NULL);
#endif
    for (int i = 0; i < t->capacity; i++) pri_data_free(&t->slots[i]);
    pri_data_free(&t->spare);
    mutex_destroy(&t->lock);
    cond_destroy(&t->notEmpty);
    cond_destroy(&t->notFull);
    free(t->slots);
    free(t);
}

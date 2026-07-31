// MIT License

/**
 * @file progress.c
 * @brief Implementation of the stage-by-stage progress log.
 */

#include <pthread.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <pari/pari.h>
#include "../headers/progress.h"

static long progress_level = 0;
static double progress_started = 0.0;
static long heartbeat_seconds = 300;

/* The most recent stage, kept for the heartbeat thread to quote.  Both
 * threads write whole lines, so one mutex serializes the stage and the
 * output together. */
static pthread_mutex_t progress_lock = PTHREAD_MUTEX_INITIALIZER;
static char progress_stage[256] = "starting up";
static pthread_t heartbeat_thread;
static volatile int heartbeat_running = 0;

static double
progress_now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Resident set size in GiB.  This is deliberately not the PARI stack: what
 * ends a run is the memory the operating system sees, and unlike `avma` it
 * can be read from a second thread without touching PARI's state. */
static double
progress_rss(void)
{
    FILE *fh = fopen("/proc/self/statm", "r");
    unsigned long total, resident = 0;
    if (!fh) return 0.0;
    if (fscanf(fh, "%lu %lu", &total, &resident) != 2) resident = 0;
    fclose(fh);
    return (double)resident * (double)sysconf(_SC_PAGESIZE) / 1073741824.0;
}

static void
progress_stamp(double elapsed, double rss)
{
    printf("[+%02d:%02d:%02d  RSS %.2f GiB]  ",
           (int)(elapsed / 3600.0), ((int)elapsed / 60) % 60,
           (int)elapsed % 60, rss);
}

static void *
heartbeat_main(void *unused)
{
    (void)unused;
    while (heartbeat_running)
    {
        long slept = 0;
        /* Sleep in one-second slices so the thread leaves promptly. */
        while (heartbeat_running && slept < heartbeat_seconds)
        { sleep(1); ++slept; }
        if (!heartbeat_running) break;
        pthread_mutex_lock(&progress_lock);
        progress_stamp(progress_now() - progress_started, progress_rss());
        printf("still running: %s\n", progress_stage);
        fflush(stdout);
        pthread_mutex_unlock(&progress_lock);
    }
    return NULL;
}

void
my_progress_init(void)
{
    const char *level = getenv("MASSEY_LOG_LEVEL");
    const char *beat = getenv("MASSEY_LOG_HEARTBEAT");
    progress_level = (level && *level) ? atol(level) : 0;
    if (beat && *beat) heartbeat_seconds = atol(beat);
    if (heartbeat_seconds < 1) heartbeat_seconds = 1;
    progress_started = progress_now();

    /* PARI's diagnostics are selected per domain, not by one global level.
     * Level 2 takes the class field construction alone, and only its first
     * level: that reports how long the auxiliary field K(zeta_p) took, which
     * is the measurement the cost of a hard field turns on.  Its higher
     * levels would name the steps of the Kummer algorithm, but not without
     * also dumping the intermediate elements, which run to kilobytes each;
     * the heartbeat says what is running at far less cost.  Relation
     * collection (domain "bnf") is the bulk of the output and almost none of
     * the meaning, so it waits for level 3, which turns everything on. */
    if (progress_level >= 3)
        setalldebug(1);
    else if (progress_level >= 2)
    {
        setalldebug(0);
        setdebug("bnrclassfield", 1);
    }
    else
        setalldebug(0);

    if (progress_level >= 2)
    {
        heartbeat_running = 1;
        if (pthread_create(&heartbeat_thread, NULL, heartbeat_main, NULL))
            heartbeat_running = 0;
    }
}

void
my_progress_stop(void)
{
    if (!heartbeat_running) return;
    heartbeat_running = 0;
    pthread_join(heartbeat_thread, NULL);
}

long
my_progress_level(void)
{
    return progress_level;
}

void
my_progress(const char *fmt, ...)
{
    va_list ap;

    if (progress_level < 1) return;

    pthread_mutex_lock(&progress_lock);
    va_start(ap, fmt);
    vsnprintf(progress_stage, sizeof(progress_stage), fmt, ap);
    va_end(ap);
    progress_stamp(progress_now() - progress_started, progress_rss());
    printf("%s\n", progress_stage);
    fflush(stdout);
    pthread_mutex_unlock(&progress_lock);
}

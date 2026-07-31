// MIT License

/**
 * @file progress.h
 * @brief Stage-by-stage progress log for long-running field computations.
 *
 * One character of a hard field can run for hours, and until now its log
 * recorded only the beginning and the end.  When such a run dies -- as the
 * runs for D = -109909943 did -- the log says nothing about which part of
 * the computation was in progress or how the memory got there.  This module
 * names the stage that is running, and above level 1 it also keeps a pulse
 * so that a silent hour can be told from a stalled one.
 *
 * The level is read once from the environment variable MASSEY_LOG_LEVEL:
 *
 *   0  silent (the default): no stage lines, and the arithmetic and its
 *      output are those of every run recorded so far; the driver banner
 *      carries the level, so a log always states its own verbosity
 *   1  one line per stage of the computation
 *   2  additionally a heartbeat every MASSEY_LOG_HEARTBEAT seconds (default
 *      300) naming the stage still in progress, and PARI's own account of
 *      the class field construction: the steps of the Kummer algorithm and
 *      the cost of the auxiliary field K(zeta_p)
 *   3  additionally every other PARI diagnostic, which is voluminous:
 *      relation collection, factor base, Bach constant, internal resultant
 *      bounds
 *
 * Every line carries the elapsed time and the resident set size, the two
 * numbers that decided the fate of the expensive computations so far.  Lines
 * go to stdout and are flushed immediately; the parallel runner pipes them
 * through a stamper that prefixes the wall-clock time.
 */

#ifndef MASSEY_PROGRESS_H
#define MASSEY_PROGRESS_H

/** Read the environment, start the clock, select PARI's debug domains and,
 *  above level 1, start the heartbeat. */
void my_progress_init(void);

/** Stop the heartbeat and join it; safe to call when none was started. */
void my_progress_stop(void);

/** The active level, for callers that want to skip expensive log arguments. */
long my_progress_level(void);

/** Log one stage; a no-op below level 1.  No trailing newline in @p fmt. */
void my_progress(const char *fmt, ...);

#endif

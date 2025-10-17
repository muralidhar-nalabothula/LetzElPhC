/**
 * @file ELPH_timers.c
 * @brief Wall-clock timing utilities for performance profiling
 * 
 * Implements a timer system using MPI_Wtime() for measuring execution times
 * of different code sections. Timers are identified by string keys and can
 * be started/stopped multiple times (accumulating total time and count).
 */

#include "ELPH_timers.h"
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include "common/ELPH_hash_map/ELPH_hmap.h"

/**
 * @struct ELPH_timer
 * @brief Internal structure to store timing information for a named timer
 */
struct ELPH_timer
{
    double Wtime;  /**< Accumulated wall-clock time in seconds */
    size_t count;  /**< Number of times this timer has been stopped */
};

/**
 * @typedef map_timer_t
 * @brief Hash map type for storing timers with string keys
 */
typedef map_t(struct ELPH_timer) map_timer_t;

/**
 * @brief Global hash map storing all active timers
 */
static map_timer_t timer_map;

/**
 * @brief Initializes the global timer system
 * 
 * Must be called before using any other timer functions. Also automatically
 * starts the "Total time" timer to track overall program execution.
 * 
 * @note This function should be called once at program initialization
 */
void init_ELPH_clocks(void)
{
    map_init(&timer_map);
    ELPH_start_clock("Total time");
}

/**
 * @brief Starts or restarts a named timer
 * 
 * If the timer doesn't exist, creates it. If it already exists, records
 * the current time to resume accumulating from the previous total.
 * The implementation stores a negative timestamp to enable accumulation:
 * \f$ t_{stored} = -(t_{current} - t_{previous}) \f$
 * 
 * @param str Name/identifier for the timer (must not be NULL)
 * 
 * @note Safe to call multiple times for the same timer
 * @note Does nothing if str is NULL
 */
void ELPH_start_clock(const char *str)
{
    if (!str)
    {
        return;
    }
    double tic = MPI_Wtime();
    struct ELPH_timer *etime = map_get(&timer_map, str);
    size_t count_tmp = 0;
    if (etime)
    {
        tic -= etime->Wtime;
        count_tmp = etime->count;
    }
    tic = -tic;
    struct ELPH_timer time_set;
    time_set.count = count_tmp;
    time_set.Wtime = tic;
    map_set(&timer_map, str, time_set);
}

/**
 * @brief Stops a named timer and accumulates elapsed time
 * 
 * Computes elapsed time since the corresponding start call and adds it to
 * the timer's accumulated total. Also increments the call count.
 * The elapsed time is computed as:
 * \f$ t_{elapsed} = t_{stop} + t_{stored} \f$
 * where \f$ t_{stored} \f$ is the negative value from start_clock.
 * 
 * @param str Name/identifier for the timer (must not be NULL)
 * 
 * @note Does nothing if str is NULL or timer doesn't exist
 * @note Must be preceded by a corresponding ELPH_start_clock() call
 */
void ELPH_stop_clock(const char *str)
{
    if (!str)
    {
        return;
    }
    double tok = MPI_Wtime();
    struct ELPH_timer *etime = map_get(&timer_map, str);
    size_t count_tmp = 0;
    if (etime)
    {
        tok += etime->Wtime;
        count_tmp = etime->count + 1;
    }
    else
    {
        return;
    }
    struct ELPH_timer time_set;
    time_set.count = count_tmp;
    time_set.Wtime = tok;
    map_set(&timer_map, str, time_set);
}

/**
 * @brief Cleans up and deallocates all timer resources
 * 
 * Frees all memory used by the timer system. No timer functions should
 * be called after this.
 * 
 * @note Should be called at program termination
 */
void cleanup_ELPH_clocks(void)
{
    map_deinit(&timer_map);
}

/**
 * @brief Prints formatted summary of all timers to stdout
 * 
 * Automatically stops the "Total time" timer before printing.
 * Displays all timers except "Total time" in a formatted table showing
 * accumulated wall time and call counts, followed by the total execution time.
 * 
 * Output format:
 * - Function name (left-aligned, 20 chars)
 * - Wall time in seconds (12.4f format)
 * - Call count (8 digits)
 * 
 * @note Should typically be called once before program termination
 */
void print_ELPH_clock_summary(void)
{
    ELPH_stop_clock("Total time");
    const char *key;
    map_iter_t iter = map_iter(&timer_map);
    fputs("\n", stdout);
    fputs("===================== Wall times =====================\n", stdout);
    fprintf(stdout, "%-20s  :    Wtime (s)     (  counts  )\n",
            "Function name");
    fputs("------------------------------------------------------\n", stdout);
    while ((key = map_next(&timer_map, &iter)))
    {
        if (0 == strcmp(key, "Total time"))
        {
            continue;
        }
        struct ELPH_timer *etime = map_get(&timer_map, key);
        fprintf(stdout, "%-20s  : %12.4f     ( %8zu )\n", key, etime->Wtime,
                etime->count);
    }
    fputs("------------------------------------------------------\n", stdout);
    key = "Total time";
    fprintf(stdout, "%-20s  : %12.4f s.\n", key,
            map_get(&timer_map, key)->Wtime);
    fputs("\n", stdout);
}

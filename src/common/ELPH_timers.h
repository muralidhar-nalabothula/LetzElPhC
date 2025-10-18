/**
 * @file ELPH_timers.h
 * @brief Wall-clock timing utilities for performance profiling
 *
 * Provides functions for measuring and reporting execution times of code
 * sections. Timers are identified by string names and support multiple
 * start/stop cycles with automatic accumulation and counting.
 *
 * Typical usage:
 * @code
 * init_ELPH_clocks();
 * ELPH_start_clock("my_function");
 * // ... code to time ...
 * ELPH_stop_clock("my_function");
 * print_ELPH_clock_summary();
 * cleanup_ELPH_clocks();
 * @endcode
 */

#pragma once

/**
 * @brief Initializes the global timer system
 *
 * Must be called once before using any timer functions.
 */
void init_ELPH_clocks(void);

/**
 * @brief Starts or restarts a named timer
 *
 * @param str Timer identifier string
 */
void ELPH_start_clock(const char* str);

/**
 * @brief Stops a named timer and accumulates elapsed time
 *
 * @param str Timer identifier string
 */
void ELPH_stop_clock(const char* str);

/**
 * @brief Cleans up all timer resources
 *
 * Should be called at program termination.
 */
void cleanup_ELPH_clocks(void);

/**
 * @brief Prints formatted summary of all timers to stdout
 *
 * Displays accumulated times and call counts for all timers.
 */
void print_ELPH_clock_summary(void);

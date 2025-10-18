/**
 * @file
 * @brief Progress bar structure and function declarations
 */

#pragma once
#include "elphC.h"

/**
 * @def PROGRESS_BAR_LEN
 * @brief Number of segments in progress bar display
 */
#define PROGRESS_BAR_LEN 5

/**
 * @def PROGRESS_BAR_SIGN
 * @brief Filled progress indicator string
 */
#define PROGRESS_BAR_SIGN "#####"

/**
 * @def PROGRESS_BAR_SPACE
 * @brief Empty progress indicator string
 */
#define PROGRESS_BAR_SPACE "     "

/**
 * @struct progress_bar
 * @brief Progress tracking state for iterative calculations
 */
struct progress_bar
{
    int rank;             /**< MPI rank (only rank 0 prints) */
    ND_int iiter;         /**< Current iteration number */
    ND_int niter;         /**< Total number of iterations */
    ND_int isign;         /**< Number of progress segments filled */
    double prev_time;     /**< Wall time when previous iteration finished */
    double elap_time;     /**< Total elapsed time in seconds */
    double max_iter_time; /**< Maximum time taken by any single iteration */
};

/**
 * @brief Initializes progress bar and prints empty bar
 *
 * @param pbar Progress bar structure to initialize
 * @param mpi_rank MPI rank of calling process
 * @param niter Total number of iterations
 */
void start_progressbar(struct progress_bar* pbar, int mpi_rank, ND_int niter);

/**
 * @brief Updates and prints progress bar
 *
 * @param pbar Progress bar structure to update
 */
void print_progressbar(struct progress_bar* pbar);

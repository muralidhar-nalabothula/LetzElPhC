/**
 * @file mpi_parallel.h
 * @brief Functions for parallel communication and domain decomposition using
 * MPI.
 *
 * This file contains utility functions for determining local array sizes and
 * indices for parallel processing, as well as functions for setting up and
 * tearing down MPI communicators used throughout the ELPH project.
 */
#pragma once
#include <mpi.h>

#include "dtypes.h"
#include "elphC.h"

/**
 * @brief Calculates the local size and starting index for an array distributed
 * among MPI processes.
 *
 * This function is used for 1D domain decomposition, ensuring a relatively
 * even distribution of $n$ elements across the processes in the given
 * communicator.
 *
 * @param n The total number of elements to be distributed.
 * @param start_idx Pointer to an @ref ND_int where the starting index of the
 * local chunk will be stored.
 * @param Comm The MPI communicator over which the distribution is performed.
 * @return @ref ND_int The local size (number of elements) for the calling
 * process.
 */
ND_int get_mpi_local_size_idx(const ND_int n, ND_int* start_idx, MPI_Comm Comm);

/**
 * @brief Creates a set of parallel communicators for the ELPH calculation.
 *
 * This function sets up specialized MPI communicators based on the specified
 * number of q-point and k-point pools for efficient parallel processing in the
 * ELPH code. The pools typically partition the world communicator.
 *
 * @param nqpools The number of q-point pools to create.
 * @param nkpools The number of k-point pools to create.
 * @param MPI_world_comm The global MPI communicator (usually @c MPI_COMM_WORLD)
 * to be split.
 * @param Comm Pointer to the @ref ELPH_MPI_Comms structure where the newly
 * created communicators will be stored.
 */
void create_parallel_comms(const int nqpools, const int nkpools,
                           const MPI_Comm MPI_world_comm,
                           struct ELPH_MPI_Comms* Comm);

/**
 * @brief Frees the parallel communicators stored in the @ref ELPH_MPI_Comms
 * structure.
 *
 * This function calls @c MPI_Comm_free on all communicators that were created
 * by @ref create_parallel_comms.
 *
 * @param Comm Pointer to the @ref ELPH_MPI_Comms structure containing the
 * communicators to be freed.
 */
void free_parallel_comms(struct ELPH_MPI_Comms* Comm);

/**
 * @brief Distributes a total number of elements (@p n) into groups.
 *
 * This function calculates the number of elements assigned to a specific group
 * (@p igrp) out of a total number of groups (@p ngrp) and determines the
 * starting shift (index) for that group.
 *
 * The number of elements assigned to group $i$ (where $i = \mathtt{igrp}$) is
 * calculated as:
 * $$
 * \text{size}_i = \left\lfloor \frac{n}{\text{ngrp}} \right\rfloor +
 * \mathbb{I}(i < (n \bmod \text{ngrp}))
 * $$
 * where $\mathbb{I}(\cdot)$ is the indicator function. The starting shift is
 * the cumulative sum of the sizes of all preceding groups.
 *
 * @param n The total number of elements to distribute.
 * @param ngrp The total number of groups (e.g., number of MPI processes/pools).
 * @param igrp The zero-based index of the current group ( $0 \le \mathtt{igrp}
 * < \mathtt{ngrp}$).
 * @param shift Pointer to an @ref ND_int where the starting shift (index) for
 * the current group will be stored.
 * @return @ref ND_int The size (number of elements) assigned to the current
 * group (@p igrp).
 */
ND_int distribute_to_grps(const ND_int n, const ND_int ngrp, const ND_int igrp,
                          ND_int* shift);

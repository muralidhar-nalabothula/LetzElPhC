/**
 * @file
 * @brief MPI process distribution and communicator management
 *
 * Implements functions for distributing work among MPI processes and
 * creating hierarchical communicator topology for q-pool and k-pool
 * parallelization schemes.
 */

#include "parallel.h"

#include <mpi.h>
#include <stddef.h>

#include "dtypes.h"
#include "elphC.h"
#include "error.h"

/**
 * @brief Computes local size and starting index for distributed dimension
 *
 * Distributes n elements among MPI processes in a communicator using
 * a balanced distribution scheme. Assigns extra elements to lower-ranked
 * processes when n is not evenly divisible.
 *
 * Distribution formula:
 * - Base size per process: q = floor(n / size)
 * - Remainder: r = n mod size
 * - Processes 0 to r-1 get (q+1) elements
 * - Processes r to size-1 get q elements
 *
 * Starting index for rank i:
 * - If i < r: start = i * (q+1) = q*i + i
 * - If i >= r: start = q*i + r
 *
 * @param n Total number of elements to distribute
 * @param start_idx Output pointer for starting index (can be NULL)
 * @param Comm MPI communicator
 * @return Number of elements assigned to calling process
 */
ND_int get_mpi_local_size_idx(const ND_int n, ND_int* start_idx, MPI_Comm Comm)
{
    int my_rank, total_size;
    int mpi_error;
    mpi_error = MPI_Comm_rank(Comm, &my_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm, &total_size);
    MPI_error_msg(mpi_error);

    ND_int n_q = n / total_size;
    ND_int n_r = n % total_size;

    ND_int n_this_cpu = n_q;
    if (my_rank < n_r)
    {
        ++n_this_cpu;
    }

    if (start_idx != NULL)
    {
        *start_idx = n_q * my_rank;
        if (my_rank < n_r)
        {
            *start_idx += my_rank;
        }
        else
        {
            *start_idx += n_r;
        }
    }
    return n_this_cpu;
}

/**
 * @brief Distributes elements among groups and computes offset
 *
 * Similar to get_mpi_local_size_idx but for generic group distribution
 * without MPI communication. Used for distributing k/q points among pools.
 *
 * @param n Total number of elements to distribute
 * @param ngrp Total number of groups
 * @param igrp Group index (0 to ngrp-1)
 * @param shift Output pointer for starting index/shift
 * @return Number of elements assigned to group igrp
 */
ND_int distribute_to_grps(const ND_int n, const ND_int ngrp, const ND_int igrp,
                          ND_int* shift)
{
    ND_int n_q = n / ngrp;
    ND_int n_r = n % ngrp;

    ND_int n_this_grp = n_q;
    if (igrp < n_r)
    {
        ++n_this_grp;
    }

    *shift = n_q * igrp;
    if (igrp < n_r)
    {
        *shift += igrp;
    }
    else
    {
        *shift += n_r;
    }

    return n_this_grp;
}

/**
 * @brief Creates hierarchical MPI communicator topology
 *
 * Establishes a two-level parallelization hierarchy:
 * - Level 1: Q-pools (nqpools groups)
 * - Level 2: K-pools within each Q-pool (nkpools groups)
 *
 * Creates five communicators:
 * - commQ: Processes within same q-pool
 * - commK: Processes within same k-pool
 * - commR: Processes with same k-pool rank across all q-pools
 * - commRq: Processes with same k-pool rank within a q-pool
 *
 * @param nqpools Number of q-point pools
 * @param nkpools Number of k-point pools per q-pool
 * @param MPI_world_comm Global MPI communicator (typically MPI_COMM_WORLD)
 * @param Comm Output structure to store all communicators and metadata
 *
 * @note Total processes must be divisible by (nqpools * nkpools)
 * @note The key parameter in MPI_Comm_split preserves original rank ordering
 */
void create_parallel_comms(const int nqpools, const int nkpools,
                           const MPI_Comm MPI_world_comm,
                           struct ELPH_MPI_Comms* Comm)
{
    int mpi_error;
    // first set the basic things
    Comm->commW = MPI_world_comm;
    Comm->nqpools = nqpools;
    Comm->nkpools = nkpools;

    mpi_error = MPI_Comm_rank(Comm->commW, &Comm->commW_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm->commW, &Comm->commW_size);
    MPI_error_msg(mpi_error);

    // check Comm->commW_size is divisible by nqpools*nkpools
    if (nqpools < 1 || nkpools < 1)
    {
        error_msg(" number of pools must be greater than 0 ");
    }
    if (Comm->commW_size % (nqpools * nkpools) != 0)
    {
        error_msg(" product of kpools and qpools must divide total cpus.");
    }

    Comm->commQ_size = Comm->commW_size / nqpools;
    Comm->commK_size = Comm->commQ_size / nkpools;

    /* create comms */
    // commQ
    mpi_error = MPI_Comm_split(Comm->commW, Comm->commW_rank / Comm->commQ_size,
                               Comm->commW_rank, &Comm->commQ);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commQ, &Comm->commQ_rank);
    MPI_error_msg(mpi_error);

    // commK
    mpi_error = MPI_Comm_split(Comm->commQ, Comm->commQ_rank / Comm->commK_size,
                               Comm->commW_rank, &Comm->commK);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commK, &Comm->commK_rank);
    MPI_error_msg(mpi_error);
    // commR
    mpi_error = MPI_Comm_split(Comm->commW, Comm->commK_rank, Comm->commW_rank,
                               &Comm->commR);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commR, &Comm->commR_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm->commR, &Comm->commR_size);
    MPI_error_msg(mpi_error);
    // commRq
    mpi_error = MPI_Comm_split(Comm->commQ, Comm->commK_rank, Comm->commW_rank,
                               &Comm->commRq);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_rank(Comm->commRq, &Comm->commRq_rank);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_size(Comm->commRq, &Comm->commRq_size);
    MPI_error_msg(mpi_error);

    // Sanity check
    if (Comm->commK_rank != Comm->commW_rank % Comm->commK_size)
    {
        error_msg("MPI Comm 1 creation failed.");
    }
    if (Comm->commR_size != (nqpools * nkpools))
    {
        error_msg("MPI Comm 2 creation failed.");
    }
    if (Comm->commRq_size != nkpools)
    {
        error_msg("MPI Comm 3 creation failed.");
    }
}

/**
 * @brief Frees all allocated MPI communicators
 *
 * Releases communicators created by create_parallel_comms.
 * Does not free commW as it is not owned by this structure.
 *
 * @param Comm Communicator structure to clean up
 */
void free_parallel_comms(struct ELPH_MPI_Comms* Comm)
{
    int mpi_error;
    mpi_error = MPI_Comm_free(&Comm->commQ);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_free(&Comm->commK);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_free(&Comm->commR);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Comm_free(&Comm->commRq);
    MPI_error_msg(mpi_error);
}

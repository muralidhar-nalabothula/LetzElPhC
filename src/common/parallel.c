/**
 * @file mpi_parallel.c
 * @brief Implementation of parallel communication and domain decomposition
 * utilities using MPI.
 *
 * This file contains the implementation for functions related to distributing
 * array dimensions and indices across MPI processes, as well as the creation
 * and management of specialized MPI communicators for the ELPH calculation.
 */
/*
This file contains function which distributes cpus
*/
#include "parallel.h"

#include <mpi.h>
#include <stddef.h>

#include "dtypes.h"
#include "elphC.h"
#include "error.h"

/**
 * @brief Calculates the local size and starting index for an array distributed
 * among MPI processes.
 *
 * This function implements a 1D block distribution of $n$ elements across
 * processes in @p Comm, ensuring that the first $n \bmod \text{size}$ processes
 * get one extra element compared to the rest.
 *
 * The local size $\text{size}_{\text{rank}}$ for a process with rank
 * $\text{rank}$ in a communicator of size $\text{total\_size}$ is:
 * $$
 * \text{size}_{\text{rank}} = \left\lfloor \frac{n}{\text{total\_size}}
 * \right\rfloor +
 * \mathbb{I}(\text{rank} < (n \bmod \text{total\_size}))
 * $$
 *
 * The starting index $\text{start}_{\text{rank}}$ is the cumulative sum of the
 * sizes of all preceding processes.
 *
 * @param n The total number of elements to be distributed.
 * @param start_idx Pointer to an @ref ND_int where the starting index of the
 * local chunk will be stored. Can be @c NULL if the starting index is not
 * needed.
 * @param Comm The MPI communicator over which the distribution is performed.
 * @return @ref ND_int The local size (number of elements) for the calling
 * process.
 */
ND_int get_mpi_local_size_idx(const ND_int n, ND_int *start_idx, MPI_Comm Comm)
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
 * @brief Distributes a total number of elements (@p n) into groups.
 *
 * This function calculates the number of elements assigned to a specific group
 * (@p igrp) out of a total number of groups (@p ngrp) and determines the
 * starting shift (index) for that group. This implements a 1D block
 * distribution where the first $n \bmod \text{ngrp}$ groups get one extra
 * element.
 *
 * The number of elements assigned to group $i$ (where $i = \mathtt{igrp}$) is
 * calculated as:
 * $$
 * \text{size}_i = \left\lfloor \frac{n}{\text{ngrp}} \right\rfloor +
 * \mathbb{I}(i < (n \bmod \text{ngrp}))
 * $$
 *
 * @param n The total number of elements to distribute (e.g., number of k-points
 * or q-points).
 * @param ngrp The total number of groups (e.g., number of MPI processes/pools).
 * @param igrp The zero-based index of the current group ($0 \le \mathtt{igrp} <
 * \mathtt{ngrp}$).
 * @param shift Pointer to an @ref ND_int where the starting shift (index) for
 * the current group will be stored.
 * @return @ref ND_int The size (number of elements) assigned to the current
 * group (@p igrp).
 */
ND_int distribute_to_grps(const ND_int n, const ND_int ngrp, const ND_int igrp,
                          ND_int *shift)
{
        /*
        This function is used to distribute k/q points between  n groups
        return k/q for each group, and the shift
        */

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
 * @brief Creates a set of parallel communicators for the ELPH calculation.
 *
 * This function initializes and splits the @c MPI_COMM_WORLD (or equivalent)
 * into several specialized communicators (@c commQ, @c commK, @c commR, @c
 * commRq) based on the specified number of q-point and k-point pools. It
 * performs a Cartesian-like decomposition of the total MPI processes.
 *
 * **Note on Splitting:** The splitting logic ensures a specific rank ordering
 * in the new communicators, which is crucial for the ELPH parallelization
 * strategy. The key given to @c MPI_Comm_split is deliberately chosen to be the
 * world rank in some cases to maintain this ordering.
 *
 * **Constraints:** It is required that $\mathtt{nqpools} \times
 * \mathtt{nkpools}$ must exactly divide the total number of processes in @p
 * MPI_world_comm.
 *
 * @param nqpools The number of q-point pools (groups of processes that share
 * k-points).
 * @param nkpools The number of k-point pools (groups of processes that share
 * q-points).
 * @param MPI_world_comm The global MPI communicator (usually @c MPI_COMM_WORLD)
 * to be split.
 * @param Comm Pointer to the @ref ELPH_MPI_Comms structure where the newly
 * created communicators will be stored.
 */
void create_parallel_comms(const int nqpools, const int int nkpools,
                           const MPI_Comm MPI_world_comm,
                           struct ELPH_MPI_Comms *Comm)
{
        /*
        !! Warning : Order of arrangement of cpus is very important. the key
    given     to comm_split must be the rank of world comm     */
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
                error_msg(
            " product of kpools and qpools must divide total cpus.");
           
    }

        Comm->commQ_size = Comm->commW_size / nqpools;
        Comm->commK_size = Comm->commQ_size / nkpools;

        /* create comms */
     // commQ
    mpi_error =
        MPI_Comm_split(Comm->commW, Comm->commW_rank / Comm->commQ_size,
                               Comm->commW_rank, &Comm->commQ);
        MPI_error_msg(mpi_error);

        mpi_error = MPI_Comm_rank(Comm->commQ, &Comm->commQ_rank);
        MPI_error_msg(mpi_error);

         // commK
    mpi_error =
        MPI_Comm_split(Comm->commQ, Comm->commQ_rank / Comm->commK_size,
                               Comm->commW_rank, &Comm->commK);
        MPI_error_msg(mpi_error);

        mpi_error = MPI_Comm_rank(Comm->commK, &Comm->commK_rank);
        MPI_error_msg(mpi_error);
         // commR
    mpi_error = MPI_Comm_split(Comm->commW, Comm->commK_rank, Comm->commW_rank,
                               & Comm->commR);
        MPI_error_msg(mpi_error);

        mpi_error = MPI_Comm_rank(Comm->commR, &Comm->commR_rank);
        MPI_error_msg(mpi_error);

        mpi_error = MPI_Comm_size(Comm->commR, &Comm->commR_size);
        MPI_error_msg(mpi_error);
         // commRq
    mpi_error = MPI_Comm_split(Comm->commQ, Comm->commK_rank, Comm->commW_rank,
                               & Comm->commRq);
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
 * @brief Frees the parallel communicators stored in the @ref ELPH_MPI_Comms
 * structure.
 *
 * This function calls @c MPI_Comm_free on the four specialized communicators:
 * @c commQ, @c commK, @c commR, and @c commRq. The world communicator
 * (@c commW) is not freed as it is typically @c MPI_COMM_WORLD.
 *
 * @param Comm Pointer to the @ref ELPH_MPI_Comms structure containing the
 * communicators to be freed.
 */
void free_parallel_comms(struct ELPH_MPI_Comms *Comm)
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

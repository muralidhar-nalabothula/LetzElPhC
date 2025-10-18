/**
 * @file
 * @brief Transpose routines for parallel FFT operations
 *
 * This file implements MPI-based transpose operations required for parallel
 * 3D FFT transforms. Data is redistributed across MPI ranks to enable
 * FFTs along different dimensions.
 */

#include <complex.h>
#include <fftw3.h>
#include <mpi.h>
#include <string.h>

#include "common/error.h"
#include "elphC.h"
#include "fft.h"

/**
 * @brief Forward transpose: (nGxy, Nz_loc) -> (nGxy_loc, Nz)
 *
 * Performs MPI transpose to redistribute data from Z-distributed layout
 * to (Gx,Gy)-distributed layout. This is needed before FFT along Z direction
 * in the forward transform.
 *
 * Data flow:
 * - Input: plan->nz_buf contains (nGxy, Nz_loc) - all (Gx,Gy) pairs, local Z
 * - Intermediate: plan->fft_data used as scratch space
 * - Output: plan->nz_buf contains (nGxy_loc, Nz) - local (Gx,Gy) pairs, all Z
 *
 * @param[in,out] plan FFT plan containing communication buffers and data
 *                     Input and output both stored in plan->nz_buf
 *
 * @note Uses MPI_Alltoallv for data redistribution
 * @note plan->fft_data is used as temporary scratch space
 */
void fwd_transpose(struct ELPH_fft_plan* plan)
{
    int mpi_error;
    int ncpus;
    mpi_error = MPI_Comm_size(plan->comm, &ncpus);
    MPI_error_msg(mpi_error);
    int* gxy_counts = plan->comm_bufs;
    int* xy_disp = gxy_counts + ncpus;
    int* z_counts = gxy_counts + 2 * ncpus;
    int* z_disp = gxy_counts + 3 * ncpus;
    mpi_error = MPI_Alltoallv(plan->nz_buf, gxy_counts, xy_disp, ELPH_MPI_cmplx,
                              plan->fft_data, z_counts, z_disp, ELPH_MPI_cmplx,
                              plan->comm);
    MPI_error_msg(mpi_error);
    /* Now plan->fft_data has the data.
    we need to redistribute and store to plan->nz_buf */
    for (ND_int ig = 0; ig < plan->nGxyloc; ++ig)
    {
        ELPH_cmplx* dest_ptr = plan->nz_buf + ig * plan->fft_dims[2];
        for (ND_int i = 0; i < ncpus; ++i)
        {
            ND_int zcount_tmp = z_counts[i] / plan->nGxyloc;
            ND_int z_disp_tmp = z_disp[i] / plan->nGxyloc;
            ELPH_cmplx* src_ptr = plan->fft_data + z_disp[i] + ig * zcount_tmp;
            memcpy(dest_ptr + z_disp_tmp, src_ptr,
                   sizeof(ELPH_cmplx) * zcount_tmp);
        }
    }
}

/**
 * @brief Backward transpose: (nGxy_loc, Nz) -> (nGxy, Nz_loc)
 *
 * Performs MPI transpose to redistribute data from (Gx,Gy)-distributed layout
 * to Z-distributed layout. This is needed after inverse FFT along Z direction
 * in the backward transform.
 *
 * Data flow:
 * - Input: plan->nz_buf contains (nGxy_loc, Nz) - local (Gx,Gy) pairs, all Z
 * - Intermediate: plan->fft_data used as scratch space
 * - Output: plan->nz_buf contains (nGxy, Nz_loc) - all (Gx,Gy) pairs, local Z
 *
 * @param[in,out] plan FFT plan containing communication buffers and data
 *                     Input and output both stored in plan->nz_buf
 *
 * @note Uses MPI_Alltoallv for data redistribution
 * @note plan->fft_data is used as temporary scratch space
 */
void bwd_transpose(struct ELPH_fft_plan* plan)
{
    int mpi_error;
    int ncpus;
    mpi_error = MPI_Comm_size(plan->comm, &ncpus);
    MPI_error_msg(mpi_error);
    int* gxy_counts = plan->comm_bufs;
    int* xy_disp = gxy_counts + ncpus;
    int* z_counts = gxy_counts + 2 * ncpus;
    int* z_disp = gxy_counts + 3 * ncpus;
    /*
    we need to redistribute and store to plan->fft_data */
    for (ND_int ig = 0; ig < plan->nGxyloc; ++ig)
    {
        ELPH_cmplx* src_ptr = plan->nz_buf + ig * plan->fft_dims[2];
        for (ND_int i = 0; i < ncpus; ++i)
        {
            ND_int zcount_tmp = z_counts[i] / plan->nGxyloc;
            ND_int z_disp_tmp = z_disp[i] / plan->nGxyloc;
            ELPH_cmplx* dest_ptr = plan->fft_data + z_disp[i] + ig * zcount_tmp;
            memcpy(dest_ptr, src_ptr + z_disp_tmp,
                   sizeof(ELPH_cmplx) * zcount_tmp);
        }
    }
    // Now plan->fft_data has the data. scatter back to plan->nz_buf
    mpi_error = MPI_Alltoallv(plan->fft_data, z_counts, z_disp, ELPH_MPI_cmplx,
                              plan->nz_buf, gxy_counts, xy_disp, ELPH_MPI_cmplx,
                              plan->comm);
    MPI_error_msg(mpi_error);
}

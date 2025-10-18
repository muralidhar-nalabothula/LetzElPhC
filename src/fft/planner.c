/**
 * @file
 * @brief FFT plan creation and destruction for wavefunction transforms
 *
 * This file implements the creation and destruction of FFTW plans for
 * parallel 3D FFT transforms of wavefunctions distributed across MPI ranks.
 */

#include <complex.h>
#include <fftw3.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "common/error.h"
#include "common/parallel.h"
#include "elphC.h"
#include "fft.h"

/**
 * @brief Creates FFT plans for parallel 3D wavefunction transforms
 *
 * Initializes all FFTW plans and communication buffers needed for parallel
 * 3D FFT transforms. The wavefunctions are distributed in G-space by
 * (Gx, Gy) pairs, and in real space by Z slices across MPI ranks.
 *
 * Creates multiple aligned plans for X and Y directions to handle different
 * memory alignments, and a single plan for Z direction transforms.
 *
 * @param[out] plan FFT plan structure to initialize
 * @param[in] ngvecs_loc Number of local G-vectors on this MPI rank
 * @param[in] nzloc Number of local Z grid points on this MPI rank
 * @param[in] nGxyloc Number of local (Gx,Gy) pairs on this MPI rank
 * @param[in] gvecs Array of G-vectors, shape: (ngvecs_loc,3)
 * @param[in] fft_dims FFT grid dimensions [Nx, Ny, Nz]
 * @param[in] fft_flags FFTW planner flags (e.g., FFTW_MEASURE, FFTW_ESTIMATE)
 * @param[in] comm MPI communicator for parallel FFT operations
 *
 * @note G-vectors must be sorted by (Gx, Gy) to group Z components together
 * @note Negative G-vector components are wrapped: G < 0 becomes G + N
 * @note Multiple aligned plans are created to handle SIMD alignment
 * requirements
 */
void wfc_plan(struct ELPH_fft_plan* plan, const ND_int ngvecs_loc,
              const ND_int nzloc, const ND_int nGxyloc, const int* gvecs,
              const ND_int* fft_dims, unsigned fft_flags, MPI_Comm comm)
{
    int ncpus, mpi_error;

    mpi_error = MPI_Comm_size(comm, &ncpus);
    MPI_error_msg(mpi_error);

    plan->comm_bufs = malloc(sizeof(int) * 4 * ncpus);
    CHECK_ALLOC(plan->comm_bufs);

    plan->nzloc = get_mpi_local_size_idx(fft_dims[2], NULL, comm);
    if (plan->nzloc != nzloc)
    {
        error_msg("Wrong z local dimensions to planner.");
    }

    ND_int size_fft_data = fft_dims[0] * fft_dims[1] * plan->nzloc;
    if (size_fft_data < (fft_dims[2] * nGxyloc))
    {
        size_fft_data = fft_dims[2] * nGxyloc;
    }

    plan->align_len = alignment_len();
    size_fft_data += plan->align_len;

    // alloc memory for fft_data and nz_buf;
    plan->fft_data = fftw_fun(malloc)(size_fft_data * sizeof(ELPH_cmplx));
    CHECK_ALLOC(plan->fft_data);

    plan->nz_buf = fftw_fun(malloc)(size_fft_data * sizeof(ELPH_cmplx));
    CHECK_ALLOC(plan->nz_buf);

    memset(plan->fft_data, 0, size_fft_data * sizeof(ELPH_cmplx));
    memset(plan->nz_buf, 0, size_fft_data * sizeof(ELPH_cmplx));

    plan->comm = comm;
    plan->gvecs = gvecs;
    plan->nGxyloc = nGxyloc;
    plan->ngvecs_loc = ngvecs_loc;

    memcpy(plan->fft_dims, fft_dims, sizeof(ND_int) * 3);

    mpi_error = MPI_Allreduce(&nGxyloc, &(plan->nGxy), 1, ELPH_MPI_ND_INT,
                              MPI_SUM, comm);
    MPI_error_msg(mpi_error);

    if (nGxyloc != get_mpi_local_size_idx(plan->nGxy, NULL, comm))
    {
        error_msg("Wrong xy local dimensions to planner.");
    }

    int* Gxy_loc = malloc(sizeof(int) * 2 * nGxyloc);
    CHECK_ALLOC(Gxy_loc);

    plan->Gxy_total = malloc(sizeof(int) * 2 * plan->nGxy);
    CHECK_ALLOC(plan->Gxy_total);

    plan->gx_inGvecs = malloc(sizeof(bool) * fft_dims[0]);
    CHECK_ALLOC(plan->gx_inGvecs);
    // initiate to false
    for (ND_int igx = 0; igx < fft_dims[0]; ++igx)
    {
        plan->gx_inGvecs[igx] = false;
    }

    int* gxy_counts = plan->comm_bufs;
    int* xy_disp = gxy_counts + ncpus;
    int* z_counts = gxy_counts + 2 * ncpus;
    int* z_disp = gxy_counts + 3 * ncpus;

    int xycount = 0;
    // initialize with gvecs which are out of the box
    int Gx_prev = fft_dims[0] + 10;
    int Gy_prev = fft_dims[1] + 10;

    plan->ngxy_z = malloc(sizeof(int) * nGxyloc);
    CHECK_ALLOC(plan->ngxy_z);
    // zero the buffer
    for (ND_int i = 0; i < nGxyloc; ++i)
    {
        plan->ngxy_z[i] = 0;
    }

    for (ND_int ig = 0; ig < plan->ngvecs_loc; ++ig)
    {
        if (gvecs[3 * ig + 1] != Gy_prev || gvecs[3 * ig] != Gx_prev)
        {
            Gx_prev = gvecs[3 * ig];
            Gy_prev = gvecs[3 * ig + 1];

            Gxy_loc[2 * xycount] = Gx_prev;
            Gxy_loc[2 * xycount + 1] = Gy_prev;

            if (Gx_prev < 0)
            {
                Gxy_loc[2 * xycount] += fft_dims[0];
            }
            if (Gy_prev < 0)
            {
                Gxy_loc[2 * xycount + 1] += fft_dims[1];
            }

            // sanity checks
            if (Gxy_loc[2 * xycount] >= fft_dims[0] ||
                Gxy_loc[2 * xycount + 1] >= fft_dims[1])
            {
                error_msg("gvec > Nfft dim");
            }
            if (Gxy_loc[2 * xycount] < 0 || Gxy_loc[2 * xycount + 1] < 0)
            {
                error_msg("Bad gvectors");
            }

            plan->gx_inGvecs[Gxy_loc[2 * xycount]] = true;

            ++xycount;
        }

        plan->ngxy_z[xycount - 1] += 1;
    }

    if (xycount != nGxyloc)
    {
        error_msg("Wrong number of xy components in each cpu");
    }

    // get global plan->gx_inGvecs
    mpi_error = MPI_Allreduce(MPI_IN_PLACE, plan->gx_inGvecs, fft_dims[0],
                              MPI_C_BOOL, MPI_LOR, comm);
    MPI_error_msg(mpi_error);

    mpi_error =
        MPI_Allgather(&xycount, 1, MPI_INT, gxy_counts, 1, MPI_INT, comm);
    MPI_error_msg(mpi_error);

    int nzloc_temp = nzloc;
    mpi_error =
        MPI_Allgather(&nzloc_temp, 1, MPI_INT, z_counts, 1, MPI_INT, comm);
    MPI_error_msg(mpi_error);

    int dispdisp_temp = 0;
    for (int i = 0; i < ncpus; ++i)
    {
        gxy_counts[i] *= 2;
        xy_disp[i] = dispdisp_temp;
        dispdisp_temp += gxy_counts[i];
    }

    mpi_error = MPI_Allgatherv(Gxy_loc, 2 * xycount, MPI_INT, plan->Gxy_total,
                               gxy_counts, xy_disp, MPI_INT, comm);
    MPI_error_msg(mpi_error);

    free(Gxy_loc);

    /* fill the comm buffer */
    dispdisp_temp = 0;
    for (int i = 0; i < ncpus; ++i)
    {
        gxy_counts[i] = nzloc * (gxy_counts[i] / 2);
        xy_disp[i] = nzloc * (xy_disp[i] / 2);

        z_counts[i] *= xycount;
        z_disp[i] = dispdisp_temp;
        dispdisp_temp += z_counts[i];
    }

    /* Now time to create plans */
    // i) create plan along entire X direction
    // ii) create two plans along y for -Gx_min < Gx < Gx_max
    // iii) create plan for along z only for set of (Gx,Gy) pairs

    // create plan buffers
    plan->fplan_x = malloc(sizeof(fftw_generic_plan) * 5 * plan->align_len);
    CHECK_ALLOC(plan->fplan_x);

    plan->fplan_y = plan->fplan_x + plan->align_len;

    /* backward plans G->r */
    plan->bplan_x = plan->fplan_x + 2 * plan->align_len;
    plan->bplan_y = plan->fplan_x + 3 * plan->align_len;

    /* convolution plan x */
    plan->cplan_x = plan->fplan_x + 4 * plan->align_len;

    // i) create forward plan and bwd plan along X
    for (ND_int i = 0; i < plan->align_len; ++i)
    {
        ELPH_cmplx* tmp_fft_plan_ptr = plan->fft_data + i;

        ND_int ia = fftw_fun(alignment_of)((void*)tmp_fft_plan_ptr);
        ia /= sizeof(ELPH_cmplx);

        int nff_dimx = fft_dims[0];
        // (Nx, k, Nz_loc)
        plan->fplan_x[ia] = fftw_fun(plan_many_dft)(
            1, &nff_dimx, fft_dims[1] * nzloc, tmp_fft_plan_ptr, NULL,
            fft_dims[1] * nzloc, 1, tmp_fft_plan_ptr, NULL, fft_dims[1] * nzloc,
            1, FFTW_FORWARD, fft_flags);
        if (plan->fplan_x[ia] == NULL)
        {
            error_msg("X forward plan failed");
        }

        // backward plan G->r
        plan->bplan_x[ia] = fftw_fun(plan_many_dft)(
            1, &nff_dimx, fft_dims[1] * nzloc, tmp_fft_plan_ptr, NULL,
            fft_dims[1] * nzloc, 1, tmp_fft_plan_ptr, NULL, fft_dims[1] * nzloc,
            1, FFTW_BACKWARD, fft_flags);
        if (plan->bplan_x[ia] == NULL)
        {
            error_msg("X backward plan failed");
        }
    }

    // ii) create Y ffts plans
    for (ND_int i = 0; i < plan->align_len; ++i)
    {
        ELPH_cmplx* tmp_fft_plan_ptr = plan->fft_data + i;

        ND_int ia = fftw_fun(alignment_of)((void*)tmp_fft_plan_ptr);
        ia /= sizeof(ELPH_cmplx);

        int nff_dimy = fft_dims[1];
        // (k, Nz_loc)
        plan->fplan_y[ia] = fftw_fun(plan_many_dft)(
            1, &nff_dimy, nzloc, tmp_fft_plan_ptr, NULL, nzloc, 1,
            tmp_fft_plan_ptr, NULL, nzloc, 1, FFTW_FORWARD, fft_flags);
        if (plan->fplan_y[ia] == NULL)
        {
            error_msg("Y forward plan failed");
        }

        // backward plan G->r
        plan->bplan_y[ia] = fftw_fun(plan_many_dft)(
            1, &nff_dimy, nzloc, tmp_fft_plan_ptr, NULL, nzloc, 1,
            tmp_fft_plan_ptr, NULL, nzloc, 1, FFTW_BACKWARD, fft_flags);
        if (plan->bplan_y[ia] == NULL)
        {
            error_msg("Y backward plan failed");
        }
    }

    // iii) create a single z plan.
    // forward plan
    plan->fplan_z = fftw_fun(plan_many_dft)(
        1, (int[1]){fft_dims[2]}, nGxyloc, plan->nz_buf, NULL, 1, fft_dims[2],
        plan->nz_buf, NULL, 1, fft_dims[2], FFTW_FORWARD, fft_flags);
    if (plan->fplan_z == NULL)
    {
        error_msg("Z forward plan failed");
    }

    // backward plan
    plan->bplan_z = fftw_fun(plan_many_dft)(
        1, (int[1]){fft_dims[2]}, nGxyloc, plan->nz_buf, NULL, 1, fft_dims[2],
        plan->nz_buf, NULL, 1, fft_dims[2], FFTW_BACKWARD, fft_flags);
    if (plan->bplan_z == NULL)
    {
        error_msg("Z backward plan failed");
    }

    /*
    Finally create convolution plan for FFT along x
    In convolutions we use these plans instead of fplan_x
    */
    for (ND_int i = 0; i < plan->align_len; ++i)
    {
        ELPH_cmplx* tmp_pln_ptr = plan->fft_data + i;
        ND_int ia = fftw_fun(alignment_of)((void*)tmp_pln_ptr);
        ia /= sizeof(ELPH_cmplx);
        // (Nx, k, Nz_loc)
        plan->cplan_x[ia] = fftw_fun(plan_many_dft)(
            1, (int[1]){fft_dims[0]}, nzloc, tmp_pln_ptr, NULL,
            fft_dims[1] * nzloc, 1, tmp_pln_ptr, NULL, fft_dims[1] * nzloc, 1,
            FFTW_FORWARD, fft_flags);

        if (plan->cplan_x[ia] == NULL)
        {
            error_msg("convolution forward X plan failed");
        }
    }
}

/**
 * @brief Destroys FFT plan and frees all allocated resources
 *
 * Frees all FFTW plans, work buffers, and communication buffers allocated
 * during plan creation.
 *
 * @param[in,out] plan FFT plan to destroy, can be NULL (no-op)
 */
void wfc_destroy_plan(struct ELPH_fft_plan* plan)
{
    if (NULL == plan)
    {
        return;
    }
    // destroy x,y,z buffer
    fftw_fun(free)(plan->fft_data);
    fftw_fun(free)(plan->nz_buf);

    // destroy plans

    for (ND_int i = 0; i < plan->align_len; ++i)
    {
        fftw_fun(destroy_plan)(plan->fplan_x[i]);
        fftw_fun(destroy_plan)(plan->bplan_x[i]);

        fftw_fun(destroy_plan)(plan->fplan_y[i]);
        fftw_fun(destroy_plan)(plan->bplan_y[i]);
        fftw_fun(destroy_plan)(plan->cplan_x[i]);
    }

    fftw_fun(destroy_plan)(plan->fplan_z);
    fftw_fun(destroy_plan)(plan->bplan_z);

    // free remaining buffers
    free(plan->gx_inGvecs);
    free(plan->fplan_x);
    free(plan->comm_bufs);
    free(plan->Gxy_total);
    free(plan->ngxy_z);
}

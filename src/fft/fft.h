/**
 * @file
 * @brief FFT operations and plan management using FFTW3
 *
 * Provides 3D FFT operations for wavefunctions and potentials with MPI
 * parallelization. Supports both single and double precision through
 * compile-time selection.
 */

#pragma once
#include <complex.h>
#include <fftw3.h>
#include <mpi.h>
#include <stdbool.h>

#include "elphC.h"

#if defined(COMPILE_ELPH_DOUBLE)
/**
 * @def fftw_fun(FUN_NAME)
 * @brief Macro to select double precision FFTW function
 *
 * Expands to fftw_FUN_NAME for double precision
 */
#define fftw_fun(FUN_NAME) fftw3_fun_HIDDEN(FUN_NAME)
#define fftw3_fun_HIDDEN(FUN_NAME) fftw_##FUN_NAME
#else
/**
 * @def fftw_fun(FUN_NAME)
 * @brief Macro to select single precision FFTW function
 *
 * Expands to fftwf_FUN_NAME for single precision
 */
#define fftw_fun(FUN_NAME) fftw3_fun_HIDDEN(FUN_NAME)
#define fftw3_fun_HIDDEN(FUN_NAME) fftwf_##FUN_NAME
#endif

/**
 * @typedef fftw_generic_plan
 * @brief Generic FFTW plan type (precision-dependent)
 */
typedef fftw_fun(plan) fftw_generic_plan;

/**
 * @struct ELPH_fft_plan
 * @brief FFT plan structure for 3D transforms with MPI parallelization
 *
 * Stores FFTW plans for forward/backward transforms along with buffers
 * and metadata for distributed FFT operations. Parallelization uses
 * pencil decomposition with z-direction distributed across MPI processes.
 */
struct ELPH_fft_plan
{
    ELPH_cmplx* fft_data; /**< Main FFT data buffer (created by planner, must be
                             freed) */
    ELPH_cmplx*
        nz_buf; /**< Buffer for MPI_Alltoallv operations (created by planner) */
    ND_int fft_dims[3]; /**< FFT grid dimensions (Nx, Ny, Nz) */
    ND_int nzloc;       /**< Number of local z elements on this MPI rank */
    ND_int align_len;   /**< SIMD alignment length */
    const int* gvecs;   /**< Local G-vectors  */
    ND_int ngvecs_loc;  /**< Number of local G-vectors */
    ND_int nGxyloc;     /**< Number of local G-vectors in xy plane */
    ND_int nGxy;        /**< Total number of G-vectors in xy plane */
    bool* gx_inGvecs; /**< Boolean array (Nx): true if gx is present in gvecs */
    int* Gxy_total;   /**< All xy G-vector components (nGxy,2) */
    int* ngxy_z;      /**< Number of z components for each (gx,gy) pair */
    int* comm_bufs;   /**< Communication buffers for MPI_Alltoallv (4,comm_size)
                         arrays: 0: gxy_counts - G-vectors per rank times nzloc 1:
                         xy_disp - Starting index of gxy for each rank times nzloc
                           2: z_counts - Nz elements per rank times nGxy_loc
                           3: z_disp - Starting z index for each rank times
                         nGxy_loc */

    /* Forward transform plans (real-space to reciprocal space) */
    fftw_generic_plan*
        fplan_x; /**< Array of FFT plans for x direction (align_len plans) */
    fftw_generic_plan*
        fplan_y; /**< Array of FFT plans for y direction (align_len plans) */
    fftw_generic_plan fplan_z; /**< FFT plan for z direction */

    /* Backward transform plans (reciprocal space to real-space) */
    fftw_generic_plan*
        bplan_x; /**< Array of invFFT plans for x direction (align_len plans) */
    fftw_generic_plan*
        bplan_y; /**< Array of invFFT plans for y direction (align_len plans) */
    fftw_generic_plan bplan_z; /**< invFFT plan for z direction */

    /* Convolution plans */
    fftw_generic_plan* cplan_x; /**< Array of FFT plans for x in convolutions
                                   (align_len plans) */

    MPI_Comm comm; /**< MPI communicator for parallel FFT operations */
};

/**
 * @brief Determines SIMD alignment length for FFTW
 *
 * @return Alignment length in units of sizeof(ELPH_cmplx)
 */
ND_int alignment_len(void);

/**
 * @brief Creates FFT plan for wavefunction transforms
 *
 * Initializes FFTW plans for 3D FFT operations with MPI parallelization.
 * Allocates buffers and creates forward/backward/convolution plans.
 *
 * @param plan Output FFT plan structure (allocated by caller)
 * @param ngvecs_loc Number of local G-vectors
 * @param nzloc Number of local z grid points
 * @param nGxyloc Number of local xy G-vectors
 * @param gvecs Array of local G-vectors (must remain valid during plan
 * lifetime)
 * @param fft_dims FFT grid dimensions (3 elements: Nx, Ny, Nz)
 * @param fft_flags FFTW planner flags (e.g., FFTW_MEASURE, FFTW_ESTIMATE)
 * @param comm MPI communicator
 */
void wfc_plan(struct ELPH_fft_plan* plan, const ND_int ngvecs_loc,
              const ND_int nzloc, const ND_int nGxyloc, const int* gvecs,
              const ND_int* fft_dims, unsigned fft_flags, MPI_Comm comm);

/**
 * @brief Destroys FFT plan and frees associated resources
 *
 * Frees all FFTW plans and buffers allocated during plan creation.
 *
 * @param plan FFT plan to destroy
 */
void wfc_destroy_plan(struct ELPH_fft_plan* plan);

/**
 * @brief Performs 3D forward FFT (real-space to reciprocal space)
 *
 * Transforms wavefunction from real-space to reciprocal space:
 * wfcG(G) = FFT[wfcr(r)]
 *
 * @param plan FFT plan
 * @param nsets Number of wavefunctions to transform
 * @param wfcr Input: real-space wavefunction (nsets,fft_grid)
 * @param wfcG Output: reciprocal-space wavefunction (nsets,ngvecs_loc)
 * @param conjugate If true, apply complex conjugation
 */
void fft3D(struct ELPH_fft_plan* plan, const ND_int nsets, ELPH_cmplx* wfcr,
           ELPH_cmplx* wfcG, const bool conjugate);

/**
 * @brief Performs 3D inverse FFT (reciprocal space to real-space)
 *
 * Transforms wavefunction from reciprocal space to real-space:
 * wfcr(r) = IFFT[wfcG(G)]
 *
 * @param plan FFT plan
 * @param nsets Number of wavefunctions to transform
 * @param wfcG Input: reciprocal-space wavefunction (nsets,ngvecs_loc)
 * @param wfcr Output: real-space wavefunction (nsets,fft_grid)
 * @param conjugate If true, apply complex conjugation
 */
void invfft3D(struct ELPH_fft_plan* plan, const ND_int nsets, ELPH_cmplx* wfcG,
              ELPH_cmplx* wfcr, const bool conjugate);

/**
 * @brief Performs forward transpose for distributed FFT
 *
 * Redistributes data for parallel FFT computation.
 *
 * @param plan FFT plan
 */
void fwd_transpose(struct ELPH_fft_plan* plan);

/**
 * @brief Performs backward transpose for distributed FFT
 *
 * Redistributes data back from parallel FFT computation.
 *
 * @param plan FFT plan
 */
void bwd_transpose(struct ELPH_fft_plan* plan);

/**
 * @brief Performs 3D convolution using FFT
 *
 * Computes convolution of potential and wavefunction in real space
 * via multiplication in Fourier space:
 * wfcG_out = FFT[Vpot(r) * psi(r)]
 *
 * @param plan FFT plan
 * @param nspinor Number of spinor components
 * @param nmag Number of magnetic components
 * @param Vpotr Potential in real space (fft_grid)
 * @param psir Wavefunction in real space (nspinor,fft_grid)
 * @param wfcG Output in reciprocal space (nspinor,ngvecs_loc)
 * @param conjugate If true, apply complex conjugation
 */
void fft_convolution3D(struct ELPH_fft_plan* plan, const ND_int nspinor,
                       ND_int nmag, const ELPH_cmplx* Vpotr,
                       const ELPH_cmplx* psir, ELPH_cmplx* wfcG,
                       const bool conjugate);

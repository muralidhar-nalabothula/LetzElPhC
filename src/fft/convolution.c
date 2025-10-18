/**
 * @file
 * @brief 3D FFT convolution of potential and wavefunction in real space
 *
 * This file implements 3D FFT convolution computing FFT(V(r)*psi(r)) where
 * V(r) is a potential and psi(r) is a wavefunction in real space.
 */

#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#include <string.h>

#include "common/error.h"
#include "elphC.h"
#include "fft.h"

/**
 * @brief Performs 3D FFT convolution of potential and wavefunction
 *
 * Computes the Fourier transform of the product V(r)*psi(r) where V(r) is
 * a potential and psi(r) is a wavefunction. The algorithm proceeds as:
 *
 * a) Compute V(r)*psi(r) and perform FFT along X and Y directions
 *    - For nmag = 1 or 2: element-wise multiplication "i,si=>si"
 *    - For nmag = 4: spinor contraction "rsxyz,sxyz->rxyz"
 *
 * b) Scatter (Gx,Gy) components to different MPI ranks for parallel FFT along Z
 *
 * c) Perform FFT along Z direction
 *
 * d) Extract sphere of G-vectors from FFT box
 *
 * @param[in,out] plan FFT plan containing buffers, dimensions, and FFTW plans
 * @param[in] nspinor Number of spinor components (1 or 2)
 * @param[in] nmag Number of magnetic/spinor components (1, 2, or 4)
 *                 - nmag=1: non-magnetic or collinear with nspinor=1
 *                 - nmag=2: LSDA (two spin channels)
 *                 - nmag=4: non-collinear with nspinor=2
 * @param[in] Vpotr Potential in real space
 *                  - Shape: (nmag, Nx, Ny, nzloc) for nmag=1,2
 *                  - Shape: (2, 2, Nx, Ny, nzloc) for nmag=4
 * @param[in] psir Wavefunction in real space, shape: (nspinor, Nx, Ny, nzloc)
 * @param[out] wfcG Output in G-space, shape: (nspinor, ngvecs_loc)
 * @param[in] conjugate If true, conjugate the output
 *
 * @note For nspinor=1, nmag is forced to 1
 * @note For nmag=2 (LSDA), two separate V(r) exist for each spin, treated as
 * nmag=1
 * @note For nmag=4, requires nspinor=2 (non-collinear magnetism)
 */
void fft_convolution3D(struct ELPH_fft_plan* plan, const ND_int nspinor,
                       ND_int nmag, const ELPH_cmplx* Vpotr,
                       const ELPH_cmplx* psir, ELPH_cmplx* wfcG,
                       const bool conjugate)
{
    // basic checks
    if (nspinor == 1)
    {
        nmag = 1;
    }
    if (nmag != 1 && nmag != 4)
    {
        error_msg("Error wrong nmag value.");
    }
    if (nmag == 4 && nspinor != 2)
    {
        error_msg("Error wrong nmag and nspinor combination.");
    }
    if (nspinor > 2)
    {
        error_msg("nspinor must be <= 2.");
    }

    // first some basic stuff
    ND_int Nx = plan->fft_dims[0];
    ND_int Ny = plan->fft_dims[1];
    ND_int Nz = plan->fft_dims[2];

    ELPH_float norm = Nx * Ny * Nz;
    norm = 1.0 / norm;

    ND_int Ny_stride = plan->nzloc * Ny;

    ND_int fft_buf_size = Nx * Ny * plan->nzloc;

    for (ND_int ispinor = 0; ispinor < nspinor; ++ispinor)
    {
        ELPH_cmplx* wfcG_tmp = wfcG + ispinor * plan->ngvecs_loc;
        const ELPH_cmplx* dV_r = Vpotr;

        // we store V(r)*psi(r) in plan->fft_data and perform FFT
        /* now we do blocking to reuse cache. we divide into Ny blocks and
            perform convolution and FFT */
        if (nmag == 1)
        {
            const ELPH_cmplx* psi_spinor = psir + ispinor * fft_buf_size;
            for (ND_int iy = 0; iy < Ny; ++iy)
            {
                ND_int iyshift = iy * plan->nzloc;
                for (ND_int ix = 0; ix < Nx; ++ix)
                {
                    ND_int ixshift = ix * Ny * plan->nzloc + iyshift;
                    ELPH_cmplx* dvpsi_xyz = plan->fft_data + ixshift;
                    const ELPH_cmplx* psi_xyz = psi_spinor + ixshift;
                    const ELPH_cmplx* dV_xyz = dV_r + ixshift;
                    for (ND_int iz = 0; iz < plan->nzloc; ++iz)
                    {
                        dvpsi_xyz[iz] = psi_xyz[iz] * dV_xyz[iz];
                    }
                }
                // perform the fft along X
                ELPH_cmplx* dvpsi_ptr = plan->fft_data + iyshift;
                ND_int iax = fftw_fun(alignment_of)((void*)(dvpsi_ptr));
                iax /= sizeof(ELPH_cmplx);
                fftw_fun(execute_dft)(plan->cplan_x[iax], dvpsi_ptr, dvpsi_ptr);
            }
        }
        else
        {
            if (ispinor == 1)
            {
                dV_r = Vpotr + 2 * fft_buf_size;
            }

            for (ND_int iy = 0; iy < Ny; ++iy)
            {
                ND_int iyshift = iy * plan->nzloc;
                for (ND_int ix = 0; ix < Nx; ++ix)
                {
                    ND_int ixshift = ix * Ny * plan->nzloc + iyshift;
                    ELPH_cmplx* dvpsi_out = plan->fft_data + ixshift;
                    const ELPH_cmplx* psi0 = psir + ixshift;
                    const ELPH_cmplx* psi1 = psir + ixshift + fft_buf_size;

                    const ELPH_cmplx* dV0 = dV_r + ixshift;
                    const ELPH_cmplx* dV1 = dV_r + ixshift + fft_buf_size;
                    for (ND_int iz = 0; iz < plan->nzloc; ++iz)
                    {
                        dvpsi_out[iz] = psi0[iz] * dV0[iz] + psi1[iz] * dV1[iz];
                    }
                }
                // perform the fft along X
                ELPH_cmplx* dvpsi_ptr = plan->fft_data + iyshift;
                ND_int iax = fftw_fun(alignment_of)((void*)(dvpsi_ptr));
                iax /= sizeof(ELPH_cmplx);
                fftw_fun(execute_dft)(plan->cplan_x[iax], dvpsi_ptr, dvpsi_ptr);
            }
        }

        // a) (ii) FFT along Y
        for (ND_int ix = 0; ix < plan->fft_dims[0]; ++ix)
        {
            if (!(plan->gx_inGvecs[ix]))
            {
                continue;
            }
            ELPH_cmplx* wfcr_tmp_y = plan->fft_data + ix * Ny_stride;
            ND_int iax = fftw_fun(alignment_of)((void*)wfcr_tmp_y);
            iax /= sizeof(ELPH_cmplx);
            fftw_fun(execute_dft)(plan->fplan_y[iax], wfcr_tmp_y, wfcr_tmp_y);
        }

        // b) (i) pack the data for transpose
        for (ND_int ixy = 0; ixy < plan->nGxy; ++ixy)
        {
            ND_int Gx = plan->Gxy_total[2 * ixy];
            ND_int Gy = plan->Gxy_total[2 * ixy + 1];
            // [Gx,Gy,0] =  Gx*Ny*Nzloc + Gy*Nzloc
            ELPH_cmplx* xy_buf = plan->fft_data + plan->nzloc * (Gy + Gx * Ny);
            // (Nxy,Zloc)
            memcpy(plan->nz_buf + ixy * plan->nzloc, xy_buf,
                   sizeof(ELPH_cmplx) * plan->nzloc);
        }
        // b) (ii)  transpose the data
        fwd_transpose(plan);

        // c) perform fft along z
        fftw_fun(execute)(plan->fplan_z);

        // d) box to sphere
        ND_int igvec = 0;
        for (ND_int ixy = 0; ixy < plan->nGxyloc; ++ixy)
        {
            ELPH_cmplx* zfft_ptr = plan->nz_buf + ixy * Nz;
            for (ND_int ig = 0; ig < plan->ngxy_z[ixy]; ++ig)
            {
                ND_int Gz = plan->gvecs[3 * igvec + 2];
                if (Gz < 0)
                {
                    Gz += Nz;
                }
                wfcG_tmp[igvec] = zfft_ptr[Gz] * norm;
                if (conjugate)
                {
                    wfcG_tmp[igvec] = conj(wfcG_tmp[igvec]);
                }
                ++igvec;
            }
        }
        if (igvec != plan->ngvecs_loc)
        {
            error_msg("Gvec mismatch in fwd execute. ");
        }
    }
}

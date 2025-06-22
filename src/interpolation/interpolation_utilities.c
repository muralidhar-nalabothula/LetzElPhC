// This file contains some helper functions used in interpolation

#include "interpolation_utilities.h"

#include <math.h>
#include <stdlib.h>

#include "common/constants.h"
#include "common/error.h"
#include "elphC.h"
#include "fft/fft.h"

static int qpt_sort_cmp(const void* a, const void* b);

void fft_q2R(ELPH_cmplx* data, ND_int* qgrid, ND_int nsets)
{
    // given a data (qx,qy,qz,nsets) -> (nsets, qx,qy,qz).
    // The reasone we want to do transpose is because
    // we we will later do matmul where we may encoutor
    // large strides > INT_MAX
    // performs inplace fourier transform q->R
    // Total number of complex elements
    //
    ND_int qx = qgrid[0];
    ND_int qy = qgrid[1];
    ND_int qz = qgrid[2];

    fftw_fun(iodim64) dims[3];
    dims[0].n = qx;
    dims[0].is = qy * qz * nsets;
    dims[0].os = qy * qz;

    dims[1].n = qy;
    dims[1].is = qz * nsets;
    dims[1].os = qz;

    dims[2].n = qz;
    dims[2].is = nsets;
    dims[2].os = 1;

    // How the transforms are repeated over nsets
    fftw_fun(iodim64) howmany_dims[1];
    howmany_dims[0].n = nsets;
    howmany_dims[0].is = 1;
    howmany_dims[0].os = qx * qy * qz;

    // DOnot overwrite the data, so always set to
    // FFTW_ESTIMATE.
    fftw_generic_plan plan = fftw_fun(plan_guru64_dft)(
        3, dims,          // 3D transform
        1, howmany_dims,  // number of transforms = nsets
        data, data, FFTW_FORWARD, FFTW_ESTIMATE);

    if (NULL == plan)
    {
        error_msg("FFT q->R plan failed.");
    }

    fftw_fun(execute)(plan);
    fftw_fun(destroy_plan)(plan);
}

void Sorted_qpts_idxs(const ND_int nqpts, ELPH_float* qpts, ND_int* indices)
{
    /*
    given list of qpts, sort then in fft_grid order.
    Return : sorted indices

    // qpts in crystal coordinates
    */
    if (NULL == qpts)
    {
        return;
    }

    ELPH_float** qvec_ptrs = malloc(sizeof(ELPH_float*) * nqpts);
    CHECK_ALLOC(qvec_ptrs);

    /* fill the struct */
    for (ND_int i = 0; i < nqpts; ++i)
    {
        qvec_ptrs[i] = qpts + 3 * i;
    }

    qsort(qvec_ptrs, nqpts, sizeof(ELPH_float*), qpt_sort_cmp);

    // store the sorted indices
    for (ND_int i = 0; i < nqpts; ++i)
    {
        indices[i] = (qvec_ptrs[i] - qpts) / 3;
    }

    free(qvec_ptrs);
}

void rearrange_qpt_grid(const ND_int nqpts, const ELPH_cmplx* in_buf,
                        const ND_int* idx, ELPH_cmplx* out_buf)
{
    // out[i] = in[idx[i]]
    for (ND_int i = 0; i < nqpts; ++i)
    {
        out_buf[i] = in_buf[idx[i]];
    }
}

void find_qpt_grid(const ND_int nqpts, const ELPH_float* qpts, ND_int* q_grid)
{
    // given list of qpoints, it finds the phonon grid.
    // q_grid is 3 ints

    for (ND_int i = 0; i < 3; ++i)
    {
        q_grid[i] = 1;
    }

    for (ND_int i = 0; i < nqpts; ++i)
    {
        for (ND_int ix = 0; ix < 3; ++ix)
        {
            ELPH_float qx = qpts[3 * i + ix] - floor(qpts[3 * i + ix]);

            if (fabs(qx) > ELPH_EPS && fabs(qx - 1) > ELPH_EPS)
            {
                qx = 1.0 / qx;
                ND_int qx_int = rint(qx);
                qx -= qx_int;

                if (fabs(qx) < ELPH_EPS)
                {
                    q_grid[ix] = (qx_int > q_grid[ix]) ? qx_int : q_grid[ix];
                }
            }
        }
    }

    // finally do a sanity check i.e the product of three dimensions are equal
    // to nqpts
    ND_int nqpt_prod = q_grid[0] * q_grid[1] * q_grid[2];

    if (nqpts != nqpt_prod)
    {
        error_msg("Number of qpts != product of found qpt grid.");
    }
}

/*  Static functions */
static int qpt_sort_cmp(const void* a, const void* b)
{
    /*
    First sorts along x, then y and finally z
    This will ensure that the fastest will be z, followed by y and then x
    */
    ELPH_float* v1 = *(ELPH_float**)a;
    ELPH_float* v2 = *(ELPH_float**)b;

    ELPH_float kdiff[3];

    // first bring them to [0,1)
    for (int i = 0; i < 3; ++i)
    {
        ELPH_float k1 = v1[i];
        k1 -= floor(k1);
        k1 += ELPH_EPS;
        k1 -= floor(k1);

        ELPH_float k2 = v2[i];
        k2 -= floor(k2);
        k2 += ELPH_EPS;
        k2 -= floor(k2);

        kdiff[i] = k1 - k2;
    }

    if (fabs(kdiff[0]) < ELPH_EPS)
    {
        if (fabs(kdiff[1]) < ELPH_EPS)
        {
            if (fabs(kdiff[2]) < ELPH_EPS)
            {
                return 0;
            }
            else
            {
                return kdiff[2] > 0.0 ? 1 : -1;
            }
        }
        else
        {
            return kdiff[1] > 0.0 ? 1 : -1;
        }
    }
    else
    {
        return kdiff[0] > 0.0 ? 1 : -1;
    }
}

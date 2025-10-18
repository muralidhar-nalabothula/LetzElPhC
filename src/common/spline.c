/**
 * @file
 * @brief Cubic spline interpolation functions
 *
 * Implements cubic spline interpolation with not-a-knot boundary conditions.
 * Uses tridiagonal matrix solver (Thomas algorithm) for O(N) computation.
 *
 * Reference: https://en.wikipedia.org/wiki/Spline_interpolation
 */

#include <stdlib.h>

#include "elphC.h"
#include "error.h"
#include "numerical_func.h"

/**
 * @brief Performs cubic spline interpolation at given point
 *
 * Evaluates the cubic spline at position x using precomputed second
 * derivatives. Uses Hermite interpolation formula:
 *
 * S(x) = (1-t)*y_i + t*y_{i+1} + t(1-t)*[(1-t)*a + t*b]
 *
 * where:
 * - t = (x - x_i) / (x_{i+1} - x_i)
 * - a = y''_i * h - (y_{i+1} - y_i)
 * - b = -y''_{i+1} * h + (y_{i+1} - y_i)
 * - h = x_{i+1} - x_i
 *
 * @param x Point at which to interpolate
 * @param inear Index of nearest grid point (should satisfy x_i <= x <= x_{i+1})
 * @param xi Array of x coordinates (grid points)
 * @param yi Array of y values at grid points
 * @param dy Array of second derivatives (from prepare_spline)
 * @return Interpolated value at x
 */
ELPH_float spline_interpolate(const ELPH_float x, ND_int inear,
                              const ELPH_float* xi, const ELPH_float* yi,
                              const ELPH_float* dy)
{
    ELPH_float tx = (x - xi[inear]) / (xi[inear + 1] - xi[inear]);
    ELPH_float ax =
        dy[inear] * (xi[inear + 1] - xi[inear]) - (yi[inear + 1] - yi[inear]);
    ELPH_float bx = -dy[inear + 1] * (xi[inear + 1] - xi[inear]) +
                    (yi[inear + 1] - yi[inear]);
    return (1 - tx) * yi[inear] + tx * yi[inear + 1] +
           tx * (1 - tx) * ((1 - tx) * ax + tx * bx);
}

/**
 * @brief Prepares cubic spline by computing second derivatives
 *
 * Solves tridiagonal system to find second derivatives at each point
 * using not-a-knot boundary conditions. These boundary conditions enforce
 * continuity of the third derivative at the second and second-to-last points.
 *
 * The natural cubic spline satisfies:
 * S''_i * h_{i-1} + 2*S''_{i-1}*(h_i + h_{i-1}) + S''_{i+1}*h_i =
 *     3*(s_{i-1}*h_i + s_i*h_{i-1})
 *
 * where h_i = x_{i+1} - x_i and s_i = (y_{i+1} - y_i) / h_i
 *
 * Not-a-knot boundary conditions:
 * - At left: S'''(x_1-) = S'''(x_1+)
 * - At right: S'''(x_{n-2}-) = S'''(x_{n-2}+)
 *
 * @param nvals Number of data points
 * @param xin Array of x coordinates (nvals, should be sorted)
 * @param yin Array of y values at grid points (nvals)
 * @param dy Output array for second derivatives (nvals)
 *
 * @note Uses Thomas algorithm (forward sweep + back substitution) for O(N)
 * solve
 */
void prepare_spline(const ND_int nvals, ELPH_float* xin, ELPH_float* yin,
                    ELPH_float* dy)
{
    /* first allocate a buffer required for scratch space*/
    ELPH_float* buf = malloc(sizeof(ELPH_float) * 4 * nvals);
    CHECK_ALLOC(buf);

    // // sort xin
    // for (ND_int i = 0 ; i<nvals; ++i)
    // {
    //     buf[2*i]   = xin[i];
    //     buf[2*i+1] = yin[i];
    // }
    // // sort
    // qsort(buf, nvals, 2*sizeof(ELPH_float), xin_xmp);
    // for (ND_int i = 0 ; i<nvals; ++i)
    // {
    //     xin[i] = buf[2*i];
    //     yin[i] = buf[2*i+1];
    // }

    /* Now Compute the derivates */
    ELPH_float* ai = buf;
    ELPH_float* bi = buf + nvals;
    ELPH_float* ci = buf + 2 * nvals;
    ELPH_float* di = buf + 3 * nvals;

    // not a knot BCs
    /* create a small scope */
    {
        /* start point BC */
        ELPH_float dx1 = xin[1] - xin[0];
        ELPH_float s1 = (yin[1] - yin[0]) / dx1;
        ELPH_float dx2 = xin[2] - xin[1];
        ELPH_float s2 = (yin[2] - yin[1]) / dx2;

        ai[0] = dx2 * dx2;
        bi[0] = dx2 * dx2 - dx1 * dx1;
        ci[0] = -dx1 * dx1;
        di[0] = 2 * (s1 * dx2 * dx2 - s2 * dx1 * dx1);

        /* end point BC */
        dx1 = xin[nvals - 2] - xin[nvals - 3];
        s1 = (yin[nvals - 2] - yin[nvals - 3]) / dx1;
        dx2 = xin[nvals - 1] - xin[nvals - 2];
        s2 = (yin[nvals - 1] - yin[nvals - 2]) / dx2;

        ai[nvals - 1] = dx2 * dx2;
        bi[nvals - 1] = dx2 * dx2 - dx1 * dx1;
        ci[nvals - 1] = -dx1 * dx1;
        di[nvals - 1] = 2 * (s1 * dx2 * dx2 - s2 * dx1 * dx1);
    }

    for (ND_int i = 1; i < (nvals - 1); ++i)
    {
        ELPH_float dxi = xin[i + 1] - xin[i];
        ELPH_float dxi_n = xin[i] - xin[i - 1];
        ELPH_float si = (yin[i + 1] - yin[i]) / dxi;
        ELPH_float si_n = (yin[i] - yin[i - 1]) / dxi_n;

        ai[i] = dxi;
        bi[i] = 2.0 * (dxi + dxi_n);
        ci[i] = dxi_n;
        di[i] = 3 * (si_n * dxi + si * dxi_n);
    }
    // note we must remove c[0] and a[n-1] terms
    // Gauss Elimination on top and bottom row
    ELPH_float alpha;
    // alpha = -c0/c1
    alpha = -ci[0] / ci[1];
    ai[0] += alpha * ai[1];
    bi[0] += alpha * bi[1];
    ci[0] += alpha * ci[1];
    di[0] += alpha * di[1];

    // now rearrange
    ci[0] = bi[0];
    bi[0] = ai[0];
    ai[0] = 0;

    // alpha = -an-1/an-2
    alpha = -ai[nvals - 1] / ai[nvals - 2];
    ai[nvals - 1] += alpha * ai[nvals - 2];
    bi[nvals - 1] += alpha * bi[nvals - 2];
    ci[nvals - 1] += alpha * ci[nvals - 2];
    di[nvals - 1] += alpha * di[nvals - 2];

    // now rearrange the last row
    ai[nvals - 1] = bi[nvals - 1];
    bi[nvals - 1] = ci[nvals - 1];
    ci[nvals - 1] = 0;
    /*solve Thomas algorithm*/
    // a) Inplace forward sweep
    for (ND_int i = 1; i < nvals; ++i)
    {
        ELPH_float ww = ai[i] / bi[i - 1];
        bi[i] -= ww * ci[i - 1];
        di[i] -= ww * di[i - 1];
    }

    // b) back substitution
    dy[nvals - 1] = di[nvals - 1] / bi[nvals - 1];
    for (ND_int i = (nvals - 2); i >= 0; --i)
    {
        dy[i] = (di[i] - ci[i] * dy[i + 1]) / bi[i];
    }

    free(buf);
}

// /**** static functions *****/
// int xin_xmp(const void* a, const void* b)
// {
//     const ELPH_float arg1 = *(const ELPH_float*)a;
//     const ELPH_float arg2 = *(const ELPH_float*)b;

// return (arg1 > arg2) - (arg1 < arg2);
// }

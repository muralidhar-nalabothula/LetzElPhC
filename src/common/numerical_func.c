/**
 * @file
 * @brief High-frequency numerical functions with optimization focus
 *
 * Contains commonly-used numerical routines including special functions,
 * linear algebra operations, FFT utilities, and vector operations.
 * These functions are performance-critical and rely on compiler optimization
 * (compile with -O3 -march=native for auto-vectorization).
 */

#include "numerical_func.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "constants.h"
#include "elphC.h"
#include "error.h"
#include "omp_pragma_def.h"

/* Static function declarations */
static ELPH_float factorial(ND_int n);

/* Function bodies */

/**
 * @brief Computes associated Legendre polynomial P_l^m(x)
 *
 * Calculates the associated Legendre polynomial without the Condon-Shortley
 * phase. Uses the recurrence relations from Eqs. 9-12 in: M.A. Blanco et al.,
 * J. Mol. Struct. THEOCHEM 419, 19-27 (1997)
 *
 * Recurrence relations:
 * - P_l^l(x) = (2l-1) * sqrt(1-x^2) * P_{l-1}^{l-1}(x)
 * - P_l^{l+1}(x) = (2l+1) * x * P_l^l(x)
 * - P_l^m(x) = [(2l-1)*x*P_{l-1}^m(x) - (l+m-1)*P_{l-2}^m(x)] / (l-m)
 *
 * Also uses the relations:
 * - P_{-l-1}^m(x) = P_l^m(x)
 * - P_l^{-m}(x) = (-1)^m * [(l-m)!/(l+m)!] * P_l^m(x)
 *
 * @param l_val Degree l (l >= 0)
 * @param m_val Order m (|m| <= l)
 * @param x_in Argument x in range [-1, 1]
 * @return Associated Legendre polynomial value P_l^m(x)
 *
 * @note Returns 0 if |m| > l
 */
ELPH_float legendre(int l_val, int m_val, ELPH_float x_in)
{
    /* Order of below if-else conditions is important */
    if (l_val == 0 && m_val == 0)
    {
        return 1.0;
    }

    else if (l_val < 0)
    {
        return legendre(-l_val - 1, m_val, x_in);
    }

    else if (abs(m_val) > l_val)
    {
        return 0.0;
    }

    else if (m_val < 0)
    {
        return legendre(l_val, -m_val, x_in) * pow(-1, -m_val) *
               factorial(l_val + m_val) / factorial(l_val - m_val);
    }
    else
    {
        ELPH_float legendre_val = 1.0;
        ELPH_float sqrt_xin = sqrt(1 - (x_in * x_in));

        /* P^l_l = (2l-1)*sqrt(1-x^2)*P^{l-1}_{l-1}(x)  */
        for (int im = 1; im <= m_val; ++im)
        {
            legendre_val *= (2 * im - 1) * sqrt_xin;
        }

        if (l_val == m_val)
        {
            return legendre_val;
        }

        ELPH_float legendre_val_mm = legendre_val;
        /* P^l_l+1 = (2l+1)*x*P^l_l*/
        legendre_val *= (2 * m_val + 1) * x_in;

        if (l_val == (m_val + 1))
        {
            return legendre_val;
        }

        /* For all other l,m combination */
        for (int im = m_val + 2; im <= l_val; ++im)
        {
            ELPH_float legendre_temp = legendre_val;
            legendre_val = (2 * im - 1) * x_in * legendre_val -
                           (im + m_val - 1) * legendre_val_mm;
            legendre_val = legendre_val / (im - m_val);
            legendre_val_mm = legendre_temp;
        }
        return legendre_val;
    }
}

/**
 * @brief Computes real spherical harmonic Y_l^m for a given direction vector
 *
 * Calculates real spherical harmonics using the convention from:
 * M.A. Blanco et al., J. Mol. Struct. THEOCHEM 419, 19-27 (1997)
 *
 * The real spherical harmonic is defined as:
 * - For m > 0: Y_l^m = sqrt(2) * N_l^m * P_l^m(cos theta) * cos(m*phi)
 * - For m < 0: Y_l^m = sqrt(2) * N_l^|m| * P_l^|m|(cos theta) * sin(|m|*phi)
 * - For m = 0: Y_l^0 = N_l^0 * P_l^0(cos theta)
 *
 * where N_l^m = sqrt[(2l+1)*(l-m)! / (4*pi*(l+m)!)] is the normalization.
 *
 * @param l_val Degree l (l >= 0)
 * @param m_val Order m (-l <= m <= l)
 * @param vec Direction vector (3 components, Cartesian)
 * @return Real spherical harmonic value Y_l^m(theta, phi)
 *
 * @note Returns 0 for l < 0
 * @note For zero vector, returns 1/sqrt(4*pi) if l=0, otherwise 0
 */
ELPH_float Ylm(int l_val, int m_val, ELPH_float* vec)
{
    if (l_val < 0)
    {
        return 0.0;  // error
    }

    ELPH_float cost, phi;
    ELPH_float norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);

    if (norm < ELPH_EPS)
    {
        if (l_val != 0)
        {
            return 0;
        }
        else
        {
            return 1.0 / sqrt(4 * ELPH_PI);
        }
    }

    cost = vec[2] / norm;  // z/r
    phi = atan2(vec[1], vec[0]);

    ELPH_float sh = 0.0;

    /* phi part*/
    if (m_val > 0)
    {
        sh = sqrt(2) * cos(m_val * phi);
    }
    else if (m_val < 0)
    {
        m_val = -m_val;
        sh = sqrt(2) * sin(m_val * phi);
    }
    else
    {
        sh = 1;  // 1 for m =0
    }

    /* Multiply with legendre polynomial */
    sh *= legendre(l_val, m_val, cost);

    /* Multiply with normalization factor*/
    sh *= sqrt(factorial(l_val - m_val) * (2 * l_val + 1) /
               (4 * ELPH_PI * factorial(l_val + m_val)));

    return sh;
}

/**
 * @brief Performs numerical integration using Simpson's 1/3 rule
 *
 * Computes the integral: integral f(x) dx using composite Simpson's rule.
 *
 * For odd number of points (n = 2k+1):
 * I = (dx/3) * [f_0 + 4*(f_1+f_3+...+f_{n-2}) + 2*(f_2+f_4+...+f_{n-3}) + f_n]
 *
 * For even number of points, uses Simpson's 3/8 rule for last segment.
 *
 * @param func_vals Function values at grid points (npts)
 * @param dx Grid spacing at each point (npts)
 * @param npts Number of grid points
 * @return Approximate integral value
 *
 * @note Returns 0 if npts < 3 (insufficient points)
 */
ELPH_float simpson(const ELPH_float* func_vals, const ELPH_float* dx,
                   ND_int npts)
{
    if (npts < 3)
    {
        return 0;  // Need atleast 3 points . Return an error instead
    }

    ELPH_float sum = func_vals[0] * dx[0];

    /* Odd terms*/
    for (ND_int i = 1; i < npts - 1; i += 2)
    {
        sum += 4.0 * func_vals[i] * dx[i];
    }

    /* even terms*/
    for (ND_int i = 2; i < npts - 1; i += 2)
    {
        sum += 2.0 * func_vals[i] * dx[i];
    }

    if (npts % 2 != 0)
    {
        sum += func_vals[npts - 1] * dx[npts - 1];
    }
    else
    {
        sum -= func_vals[npts - 3] * dx[npts - 3] * 0.25;
        sum += func_vals[npts - 2] * dx[npts - 2];
        sum += func_vals[npts - 1] * dx[npts - 1] * 1.25;
    }

    return sum / 3.0;
}

/**
 * @brief Computes cosine of angle between two 3D vectors
 *
 * Calculates: cos(omega) = (vec1 . vec2) / (|vec1| * |vec2|)
 *
 * @param vec1 First vector (3 components)
 * @param vec2 Second vector (3 components)
 * @return Cosine of angle between vectors
 *
 * @note Returns 0 if either vector has norm < ELPH_EPS
 */
ELPH_float cos_angle_bw_Vec(const ELPH_float* vec1, const ELPH_float* vec2)
{
    ELPH_float norm1 =
        vec1[0] * vec1[0] + vec1[1] * vec1[1] + vec1[2] * vec1[2];
    ELPH_float norm2 =
        vec2[0] * vec2[0] + vec2[1] * vec2[1] + vec2[2] * vec2[2];
    ELPH_float norm = sqrt(norm1 * norm2);
    if (norm < ELPH_EPS)
    {
        return 0;
    }
    ELPH_float dot_12 =
        vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
    return dot_12 / norm;
}

/** Static functions **/

/**
 * @brief Computes factorial n! using gamma function
 *
 * Calculates: n! = Gamma(n+1)
 *
 * @param n Non-negative integer
 * @return Factorial value n!
 *
 * @note Returns 0 for n < 0 (error case)
 */
static ELPH_float factorial(ND_int n)
{
    if (n < 0)
    {
        return 0.0;  // error
    }
    return tgamma(n + 1);
}

/**
 * @brief Matrix-vector multiplication for 3x3 matrix and 3D vector
 *
 * Computes:
 * - If trans = false: out = Mat * vec
 * - If trans = true:  out = Mat^T * vec
 *
 * Matrix is stored in row-major order (C-style).
 *
 * @param Mat 3x3 matrix (9 elements, row-major)
 * @param vec Input vector (3 elements)
 * @param trans If true, use transpose of Mat
 * @param out Output vector (3 elements)
 */
void MatVec3f(const ELPH_float* Mat, const ELPH_float* vec, const bool trans,
              ELPH_float* restrict out)
{
    if (trans)
    {
        out[0] = Mat[0] * vec[0] + Mat[3] * vec[1] + Mat[6] * vec[2];
        out[1] = Mat[1] * vec[0] + Mat[4] * vec[1] + Mat[7] * vec[2];
        out[2] = Mat[2] * vec[0] + Mat[5] * vec[1] + Mat[8] * vec[2];
    }
    else
    {
        out[0] = Mat[0] * vec[0] + Mat[1] * vec[1] + Mat[2] * vec[2];
        out[1] = Mat[3] * vec[0] + Mat[4] * vec[1] + Mat[5] * vec[2];
        out[2] = Mat[6] * vec[0] + Mat[7] * vec[1] + Mat[8] * vec[2];
    }
}

/**
 * @brief Computes complex dot product of two vectors
 *
 * Calculates: sum = conj(vec1) . vec2 = sum_i conj(vec1[i]) * vec2[i]
 *
 * @param vec1 First complex vector (n elements)
 * @param vec2 Second complex vector (n elements)
 * @param n Number of elements
 * @return Complex dot product
 */
ELPH_cmplx Cmplxdot(const ELPH_cmplx* vec1, const ELPH_cmplx* vec2,
                    const ND_int n)
{
    ELPH_cmplx sum = 0.0;

    for (ND_int i = 0; i < n; ++i)
    {
        sum += conj(vec1[i]) * vec2[i];
    }

    return sum;
}

/**
 * @brief Normalizes a complex vector in-place
 *
 * Normalizes vec such that sqrt(vec^dagger * vec) = 1.
 *
 * @param vec Complex vector to normalize (n elements, modified in-place)
 * @param n Number of elements
 *
 * @note Does nothing if norm < ELPH_EPS (zero vector)
 */
void normalize_Cmplx_vec(ELPH_cmplx* vec, const ND_int n)
{
    ELPH_float norm = sqrt(cabs(Cmplxdot(vec, vec, n)));

    if (norm < ELPH_EPS)
    {
        return;
    }

    for (ND_int i = 0; i < n; ++i)
    {
        vec[i] = vec[i] / norm;
    }
}

/**
 * @brief Computes determinant of 3x3 matrix
 *
 * Calculates: det(A) using cofactor expansion along first row:
 * det = a_00*(a_11*a_22 - a_12*a_21) - a_01*(a_10*a_22 - a_12*a_20) +
 * a_02*(a_10*a_21 - a_11*a_20)
 *
 * @param mat 3x3 matrix (9 elements, row-major)
 * @return Determinant value
 */
ELPH_float det3x3(const ELPH_float* mat)
{
    ELPH_float det = 0;
    det += mat[0] * (mat[4] * mat[8] - mat[5] * mat[7]);
    det -= mat[1] * (mat[3] * mat[8] - mat[5] * mat[6]);
    det += mat[2] * (mat[7] * mat[3] - mat[4] * mat[6]);
    return det;
}

/**
 * @brief Computes reciprocal lattice vectors from direct lattice vectors
 *
 * Calculates reciprocal lattice vectors using:
 * b_i = 2*pi * (a_j x a_k) / V
 * where V = det(a) is the unit cell volume.
 *
 * Result includes the 2*pi factor.
 *
 * @param lat_vec Direct lattice vectors (9 elements, row-major: a_1, a_2, a_3)
 * @param blat Reciprocal lattice vectors output (9 elements, row-major: b_1,
 * b_2, b_3)
 *
 * @note Calls error_msg if determinant < ELPH_EPS (singular matrix)
 */
void reciprocal_vecs(const ELPH_float* lat_vec, ELPH_float* restrict blat)
{
    ELPH_float det = det3x3(lat_vec);
    if (det < ELPH_EPS)
    {
        error_msg("Inverting singular matrix");
    }
    det = 2 * ELPH_PI / det;
    blat[0] = (lat_vec[4] * lat_vec[8] - lat_vec[7] * lat_vec[5]) * det;
    blat[1] = (lat_vec[5] * lat_vec[6] - lat_vec[3] * lat_vec[8]) * det;
    blat[2] = (lat_vec[3] * lat_vec[7] - lat_vec[6] * lat_vec[4]) * det;
    blat[3] = (lat_vec[2] * lat_vec[7] - lat_vec[1] * lat_vec[8]) * det;
    blat[4] = (lat_vec[0] * lat_vec[8] - lat_vec[2] * lat_vec[6]) * det;
    blat[5] = (lat_vec[1] * lat_vec[6] - lat_vec[0] * lat_vec[7]) * det;
    blat[6] = (lat_vec[1] * lat_vec[5] - lat_vec[2] * lat_vec[4]) * det;
    blat[7] = (lat_vec[2] * lat_vec[3] - lat_vec[0] * lat_vec[5]) * det;
    blat[8] = (lat_vec[0] * lat_vec[4] - lat_vec[1] * lat_vec[3]) * det;
}

/**
 * @brief Performs AXPY operation: Y = a*X + Y
 *
 * @param n Number of elements
 * @param a Complex scalar multiplier
 * @param X Complex input vector (n elements)
 * @param Y Complex input/output vector (n elements, modified in-place)
 */
void aXpY(const ND_int n, const ELPH_cmplx a, const ELPH_cmplx* X,
          ELPH_cmplx* Y)
{
    // ELPH_OMP_PAR_FOR_SIMD
    for (ND_int i = 0; i < n; ++i)
    {
        Y[i] += a * X[i];
    }
}

/**
 * @brief Transposes a 3x3 matrix (out-of-place)
 *
 * @param inmat Input 3x3 matrix (9 elements, row-major)
 * @param outmat Output transposed matrix (9 elements, row-major)
 */
void transpose3x3f(const ELPH_float* inmat, ELPH_float* restrict outmat)
{
    outmat[0] = inmat[0];
    outmat[1] = inmat[3];
    outmat[2] = inmat[6];

    outmat[3] = inmat[1];
    outmat[4] = inmat[4];
    outmat[5] = inmat[7];

    outmat[6] = inmat[2];
    outmat[7] = inmat[5];
    outmat[8] = inmat[8];
}

/**
 * @brief Transposes a 3x3 matrix in-place
 *
 * @param mat 3x3 matrix (9 elements, row-major, modified in-place)
 */
void transpose3x3f_inplace(ELPH_float* mat)
{
    swap_floats(mat + 1, mat + 3);
    swap_floats(mat + 2, mat + 6);
    swap_floats(mat + 5, mat + 7);
}

/**
 * @brief Finds maximum value in integer array
 *
 * @param in_arr Input integer array (nelements)
 * @param nelements Number of elements
 * @return Maximum value in array
 */
ND_int find_maxint(ND_int* in_arr, ND_int nelements)
{
    ND_int max = in_arr[0];
    for (ND_int imax = 0; imax < nelements; ++imax)
    {
        if (in_arr[imax] > max)
        {
            max = in_arr[imax];
        }
    }
    // printf("Max : %llu \n ", max);
    return max;
}

/**
 * @brief Finds maximum value in floating-point array
 *
 * @param in_arr Input floating-point array (nelements)
 * @param nelements Number of elements
 * @return Maximum value in array
 */
ELPH_float find_maxfloat(ELPH_float* in_arr, ND_int nelements)
{
    ELPH_float max = in_arr[0];
    for (ND_int imax = 0; imax < nelements; ++imax)
    {
        if (in_arr[imax] > max)
        {
            max = in_arr[imax];
        }
    }
    // printf("Max : %llu \n ", max);
    return max;
}

/**
 * @brief 3x3 matrix multiplication with optional transposes: C = op(A) * op(B)
 *
 * Performs matrix multiplication with optional transposes on A and/or B.
 *
 * @param A First 3x3 matrix (9 elements, row-major)
 * @param transA 'N' for A, 'T' for A^T
 * @param B Second 3x3 matrix (9 elements, row-major)
 * @param transB 'N' for B, 'T' for B^T
 * @param C Output 3x3 matrix (9 elements, row-major)
 *
 * @note Calls error_msg for invalid transpose flags
 */
void Gemm3x3f(const ELPH_float* A, const char transA, const ELPH_float* B,
              const char transB, ELPH_float* restrict C)
{
    if (transA == 'N' && transB == 'N')
    {
        C[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
        C[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
        C[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];
        C[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
        C[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
        C[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];
        C[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
        C[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
        C[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];
    }
    else if (transA == 'T' && transB == 'N')
    {
        C[0] = A[0] * B[0] + A[3] * B[3] + A[6] * B[6];
        C[1] = A[0] * B[1] + A[3] * B[4] + A[6] * B[7];
        C[2] = A[0] * B[2] + A[3] * B[5] + A[6] * B[8];
        C[3] = A[1] * B[0] + A[4] * B[3] + A[7] * B[6];
        C[4] = A[1] * B[1] + A[4] * B[4] + A[7] * B[7];
        C[5] = A[1] * B[2] + A[4] * B[5] + A[7] * B[8];
        C[6] = A[2] * B[0] + A[5] * B[3] + A[8] * B[6];
        C[7] = A[2] * B[1] + A[5] * B[4] + A[8] * B[7];
        C[8] = A[2] * B[2] + A[5] * B[5] + A[8] * B[8];
    }
    else if (transA == 'N' && transB == 'T')
    {
        C[0] = A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
        C[1] = A[0] * B[3] + A[1] * B[4] + A[2] * B[5];
        C[2] = A[0] * B[6] + A[1] * B[7] + A[2] * B[8];
        C[3] = A[3] * B[0] + A[4] * B[1] + A[5] * B[2];
        C[4] = A[3] * B[3] + A[4] * B[4] + A[5] * B[5];
        C[5] = A[3] * B[6] + A[4] * B[7] + A[5] * B[8];
        C[6] = A[6] * B[0] + A[7] * B[1] + A[8] * B[2];
        C[7] = A[6] * B[3] + A[7] * B[4] + A[8] * B[5];
        C[8] = A[6] * B[6] + A[7] * B[7] + A[8] * B[8];
    }
    else if (transA == 'T' && transB == 'T')
    {
        C[0] = A[0] * B[0] + A[3] * B[1] + A[6] * B[2];
        C[1] = A[0] * B[3] + A[3] * B[4] + A[6] * B[5];
        C[2] = A[0] * B[6] + A[3] * B[7] + A[6] * B[8];
        C[3] = A[1] * B[0] + A[4] * B[1] + A[7] * B[2];
        C[4] = A[1] * B[3] + A[4] * B[4] + A[7] * B[5];
        C[5] = A[1] * B[6] + A[4] * B[7] + A[7] * B[8];
        C[6] = A[2] * B[0] + A[5] * B[1] + A[8] * B[2];
        C[7] = A[2] * B[3] + A[5] * B[4] + A[8] * B[5];
        C[8] = A[2] * B[6] + A[5] * B[7] + A[8] * B[8];
    }
    else
    {
        error_msg("Wrong Trans type");
    }
    return;
}

/**
 * @brief Multiplies two 2x2 complex matrices: out = mat1 * mat2
 *
 * @param mat1 First 2x2 complex matrix (4 elements, row-major)
 * @param mat2 Second 2x2 complex matrix (4 elements, row-major)
 * @param out Output 2x2 complex matrix (4 elements, row-major)
 */
void matmul_Cmpl2x2(ELPH_cmplx* mat1, ELPH_cmplx* mat2,
                    ELPH_cmplx* restrict out)
{
    out[0] = mat1[0] * mat2[0] + mat1[1] * mat2[2];
    out[1] = mat1[0] * mat2[1] + mat1[1] * mat2[3];
    out[2] = mat1[2] * mat2[0] + mat1[3] * mat2[2];
    out[3] = mat1[2] * mat2[1] + mat1[3] * mat2[3];
}

/* functions related to fft */

/**
 * @brief Converts FFT index to Miller index
 *
 * Maps FFT indices [0, N) to Miller indices:
 * - For even N: [-N/2, N/2)
 * - For odd N: [-(N-1)/2, (N-1)/2]
 *
 * @param idx_in FFT index in range [0, FFT_dimension)
 * @param FFT_dimension FFT grid size
 * @return Miller index
 */
ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension)
{
    ND_int mid_pnt = (FFT_dimension - 1) / 2 + 1;
    if (idx_in < mid_pnt)
    {
        return idx_in;
    }
    else
    {
        return (idx_in - mid_pnt) - FFT_dimension / 2;
    }
}

/**
 * @brief Converts Miller index to FFT index
 *
 * Maps Miller indices (typically [-N/2, N/2)) to FFT indices [0, N).
 * FFT libraries assume indices run from [0, N), but Miller indices
 * conventionally run from [-N/2, N/2).
 *
 * @param idx_in Miller index (will be rounded to nearest integer)
 * @param FFT_dimension FFT grid size
 * @return FFT index in range [0, N-1]
 */
int get_fft_idx(ELPH_float idx_in, int FFT_dimension)
{
    int temp_idx = rint(idx_in);
    if (temp_idx >= 0)
    {
        return temp_idx;
    }
    else
    {
        return FFT_dimension + temp_idx;
    }
}

/**
 * @brief Finds index of k-point in list using crystal coordinates
 *
 * Searches for a k-point in the list by comparing crystal coordinates
 * with periodic boundary conditions (difference modulo reciprocal lattice
 * vector). Two k-points are considered equivalent if their difference (mod G)
 * has norm < ELPH_EPS.
 *
 * @param nkpts Number of k-points in list
 * @param kpts_list k-point list in crystal coordinates (nkpts,3)
 * @param kpt k-point to search for in crystal coordinates (3)
 * @return Index of k-point in list, or -1 if not found
 */
ND_int find_kidx_in_list(ND_int nkpts, const ELPH_float* kpts_list,
                         const ELPH_float* kpt)
{
    ND_int kidx = -1;
    for (ND_int ik = 0; ik < nkpts; ++ik)
    {
        const ELPH_float* ik_vec_tmp = kpts_list + 3 * ik;
        ELPH_float sum = 0;
        for (int i = 0; i < 3; ++i)
        {
            ELPH_float diff_tmp = ik_vec_tmp[i] - kpt[i];
            diff_tmp = diff_tmp - rint(diff_tmp);
            sum += diff_tmp * diff_tmp;
        }
        sum = sqrt(sum);
        if (sum < ELPH_EPS)
        {
            kidx = ik;
            break;
        }
    }
    return kidx;
}

/**
 * @brief Finds K+Q point indices in k-point grid
 *
 * For each k-point in the grid, finds the index of k+Q (with periodic
 * boundary conditions). Used to establish k to k+q mapping for electron-phonon
 * calculations.
 *
 * @param Nbz Number of k-points in Brillouin zone
 * @param kpoints k-point grid in crystal coordinates (Nbz,3)
 * @param Q_pt q-point in crystal coordinates (3)
 * @param KplusQidxs Output array of k+Q indices (Nbz)
 *
 * @note Calls error_msg if any k+Q point is not found (incommensurate grids)
 */
void get_KplusQ_idxs(const ND_int Nbz, const ELPH_float* kpoints,
                     const ELPH_float* Q_pt, int* KplusQidxs)
{
    for (ND_int i = 0; i < Nbz; ++i)
    {
        const ELPH_float* ktemp = kpoints + 3 * i;

        ELPH_float KplusQ[3];
        KplusQ[0] = ktemp[0] + Q_pt[0];
        KplusQ[1] = ktemp[1] + Q_pt[1];
        KplusQ[2] = ktemp[2] + Q_pt[2];

        KplusQidxs[i] = find_kidx_in_list(Nbz, kpoints, KplusQ);

        if (KplusQidxs[i] < 0)
        {
            error_msg(
                "K+Q point cannot be found in kgrid. "
                "The k-grid and q-grid are not commensurate! \n");
        }
    }
}

/**
 * @brief Swaps two integer values
 *
 * @param a Pointer to first integer
 * @param b Pointer to second integer
 */
void swap_ints(int* a, int* b)
{
    const int c = *b;
    *b = *a;
    *a = c;
}

/**
 * @brief Swaps two floating-point values
 *
 * @param a Pointer to first float
 * @param b Pointer to second float
 */
void swap_floats(ELPH_float* a, ELPH_float* b)
{
    const ELPH_float c = *b;
    *b = *a;
    *a = c;
}

/**
 * @file
 * @brief Numerical functions and linear algebra utilities
 *
 * Provides special functions (Legendre polynomials, spherical harmonics),
 * linear algebra operations (matrix multiplication, dot products, transposes),
 * FFT utilities, spline interpolation, and various helper functions.
 */

#pragma once
#include <stdbool.h>

#include "elphC.h"

/**
 * @brief Finds index of k-point in list using crystal coordinates
 *
 * @param nkpts Number of k-points in list
 * @param kpts_list k-point list in crystal coordinates (nkpts,3)
 * @param kpt k-point to search for in crystal coordinates (3)
 * @return Index of k-point in list, or -1 if not found
 */
ND_int find_kidx_in_list(ND_int nkpts, const ELPH_float* kpts_list,
                         const ELPH_float* kpt);

/**
 * @brief Finds K+Q point indices in k-point grid
 *
 * @param Nbz Number of k-points in Brillouin zone
 * @param kpoints k-point grid in crystal coordinates (Nbz,3)
 * @param Q_pt q-point in crystal coordinates (3)
 * @param KplusQidxs Output array of k+Q indices (Nbz)
 */
void get_KplusQ_idxs(const ND_int Nbz, const ELPH_float* kpoints,
                     const ELPH_float* Q_pt, int* KplusQidxs);

/**
 * @def dot3_macro(a, b)
 * @brief Computes 3D dot product as macro
 *
 * Calculates: a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
 *
 * @param a First 3D vector
 * @param b Second 3D vector
 * @return Dot product value
 *
 * @warning Arguments evaluated multiple times - avoid expressions with side
 * effects
 */
#define dot3_macro(a, b) ((a)[0] * (b)[0] + (a)[1] * (b)[1] + (a)[2] * (b)[2])

/**
 * @def MIN(a, b)
 * @brief Returns minimum of two values
 *
 * @param a First value
 * @param b Second value
 * @return Minimum value
 *
 * @warning Arguments evaluated multiple times - avoid expressions with side
 * effects
 */
#if !defined(MIN)
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

/**
 * @def MAX(a, b)
 * @brief Returns maximum of two values
 *
 * @param a First value
 * @param b Second value
 * @return Maximum value
 *
 * @warning Arguments evaluated multiple times - avoid expressions with side
 * effects
 */
#if !defined(MAX)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

/* numerical_func.c */

/**
 * @brief Computes associated Legendre polynomial P_l^m(x)
 *
 * @param l_val Degree l
 * @param m_val Order m
 * @param x_in Argument x in range [-1, 1]
 * @return Associated Legendre polynomial value
 */
ELPH_float legendre(int l_val, int m_val, ELPH_float x_in);

/**
 * @brief Computes real spherical harmonic Y_l^m for a direction vector
 *
 * @param l_val Degree l
 * @param m_val Order m
 * @param vec Direction vector (3 components, Cartesian)
 * @return Real spherical harmonic value
 */
ELPH_float Ylm(int l_val, int m_val, ELPH_float* vec);

/**
 * @brief Performs numerical integration using Simpson's 1/3 rule
 *
 * @param func_vals Function values at grid points (npts)
 * @param dx Grid spacing at each point (npts)
 * @param npts Number of grid points
 * @return Approximate integral value
 */
ELPH_float simpson(const ELPH_float* func_vals, const ELPH_float* dx,
                   ND_int npts);

/**
 * @brief Computes cosine of angle between two 3D vectors
 *
 * @param vec1 First vector (3 components)
 * @param vec2 Second vector (3 components)
 * @return Cosine of angle between vectors
 */
ELPH_float cos_angle_bw_Vec(const ELPH_float* vec1, const ELPH_float* vec2);

/**
 * @brief Matrix-vector multiplication for 3x3 matrix and 3D vector
 *
 * @param Mat 3x3 matrix (9 elements, row-major)
 * @param vec Input vector (3 elements)
 * @param trans If true, use transpose of Mat
 * @param out Output vector (3 elements)
 */
void MatVec3f(const ELPH_float* Mat, const ELPH_float* vec, const bool trans,
              ELPH_float* restrict out);

/**
 * @brief Computes complex dot product: conj(vec1) . vec2
 *
 * @param vec1 First complex vector (n elements)
 * @param vec2 Second complex vector (n elements)
 * @param n Number of elements
 * @return Complex dot product
 */
ELPH_cmplx Cmplxdot(const ELPH_cmplx* vec1, const ELPH_cmplx* vec2,
                    const ND_int n);

/**
 * @brief Normalizes a complex vector in-place
 *
 * @param vec Complex vector to normalize (n elements)
 * @param n Number of elements
 */
void normalize_Cmplx_vec(ELPH_cmplx* vec, const ND_int n);

/**
 * @brief Computes determinant of 3x3 matrix
 *
 * @param mat 3x3 matrix (9 elements, row-major)
 * @return Determinant value
 */
ELPH_float det3x3(const ELPH_float* mat);

/**
 * @brief Computes reciprocal lattice vectors from direct lattice vectors
 *
 * Result includes the 2*pi factor.
 *
 * @param lat_vec Direct lattice vectors (9 elements, row-major)
 * @param blat Reciprocal lattice vectors output (9 elements, row-major)
 */
void reciprocal_vecs(const ELPH_float* lat_vec, ELPH_float* restrict blat);

/**
 * @brief Performs AXPY operation: Y = a*X + Y
 *
 * @param n Number of elements
 * @param a Complex scalar multiplier
 * @param X Complex input vector (n elements)
 * @param Y Complex input/output vector (n elements)
 */
void aXpY(const ND_int n, const ELPH_cmplx a, const ELPH_cmplx* X,
          ELPH_cmplx* Y);

/**
 * @brief Transposes a 3x3 matrix (out-of-place)
 *
 * @param inmat Input 3x3 matrix (9 elements, row-major)
 * @param outmat Output transposed matrix (9 elements, row-major)
 */
void transpose3x3f(const ELPH_float* inmat, ELPH_float* restrict outmat);

/**
 * @brief Transposes a 3x3 matrix in-place
 *
 * @param mat 3x3 matrix (9 elements, row-major)
 */
void transpose3x3f_inplace(ELPH_float* mat);

/**
 * @brief Finds maximum value in integer array
 *
 * @param in_arr Input integer array (nelements)
 * @param nelements Number of elements
 * @return Maximum value in array
 */
ND_int find_maxint(ND_int* in_arr, ND_int nelements);

/**
 * @brief Finds maximum value in floating-point array
 *
 * @param in_arr Input floating-point array (nelements)
 * @param nelements Number of elements
 * @return Maximum value in array
 */
ELPH_float find_maxfloat(ELPH_float* in_arr, ND_int nelements);

/**
 * @brief 3x3 matrix multiplication: C = op(A) * op(B)
 *
 * @param A First 3x3 matrix (9 elements, row-major)
 * @param transA 'N' for A, 'T' for A^T
 * @param B Second 3x3 matrix (9 elements, row-major)
 * @param transB 'N' for B, 'T' for B^T
 * @param C Output 3x3 matrix (9 elements, row-major)
 */
void Gemm3x3f(const ELPH_float* A, const char transA, const ELPH_float* B,
              const char transB, ELPH_float* restrict C);

/**
 * @brief Multiplies two 2x2 complex matrices: out = mat1 * mat2
 *
 * @param mat1 First 2x2 complex matrix (4 elements, row-major)
 * @param mat2 Second 2x2 complex matrix (4 elements, row-major)
 * @param out Output 2x2 complex matrix (4 elements, row-major)
 */
void matmul_Cmpl2x2(ELPH_cmplx* mat1, ELPH_cmplx* mat2,
                    ELPH_cmplx* restrict out);

/**
 * @brief Converts Miller index to FFT index
 *
 * @param idx_in Miller index (will be rounded to nearest integer)
 * @param FFT_dimension FFT grid size
 * @return FFT index in range [0, N-1]
 */
int get_fft_idx(ELPH_float idx_in, int FFT_dimension);

/**
 * @brief Converts FFT index to Miller index
 *
 * @param idx_in FFT index in range [0, FFT_dimension)
 * @param FFT_dimension FFT grid size
 * @return Miller index
 */
ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension);

/**
 * @brief Swaps two integer values
 *
 * @param a Pointer to first integer
 * @param b Pointer to second integer
 */
void swap_ints(int* a, int* b);

/**
 * @brief Swaps two floating-point values
 *
 * @param a Pointer to first float
 * @param b Pointer to second float
 */
void swap_floats(ELPH_float* a, ELPH_float* b);

/* spline.c*/

/**
 * @brief Performs cubic spline interpolation at given point
 *
 * @param x Point at which to interpolate
 * @param inear Index hint for nearest point (used for efficiency)
 * @param xi Array of x coordinates (grid points)
 * @param yi Array of y values at grid points
 * @param dy Array of second derivatives for spline
 * @return Interpolated value at x
 */
ELPH_float spline_interpolate(const ELPH_float x, ND_int inear,
                              const ELPH_float* xi, const ELPH_float* yi,
                              const ELPH_float* dy);

/**
 * @brief Prepares cubic spline coefficients for interpolation
 *
 * @param nvals Number of data points
 * @param xin Array of x coordinates (grid points, nvals)
 * @param yin Array of y values at grid points (nvals)
 * @param dy Output array for second derivatives (nvals)
 */
void prepare_spline(const ND_int nvals, ELPH_float* xin, ELPH_float* yin,
                    ELPH_float* dy);

/**
 * @brief Complex matrix multiplication wrapper (ZGEMM/CGEMM)
 *
 * Computes: C = alpha * op(A) * op(B) + beta * C
 *
 * @param TransA Operation on matrix A: 'N', 'T', or 'C'
 * @param TransB Operation on matrix B: 'N', 'T', or 'C'
 * @param arr_A Pointer to matrix A
 * @param arr_B Pointer to matrix B
 * @param arr_C Pointer to matrix C (input/output)
 * @param alpha Scalar multiplier for product
 * @param beta Scalar multiplier for C
 * @param ldA Leading dimension of A
 * @param ldB Leading dimension of B
 * @param ldC Leading dimension of C
 * @param m Number of rows in op(A) and C
 * @param n Number of columns in op(B) and C
 * @param k Number of columns in op(A) and rows in op(B)
 */
void matmul_cmplx(const char TransA, const char TransB, const ELPH_cmplx* arr_A,
                  const ELPH_cmplx* arr_B, ELPH_cmplx* arr_C,
                  const ELPH_cmplx alpha, const ELPH_cmplx beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k);

/**
 * @brief Real matrix multiplication wrapper (DGEMM/SGEMM)
 *
 * Computes: C = alpha * op(A) * op(B) + beta * C
 *
 * @param TransA Operation on matrix A: 'N' or 'T'
 * @param TransB Operation on matrix B: 'N' or 'T'
 * @param arr_A Pointer to matrix A
 * @param arr_B Pointer to matrix B
 * @param arr_C Pointer to matrix C (input/output)
 * @param alpha Scalar multiplier for product
 * @param beta Scalar multiplier for C
 * @param ldA Leading dimension of A
 * @param ldB Leading dimension of B
 * @param ldC Leading dimension of C
 * @param m Number of rows in op(A) and C
 * @param n Number of columns in op(B) and C
 * @param k Number of columns in op(A) and rows in op(B)
 */
void matmul_float(const char TransA, const char TransB, const ELPH_float* arr_A,
                  const ELPH_float* arr_B, ELPH_float* arr_C,
                  const ELPH_float alpha, const ELPH_float beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k);

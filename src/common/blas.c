/**
 * @file
 * @brief BLAS wrapper functions for matrix operations
 *
 * Provides simplified wrappers around CBLAS matrix multiplication routines
 * (GEMM) with automatic precision selection and overflow checking.
 * Supports both real and complex matrix operations.
 */

#define CBLAS_INT int
#include <limits.h>

#include "cblas.h"
#include "elphC.h"
#include "error.h"
#include "numerical_func.h"

/**
 * @def BLAS_INT_MAX
 * @brief Maximum value for BLAS integer parameters
 */
#define BLAS_INT_MAX INT_MAX

#if defined(COMPILE_ELPH_DOUBLE)
/**
 * @def CBLAS_cmplx(FUN_NAME)
 * @brief Macro to select double precision complex BLAS function (zgemm, etc.)
 */
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_z##FUN_NAME

/**
 * @def CBLAS_float(FUN_NAME)
 * @brief Macro to select double precision real BLAS function (dgemm, etc.)
 */
#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_d##FUN_NAME
#else
/**
 * @def CBLAS_cmplx(FUN_NAME)
 * @brief Macro to select single precision complex BLAS function (cgemm, etc.)
 */
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_c##FUN_NAME

/**
 * @def CBLAS_float(FUN_NAME)
 * @brief Macro to select single precision real BLAS function (sgemm, etc.)
 */
#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_s##FUN_NAME
#endif

static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans);

/**
 * @brief Complex matrix multiplication wrapper (ZGEMM/CGEMM)
 *
 * Computes the matrix operation:
 * \f$ C = \alpha \cdot \text{op}(A) \cdot \text{op}(B) + \beta \cdot C \f$
 *
 * where \f$ \text{op}(X) \f$ can be \f$ X \f$, \f$ X^T \f$, or \f$ X^H \f$
 * (identity, transpose, or conjugate transpose).
 *
 * Matrix dimensions after applying operations:
 * - \f$ \text{op}(A) \f$: \f$ m \times k \f$
 * - \f$ \text{op}(B) \f$: \f$ k \times n \f$
 * - \f$ C \f$: \f$ m \times n \f$
 *
 * @param TransA Operation on matrix A: 'N' (none), 'T' (transpose), 'C'
 * (conjugate transpose)
 * @param TransB Operation on matrix B: 'N' (none), 'T' (transpose), 'C'
 * (conjugate transpose)
 * @param arr_A Pointer to matrix A (row-major order)
 * @param arr_B Pointer to matrix B (row-major order)
 * @param arr_C Pointer to matrix C (row-major order, input/output)
 * @param alpha Scalar multiplier for the product
 * @param beta Scalar multiplier for C
 * @param ldA Leading dimension of A (number of columns in storage)
 * @param ldB Leading dimension of B (number of columns in storage)
 * @param ldC Leading dimension of C (number of columns in storage)
 * @param m Number of rows in op(A) and C
 * @param n Number of columns in op(B) and C
 * @param k Number of columns in op(A) and rows in op(B)
 *
 * @note All dimension parameters are checked for overflow before passing to
 * BLAS
 * @note Uses row-major storage order (C-style arrays)
 */
void matmul_cmplx(const char TransA, const char TransB, const ELPH_cmplx* arr_A,
                  const ELPH_cmplx* arr_B, ELPH_cmplx* arr_C,
                  const ELPH_cmplx alpha, const ELPH_cmplx beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k)
{
    CHECK_OVERFLOW_ERROR(ldA, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldB, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldC, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(m, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(n, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(k, BLAS_INT_MAX);
    CBLAS_cmplx(gemm)(CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB),
                      m, n, k, &alpha, arr_A, ldA, arr_B, ldB, &beta, arr_C,
                      ldC);
}

/**
 * @brief Real matrix multiplication wrapper (DGEMM/SGEMM)
 *
 * Computes the matrix operation:
 * \f$ C = \alpha \cdot \text{op}(A) \cdot \text{op}(B) + \beta \cdot C \f$
 *
 * where \f$ \text{op}(X) \f$ can be \f$ X \f$ or \f$ X^T \f$
 * (identity or transpose).
 *
 * Matrix dimensions after applying operations:
 * - \f$ \text{op}(A) \f$: \f$ m \times k \f$
 * - \f$ \text{op}(B) \f$: \f$ k \times n \f$
 * - \f$ C \f$: \f$ m \times n \f$
 *
 * @param TransA Operation on matrix A: 'N' (none) or 'T' (transpose)
 * @param TransB Operation on matrix B: 'N' (none) or 'T' (transpose)
 * @param arr_A Pointer to matrix A (row-major order)
 * @param arr_B Pointer to matrix B (row-major order)
 * @param arr_C Pointer to matrix C (row-major order, input/output)
 * @param alpha Scalar multiplier for the product
 * @param beta Scalar multiplier for C
 * @param ldA Leading dimension of A (number of columns in storage)
 * @param ldB Leading dimension of B (number of columns in storage)
 * @param ldC Leading dimension of C (number of columns in storage)
 * @param m Number of rows in op(A) and C
 * @param n Number of columns in op(B) and C
 * @param k Number of columns in op(A) and rows in op(B)
 *
 * @note All dimension parameters are checked for overflow before passing to
 * BLAS
 * @note Uses row-major storage order (C-style arrays)
 */
void matmul_float(const char TransA, const char TransB, const ELPH_float* arr_A,
                  const ELPH_float* arr_B, ELPH_float* arr_C,
                  const ELPH_float alpha, const ELPH_float beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k)
{
    CHECK_OVERFLOW_ERROR(ldA, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldB, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldC, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(m, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(n, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(k, BLAS_INT_MAX);
    CBLAS_float(gemm)(CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB),
                      m, n, k, alpha, arr_A, ldA, arr_B, ldB, beta, arr_C, ldC);
}

/**
 * @brief Converts character transpose flag to CBLAS enumeration
 *
 * @param Trans Character flag: 'N' (no transpose), 'T' (transpose), 'C'
 * (conjugate transpose)
 * @return enum CBLAS_TRANSPOSE Corresponding CBLAS transpose operation
 *
 * @note Calls error_msg() and returns CblasNoTrans for invalid input
 */
static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans)
{
    if (Trans == 'N')
    {
        return CblasNoTrans;
    }
    else if (Trans == 'T')
    {
        return CblasTrans;
    }
    else if (Trans == 'C')
    {
        return CblasConjTrans;
    }
    else
    {
        error_msg("Can only take 'C' or 'T' or 'N' for Trans input");
        return CblasNoTrans;
    }
}

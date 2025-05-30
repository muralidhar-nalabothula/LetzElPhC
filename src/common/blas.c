/*
This file contails blas wrappers
*/
#define CBLAS_INT int
// note the above header must be called before #include "cblas.h"
#include <limits.h>

#include "cblas.h"
#include "elphC.h"
#include "error.h"
#include "numerical_func.h"

#define BLAS_INT_MAX INT_MAX

#if defined(COMPILE_ELPH_DOUBLE)
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_z##FUN_NAME

#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_d##FUN_NAME

#else
#define CBLAS_cmplx(FUN_NAME) CBLAS_cmplx_hidden(FUN_NAME)
#define CBLAS_cmplx_hidden(FUN_NAME) cblas_c##FUN_NAME

#define CBLAS_float(FUN_NAME) CBLAS_float_hidden(FUN_NAME)
#define CBLAS_float_hidden(FUN_NAME) cblas_s##FUN_NAME
#endif

static enum CBLAS_TRANSPOSE get_gemmn_T(char Trans);

// =================================================================
// NM : Fix me : Check for overflows

void matmul_cmplx(const char TransA, const char TransB, const ELPH_cmplx* arr_A,
                  const ELPH_cmplx* arr_B, ELPH_cmplx* arr_C,
                  const ELPH_cmplx alpha, const ELPH_cmplx beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k)
{
    // first check for interger overflows before passing to blas library
    //
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

void matmul_float(const char TransA, const char TransB, const ELPH_float* arr_A,
                  const ELPH_float* arr_B, ELPH_float* arr_C,
                  const ELPH_float alpha, const ELPH_float beta,
                  const ND_int ldA, const ND_int ldB, const ND_int ldC,
                  const ND_int m, const ND_int n, const ND_int k)
{
    // first check for interger overflows before passing to blas library
    //
    CHECK_OVERFLOW_ERROR(ldA, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldB, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(ldC, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(m, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(n, BLAS_INT_MAX);
    CHECK_OVERFLOW_ERROR(k, BLAS_INT_MAX);

    CBLAS_float(gemm)(CblasRowMajor, get_gemmn_T(TransA), get_gemmn_T(TransB),
                      m, n, k, alpha, arr_A, ldA, arr_B, ldB, beta, arr_C, ldC);
}

// =================================================================
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

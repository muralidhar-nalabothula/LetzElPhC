/**
 * @file
 * @brief OpenMP pragma definitions for conditional compilation
 *
 * Provides macro definitions for OpenMP directives that expand to actual
 * pragmas when ELPH_OMP_PARALLEL_BUILD is defined, or to nothing otherwise.
 * This allows the same source code to be compiled with or without OpenMP
 * support without #ifdef cluttering the implementation files.
 */

#pragma once
#include "elphC.h"

#if defined(ELPH_OMP_PARALLEL_BUILD)
#include <omp.h>

/**
 * @def ELPH_OMP_PAR_FOR_SIMD
 * @brief OpenMP parallel for loop with SIMD vectorization hint
 *
 * Expands to: #pragma omp parallel for simd
 */
#define ELPH_OMP_PAR_FOR_SIMD _Pragma("omp parallel for simd")

/**
 * @def ELPH_OMP_PAR_FOR
 * @brief OpenMP parallel for loop
 *
 * Expands to: #pragma omp parallel for
 */
#define ELPH_OMP_PAR_FOR _Pragma("omp parallel for")

/**
 * @def ELPH_OMP_FOR
 * @brief OpenMP work-sharing for loop (within existing parallel region)
 *
 * Expands to: #pragma omp for
 */
#define ELPH_OMP_FOR _Pragma("omp for")

/**
 * @def ELPH_OMP_PAR_SIMD
 * @brief OpenMP SIMD vectorization hint
 *
 * Expands to: #pragma omp simd
 */
#define ELPH_OMP_PAR_SIMD _Pragma("omp simd")

/**
 * @def ELPH_OMP_PAR_CRITICAL
 * @brief OpenMP critical section (mutual exclusion)
 *
 * Expands to: #pragma omp critical
 */
#define ELPH_OMP_PAR_CRITICAL _Pragma("omp critical")

/**
 * @def ELPH_OMP_PAR
 * @brief OpenMP parallel region
 *
 * Expands to: #pragma omp parallel
 */
#define ELPH_OMP_PAR _Pragma("omp parallel")

/**
 * @def ELPH_OMP_ATOMIC
 * @brief OpenMP atomic operation
 *
 * Expands to: #pragma omp atomic
 */
#define ELPH_OMP_ATOMIC _Pragma("omp atomic")

/**
 * @def ELPH_OMP_PAR_COLLAPSE_3
 * @brief OpenMP parallel for with 3 nested loops collapsed
 *
 * Expands to: #pragma omp parallel for collapse(3)
 */
#define ELPH_OMP_PAR_COLLAPSE_3 _Pragma("omp parallel for collapse(3)")

/**
 * @def ELPH_OMP_PAR_COLLAPSE_2
 * @brief OpenMP parallel for with 2 nested loops collapsed
 *
 * Expands to: #pragma omp parallel for collapse(2)
 */
#define ELPH_OMP_PAR_COLLAPSE_2 _Pragma("omp parallel for collapse(2)")

/**
 * @def ELPH_OMP_SINGLE
 * @brief OpenMP single construct (executed by only one thread)
 *
 * Expands to: #pragma omp single
 */
#define ELPH_OMP_SINGLE _Pragma("omp single")

#else
/* When OpenMP is disabled, all macros expand to nothing */

/**
 * @def ELPH_OMP_PAR_FOR_SIMD
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR_FOR_SIMD

/**
 * @def ELPH_OMP_PAR_FOR
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR_FOR

/**
 * @def ELPH_OMP_FOR
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_FOR

/**
 * @def ELPH_OMP_PAR_SIMD
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR_SIMD

/**
 * @def ELPH_OMP_PAR_CRITICAL
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR_CRITICAL

/**
 * @def ELPH_OMP_PAR
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR

/**
 * @def ELPH_OMP_ATOMIC
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_ATOMIC

/**
 * @def ELPH_OMP_PAR_COLLAPSE_3
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR_COLLAPSE_3

/**
 * @def ELPH_OMP_PAR_COLLAPSE_2
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_PAR_COLLAPSE_2

/**
 * @def ELPH_OMP_SINGLE
 * @brief No-op when OpenMP is disabled
 */
#define ELPH_OMP_SINGLE

#endif

/**
 * @file elphC.h
 * @brief Core type definitions and configuration for the ELPH library
 * 
 * This header defines the fundamental data types used throughout the ELPH
 * library, including floating-point precision selection and MPI type mappings.
 * It enforces C99 standard compliance and complex number support.
 */

#pragma once

/* Ensure C99 standard compliance */
#if defined __STDC_VERSION__ && __STDC_VERSION__ >= 199901L
#ifdef __STDC_NO_COMPLEX__
#error Your compiler does not support C99 complex numbers, Please use a supported compiler.
#endif
#else
#error Your compiler does not support C99 standard.
#endif

/**
 * @typedef ND_int
 * @brief Standard integer type for array dimensions and indices
 *
 * Defined as a 64-bit signed integer to support large array dimensions
 */
typedef long long int ND_int;

/**
 * @def ELPH_MPI_ND_INT
 * @brief MPI datatype corresponding to ND_int
 */
#define ELPH_MPI_ND_INT MPI_LONG_LONG_INT

/* Precision-dependent type definitions */
#if defined(COMPILE_ELPH_DOUBLE)
/**
 * @typedef ELPH_float
 * @brief Floating-point type (double precision when COMPILE_ELPH_DOUBLE is defined)
 */
typedef double ELPH_float;
typedef double _Complex ELPH_cmplx;

#define ELPH_MPI_float MPI_DOUBLE
#define ELPH_MPI_cmplx MPI_C_DOUBLE_COMPLEX

#define ELPH_NC4_IO_FLOAT NC_DOUBLE
#else
typedef float ELPH_float;
typedef float _Complex ELPH_cmplx;

#define ELPH_MPI_float MPI_FLOAT
#define ELPH_MPI_cmplx MPI_C_FLOAT_COMPLEX

#define ELPH_NC4_IO_FLOAT NC_FLOAT
#endif


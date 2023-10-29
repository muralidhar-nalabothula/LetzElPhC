#pragma once

#ifdef __STDC_NO_COMPLEX__
#error Your compiler does not C99 complex numbers, Please use a supported compiler.
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "nd_array/src/nd_array.h"
#include <mpi.h>
#include "common/omp_pragma_def.h"
#include <fftw3.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <limits.h>
#if defined(ELPH_OMP_PARALLEL_BUILD)
    #include <omp.h>
#endif



#define ELPH_MPI_ND_INT MPI_LONG_LONG_INT

#if defined(COMPILE_ELPH_DOUBLE)

    typedef double ELPH_float;
    typedef double complex ELPH_cmplx;

    /* Constants */
    #define ELPH_EPS 1e-6
    #define Nd_floatS d
    #define Nd_cmplxS z


    #define ELPH_PI 3.1415927
    #define ELPH_SQRT2 1.4142136
    #define ELPH_e2 2.0 // in Ry units

    #define ELPH_MPI_float MPI_DOUBLE
    #define ELPH_MPI_cmplx MPI_C_DOUBLE_COMPLEX

    // fftw function
    #define fftw_fun(FUN_NAME) fftw3_fun_HIDDEN(FUN_NAME)
    #define fftw3_fun_HIDDEN(FUN_NAME) fftw_##FUN_NAME


#else
    typedef float ELPH_float;
    typedef float complex ELPH_cmplx;

    /* Constants */
    #define ELPH_EPS 1e-6
    #define Nd_floatS s
    #define Nd_cmplxS c


    #define ELPH_PI 3.1415927f
    #define ELPH_SQRT2 1.4142136f
    #define ELPH_e2 2.0f // in Ry units

    #define ELPH_MPI_float MPI_FLOAT
    #define ELPH_MPI_cmplx MPI_C_FLOAT_COMPLEX

    // fftw function
    #define fftw_fun(FUN_NAME) fftw3_fun_HIDDEN(FUN_NAME)
    #define fftw3_fun_HIDDEN(FUN_NAME) fftwf_##FUN_NAME

#endif






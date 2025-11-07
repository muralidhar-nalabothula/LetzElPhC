# Installing the Code

## Mandatory Requirements

- GNU Make
- C99 compiler with complex number support (GCC, Clang, ICC, AMD C-Compiler, MinGW, PGI, Arm C compilers)
- MPI implementation (Open-MPI, MPICH, Intel MPI, Microsoft MPI)
- FFTW-3 or Intel MKL
- HDF5 and NETCDF-4 libraries with Parallel IO support
- BLAS library (OpenBLAS, BLIS, Intel-MKL, Atlas)

## Installation Process

LetzElPhC uses a standard **Make** build system. Sample make files are in the `sample_config` directory.  

1. Copy a sample to the `src` directory and rename it `make.inc`.  
2. Edit `make.inc` according to your requirements.  
3. In `src/` execute:

```bash
$ make
$ make -j n   # where n is the number of processes
```

Upon success, the `lelphc` executable will be in the `src/` directory.

### Key Variables in `make.inc`

```make
CC        := mpicc          # MPI C compiler
CFLAGS    := -O3            # Optimization
LD_FLAGS  :=                # Additional linker flags

# Optional OpenMP support
# OPENMP_FLAGS := -DELPH_OMP_PARALLEL_BUILD
# CFLAGS      += -fopenmp   # or -qopenmp for Intel
# LD_FLAGS    += -fopenmp   # or -qopenmp for Intel

# FFTW3
FFTW_INC  := -I/opt/homebrew/include
FFTW3_LIB := -L/opt/homebrew/lib -lfftw3_threads -lfftw3f -lfftw3f_omp -lfftw3_omp -lfftw3

# BLAS
BLAS_LIB  := -L/opt/homebrew/opt/openblas/lib -lopenblas

# NetCDF
NETCDF_INC := -I/Users/murali/softwares/core/include
NETCDF_LIB := -L/Users/murali/softwares/core/lib -lnetcdf

# HDF5
HDF5_LIB := -L/opt/homebrew/lib -lhdf5

# Extra include dirs / libs
INC_DIRS := 
LIBS     := 
```

> Notes: Add `-DCOMPILE_ELPH_DOUBLE` to `CFLAGS` for double precision.  
> Use `-DELPH_OMP_PARALLEL_BUILD` and enable OpenMP flags for parallel builds.
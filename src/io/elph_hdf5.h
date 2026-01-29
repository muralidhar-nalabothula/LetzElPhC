#pragma once

#include <hdf5.h>
#include <mpi.h>

#include "common/dtypes.h"
#include "common/error.h"

// Define HDF5 types based on compilation
#if defined(COMPILE_ELPH_DOUBLE)
#define ELPH_H5_IO_FLOAT H5T_NATIVE_DOUBLE
#else
#define ELPH_H5_IO_FLOAT H5T_NATIVE_FLOAT
#endif

// Error handling macro for HDF5
#define H5_ERR(e)                                                             \
    if ((e) < 0)                                                              \
    {                                                                         \
        fprintf(stderr, "HDF5 Error in %s at line %d\n", __func__, __LINE__); \
        H5Eprint2(H5E_DEFAULT, stderr);                                       \
        error_msg("hdf5_error");                                              \
    }

// Wrapper functions

// File operations
hid_t elph_h5_create_file_par(const char* filename, MPI_Comm comm,
                              MPI_Info info);
hid_t elph_h5_open_file_par(const char* filename, int read_only, MPI_Comm comm,
                            MPI_Info info);
void elph_h5_close_file(hid_t file_id);

// Dataset operations
// Returns dataset_id in *dset_id.
void elph_h5_def_var(hid_t file_id, const char* name, hid_t type_id, int rank,
                     const ND_int* dims, const char** dim_names,
                     const size_t* chunk_sizes, hid_t* dset_id);

void elph_h5_write_var(hid_t dset_id, hid_t mem_type_id, const void* data);
void elph_h5_write_vara(hid_t dset_id, hid_t mem_type_id, const hsize_t* start,
                        const hsize_t* count, const void* data);

void elph_h5_read_var(hid_t dset_id, hid_t mem_type_id, void* data);
void elph_h5_read_vara(hid_t dset_id, hid_t mem_type_id, const hsize_t* start,
                       const hsize_t* count, void* data);

// Helper to close dataset
void elph_h5_close_var(hid_t dset_id);

// Open variable (dataset) by name
hid_t elph_h5_open_var(hid_t file_id, const char* name);

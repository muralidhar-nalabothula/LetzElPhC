#include "elph_hdf5.h"

#include <stdlib.h>
#include <string.h>

hid_t elph_h5_create_file_par(const char* filename, MPI_Comm comm,
                              MPI_Info info)
{
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5_ERR(plist_id);

    herr_t ret = H5Pset_fapl_mpio(plist_id, comm, info);
    H5_ERR(ret);

    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5_ERR(file_id);

    ret = H5Pclose(plist_id);
    H5_ERR(ret);

    return file_id;
}

hid_t elph_h5_open_file_par(const char* filename, int read_only, MPI_Comm comm,
                            MPI_Info info)
{
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5_ERR(plist_id);

    herr_t ret = H5Pset_fapl_mpio(plist_id, comm, info);
    H5_ERR(ret);

    unsigned flags = read_only ? H5F_ACC_RDONLY : H5F_ACC_RDWR;
    hid_t file_id = H5Fopen(filename, flags, plist_id);
    // If open fails, H5Fopen returns negative. We check this by H5_ERR.
    // However, some code might check for return value to determine existence.
    // For now we assume strict error handling.
    H5_ERR(file_id);

    ret = H5Pclose(plist_id);
    H5_ERR(ret);

    return file_id;
}

void elph_h5_close_file(hid_t file_id)
{
    herr_t ret = H5Fclose(file_id);
    H5_ERR(ret);
}

void elph_h5_def_var(hid_t file_id, const char* name, hid_t type_id, int rank,
                     const ND_int* dims, const char** dim_names,
                     const size_t* chunk_sizes, hid_t* dset_id)
{
    // 1. Create Dataspace
    hsize_t* h5_dims = malloc(rank * sizeof(hsize_t));
    CHECK_ALLOC(h5_dims);
    for (int i = 0; i < rank; ++i)
    {
        h5_dims[i] = (hsize_t)dims[i];
    }

    hid_t space_id = H5Screate_simple(rank, h5_dims, NULL);
    H5_ERR(space_id);

    // 2. Create Property List (Chunking)
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    H5_ERR(plist_id);

    if (chunk_sizes)
    {
        hsize_t* h5_chunk = malloc(rank * sizeof(hsize_t));
        CHECK_ALLOC(h5_chunk);
        for (int i = 0; i < rank; ++i)
        {
            h5_chunk[i] = (hsize_t)chunk_sizes[i];
        }
        herr_t ret = H5Pset_chunk(plist_id, rank, h5_chunk);
        H5_ERR(ret);
        free(h5_chunk);
    }

    // 3. Create Dataset
    *dset_id = H5Dcreate2(file_id, name, type_id, space_id, H5P_DEFAULT,
                          plist_id, H5P_DEFAULT);
    H5_ERR(*dset_id);

    // 4. Cleanup
    H5_ERR(H5Sclose(space_id));
    H5_ERR(H5Pclose(plist_id));
    free(h5_dims);

    // Write dimension names as attribute if provided
    if (dim_names)
    {
        hid_t attr_type = H5Tcopy(H5T_C_S1);
        H5_ERR(attr_type);
        H5_ERR(H5Tset_size(attr_type, H5T_VARIABLE));

        hsize_t dims_attr = rank;
        hid_t attr_space = H5Screate_simple(1, &dims_attr, NULL);
        H5_ERR(attr_space);

        hid_t attr_id = H5Acreate2(*dset_id, "DIM_NAMES", attr_type, attr_space,
                                   H5P_DEFAULT, H5P_DEFAULT);
        H5_ERR(attr_id);

        H5_ERR(H5Awrite(attr_id, attr_type, dim_names));

        H5_ERR(H5Aclose(attr_id));
        H5_ERR(H5Sclose(attr_space));
        H5_ERR(H5Tclose(attr_type));
    }
}

void elph_h5_write_var(hid_t dset_id, hid_t mem_type_id, const void* data)
{
    // Write the entire dataset
    herr_t ret =
        H5Dwrite(dset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5_ERR(ret);
}

void elph_h5_write_vara(hid_t dset_id, hid_t mem_type_id, const hsize_t* start,
                        const hsize_t* count, const void* data)
{
    hid_t filespace = H5Dget_space(dset_id);
    H5_ERR(filespace);

    int rank = H5Sget_simple_extent_ndims(filespace);

    herr_t ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL,
                                     count, NULL);
    H5_ERR(ret);

    hid_t memspace = H5Screate_simple(rank, count, NULL);
    H5_ERR(memspace);

    // Create property list for collective IO if needed
    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5_ERR(plist_id);
    // ret = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    // Defaulting to INDEPENDENT for now unless we know it's collective.
    // Code analysis showed NC_INDEPENDENT in many places, but also
    // NC_COLLECTIVE. Ideally we should pass this as a flag. For safety with
    // parallel HDF5, collective is often preferred if all ranks participate.
    // Let's stick to H5FD_MPIO_INDEPENDENT default or allow it to be passed?
    // The user code uses NC_INDEPENDENT and NC_COLLECTIVE explicitly.
    // To keep the wrapper simple, I'll use H5P_DEFAULT (Independent).
    // If collective is needed, we might need an extra argument.

    ret = H5Dwrite(dset_id, mem_type_id, memspace, filespace, plist_id, data);
    H5_ERR(ret);

    H5_ERR(H5Sclose(filespace));
    H5_ERR(H5Sclose(memspace));
    H5_ERR(H5Pclose(plist_id));
}

void elph_h5_read_var(hid_t dset_id, hid_t mem_type_id, void* data)
{
    herr_t ret =
        H5Dread(dset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5_ERR(ret);
}

void elph_h5_read_vara(hid_t dset_id, hid_t mem_type_id, const hsize_t* start,
                       const hsize_t* count, void* data)
{
    hid_t filespace = H5Dget_space(dset_id);
    H5_ERR(filespace);

    int rank = H5Sget_simple_extent_ndims(filespace);

    herr_t ret = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL,
                                     count, NULL);
    H5_ERR(ret);

    hid_t memspace = H5Screate_simple(rank, count, NULL);
    H5_ERR(memspace);

    herr_t ret_read =
        H5Dread(dset_id, mem_type_id, memspace, filespace, H5P_DEFAULT, data);
    H5_ERR(ret_read);

    H5_ERR(H5Sclose(filespace));
    H5_ERR(H5Sclose(memspace));
}

void elph_h5_close_var(hid_t dset_id)
{
    herr_t ret = H5Dclose(dset_id);
    H5_ERR(ret);
}

hid_t elph_h5_open_var(hid_t file_id, const char* name)
{
    hid_t dset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
    H5_ERR(dset_id);
    return dset_id;
}

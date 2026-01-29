#include <hdf5.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "elphC.h"
#include "elph_hdf5.h"
#include "io.h"

// write basic data related to the dft and ph
//
void write_basic_data(const hid_t file_id, struct Lattice* lattice,
                      struct Phonon* phonon, const char* kernel_str,
                      const char* convention_str)
{
    // only single cpu must call this which implies that the file must be opened
    // by only single cpu
    //
    //  writes kpoints in (crystal units)
    //  qpoints in (crystal units)
    //  kernel type
    //  convention
    //  start and end band
    //  nsym, symmetries matrices for phonons
    hid_t dset_id;

    // kpoints
    elph_h5_def_var(file_id, "kpoints", ELPH_H5_IO_FLOAT, 2,
                    (ND_int[]){lattice->nkpts_BZ, 3},
                    (const char*[]){"nk", "pol"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, lattice->kpt_fullBZ_crys);
    elph_h5_close_var(dset_id);

    // qpoints
    elph_h5_def_var(file_id, "qpoints", ELPH_H5_IO_FLOAT, 2,
                    (ND_int[]){phonon->nq_BZ, 3}, (const char*[]){"nq", "pol"},
                    NULL, &dset_id);
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, phonon->qpts_BZ);
    elph_h5_close_var(dset_id);

    // qpoints in iBZ
    elph_h5_def_var(file_id, "qpoints_iBZ", ELPH_H5_IO_FLOAT, 2,
                    (ND_int[]){phonon->nq_iBZ, 3},
                    (const char*[]){"nq_iBZ", "pol"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, phonon->qpts_iBZ);
    elph_h5_close_var(dset_id);

    // write qmap for qpoints (for each qpt in BZ, it gives corresponding iBZ
    // and symm used)
    elph_h5_def_var(file_id, "qmap", H5T_NATIVE_INT, 2,
                    (ND_int[]){phonon->nq_BZ, 2},
                    (const char*[]){"nq", "dim_two"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_INT, phonon->qmap);
    elph_h5_close_var(dset_id);

    // write kmap for kpoints (for each kpt in BZ, it gives corresponding iBZ
    // and symm used) This is internally available in yambo, it can be used to
    // cross check if rotation is done in same way

    elph_h5_def_var(file_id, "kmap", H5T_NATIVE_INT, 2,
                    (ND_int[]){lattice->nkpts_BZ, 2},
                    (const char*[]){"nk", "dim_two"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_INT, lattice->kmap);
    elph_h5_close_var(dset_id);

    // write start and end band indices
    int bands_tmp[2] = {lattice->start_band, lattice->end_band};
    elph_h5_def_var(file_id, "bands", H5T_NATIVE_INT, 1, (ND_int[]){2},
                    (const char*[]){"two_scalars"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_INT, bands_tmp);
    elph_h5_close_var(dset_id);

    int nph_sym = phonon->nph_sym;
    // write number of phonon symmetries
    elph_h5_def_var(file_id, "number_of_phonon_symmetries", H5T_NATIVE_INT, 1,
                    (ND_int[]){1}, (const char*[]){"scalar"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_INT, &nph_sym);
    elph_h5_close_var(dset_id);

    int time_rev_present = 0;
    ELPH_float* symm_mats = malloc(sizeof(ELPH_float) * 9 * nph_sym);
    CHECK_ALLOC(symm_mats);

    ELPH_float* tau_vecs = malloc(sizeof(ELPH_float) * 3 * nph_sym);
    CHECK_ALLOC(tau_vecs);

    for (int isym = 0; isym < nph_sym; ++isym)
    {
        memcpy(symm_mats + isym * 9, phonon->ph_syms[isym].Rmat,
               sizeof(ELPH_float) * 9);
        memcpy(tau_vecs + isym * 3, phonon->ph_syms[isym].tau,
               sizeof(ELPH_float) * 3);
        if (phonon->ph_syms[isym].time_rev)
        {
            time_rev_present = 1;
        }
    }

    // write information about the time reversal symmetry of the phonon
    elph_h5_def_var(file_id, "time_reversal_phonon", H5T_NATIVE_INT, 1,
                    (ND_int[]){1}, (const char*[]){"scalar"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_INT, &time_rev_present);
    elph_h5_close_var(dset_id);

    // write info about the symmetry_matrices in cart units
    elph_h5_def_var(file_id, "symmetry_matrices", ELPH_H5_IO_FLOAT, 3,
                    (ND_int[]){nph_sym, 3, 3},
                    (const char*[]){"nsym_ph", "pol", "pol"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, symm_mats);
    elph_h5_close_var(dset_id);

    // write info about the fractional_translation in cart units
    elph_h5_def_var(file_id, "fractional_translation", ELPH_H5_IO_FLOAT, 2,
                    (ND_int[]){nph_sym, 3}, (const char*[]){"nsym_ph", "pol"},
                    NULL, &dset_id);
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, tau_vecs);
    elph_h5_close_var(dset_id);

    // write what type of kernel is used in the calculation
    size_t str_size_tmp = strlen(kernel_str) + 1;
    // for strings we also write the null terminator to the netcdf variable
    elph_h5_def_var(file_id, "kernel", H5T_NATIVE_CHAR, 1,
                    (ND_int[]){str_size_tmp},
                    (const char*[]){"kernel_str_size"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_CHAR, kernel_str);
    elph_h5_close_var(dset_id);

    // write the information about the convention used in the code
    str_size_tmp = strlen(convention_str) + 1;
    elph_h5_def_var(file_id, "convention", H5T_NATIVE_CHAR, 1,
                    (ND_int[]){str_size_tmp},
                    (const char*[]){"convention_str_size"}, NULL, &dset_id);
    elph_h5_write_var(dset_id, H5T_NATIVE_CHAR, convention_str);
    elph_h5_close_var(dset_id);

    free(symm_mats);
    free(tau_vecs);

    //
    ELPH_float* tmp_buf = calloc(27 * lattice->natom, sizeof(*tmp_buf));
    CHECK_ALLOC(tmp_buf);
    //
    for (ND_int i = 0; i < (27 * lattice->natom); ++i)
    {
        tmp_buf[i] = 0.0;
    }

    // write dielectric tensor
    elph_h5_def_var(file_id, "epsilon", ELPH_H5_IO_FLOAT, 2, (ND_int[]){3, 3},
                    (const char*[]){"pol", "pol"}, NULL, &dset_id);
    //
    ELPH_float* write_tensor_buf = tmp_buf;
    if (phonon->epsilon)
    {
        write_tensor_buf = phonon->epsilon;
    }
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, write_tensor_buf);
    elph_h5_close_var(dset_id);

    // write born charges
    elph_h5_def_var(file_id, "Born_charges", ELPH_H5_IO_FLOAT, 3,
                    (ND_int[]){lattice->natom, 3, 3},
                    (const char*[]){"natoms", "pol", "pol"}, NULL, &dset_id);
    //
    write_tensor_buf = tmp_buf;
    if (phonon->Zborn)
    {
        write_tensor_buf = phonon->Zborn;
    }
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, write_tensor_buf);
    elph_h5_close_var(dset_id);

    // write Quadrupole tensor
    elph_h5_def_var(file_id, "Quadrupole_tensor", ELPH_H5_IO_FLOAT, 4,
                    (ND_int[]){lattice->natom, 3, 3, 3},
                    (const char*[]){"natoms", "pol", "pol", "pol"}, NULL,
                    &dset_id);
    //
    write_tensor_buf = tmp_buf;
    if (phonon->Qpole)
    {
        write_tensor_buf = phonon->Qpole;
    }
    elph_h5_write_var(dset_id, ELPH_H5_IO_FLOAT, write_tensor_buf);
    elph_h5_close_var(dset_id);

    free(tmp_buf);
}

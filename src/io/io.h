#pragma once

#include <mpi.h>
#include <stddef.h>

#include "common/dtypes.h"
#include "elphC.h"
#include "elph_hdf5.h"

#define H5_DEFAULT_CHUNK_KB 2048
// default chunking for large hdf5 varaibles (in Kilobytes)

void read_and_alloc_save_data(char* SAVEdir, const struct ELPH_MPI_Comms* Comm,
                              ND_int start_band, ND_int end_band,
                              struct WFC** wfcs, char* ph_save_dir,
                              struct Lattice* lattice, struct Pseudo* pseudo,
                              struct Phonon* phonon,
                              enum ELPH_dft_code dft_code);

void free_phonon_data(struct Phonon* phonon);

void free_save_data(struct WFC* wfcs, struct Lattice* lattice,
                    struct Pseudo* pseudo, struct Phonon* phonon);

void write_basic_data(const hid_t file_id, struct Lattice* lattice,
                      struct Phonon* phonon, const char* kernel_str,
                      const char* convention_str);

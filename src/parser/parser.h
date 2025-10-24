#pragma once

#include <mpi.h>

#include "common/dtypes.h"
#include "elphC.h"

// electron-phonon input file parser
void init_elph_usr_input(struct elph_usr_input** input);
void free_elph_usr_input(struct elph_usr_input* input);
void read_elph_input_file(const char* input_file,
                          struct elph_usr_input** input_data,
                          MPI_Comm MPI_world_comm);

// interpolation input io
void init_interpolation_usr_input(struct interpolation_usr_input** input);
void free_interpolation_usr_input(struct interpolation_usr_input* input);
void read_interpolation_input_file(const char* input_file,
                                   struct interpolation_usr_input** input_data,
                                   MPI_Comm MPI_world_comm);

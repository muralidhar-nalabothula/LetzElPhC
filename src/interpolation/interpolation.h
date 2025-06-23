#pragma once

#include <mpi.h>

#include "common/dtypes.h"
#include "elphC.h"

void interpolation_driver(const char* ph_save, const char* ph_save_interpolated,
                          enum ELPH_dft_code dft_code, const ND_int* qgrid_new,
                          MPI_Comm comm_world);

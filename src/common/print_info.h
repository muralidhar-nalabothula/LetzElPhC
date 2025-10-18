/**
 * @file
 * @brief Information printing function declarations
 */

#pragma once
#include <stdio.h>

#include "dtypes.h"
#include "elphC.h"

/**
 * @brief Prints ELPH ASCII art logo
 *
 * @param mpi_rank MPI rank of calling process
 * @param output Output file stream
 */
void print_ELPH_logo(int mpi_rank, FILE* output);

/**
 * @brief Prints formatted message from rank 0 only
 *
 * @param mpi_rank MPI rank of calling process
 * @param fmt Printf-style format string
 * @param ... Variable arguments for format string
 */
void print_info_msg(int mpi_rank, const char* fmt, ...);

/**
 * @brief Prints input configuration and parallelization information
 *
 * @param save_dir Path to DFT save directory
 * @param ph_save_dir Path to phonon save directory
 * @param kernel Screening kernel description
 * @param kminusq Scattering convention flag
 * @param dft_code DFT code identifier
 * @param Comm MPI communicator structure
 */
void print_input_info(const char* save_dir, const char* ph_save_dir,
                      const char* kernel, const bool kminusq,
                      const enum ELPH_dft_code dft_code,
                      const struct ELPH_MPI_Comms* Comm);

/**
 * @brief Prints lattice and electronic structure information
 *
 * @param Comm MPI communicator structure
 * @param lattice Lattice structure
 */
void print_lattice_info(const struct ELPH_MPI_Comms* Comm,
                        const struct Lattice* lattice);

/**
 * @brief Prints phonon structure information
 *
 * @param Comm MPI communicator structure
 * @param phonon Phonon structure
 */
void print_phonon_info(const struct ELPH_MPI_Comms* Comm,
                       const struct Phonon* phonon);

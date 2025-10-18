/**
 * @file
 * @brief Initialization function declarations for data structures
 *
 * Declares initialization routines for major data structure types used
 * throughout the ELPH library.
 */

#pragma once
#include "dtypes.h"
#include "elphC.h"

/**
 * @brief Initializes a Lattice structure to safe default values
 *
 * @param lattice Pointer to Lattice structure to initialize
 */
void init_lattice_type(struct Lattice* lattice);

/**
 * @brief Initializes a Phonon structure to safe default values
 *
 * @param phonon Pointer to Phonon structure to initialize
 */
void init_phonon_type(struct Phonon* phonon);

/**
 * @brief Initializes a Vloc_table structure to safe default values
 *
 * @param vloc_tab Pointer to Vloc_table structure to initialize
 */
void init_Vloc_table_type(struct Vloc_table* vloc_tab);

/**
 * @brief Initializes a local_pseudo structure to safe default values
 *
 * @param lpseudo Pointer to local_pseudo structure to initialize
 */
void init_local_pseudo_type(struct local_pseudo* lpseudo);

/**
 * @brief Initializes a Pseudo structure to safe default values
 *
 * @param pseudo Pointer to Pseudo structure to initialize
 */
void init_Pseudo_type(struct Pseudo* pseudo);

/**
 * @file
 * @brief Initialization functions for data structure types
 *
 * Provides initialization routines that zero-initialize structures and
 * explicitly set all pointer members to NULL for safety and portability.
 *
 * @note All pointer members are explicitly set to NULL after memset, as the
 *       C standard does not guarantee that all-zero-bits represents a NULL
 *       pointer. While most compilers optimize away redundant NULL assignments,
 *       this approach ensures strict standard compliance and portability.
 */

#include "init_dtypes.h"

#include <mpi.h>
#include <string.h>

#include "dtypes.h"
#include "elphC.h"

/**
 * @brief Initializes a Lattice structure to safe default values
 *
 * Zeros all members and explicitly sets all pointer members to NULL.
 *
 * @param lattice Pointer to Lattice structure to initialize
 */
void init_lattice_type(struct Lattice* lattice)
{
    memset(lattice, 0, sizeof(*lattice));
    lattice->atomic_pos = NULL;
    lattice->kpt_iredBZ = NULL;
    lattice->kpt_fullBZ = NULL;
    lattice->kpt_fullBZ_crys = NULL;
    lattice->kmap = NULL;
    lattice->syms = NULL;
    lattice->atom_type = NULL;
    return;
}

/**
 * @brief Initializes a Phonon structure to safe default values
 *
 * Zeros all members and explicitly sets all pointer members to NULL.
 *
 * @param phonon Pointer to Phonon structure to initialize
 */
void init_phonon_type(struct Phonon* phonon)
{
    memset(phonon, 0, sizeof(*phonon));
    phonon->qpts_iBZ = NULL;
    phonon->qpts_BZ = NULL;
    phonon->ph_syms = NULL;
    phonon->qmap = NULL;
    phonon->nqstar = NULL;
    phonon->epsilon = NULL;
    phonon->Zborn = NULL;
    phonon->Qpole = NULL;
    return;
}

/**
 * @brief Initializes a Vloc_table structure to safe default values
 *
 * Zeros all members and explicitly sets all pointer members to NULL.
 *
 * @param vloc_tab Pointer to Vloc_table structure to initialize
 */
void init_Vloc_table_type(struct Vloc_table* vloc_tab)
{
    memset(vloc_tab, 0, sizeof(*vloc_tab));
    vloc_tab->g_co = NULL;
    vloc_tab->vlocg = NULL;
    vloc_tab->vploc_co = NULL;
    return;
}

/**
 * @brief Initializes a local_pseudo structure to safe default values
 *
 * Zeros all members and explicitly sets all pointer members to NULL.
 *
 * @param lpseudo Pointer to local_pseudo structure to initialize
 */
void init_local_pseudo_type(struct local_pseudo* lpseudo)
{
    memset(lpseudo, 0, sizeof(*lpseudo));
    lpseudo->Vloc_atomic = NULL;
    lpseudo->r_grid = NULL;
    lpseudo->rab_grid = NULL;
    return;
}

/**
 * @brief Initializes a Pseudo structure to safe default values
 *
 * Zeros all members, explicitly sets all pointer members to NULL,
 * and initializes the embedded vloc_table member.
 *
 * @param pseudo Pointer to Pseudo structure to initialize
 */
void init_Pseudo_type(struct Pseudo* pseudo)
{
    memset(pseudo, 0, sizeof(*pseudo));
    pseudo->loc_pseudo = NULL;
    pseudo->PP_table = NULL;
    pseudo->Fsign = NULL;
    pseudo->fCoeff = NULL;
    init_Vloc_table_type(pseudo->vloc_table);
    return;
}

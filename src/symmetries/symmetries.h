#pragma once
#include <stdbool.h>

#include "common/dtypes.h"
#include "elphC.h"

ND_int bz_expand(const ND_int Nibz, const ND_int Nsym,
                 const ELPH_float* ibz_kpts, const struct symmetry* symms,
                 const ELPH_float* lat_vec, ELPH_float* kpoints, ND_int* kstar,
                 int* kmap);

void electronic_reps(const struct WFC* wfcs, const struct Lattice* lattice,
                     const ELPH_float* Rsym_mat, const ELPH_float* tauR,
                     const bool tim_revR, const ND_int ikBZ,
                     ELPH_cmplx* Dkmn_rep, const struct ELPH_MPI_Comms* Comm);

void elph_q_rotate(const ELPH_cmplx* Dmats_l, const ELPH_cmplx* elph_mat_q,
                   const ELPH_cmplx* Dmats_r, const struct Lattice* lattice,
                   const bool tim_rev, ELPH_cmplx* elph_mat_Sq);

void rotate_eig_vecs(struct symmetry* sym, const struct Lattice* lattice,
                     const ELPH_float* qpt, const ELPH_cmplx* eig_q,
                     ELPH_cmplx* eig_Sq);

void rotate_dvscf(const ELPH_cmplx* dvscf_in, struct symmetry* sym,
                  const struct Lattice* lattice, const bool composite_form,
                  ELPH_cmplx* restrict dvscf_out, MPI_Comm commK);

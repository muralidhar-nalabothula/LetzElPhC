#pragma once
#include "elphC.h"

void mass_normalize_pol_vecs(const ELPH_float* atomic_masses,
                             const ND_int nsets, const ND_int natoms,
                             const ELPH_float power, ELPH_cmplx* pol_vecs);

void pol_vecs_to_dyn(const ELPH_float* omega, const ND_int natom,
                     const ELPH_float* atomic_masses, ELPH_cmplx* pol_vecs);

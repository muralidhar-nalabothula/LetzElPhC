#pragma once
#include "elphC.h"

void mass_normalize_pol_vecs(const ELPH_float* atomic_masses,
                             const ND_int nsets, const ND_int natoms,
                             const ELPH_float power, ELPH_cmplx* pol_vecs);

#pragma once
#include "../elphC.h"

void bare_lr_vertex(const ELPH_float* qpt, const ELPH_float* gvec,
                    const ND_int npw_loc, const ELPH_float* Zatoms,
                    const ND_int natom, const ELPH_float* atom_pos,
                    const char diminsion, const ELPH_float volume,
                    const ELPH_float zlat, ELPH_cmplx* elph_lr_out);

void frohlich_lr_vertex(const ELPH_float* qpt, const ELPH_float* gvec,
                        const ND_int npw_loc, const ELPH_float* epslion,
                        const ELPH_float* Zeu, const ELPH_float* Qpole,
                        const ND_int natom, const ELPH_float* atom_pos,
                        const char diminsion, const ELPH_float volume,
                        const ELPH_float zlat, ELPH_cmplx* elph_lr_out);

//
// Multiplies dvscf (natom, 3, ....) with structural factor
// i.e dvscf = dvscf(a,...) e^{-+i q.tau_a}
// where tau_a is atomic positon position
#include <complex.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/dtypes.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"
#include "dvloc.h"
#include "elphC.h"

void mul_dvscf_struct_fac(const ELPH_float* qpt_cart, struct Lattice* lattice,
                          const ND_int nsets, const ND_int sign,
                          ELPH_cmplx* dVscf)
{
    // add or remove structural factor to dvscf
    // exp(i q.tau) if sign < 0 (i,.e remove phase)
    // else : exp(-i q.tau) } (i.e add phase)
    // THis is inplace, operation.
    // qpt is in cart units with 2*pi included.
    // dvscf (natom, 3, nsets)
    //
    ELPH_start_clock("mul_dvscf_struct_fac");
    //
    ELPH_float factor = 1.0;
    if (sign < 0)
    {
        factor = -factor;
    }
    for (ND_int ia = 0; ia < lattice->natom; ++ia)
    {
        const ELPH_float* atom_pos_ia = lattice->atomic_pos + 3 * ia;
        //
        ELPH_cmplx qdot_tau =
            cexp(-I * factor * dot3_macro(qpt_cart, atom_pos_ia));
        //
        ELPH_cmplx* tmp_ptr = dVscf + ia * 3 * nsets;
        //
        ELPH_OMP_PAR_FOR_SIMD
        for (ND_int i = 0; i < 3 * nsets; ++i)
        {
            tmp_ptr[i] *= qdot_tau;
        }
    }
    ELPH_stop_clock("mul_dvscf_struct_fac");
}

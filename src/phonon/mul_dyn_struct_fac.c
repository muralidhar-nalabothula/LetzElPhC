//
// Multiplies dyn (natom, 3, natom, 3) with structural factor
// i.e dyn = dyn(a,x b,y) e^{+-i q.tau_a} * e^{-+i q.tau_b}
// where tau_a is atomic positon position
#include <complex.h>

#include "common/ELPH_timers.h"
#include "common/dtypes.h"
#include "common/numerical_func.h"
#include "elphC.h"
#include "phonon.h"

void mul_dyn_struct_fac(const ELPH_float* qpt_cart, struct Lattice* lattice,
                        const ND_int sign, ELPH_cmplx* dyn)
{
    // add or remove structural factor to dyn (ia,3,ib,3)
    // exp(i q.tau_b)*exp(-i q.tau_a) is sign < 0 (i,.e remove phase)
    // else : exp(-i q.tau_b)*exp(i q.tau_a) } (i.e add phase)
    // THis is inplace, operation.
    // qpt is in cart units with 2*pi included.
    //
    ELPH_start_clock("mul_dyn_struct_fac");
    //
    ELPH_float factor = 1.0;
    if (sign < 0)
    {
        factor = -factor;
    }

    ND_int nmodes = 3 * lattice->natom;
    for (ND_int ia = 0; ia < lattice->natom; ++ia)
    {
        const ELPH_float* atom_pos_ia = lattice->atomic_pos + 3 * ia;
        //
        const ELPH_cmplx qdot_tau_a =
            cexp(I * factor * dot3_macro(qpt_cart, atom_pos_ia));
        //
        ELPH_cmplx* tmp_ptr_ia = dyn + ia * 3 * nmodes;
        //
        for (ND_int ib = 0; ib < lattice->natom; ++ib)
        {
            const ELPH_float* atom_pos_ib = lattice->atomic_pos + 3 * ib;
            //
            const ELPH_cmplx qdot_tau_b =
                cexp(-I * factor * dot3_macro(qpt_cart, atom_pos_ib));
            //
            for (ND_int xa = 0; xa < 3; ++xa)
            {
                for (ND_int xb = 0; xb < 3; ++xb)
                {
                    tmp_ptr_ia[xb + ib * 3 + xa * nmodes] *=
                        (qdot_tau_a * qdot_tau_b);
                }
            }
        }
    }
    ELPH_stop_clock("mul_dyn_struct_fac");
}

// This file contains functions to compute long range part
// of the bare electron-phonon matrix elements.
//
//

#include <complex.h>
#include <math.h>
#include <string.h>

#include "../common/ELPH_timers.h"
#include "../common/constants.h"
#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../elphC.h"
#include "elph_lr_kernels.h"

void bare_lr_vertex(const ELPH_float* qpt, const ELPH_float* gvec,
                    const ND_int npw_loc, const ELPH_float* Zatoms,
                    const ND_int natom, const ELPH_float* atom_pos,
                    const char diminsion, const ELPH_float volume,
                    const ELPH_float zlat, ELPH_cmplx* elph_lr_out)
{
    // gvec in cartisian coordinates  ( no 2*pi)
    // qpt in cart units             ( no 2*pi)
    // npw_loc, number of gvecs in this cpu
    // atomic pos in cart units
    // elph_lr_out (natom,3)
    // zlat in the dimension along the out ot plane direction.
    // only used in the 2D case

    // first zero out the elph_lr_out buffer

    // note that in ca
    //
    // NOTE: Donot forget to allreduce the result over plane waves.

    ELPH_start_clock("Bare lr part");

    for (ND_int i = 0; i < (3 * natom); ++i)
    {
        elph_lr_out[i] = 0.0;
    }

    ELPH_cmplx factor = 4.0 * ELPH_PI * I * ELPH_e2 / volume;  // prefactor

    for (ND_int ia = 0; ia < natom; ++ia)
    {
        const ELPH_float Z_k = Zatoms[ia];
        const ELPH_float* tau_k = atom_pos + 3 * ia;

        ELPH_cmplx* out_tmp_buf = elph_lr_out + ia * 3;

        // No Openmp here. not thread safe !!
        for (ND_int ig = 0; ig < npw_loc; ++ig)
        {
            const ELPH_float* gtmp = gvec + 3 * ig;

            ELPH_float qplusG[3];
            for (int i = 0; i < 3; ++i)
            {
                qplusG[i] = 2.0 * ELPH_PI * (qpt[i] + gtmp[i]);
            }

            ELPH_float qplusG_abs2 = dot3_macro(qplusG, qplusG);
            if (sqrt(qplusG_abs2) < ELPH_EPS)
            {
                // skip q+G = 0 term
                continue;
            }

            ELPH_float cutoff_fac = 1;

            if (diminsion == '2')
            {
                ELPH_float qGp =
                    zlat * sqrt(qplusG[0] * qplusG[0] + qplusG[1] * qplusG[1]);
                cutoff_fac -= exp(-qGp * 0.5) * cos(qplusG[2] * zlat * 0.5);
            }

            // struture factor and multiply with cutoff factor
            ELPH_cmplx qdot_tau = cexp(-I * (dot3_macro(qplusG, tau_k)));

            // compute and multiply with a decay factor
            qdot_tau *= exp(-qplusG_abs2 * 0.25);

            for (int i = 0; i < 3; ++i)
            {
                out_tmp_buf[i] = out_tmp_buf[i] + qdot_tau * cutoff_fac *
                                                      qplusG[i] / qplusG_abs2;
            }
        }
        // multiply with prefactor
        for (int i = 0; i < 3; ++i)
        {
            out_tmp_buf[i] *= (factor * Z_k);
        }
    }

    ELPH_stop_clock("Bare lr part");
    return;
}

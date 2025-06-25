// compute non analytical term for phonon
//
// 1) Only the dipole-dipole term is include//
//

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"
#include "elphC.h"
#include "phonon.h"

void add_ph_dyn_long_range(const ELPH_float* qpt, struct Lattice* lattice,
                           struct Phonon* phonon, const ND_int* Ggrid,
                           const ND_int sign, const ELPH_float* atomic_masses,
                           ELPH_cmplx* dyn_mat)
{
    // adds or subtracts non-analytical term to the dynamical matrix.
    // sign < 0 : subtract else add
    // qpt in crystal coordinater
    //
    if (!phonon->epsilon || !phonon->Zborn)
    {
        return;
    }
    //
    ELPH_start_clock("dyn lr part");

    if (lattice->dimension == '2' && fabs(qpt[2]) > ELPH_EPS)
    {
        error_msg(
            "In 2D, only qz == 0 points are accepted when interpolating.");
    }

    ELPH_cmplx factor = 4.0 * ELPH_PI * ELPH_e2 / lattice->volume;  // prefactor
    if (sign < 0)
    {
        factor = -factor;
    }
    //
    ELPH_float eps_alpha[9];

    memcpy(eps_alpha, phonon->epsilon, sizeof(eps_alpha));

    ELPH_float zlat = lattice->alat_vec[8];
    //
    if (lattice->dimension == '2')
    {
        factor = factor * zlat / 2.0;
        // compute alpha = c/2 * (eps-1)
        for (int i = 0; i < 9; ++i)
        {
            eps_alpha[i] = 0.5 * zlat * (eps_alpha[i] - 1);
        }
    }

    ND_int Gvec_size = Ggrid[0] * Ggrid[1] * Ggrid[2];

    ND_int natom = lattice->natom;

    ELPH_cmplx* Zdotq_tau =
        malloc(2 * Gvec_size * 3 * natom * sizeof(*Zdotq_tau));
    CHECK_ALLOC(Zdotq_tau);

    //
    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ig = 0; ig < Gvec_size; ++ig)
    {
        ELPH_cmplx* out_tmp1 = Zdotq_tau + ig * 3 * natom;
        ELPH_cmplx* out_tmp2 =
            Zdotq_tau + ig * 3 * natom + Gvec_size * 3 * natom;
        const ND_int Gx = ig / (Ggrid[1] * Ggrid[2]);
        const ND_int Gy = (ig % (Ggrid[1] * Ggrid[2])) / Ggrid[2];
        const ND_int Gz = (ig % (Ggrid[1] * Ggrid[2])) % Ggrid[2];
        //
        ELPH_float qplusG[3], tmp_buf[3];
        tmp_buf[0] = (qpt[0] + get_miller_idx(Gx, Ggrid[0]));
        tmp_buf[1] = (qpt[1] + get_miller_idx(Gy, Ggrid[1]));
        tmp_buf[2] = (qpt[2] + get_miller_idx(Gz, Ggrid[2]));
        // COnvert to cart units (2*pi) is included
        MatVec3f(lattice->blat_vec, tmp_buf, false, qplusG);
        //
        ELPH_float qplusG_norm = sqrt(dot3_macro(qplusG, qplusG));
        if (qplusG_norm < ELPH_EPS)
        {
            continue;
        }
        //
        MatVec3f(eps_alpha, qplusG, false, tmp_buf);
        ELPH_float q_eps_q = dot3_macro(tmp_buf, qplusG);
        if (lattice->dimension == '2')
        {
            q_eps_q += qplusG_norm;
        }
        ELPH_float decay_fac = exp(-qplusG_norm * qplusG_norm * 0.125);
        //
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            const ELPH_float* Z_k =
                phonon->Zborn ? (phonon->Zborn + 9 * ia) : NULL;
            const ELPH_float* Q_k =
                phonon->Qpole ? (phonon->Qpole + 27 * ia) : NULL;
            const ELPH_float* tau_k = lattice->atomic_pos + 3 * ia;

            MatVec3f(Z_k, qplusG, true, tmp_buf);
            ELPH_cmplx* out_tmp_buf = out_tmp1 + 3 * ia;
            ELPH_cmplx* out_tmp_buf2 = out_tmp2 + 3 * ia;

            // compute (q+G)_x Q_xyz * (q+G)_y
            ELPH_float Qpole_buf[3] = {0.0, 0.0, 0.0};
            if (Q_k)
            {
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        for (int k = 0; k < 3; ++k)
                        {
                            Qpole_buf[k] =
                                Qpole_buf[k] +
                                qplusG[i] * qplusG[j] * Q_k[k + 3 * j + 9 * i];
                        }
                    }
                }
            }
            // e^{-iq.tau}
            ELPH_cmplx qdot_tau = cexp(-I * (dot3_macro(qplusG, tau_k)));
            qdot_tau /= q_eps_q;
            qdot_tau *= decay_fac;
            qdot_tau /= sqrt(atomic_masses[ia]);
            //
            for (ND_int i = 0; i < 3; ++i)
            {
                out_tmp_buf[i] =
                    (tmp_buf[i] - I * 0.5 * Qpole_buf[i]) * qdot_tau;
                out_tmp_buf2[i] = out_tmp_buf[i] * q_eps_q;
            }
        }
    }
    //
    ND_int nmodes = 3 * natom;
    matmul_cmplx('C', 'N', Zdotq_tau + Gvec_size * 3 * natom, Zdotq_tau,
                 dyn_mat, factor, 1.0, nmodes, nmodes, nmodes, nmodes, nmodes,
                 Gvec_size);
    //
    free(Zdotq_tau);
    //
    //
    ELPH_stop_clock("dyn lr part");
    //
    return;
}

// This file contains functions to compute long range part
// of the electron-phonon matrix elements.
//
//

#include <complex.h>
#include <math.h>
#include <string.h>

#include "../common/constants.h"
#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/numerical_func.h"
#include "../elphC.h"

static void frohlich_dip2D_kernel(const ELPH_float* qplusG,
                                  const ELPH_float* Zborn_k,
                                  const ELPH_float* alpha,
                                  const ELPH_float* tau_k, ELPH_cmplx* out_buf);

static void frohlich_dip3D_kernel(const ELPH_float* qplusG,
                                  const ELPH_float* Zborn_k,
                                  const ELPH_float* epslion,
                                  const ELPH_float* tau_k, ELPH_cmplx* out_buf);

void frohlich_lr_vertex(const ELPH_float* qpt, const ELPH_float* gvec,
                        const ND_int npw_loc, const ELPH_float* epslion,
                        const ELPH_float* Zeu, const ELPH_float* Qpole,
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
    // Donot forget to allreduce the result.

    // do a basic check
    if (diminsion == '2' && fabs(qpt[2]) > ELPH_EPS)
    {
        error_msg(
            "In 2D, only qz == 0 points are accepted when interpolating.");
    }

    for (ND_int i = 0; i < (3 * natom); ++i)
    {
        elph_lr_out[i] = 0.0;
    }

    if (!epslion || !Zeu)
    {
        return;
    }

    ELPH_cmplx factor = 2.0 * I * ELPH_e2 / volume;  // prefactor
    ELPH_float eps_alpha[9];

    memcpy(eps_alpha, epslion, sizeof(eps_alpha));

    if (diminsion == '2')
    {
        factor = factor * zlat / 2.0;
        // compute alpha = c/2 * (eps-1)
        for (int i = 0; i < 9; ++i)
        {
            eps_alpha[i] = 0.5 * zlat * (eps_alpha[i] - 1);
        }
    }

    for (ND_int ia = 0; ia < natom; ++ia)
    {
        const ELPH_float* Z_k = Zeu + 9 * ia;
        const ELPH_float* tau_k = atom_pos + 3 * ia;

        ELPH_cmplx* out_tmp_buf = elph_lr_out + ia * 3;

        // No Openmp here. not thread safe !!
        for (ND_int ig = 0; ig < npw_loc; ++ig)
        {
            const ELPH_float* gtmp = gvec + 3 * ig;

            ELPH_float qplusG[3];
            for (int i = 0; i < 3; ++i)
            {
                qplusG[i] = qpt[i] + gtmp[i];
            }

            if (diminsion == '3')
            {
                frohlich_dip3D_kernel(qplusG, Z_k, eps_alpha, tau_k,
                                      out_tmp_buf);
            }
            else if (diminsion == '2')
            {
                frohlich_dip2D_kernel(qplusG, Z_k, eps_alpha, tau_k,
                                      out_tmp_buf);
            }
        }
        // multiply with prefactor
        for (int i = 0; i < 3; ++i)
        {
            out_tmp_buf[i] *= factor;
        }
    }

    return;
}

static void frohlich_dip3D_kernel(const ELPH_float* qplusG,
                                  const ELPH_float* Zborn_k,
                                  const ELPH_float* epslion,
                                  const ELPH_float* tau_k, ELPH_cmplx* out_buf)
{
    // G space dipole term for frohlich in 3D
    // Eq :4 of PhysRevLett.115.176401 (C. Verdi, F. Giustino)
    // We are in limit (q+G)-> 0 so
    // approximate the bracket <K+q,n|e^{(q+G)ir}|k,m> = \delta_{nm}

    // tau_k are atomic coordinates of atom k in cart
    // q+G in cart ,
    // Zborn_k, born charges for atom k

    if (sqrt(dot3_macro(qplusG, qplusG)) < ELPH_EPS)
    {
        // skip q+G = 0 term
        return;
    }

    // struture factor
    ELPH_cmplx qdot_tau =
        cexp(-2.0 * ELPH_PI * I * (dot3_macro(qplusG, tau_k)));

    ELPH_float tmp_buf[3];
    // compute (q+G).eps.(q+G)
    MatVec3f(epslion, qplusG, false, tmp_buf);
    ELPH_float q_eps_q = dot3_macro(tmp_buf, qplusG);

    qdot_tau /= q_eps_q;

    // compute (q+G).Z
    MatVec3f(Zborn_k, qplusG, true, tmp_buf);

    // compute and multiply with a decay factor
    qdot_tau *= exp(-q_eps_q * 0.25);

    for (int i = 0; i < 3; ++i)
    {
        out_buf[i] = out_buf[i] + qdot_tau * tmp_buf[i];
    }

    return;
}

static void frohlich_dip2D_kernel(const ELPH_float* qplusG,
                                  const ELPH_float* Zborn_k,
                                  const ELPH_float* alpha,
                                  const ELPH_float* tau_k, ELPH_cmplx* out_buf)
{
    // G space dipole term for frohlich in 2D
    // Note that we are at (q+G) -> 0 limit.
    // approximate the bracket <K+q,n|e^{(q+G)ir}|k,m> = \delta_{nm}

    // tau_k are atomic coordinates of atom k in cart
    // q+G in cart ,
    // Zborn_k, born charges for atom k
    // alpha = c/2 * (epslion-1)
    //
    // Eq. 5 of PHYSICAL REVIEW B 103, 075410 (2021)
    // Note 1 : We donot take into \hat{z} term as we  interpolation
    // only in the inplane q direction. This avoids clumsy sgn(z-zk) term
    // computation
    //

    // skip any term with Gz > 0
    if (fabs(qplusG[2]) > ELPH_EPS)
    {
        return;
    }

    ELPH_float qplusG_abs = sqrt(dot3_macro(qplusG, qplusG));

    if (qplusG_abs < ELPH_EPS)
    {
        // skip q+G = 0 term
        return;
    }

    // struture factor
    ELPH_cmplx qdot_tau =
        cexp(-2.0 * ELPH_PI * I * (dot3_macro(qplusG, tau_k)));

    ELPH_float tmp_buf[3];
    // compute (q+G).alpha.(q+G)
    MatVec3f(alpha, qplusG, false, tmp_buf);
    ELPH_float q_eps_q = dot3_macro(tmp_buf, qplusG);

    qdot_tau = qdot_tau / (q_eps_q + qplusG_abs);

    // compute (q+G).Z
    MatVec3f(Zborn_k, qplusG, true, tmp_buf);

    // compute and multiply with a decay factor
    qdot_tau *= exp(-q_eps_q * 0.25);

    for (int i = 0; i < 3; ++i)
    {
        out_buf[i] = out_buf[i] + qdot_tau * tmp_buf[i];
    }

    return;
}

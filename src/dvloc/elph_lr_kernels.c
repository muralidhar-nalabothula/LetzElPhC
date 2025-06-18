// This file contains functions to compute long range part
// of the electron-phonon matrix elements.
// monopole + dipole + Quadrupole are coded
//
//
#include <complex.h>
#include <math.h>
#include <string.h>

#include "common/constants.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/omp_pragma_def.h"
#include "dvloc.h"
#include "elphC.h"

static void long_range_3D_kernel(const ELPH_float* qplusG,
                                 const ELPH_float* Zval,
                                 const ELPH_float* Zborn_k,
                                 const ELPH_float* Qpole_k,
                                 const ELPH_float* epslion,
                                 const ELPH_float* tau_k, ELPH_cmplx* out_buf);

static void long_range_2D_kernel(const ELPH_float* qplusG,
                                 const ELPH_float* Zval,
                                 const ELPH_float* Zborn_k,
                                 const ELPH_float* Qpole_k,
                                 const ELPH_float* epslion,
                                 const ELPH_float* tau_k, const ELPH_float zlat,
                                 const ELPH_float qz, ELPH_cmplx* out_buf);

void elph_lr_vertex(const ELPH_float* qpt, const ELPH_float* gvecs,
                    const ND_int npw_loc, const ELPH_float* Zvals,
                    const ELPH_float* epslion, const ELPH_float* Zeu,
                    const ELPH_float* Qpole, const ND_int natom,
                    const ELPH_float* atom_pos, const char diminsion,
                    const ELPH_float volume, const ELPH_float zlat,
                    const ELPH_float EcutRy, ELPH_cmplx* elph_lr_out)
{
    // gvecs in cartisian coordinates  ( no 2*pi)
    // qpt in cart units             ( no 2*pi)
    // npw_loc, number of gvecs in this cpu
    // atomic pos in cart units
    // elph_lr_out (npw_loc, natom,3)
    // zlat in the dimension along the out ot plane direction.
    // only used in the 2D case

    // first zero out the elph_lr_out buffer

    // note that in ca
    // //
    // NOTE: Donot forget to allreduce the result over plane waves.

    // EcutRy : Energy cut off in Ry. The deformation potential pws for
    // |q+G|^2 > EcutRy is set to zero.
    // //
    //
    for (ND_int i = 0; i < (3 * natom * npw_loc); ++i)
    {
        elph_lr_out[i] = 0.0;
    }

    if (!epslion && !Zvals)
    {
        return;
    }

    const ELPH_cmplx factor =
        4.0 * ELPH_PI * I * ELPH_e2 / volume;  // prefactor

    ELPH_OMP_PAR_FOR_SIMD
    for (ND_int ig = 0; ig < npw_loc; ++ig)
    {
        const ELPH_float* gtmp = gvecs + 3 * ig;
        ELPH_float qplusG[3];
        for (int i = 0; i < 3; ++i)
        {
            qplusG[i] = 2 * ELPH_PI * (qpt[i] + gtmp[i]);
        }
        const ELPH_float qplusG_norm2 = dot3_macro(qplusG, qplusG);
        if (qplusG_norm2 > EcutRy)
        {
            continue;
        }
        for (ND_int ia = 0; ia < natom; ++ia)
        {
            const ELPH_float* Zval = Zvals ? (Zvals + ia) : NULL;
            const ELPH_float* Z_k = Zeu ? (Zeu + 9 * ia) : NULL;
            const ELPH_float* Q_k = Qpole ? (Qpole + 27 * ia) : NULL;
            const ELPH_float* tau_k = atom_pos + 3 * ia;

            ELPH_cmplx* out_tmp_buf = elph_lr_out + ia * 3 + ig * natom * 3;

            if (diminsion == '3')
            {
                long_range_3D_kernel(qplusG, Zval, Z_k, Q_k, epslion, tau_k,
                                     out_tmp_buf);
            }
            else if (diminsion == '2')
            {
                long_range_2D_kernel(qplusG, Zval, Z_k, Q_k, epslion, tau_k,
                                     zlat, 2 * ELPH_PI * qpt[2], out_tmp_buf);
            }
            // multiply with prefactor
            for (ND_int i = 0; i < 3; ++i)
            {
                out_tmp_buf[i] *= factor;
            }
        }
    }

    return;
}

static void long_range_3D_kernel(const ELPH_float* qplusG,
                                 const ELPH_float* Zval,
                                 const ELPH_float* Zborn_k,
                                 const ELPH_float* Qpole_k,
                                 const ELPH_float* epslion,
                                 const ELPH_float* tau_k, ELPH_cmplx* out_buf)
{
    // Compute long range part of the electron-phonon matrix elements for 3D
    // case. Zval -> valence charge (for monopole) Z_born_k-> born charges (for
    // dipole) Qpole_k -> Quadrupoles In general we donot need monolope when we
    // considered total dVSCF i.e dV_ext + dVHa_xc. But in Q.E, only dVHa_xc is
    // considered, so sometime we need to subtract the monopole (if there is
    // external potential this will cancell the monopole long range compoment)
    //
    //
    //
    // G space dipole term for frohlich in 3D
    // Eq :4 of PhysRevLett.115.176401 (C. Verdi, F. Giustino)

    // tau_k are atomic coordinates of atom k in cart
    // q+G in cart with 2*pi included ,
    // Zborn_k, born charges for atom k

    if (!Zval && !Zborn_k && !Qpole_k)
    {
        return;
    }
    //
    ELPH_float q_G_square = dot3_macro(qplusG, qplusG);
    //
    if (sqrt(q_G_square) < ELPH_EPS)
    {
        // skip q+G = 0 term
        return;
    }
    // struture factor
    ELPH_cmplx qdot_tau = cexp(-I * (dot3_macro(qplusG, tau_k)));

    ELPH_float q_eps_q = q_G_square;
    if (epslion)
    {
        // compute (q+G).eps.(q+G)
        ELPH_float tmp_buf[3] = {0.0, 0.0, 0.0};
        MatVec3f(epslion, qplusG, false, tmp_buf);
        q_eps_q = dot3_macro(tmp_buf, qplusG);
    }

    ELPH_float Zval_buf[3] = {0.0, 0.0, 0.0};
    ELPH_float Zborn_buf[3] = {0.0, 0.0, 0.0};
    ELPH_float Qpole_buf[3] = {0.0, 0.0, 0.0};
    //
    if (Zval)
    {
        for (int i = 0; i < 3; ++i)
        {
            Zval_buf[i] = -(*Zval) * qplusG[i] / q_G_square;
        }
    }
    // compute (q+G).Z
    if (Zborn_k)
    {
        MatVec3f(Zborn_k, qplusG, true, Zborn_buf);
    }
    // compute (q+G)_x Q_xyz * (q+G)_y
    if (Qpole_k)
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    Qpole_buf[k] =
                        Qpole_buf[k] +
                        qplusG[i] * qplusG[j] * Qpole_k[k + 3 * j + 9 * i];
                }
            }
        }
    }
    // compute and multiply with a decay factor
    qdot_tau *= exp(-q_G_square * 0.25);

    for (int i = 0; i < 3; ++i)
    {
        out_buf[i] =
            qdot_tau *
            ((Zborn_buf[i] - I * 0.5 * Qpole_buf[i]) / q_eps_q + Zval_buf[i]);
    }

    return;
}

static void long_range_2D_kernel(const ELPH_float* qplusG,
                                 const ELPH_float* Zval,
                                 const ELPH_float* Zborn_k,
                                 const ELPH_float* Qpole_k,
                                 const ELPH_float* epslion,
                                 const ELPH_float* tau_k, const ELPH_float zlat,
                                 const ELPH_float qz, ELPH_cmplx* out_buf)
{
    // Compute long range part of the electron-phonon matrix elements for 2D
    // case. Zval -> valence charge (for monopole) Z_born_k-> born charges (for
    // dipole) Qpole_k -> Quadrupoles In general we donot need monolope when we
    // considered total dVSCF i.e dV_ext + dVHa_xc. But in Q.E, only dVHa_xc is
    // considered, so sometime we need to subtract the monopole (if there is
    // external potential this will cancell the monopole long range compoment)
    //
    // G space dipole term for frohlich in 2D
    // T Deng et. al Phys. Rev. B 103, 075410 (2021)
    // T Sohier et. al  Phys. Rev. B 96, 075448 (2017)

    // tau_k are atomic coordinates of atom k in cart
    // q+G, qz in cart with 2*pi included ,
    // Zborn_k, born charges for atom k
    // zlat : lattice parameter along z parameter.

    if (!Zval && !Zborn_k && !Qpole_k)
    {
        return;
    }
    //
    ELPH_float q_G_square = dot3_macro(qplusG, qplusG);
    //
    if (sqrt(q_G_square) < ELPH_EPS)
    {
        // skip q+G = 0 term
        return;
    }
    // struture factor
    ELPH_cmplx qdot_tau = cexp(-I * (dot3_macro(qplusG, tau_k)));

    ELPH_float Gp_norm = sqrt(qplusG[0] * qplusG[0] + qplusG[1] * qplusG[1]);
    ELPH_float qGp = zlat * Gp_norm;
    ELPH_float cos_sin_fac = cos(qplusG[2] * zlat * 0.5);
    if (Gp_norm > ELPH_EPS)
    {
        cos_sin_fac =
            cos_sin_fac - sin(qplusG[2] * zlat * 0.5) * qplusG[2] / Gp_norm;
    }
    ELPH_float cutoff_fac = 1 - exp(-qGp * 0.5) * cos_sin_fac;

    ELPH_float q_eps_q = q_G_square;
    if (epslion)
    {
        // compute (q+G).eps.(q+G)
        ELPH_float tmp_buf[3] = {0.0, 0.0, 0.0};
        MatVec3f(epslion, qplusG, false, tmp_buf);
        q_eps_q = dot3_macro(tmp_buf, qplusG);
        q_eps_q = 0.5 * (q_eps_q - q_G_square) / sqrt(q_G_square) + 1.0 / zlat;
    }

    ELPH_float Zval_buf[3] = {0.0, 0.0, 0.0};
    ELPH_float Zborn_buf[3] = {0.0, 0.0, 0.0};
    ELPH_float Qpole_buf[3] = {0.0, 0.0, 0.0};
    //
    if (Zval)
    {
        for (int i = 0; i < 3; ++i)
        {
            Zval_buf[i] = -(*Zval) * qplusG[i] * cutoff_fac / q_G_square;
        }
    }
    // compute (q+G).Z
    if (Zborn_k)
    {
        MatVec3f(Zborn_k, qplusG, true, Zborn_buf);
    }
    // compute (q+G)_x Q_xyz * (q+G)_y
    if (Qpole_k)
    {
        ELPH_float Q_mod_inplane = sqrt(Gp_norm * Gp_norm + qz * qz);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int k = 0; k < 3; ++k)
                {
                    Qpole_buf[k] =
                        Qpole_buf[k] +
                        qplusG[i] * qplusG[j] * Qpole_k[k + 3 * j + 9 * i];
                }
            }
        }
        // subtract Qzz*(1+|qz+Gp|)*|q+G|^2-2qz*Gz
        for (int i = 0; i < 3; ++i)
        {
            Qpole_buf[i] -= ((q_G_square - 2.0 * qz * (qplusG[2] - qz)) *
                             (1 + Q_mod_inplane) * Qpole_k[i + 24]);
        }
    }
    // compute and multiply with a decay factor
    qdot_tau *= exp(-q_G_square * 0.25);

    for (int i = 0; i < 3; ++i)
    {
        out_buf[i] =
            qdot_tau * (((Zborn_buf[i] - I * 0.5 * Qpole_buf[i]) / q_eps_q) /
                            (q_G_square - 2.0 * qz * (qplusG[2] - qz)) +
                        Zval_buf[i]);
    }

    return;
}

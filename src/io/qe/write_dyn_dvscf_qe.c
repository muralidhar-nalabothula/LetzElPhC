// THis file contains write functions
// to read fake dyn files and dvscfs.
// The file can only be read by Letzelphc
//
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/error.h"
#include "elphC.h"
#include "qe_io.h"

void write_dyn_qe(const char* file_name, ND_int natom, const ELPH_float* qpts,
                  const ELPH_cmplx* dyn_mat, const ELPH_float* atomic_masses)
{
    //
    FILE* fp = fopen(file_name, "w");
    if (!fp)
    {
        error_msg("Error creating dyn file");
    }

    // Write comment lines
    fprintf(fp, "%s\n", "Create by LetzElphC.");
    fprintf(fp, "%s\n", "This is fake dyn file only to be read by LetzElphC.");

    ND_int ntype = natom;  // set it same as natoms
    ND_int ibrav = -1;     // set to negative (should not be 0)
    //
    fprintf(fp,
            "%lld    %lld    %lld    %.8f    %.8f    %.8f    %.8f    %.8f    "
            "%.8f\n",
            (long long)ntype, (long long)natom, (long long)ibrav, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0);

    // Write atom types and masses
    for (ND_int ia = 0; ia < natom; ia++)
    {
        fprintf(fp, "%lld    'X%lld'    %.10f\n", (long long)(ia + 1),
                (long long)(ia + 1), atomic_masses[ia]);
    }

    // Write atom positions
    for (ND_int ia = 0; ia < natom; ia++)
    {
        fprintf(fp, "%lld    %lld    %.10f    %.10f    %.10f\n",
                (long long)(ia + 1), (long long)(ia + 1), 0.0, 0.0, 0.0);
    }

    // For each q-point
    ND_int nq = 1;
    for (ND_int iq = 0; iq < nq; iq++)
    {
        const ELPH_float* qpt = qpts + iq * 3;
        const ELPH_cmplx* dyn_mat_q = dyn_mat + iq * 9 * natom * natom;

        // Write q-poND_int header
        fprintf(fp, "\n");
        fprintf(fp, "Dynamical  Matrix in cartesian axes\n");
        fprintf(fp, "\n");
        fprintf(fp, "q = ( %.10f   %.10f   %.10f )\n", qpt[0], qpt[1], qpt[2]);
        fprintf(fp, "\n");

        // Write dynamical matrix blocks
        for (ND_int ia = 0; ia < natom; ia++)
        {
            for (ND_int ib = 0; ib < natom; ib++)
            {
                fprintf(fp, "%lld    %lld\n", (long long)(ia + 1),
                        (long long)(ib + 1));

                ELPH_float mass_sqrt =
                    sqrt(atomic_masses[ia] * atomic_masses[ib]);
                for (ND_int ix = 0; ix < 3; ix++)
                {
                    for (ND_int iy = 0; iy < 3; iy++)
                    {
                        // our dynamical matrix is in row major
                        // but in dynma we store in coloumn major.
                        ND_int idx =
                            iy + ib * 3 + ix * 3 * natom + ia * 3 * 3 * natom;
                        fprintf(fp, "%.8f   %.8f     ",
                                creal(dyn_mat_q[idx]) * mass_sqrt,
                                cimag(dyn_mat_q[idx]) * mass_sqrt);
                    }
                    fprintf(fp, "\n");
                }
            }
        }
    }

    fprintf(fp, "\n");
    fclose(fp);
}

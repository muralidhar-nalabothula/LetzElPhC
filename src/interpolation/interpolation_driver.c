#include <mpi.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/constants.h"
#include "common/cwalk/cwalk.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/init_dtypes.h"
#include "common/numerical_func.h"
#include "common/parallel.h"
#include "common/print_info.h"
#include "dvloc/dvloc.h"
#include "elphC.h"
#include "interpolation_utilities.h"
#include "io/io.h"
#include "io/qe/qe_io.h"
#include "phonon/phonon.h"
#include "symmetries/symmetries.h"

void interpolation_driver(const char* ph_save, const char* ph_save_interpolated,
                          enum ELPH_dft_code dft_code, const ND_int* qgrid_new,
                          MPI_Comm comm_world)
{
    //
    init_ELPH_clocks();

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    // Only plain wave plarallization is supported.
    create_parallel_comms(1, 1, comm_world, mpi_comms);

    print_ELPH_logo(mpi_comms->commW_rank, stdout);
    print_info_msg(mpi_comms->commW_rank,
                   "********** Interpolation Program started **********");

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);
    init_lattice_type(lattice);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);
    init_phonon_type(phonon);

    bool interpolate_dvscf = true;

    ELPH_float* Zvals = NULL;
    if (dft_code == DFT_CODE_QE)
    {
        get_interpolation_data_from_qe(lattice, phonon, ph_save, &Zvals,
                                       mpi_comms);
    }
    else
    {
        error_msg("Only qe is supported currently.");
    }
    //
    //
    // We need atomic masses
    ELPH_float* atomic_masses = malloc(sizeof(*atomic_masses) * lattice->natom);
    CHECK_ALLOC(atomic_masses);

    ELPH_float* dummy1 = malloc(sizeof(*dummy1) * lattice->nmodes);
    CHECK_ALLOC(dummy1);

    ELPH_cmplx* dummy2 =
        malloc(sizeof(*dummy2) * lattice->nmodes * lattice->nmodes);
    CHECK_ALLOC(dummy2);

    if (dft_code == DFT_CODE_QE)
    {
        char read_buf[1024];
        cwk_path_join(ph_save, "dyn1", read_buf, sizeof(read_buf));
        ELPH_float qpt_tmp[3];
        ND_int iq_read = read_dyn_qe(read_buf, lattice, qpt_tmp, dummy1, dummy2,
                                     atomic_masses);
        if (iq_read != 1)
        {
            error_msg("More than 1 dynmat read.");
        }

        if (interpolate_dvscf)
        {
            // read the first pattern file
            cwk_path_join(ph_save, "patterns.1.xml", read_buf,
                          sizeof(read_buf));
            read_pattern_qe(read_buf, lattice, dummy2);
        }
    }

    //
    // get the coarse q-grid
    ND_int q_grid_co[3];
    find_qpt_grid(phonon->nq_BZ, phonon->qpts_BZ, q_grid_co);
    // find qBZ to fft grid indices
    ND_int* indices_q2fft = malloc(phonon->nq_BZ * sizeof(*indices_q2fft));
    CHECK_ALLOC(indices_q2fft);
    //
    Sorted_qpts_idxs(phonon->nq_BZ, phonon->qpts_BZ, indices_q2fft);
    //
    // read dvscf and eigen_vectors
    // allocate large buffers
    ELPH_cmplx* dVscfs_co = NULL;
    ELPH_cmplx* dyns_co = NULL;

    ND_int dvscf_loc_len = lattice->nmodes * lattice->nmag *
                           lattice->nfftz_loc * lattice->fft_dims[0] *
                           lattice->fft_dims[1];
    if (interpolate_dvscf)
    {
        dVscfs_co = malloc(phonon->nq_BZ * dvscf_loc_len * sizeof(*dVscfs_co));
        CHECK_ALLOC(dVscfs_co);
    }

    dyns_co = malloc(phonon->nq_BZ * lattice->nmodes * lattice->nmodes *
                     sizeof(*dyns_co));
    CHECK_ALLOC(dyns_co);

    ELPH_float* omega_ph_co =
        malloc(phonon->nq_BZ * lattice->nmodes * sizeof(*omega_ph_co));
    CHECK_ALLOC(omega_ph_co);

    ELPH_float EcutRy = 10;
    bool dvscf_composite_form = false;
    // Always initiate to false
    bool nmags_add_long_range[4] = {false, false, false, false};
    bool only_induced_part_long_range = false;
    if (dft_code == DFT_CODE_QE)
    {
        // q.e stores dvscf in [V,Bx,By,Bz]
        dvscf_composite_form = false;
        only_induced_part_long_range = false;
        nmags_add_long_range[0] = true;
        if (lattice->nmag == 2)
        {
            nmags_add_long_range[1] = true;
        }
    }

    ND_int iqpt_tmp = 0;
    for (ND_int iqco = 0; iqco < phonon->nq_iBZ; ++iqco)
    {
        // read eigen vectors and dvscf
        // we will for now read eigenvectors in dyns_co buffer
        // latter we convert it to dynamical matrices
        ND_int iq_fft_idx = indices_q2fft[iqpt_tmp];
        //
        ELPH_cmplx* dV_co_tmp =
            dVscfs_co ? dVscfs_co + iq_fft_idx * dvscf_loc_len : NULL;
        ELPH_cmplx* eigs_co =
            dyns_co + iq_fft_idx * lattice->nmodes * lattice->nmodes;
        ELPH_float* omege_ph_co = omega_ph_co + iq_fft_idx * lattice->nmodes;
        //
        if (dft_code == DFT_CODE_QE)
        {
            get_dvscf_dyn_qe(ph_save, lattice, iqco, eigs_co, dV_co_tmp,
                             omege_ph_co, mpi_comms);
        }
        if (dV_co_tmp)
        {
            // remore long range
            dV_add_longrange(phonon->qpts_BZ + iqpt_tmp * 3, lattice, phonon,
                             Zvals, eigs_co, dV_co_tmp, -1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, mpi_comms->commK);
        }
        ++iqpt_tmp;
        // we will remove the long range part of dynmats later
        //
        for (ND_int istar = 1; istar < phonon->nqstar[iqco]; ++istar)
        {
            iq_fft_idx = indices_q2fft[iqpt_tmp];
            //
            ELPH_cmplx* dV_co_star =
                dVscfs_co ? dVscfs_co + iq_fft_idx * dvscf_loc_len : NULL;
            ELPH_cmplx* eigs_co_star =
                dyns_co + iq_fft_idx * lattice->nmodes * lattice->nmodes;
            ELPH_float* omege_ph_co_star =
                omega_ph_co + iq_fft_idx * lattice->nmodes;

            //
            ND_int iq_iBZ = phonon->qmap[2 * iqpt_tmp];
            ND_int idx_qsym = phonon->qmap[2 * iqpt_tmp + 1];
            if (iq_iBZ != iqco)
            {
                error_msg("Wrong qstar.");
            }
            // rotate dvscf

            memcpy(omege_ph_co_star, omege_ph_co,
                   lattice->nmodes * sizeof(*omega_ph_co));
            //
            struct symmetry* sym_star = phonon->ph_syms + idx_qsym;
            rotate_eig_vecs(sym_star, lattice, phonon->qpts_BZ + 3 * iqpt_tmp,
                            eigs_co, eigs_co_star);
            if (dVscfs_co)
            {
                rotate_dvscf(dV_co_tmp, sym_star, lattice, dvscf_composite_form,
                             dV_co_star, mpi_comms->commK);
            }
            ++iqpt_tmp;
        }
    }

    free(indices_q2fft);
    //
    if (dVscfs_co)
    {
        for (ND_int iq = 0; iq < phonon->nq_BZ; ++iq)
        {
            ELPH_cmplx* rot_vecs =
                dyns_co + iq * lattice->nmodes * lattice->nmodes;
            // remove mass normalization in the dynmats
            mass_normalize_pol_vecs(atomic_masses, phonon->nmodes,
                                    phonon->natom, 1.0, rot_vecs);

            dVscf_change_basis(dVscfs_co + iq * dvscf_loc_len, rot_vecs, 1,
                               lattice->nmodes, lattice->nmag,
                               lattice->fft_dims[0], lattice->fft_dims[1],
                               lattice->nfftz_loc, 'C');
            // get back to previous normalization
            mass_normalize_pol_vecs(atomic_masses, phonon->nmodes,
                                    phonon->natom, -1.0, rot_vecs);
        }
    }
    // Now perform fourier transform
    if (dVscfs_co)
    {
        fft_q2R(dVscfs_co, q_grid_co, dvscf_loc_len);
    }

    /* void fft_R2q(const ELPH_cmplx* dataR, const ELPH_float* qpt_crys, */
    /*          const ND_int* qgrid, const ND_int nsets, const ND_int Nx, */
    /*          const ND_int Ny, const ND_int Nz, ELPH_cmplx* dataq); */

    //

    int World_rank_tmp = mpi_comms->commW_rank;

    free(atomic_masses);
    free(dummy1);
    free(dummy2);

    free(omega_ph_co);
    free(dVscfs_co);
    free(dyns_co);

    free(Zvals);
    free(lattice->atomic_pos);
    free_phonon_data(phonon);
    //
    free(lattice);
    free(phonon);
    //
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    // From here nothing should happen apart from last minite things such
    // as printing clocks
    //

    if (0 == World_rank_tmp)
    {
        print_ELPH_clock_summary();
    }
    // cleanup the clocks
    cleanup_ELPH_clocks();
    // done with the calculation
    print_info_msg(World_rank_tmp,
                   "********** Interpolation Program ended **********");

    return;
}

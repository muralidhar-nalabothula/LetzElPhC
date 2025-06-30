#include <math.h>
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
#include "interpolation.h"
#include "interpolation_utilities.h"
#include "io/io.h"
#include "io/qe/qe_io.h"
#include "phonon/phonon.h"
#include "symmetries/symmetries.h"
#include "wfc/wfc.h"

void interpolation_driver(const char* ELPH_input_file,
                          enum ELPH_dft_code dft_code, MPI_Comm comm_world)
{
    //
    /* ND_int LLLLLLLL[3] = {6, 6, 1}; */
    /* const char* ph_save = "../MoS2_SOC_2D/ph_save"; */
    /* const char* ph_save_interpolated = "../ph_interpolated"; */
    ND_int LLLLLLLL[3] = {9, 9, 1};
    const char* ph_save = "../MoS2_SOC_2D/ph_save";
    const char* ph_save_interpolated = "../ph_interpolated";
    const ND_int* qgrid_new = LLLLLLLL;

    init_ELPH_clocks();

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    int mpi_error = MPI_SUCCESS;
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
    ELPH_float alat_scale[3];
    if (dft_code == DFT_CODE_QE)
    {
        get_interpolation_data_from_qe(lattice, phonon, ph_save, &Zvals,
                                       alat_scale, mpi_comms);
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

    // There are reference patern basis.
    ELPH_cmplx* ref_pat_basis =
        malloc(sizeof(*ref_pat_basis) * lattice->nmodes * lattice->nmodes);
    CHECK_ALLOC(ref_pat_basis);

    if (dft_code == DFT_CODE_QE)
    {
        if (0 == mpi_comms->commW_rank)
        {
            char read_buf[1024];
            cwk_path_join(ph_save, "dyn1", read_buf, sizeof(read_buf));
            ELPH_float qpt_tmp[3];
            ND_int iq_read = read_dyn_qe(read_buf, lattice, qpt_tmp, dummy1,
                                         ref_pat_basis, atomic_masses);
            if (iq_read != 1)
            {
                error_msg("More than 1 dynmat read.");
            }

            if (interpolate_dvscf)
            {
                // read the first pattern file
                cwk_path_join(ph_save, "patterns.1.xml", read_buf,
                              sizeof(read_buf));
                read_pattern_qe(read_buf, lattice, ref_pat_basis);
            }
        }

        //
        mpi_error = MPI_Bcast(atomic_masses, lattice->natom, ELPH_MPI_float, 0,
                              mpi_comms->commW);
        MPI_error_msg(mpi_error);

        mpi_error = MPI_Bcast(ref_pat_basis, lattice->nmodes * lattice->nmodes,
                              ELPH_MPI_cmplx, 0, mpi_comms->commW);
        MPI_error_msg(mpi_error);

        MPI_error_msg(mpi_error);
    }

    //
    // get the coarse q-grid
    ND_int q_grid_co[3];
    find_qpt_grid(phonon->nq_BZ, phonon->qpts_BZ, q_grid_co);
    // find qBZ to fft grid indices
    ND_int* indices_q2fft = malloc(2 * phonon->nq_BZ * sizeof(*indices_q2fft));
    CHECK_ALLOC(indices_q2fft);
    //
    Sorted_qpts_idxs(phonon->nq_BZ, phonon->qpts_BZ,
                     indices_q2fft + phonon->nq_BZ);
    for (ND_int iq = 0; iq < phonon->nq_BZ; ++iq)
    {
        indices_q2fft[indices_q2fft[phonon->nq_BZ + iq]] = iq;
    }
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

    ELPH_float EcutRy = 30;
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
        ELPH_float* omega_co_tmp = omega_ph_co + iq_fft_idx * lattice->nmodes;
        //
        if (dft_code == DFT_CODE_QE)
        {
            get_dvscf_dyn_qe(ph_save, lattice, iqco, eigs_co, dV_co_tmp,
                             omega_co_tmp, mpi_comms);
        }
        if (dV_co_tmp)
        {
            // remore long range
            dV_add_longrange(phonon->qpts_iBZ + iqco * 3, lattice, phonon,
                             Zvals, eigs_co, dV_co_tmp, -1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, mpi_comms->commK);
            // THe potential here is lattice periodic and not q peridoic.
            // due to e^{-iq.tau}. so we need to remove it
            // But this will be done later when dvscf in in cart basis
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
            ELPH_float* omega_co_tmp_star =
                omega_ph_co + iq_fft_idx * lattice->nmodes;

            //
            ND_int iq_iBZ = phonon->qmap[2 * iqpt_tmp];
            ND_int idx_qsym = phonon->qmap[2 * iqpt_tmp + 1];
            if (iq_iBZ != iqco)
            {
                error_msg("Wrong qstar.");
            }
            // rotate dvscf

            memcpy(omega_co_tmp_star, omega_co_tmp,
                   lattice->nmodes * sizeof(*omega_ph_co));
            //
            struct symmetry* sym_star = phonon->ph_syms + idx_qsym;
            rotate_eig_vecs(sym_star, lattice, phonon->qpts_iBZ + iqco * 3,
                            eigs_co, eigs_co_star);
            if (dVscfs_co)
            {
                rotate_dvscf(dV_co_tmp, sym_star, lattice, dvscf_composite_form,
                             dV_co_star, mpi_comms->commK);
            }
            ++iqpt_tmp;
        }
    }

    //
    if (dVscfs_co)
    {
        for (ND_int i = 0; i < phonon->nq_BZ; ++i)
        {
            ND_int iq = indices_q2fft[i];
            ELPH_cmplx* rot_vecs =
                dyns_co + iq * lattice->nmodes * lattice->nmodes;
            // remove mass normalization in the dynmats
            mass_normalize_pol_vecs(atomic_masses, lattice->nmodes,
                                    lattice->natom, 1.0, rot_vecs);

            dVscf_change_basis(dVscfs_co + iq * dvscf_loc_len, rot_vecs, 1,
                               lattice->nmodes, lattice->nmag,
                               lattice->fft_dims[0], lattice->fft_dims[1],
                               lattice->nfftz_loc, 'C');
            // get back to previous normalization
            mass_normalize_pol_vecs(atomic_masses, lattice->nmodes,
                                    lattice->natom, -1.0, rot_vecs);
            //
            // Now remove e^{-iqtau}
            ND_int iq_iBZ = phonon->qmap[2 * i];
            ND_int idx_qsym = phonon->qmap[2 * i + 1];
            //
            ELPH_float tmp_qpt[3], qpt_cart_iq[3];
            MatVec3f(lattice->blat_vec, phonon->qpts_iBZ + iq_iBZ * 3, false,
                     tmp_qpt);
            MatVec3f(phonon->ph_syms[idx_qsym].Rmat, tmp_qpt, false,
                     qpt_cart_iq);
            //
            mul_dvscf_struct_fac(qpt_cart_iq, lattice,
                                 dvscf_loc_len / lattice->nmodes, -1,
                                 dVscfs_co + iq * dvscf_loc_len);
            // Now dvscf is fully q-periodic.
            // Time for fourier transform.
        }
        //
        // Now perform fourier transform
        fft_q2R(dVscfs_co, q_grid_co, dvscf_loc_len);
    }

    // Get a grid to perform summation
    ND_int Ggrid_phonon[3];
    get_fft_box(EcutRy, lattice->blat_vec, Ggrid_phonon, mpi_comms->commK);

    // first compute the long_range asr term
    ELPH_cmplx* dyn_mat_asr_lr =
        calloc(lattice->natom * 9, sizeof(*dyn_mat_asr_lr));
    CHECK_ALLOC(dyn_mat_asr_lr);

    compute_dyn_lr_asr_correction(lattice, phonon, Ggrid_phonon, atomic_masses,
                                  dyn_mat_asr_lr);

    // FIX ME need to parallize
    for (ND_int i = 0; i < phonon->nq_BZ; ++i)
    {
        ND_int iq = indices_q2fft[i];
        // convert polarization vectors to dynamical matrix and remove long
        // range part
        ELPH_cmplx* pol_vecs_iq =
            dyns_co + iq * lattice->nmodes * lattice->nmodes;

        pol_vecs_to_dyn(omega_ph_co + iq * lattice->nmodes, lattice->natom,
                        atomic_masses, pol_vecs_iq);

        const ELPH_float* qpt_iq_tmp = phonon->qpts_BZ + 3 * i;
        // remove long range part
        add_ph_dyn_long_range(qpt_iq_tmp, lattice, phonon, Ggrid_phonon, -1,
                              atomic_masses, dyn_mat_asr_lr, pol_vecs_iq);
        //
    }

    // fourier transform phonons
    fft_q2R(dyns_co, q_grid_co, lattice->nmodes * lattice->nmodes);

    ND_int nqpts_to_interpolate = qgrid_new[0] * qgrid_new[1] * qgrid_new[2];
    // this will be over written lattern with number of qpts in iBZ

    ELPH_float* qpts_interpolation =
        malloc(sizeof(*qpts_interpolation) * 3 * nqpts_to_interpolate);
    // in crystal coordinates
    nqpts_to_interpolate = generate_iBZ_kpts(
        qgrid_new, phonon->nph_sym, phonon->ph_syms, lattice->alat_vec,
        lattice->blat_vec, qpts_interpolation, true);

    ELPH_cmplx* dvscf_interpolated = NULL;
    if (dVscfs_co)
    {
        dvscf_interpolated =
            malloc(dvscf_loc_len * sizeof(*dvscf_interpolated));
        CHECK_ALLOC(dvscf_interpolated);
    }

    ELPH_cmplx* dyn_interpolated =
        malloc(sizeof(*dyn_interpolated) * lattice->nmodes * lattice->nmodes);
    CHECK_ALLOC(dyn_interpolated);

    // now interpolate
    for (ND_int iq = 0; iq < nqpts_to_interpolate; ++iq)
    {
        char read_buf[1024];
        char dvscf_dyn_name[32];

        ELPH_float* qpt_interpolate = qpts_interpolation + 3 * iq;

        ELPH_float qpt_interpolate_cart[3];
        MatVec3f(lattice->blat_vec, qpt_interpolate, false,
                 qpt_interpolate_cart);
        //
        if (dVscfs_co)
        {
            snprintf(dvscf_dyn_name, sizeof(dvscf_dyn_name), "dvscf%lld",
                     (long long)(iq + 1));
            cwk_path_join(ph_save_interpolated, dvscf_dyn_name, read_buf,
                          sizeof(read_buf));
            fft_R2q(dVscfs_co, qpt_interpolate, q_grid_co,
                    lattice->nmodes * lattice->nmag, lattice->fft_dims[0],
                    lattice->fft_dims[1], lattice->nfftz_loc,
                    dvscf_interpolated);
            //
            //// add the phase back due to atomic coordinates (e^{-iq.tau}).
            mul_dvscf_struct_fac(qpt_interpolate_cart, lattice,
                                 dvscf_loc_len / lattice->nmodes, 1,
                                 dvscf_interpolated);
            // change to pattern basis
            dVscf_change_basis(dvscf_interpolated, ref_pat_basis, 1,
                               lattice->nmodes, lattice->nmag,
                               lattice->fft_dims[0], lattice->fft_dims[1],
                               lattice->nfftz_loc, 'N');
            //
            dV_add_longrange(qpt_interpolate, lattice, phonon, Zvals,
                             ref_pat_basis, dvscf_interpolated, 1,
                             only_induced_part_long_range, EcutRy,
                             nmags_add_long_range, mpi_comms->commK);

            // write to file
            if (dft_code == DFT_CODE_QE)
            {
                write_dvscf_qe(read_buf, lattice, dvscf_interpolated,
                               mpi_comms->commK);
            }
        }

        if (0 == mpi_comms->commW_rank)
        {
            // interpolate dyn file
            snprintf(dvscf_dyn_name, sizeof(dvscf_dyn_name), "dyn%lld",
                     (long long)(iq + 1));
            cwk_path_join(ph_save_interpolated, dvscf_dyn_name, read_buf,
                          sizeof(read_buf));

            fft_R2q(dyns_co, qpt_interpolate, q_grid_co, 1, 1, 1,
                    lattice->nmodes * lattice->nmodes, dyn_interpolated);
            // add back the long range part
            add_ph_dyn_long_range(qpt_interpolate, lattice, phonon,
                                  Ggrid_phonon, 1, atomic_masses,
                                  dyn_mat_asr_lr, dyn_interpolated);
            // write dyn file
            // FIX me write qpoint in alat units q.e
            // Convert qpoints to
            if (dft_code == DFT_CODE_QE)
            {
                write_dyn_qe(read_buf, lattice->natom, qpt_interpolate,
                             dyn_interpolated, atomic_masses);
            }
        }
        if (dft_code == DFT_CODE_QE)
        {
            // convert the interpolated qpoints in weird q.e units
            // ELPH_float* qpt_interpolate = qpts_interpolation + 3 * iq;
            for (ND_int ix = 0; ix < 3; ++ix)
            {
                qpt_interpolate[ix] =
                    alat_scale[ix] * qpt_interpolate_cart[ix] / (2 * ELPH_PI);
            }
        }
    }

    if (0 == mpi_comms->commW_rank)
    {
        char read_buf[1024];
        cwk_path_join(ph_save_interpolated, "dyn0", read_buf, sizeof(read_buf));

        write_qpts_qe(read_buf, nqpts_to_interpolate, qpts_interpolation,
                      qgrid_new);
    }

    free(dyn_mat_asr_lr);
    free(indices_q2fft);
    free(dvscf_interpolated);
    free(dyn_interpolated);
    free(qpts_interpolation);
    int World_rank_tmp = mpi_comms->commW_rank;

    free(atomic_masses);
    free(dummy1);
    free(ref_pat_basis);

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

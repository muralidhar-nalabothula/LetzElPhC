/*
THe starting point for the entire code
*/
#include <complex.h>
#include <fftw3.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/ELPH_timers.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/init_dtypes.h"
#include "common/parallel.h"
#include "common/print_info.h"
#include "dvloc/dvloc.h"
#include "elph.h"
#include "elphC.h"
#include "fft/fft.h"
#include "io/elph_hdf5.h"
#include "io/io.h"
#include "io/qe/qe_io.h"
#include "parser/parser.h"
#include "symmetries/symmetries.h"

void elph_driver(const char* ELPH_input_file, enum ELPH_dft_code dft_code,
                 MPI_Comm comm_world)
{
    struct elph_usr_input* input_data;
    // start the clocks
    init_ELPH_clocks();
    //
    // read the input file
    read_elph_input_file(ELPH_input_file, &input_data, comm_world);
    // Note input parameters are broadcasted internally
    // All the parameters in input_data must be available for all cpus in
    // comm_world

    struct kernel_info* kernel = malloc(sizeof(struct kernel_info));
    // initate the kernel with default
    init_kernel(kernel);
    // set the kernel
    set_kernel(input_data->kernel_str, kernel);

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    create_parallel_comms(input_data->nqpool, input_data->nkpool, comm_world,
                          mpi_comms);

    // print logo and stated message
    print_ELPH_logo(mpi_comms->commW_rank, stdout);
    print_info_msg(mpi_comms->commW_rank,
                   "********** Program started **********");
    print_input_info(input_data->save_dir, input_data->ph_save_dir,
                     input_data->kernel_str, input_data->kminusq, dft_code,
                     mpi_comms);

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);
    init_lattice_type(lattice);

    struct Pseudo* pseudo = malloc(sizeof(struct Pseudo));
    CHECK_ALLOC(pseudo);
    init_Pseudo_type(pseudo);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);
    init_phonon_type(phonon);

    struct WFC* wfcs;

    // read the SAVE data and phonon related data.
    read_and_alloc_save_data(input_data->save_dir, mpi_comms,
                             input_data->start_bnd, input_data->end_bnd, &wfcs,
                             input_data->ph_save_dir, lattice, pseudo, phonon,
                             dft_code);

    // print info about lattice and phonons
    print_lattice_info(mpi_comms, lattice);
    print_phonon_info(mpi_comms, phonon);
    //======= Now start the real computation =========
    // a) COmpute the D_mats and store them in the netcdf file
    // ============= Dmats =====================
    print_info_msg(mpi_comms->commW_rank, "");
    print_info_msg(mpi_comms->commW_rank,
                   "=== Computing Dmats for phonon symmetries ===");
    compute_and_write_dmats("ndb.Dmats", wfcs, lattice, phonon->nph_sym,
                            phonon->ph_syms, mpi_comms);
    // b) Compute elph
    // ============= ELPH iBZ computation =============
    ND_int nmodes = lattice->nmodes;
    ND_int nfft_loc =
        lattice->fft_dims[0] * lattice->fft_dims[1] * lattice->nfftz_loc;

    ELPH_cmplx* eigVec = malloc(sizeof(ELPH_cmplx) * nmodes * nmodes);
    //(nmodes,nmodes) // buffer to store eigen vectors
    CHECK_ALLOC(eigVec);

    ELPH_cmplx* dVscf =
        malloc(sizeof(ELPH_cmplx) * nmodes * lattice->nmag * nfft_loc);
    // (nmodes,nmag,Nx,Ny,Nz) // bufffer to store dVscf
    CHECK_ALLOC(dVscf);

    ELPH_float* omega_ph = malloc(sizeof(ELPH_float) * nmodes);
    // buffer for storing phonon freq
    CHECK_ALLOC(omega_ph);

    hid_t file_id_elph, file_id_dmat;
    hid_t dset_id_eig, dset_id_elph, dset_id_omega, dset_id_dmat;

    // Define netcdf variables
    if (mpi_comms->commK_rank == 0)
    {
        // open Dmat file
        file_id_dmat = elph_h5_open_file_par("ndb.Dmats", 1, mpi_comms->commR,
                                             MPI_INFO_NULL);

        // get dmat var id for dmats
        dset_id_dmat = elph_h5_open_var(file_id_dmat, "Dmats");

        // create elph file. Note: we overwrite any existing file
        file_id_elph = elph_h5_create_file_par("ndb.elph", mpi_comms->commR,
                                               MPI_INFO_NULL);

        elph_h5_def_var(file_id_elph, "POLARIZATION_VECTORS", ELPH_H5_IO_FLOAT,
                        5, (ND_int[]){phonon->nq_BZ, nmodes, nmodes / 3, 3, 2},
                        (const char*[]){"nq", "nmodes", "atom", "pol", "re_im"},
                        NULL, &dset_id_eig);

        elph_h5_def_var(file_id_elph, "FREQ", ELPH_H5_IO_FLOAT, 2,
                        (ND_int[]){phonon->nq_BZ, nmodes},
                        (const char*[]){"nq", "nmodes"}, NULL, &dset_id_omega);

        ND_int nk_chunk_size =
            H5_DEFAULT_CHUNK_KB * 1024;  // now this is in bytes
        // scale with complex number size to get the number of elements
        nk_chunk_size /= (sizeof(ELPH_cmplx) * nmodes * lattice->nspin *
                          lattice->nbnds * lattice->nbnds);
        // chuck the varaible elph_mat with atmost default size
        if (nk_chunk_size == 0)
        {
            nk_chunk_size = 1;
        }
        else if (nk_chunk_size > lattice->nkpts_BZ)
        {
            nk_chunk_size = lattice->nkpts_BZ;
        }

        elph_h5_def_var(
            file_id_elph, "elph_mat", ELPH_H5_IO_FLOAT, 7,
            (ND_int[]){phonon->nq_BZ, lattice->nkpts_BZ, nmodes, lattice->nspin,
                       lattice->nbnds, lattice->nbnds, 2},
            (const char*[]){"nq", "nk", "nmodes", "nspin", "initial_band",
                            "final_band_PH_abs", "re_im"},
            (size_t[]){1, (size_t)nk_chunk_size, (size_t)nmodes,
                       (size_t)lattice->nspin, (size_t)lattice->nbnds,
                       (size_t)lattice->nbnds, 2},
            &dset_id_elph);
    }

    ELPH_cmplx* eig_Sq = NULL;

    if (mpi_comms->commQ_rank == 0)
    {
        eig_Sq = calloc(nmodes * nmodes, sizeof(ELPH_cmplx));
        CHECK_ALLOC(eig_Sq);
    }

    print_info_msg(mpi_comms->commW_rank,
                   "=== Computing Electron-phonon matrix elements ===");
    print_info_msg(mpi_comms->commW_rank, "");
    for (ND_int iqpt = 0; iqpt < phonon->nq_iBZ_loc; ++iqpt)
    {
        print_info_msg(mpi_comms->commW_rank, "### q-point : %d/%d",
                       (int)(iqpt + 1), (int)phonon->nq_iBZ_loc);

        ND_int iqpt_iBZg = iqpt + phonon->nq_shift;
        // read dynamical matrix and dvscf for the iBZ qpt
        if (dft_code == DFT_CODE_QE)
        {
            ELPH_cmplx* dvscf_read = NULL;
            if (kernel->screening == ELPH_DFPT_SCREENING)
            {
                dvscf_read = dVscf;
            }
            else
            {
                // zero out the buffer. This is must !
                ND_int dvscf_num = nmodes * lattice->nmag * nfft_loc;
                for (ND_int ix = 0; ix < dvscf_num; ++ix)
                {
                    dVscf[ix] = 0.0;
                }
            }
            //
            get_dvscf_dyn_qe(input_data->ph_save_dir, lattice, iqpt_iBZg,
                             eigVec, dvscf_read, omega_ph, mpi_comms);
            // qe dvscf only contains dV_Ha + dV_xc, we need to add the local
            // part of pseudo
        }
        else
        {
            error_msg("Currently only quantum espresso supported");
        }
        // local part
        ELPH_cmplx* Vlocr = malloc(sizeof(ELPH_cmplx) * nmodes * nfft_loc);
        // buffer to store local part of the pseudo potential
        CHECK_ALLOC(Vlocr);
        //
        // compute the local part of the bare
        dVlocq(phonon->qpts_iBZ + iqpt_iBZg * 3, lattice, pseudo, eigVec, Vlocr,
               mpi_comms->commK);

        if (dft_code == DFT_CODE_QE)
        {
            add_dvscf_qe(dVscf, Vlocr,
                         lattice);  // add bare local to induce part
        }
        else
        {
            error_msg("Currently only quantum espresso supported");
        }

        free(Vlocr);
        // free the local part buffer
        ND_int qpos = 0;  // positon of this iBZ qpoint in full q point list
        for (ND_int i = 0; i < iqpt_iBZg; ++i)
        {
            qpos += phonon->nqstar[i];
        }
        // write eigen vectors and frequencies
        if (mpi_comms->commQ_rank == 0)
        {
            hsize_t startp[5] = {(hsize_t)qpos, 0, 0, 0, 0};
            hsize_t countp[5] = {1, (hsize_t)nmodes, (hsize_t)(nmodes / 3), 3,
                                 2};
            elph_h5_write_vara(dset_id_eig, ELPH_H5_IO_FLOAT, startp, countp,
                               eigVec);
            elph_h5_write_vara(dset_id_omega, ELPH_H5_IO_FLOAT, startp, countp,
                               omega_ph);

            // write down the rotate eigen vectors;
            for (ND_int istar = 1; istar < phonon->nqstar[iqpt_iBZg]; ++istar)
            {
                ND_int qpos_star = qpos + istar;

                struct symmetry* sym_rot =
                    phonon->ph_syms + phonon->qmap[2 * qpos_star + 1];

                rotate_eig_vecs(sym_rot, lattice,
                                phonon->qpts_iBZ + iqpt_iBZg * 3, eigVec,
                                eig_Sq);
                //
                startp[0] = (hsize_t)qpos_star;
                elph_h5_write_vara(dset_id_eig, ELPH_H5_IO_FLOAT, startp,
                                   countp, eig_Sq);
                elph_h5_write_vara(dset_id_omega, ELPH_H5_IO_FLOAT, startp,
                                   countp, omega_ph);
            }
        }
        // Now compute and write the electron-phonon matrix elements
        compute_and_write_elphq(wfcs, lattice, pseudo, phonon, iqpt_iBZg,
                                eigVec, dVscf, file_id_elph, dset_id_elph,
                                file_id_dmat, dset_id_dmat, kernel->non_loc,
                                input_data->kminusq, mpi_comms);
    }

    free(eig_Sq);

    if (mpi_comms->commK_rank == 0)
    {
        // Close variables/datasets
        elph_h5_close_var(dset_id_eig);
        elph_h5_close_var(dset_id_elph);
        elph_h5_close_var(dset_id_omega);
        elph_h5_close_var(dset_id_dmat);

        // close files
        elph_h5_close_file(file_id_elph);
        elph_h5_close_file(file_id_dmat);
    }

    // finally write some basic info to ndb.elph file (only master node writes
    // it)

    if (mpi_comms->commW_rank == 0)
    {
        file_id_elph =
            elph_h5_open_file_par("ndb.elph", 0, MPI_COMM_SELF, MPI_INFO_NULL);

        char convention_str[32];
        strcpy(convention_str, "standard");
        if (input_data->kminusq)
        {
            strcpy(convention_str, "yambo");
        }
        write_basic_data(file_id_elph, lattice, phonon, input_data->kernel_str,
                         convention_str);

        elph_h5_close_file(file_id_elph);
    }

    int World_rank_tmp = mpi_comms->commW_rank;
    // ELPH_cmplx Ry2Ha = pow(2,-1.5);
    free(omega_ph);
    free(eigVec);
    free(dVscf);

    // cleanup
    free(kernel);
    free_elph_usr_input(input_data);
    free_save_data(wfcs, lattice, pseudo, phonon);
    free(lattice);
    free(pseudo);
    free(phonon);
    free_parallel_comms(mpi_comms);
    free(mpi_comms);
    fftw_fun(cleanup)();

    // print the clocks
    if (0 == World_rank_tmp)
    {
        print_ELPH_clock_summary();
    }
    // cleanup the clocks
    cleanup_ELPH_clocks();
    // done with the calculation
    print_info_msg(World_rank_tmp, "********** Program ended **********");
}

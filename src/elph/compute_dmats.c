#include <stdio.h>
#include <stdlib.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "common/parallel.h"
#include "common/progess_bar.h"
#include "elph.h"
#include "elphC.h"
#include "io/elph_hdf5.h"
#include "io/io.h"
#include "symmetries/symmetries.h"

/*
 * This function contain the wrapper functions to compute and
 * write or read dmat functions to file)
 */

void compute_and_write_dmats(const char* file_name, const struct WFC* wfcs,
                             const struct Lattice* lattice,
                             const ND_int nph_sym,
                             const struct symmetry* sym_data,
                             const struct ELPH_MPI_Comms* Comm)
{
    ND_int nk_totalBZ = lattice->nkpts_BZ;
    ELPH_cmplx* Dkmn_rep_ptr = NULL;

    hid_t file_id, dset_id;

    hsize_t startp[6] = {0, 0, 0, 0, 0, 0};
    hsize_t countp[6] = {1,
                         1,
                         (hsize_t)lattice->nspin,
                         (hsize_t)lattice->nbnds,
                         (hsize_t)lattice->nbnds,
                         2};

    ND_int nk_chunk_size = H5_DEFAULT_CHUNK_KB * 1024;  // now this is in bytes
    // scale with complex number size to get the number of elements
    nk_chunk_size /=
        (sizeof(ELPH_cmplx) * lattice->nspin * lattice->nbnds * lattice->nbnds);
    // chuck the varaible elph_mat with atmost default size
    if (nk_chunk_size == 0)
    {
        nk_chunk_size = 1;
    }
    else if (nk_chunk_size > nk_totalBZ)
    {
        nk_chunk_size = nk_totalBZ;
    }

    if (Comm->commK_rank == 0)
    {
        // we overwrite any existing file
        file_id =
            elph_h5_create_file_par(file_name, Comm->commR, MPI_INFO_NULL);

        // Define variable
        elph_h5_def_var(
            file_id, "Dmats", ELPH_H5_IO_FLOAT, 6,
            (ND_int[]){nph_sym, nk_totalBZ, lattice->nspin, lattice->nbnds,
                       lattice->nbnds, 2},
            (const char*[]){"nsym_ph", "nkpts", "nspin", "Rk_band", "k_band",
                            "re_im"},
            (size_t[]){1, (size_t)nk_chunk_size, (size_t)lattice->nspin,
                       (size_t)lattice->nbnds, (size_t)lattice->nbnds, 2},
            &dset_id);

        // HDF5 parallel access is handled via property lists in write calls if
        // needed. We leave dset_id open for writing.

        Dkmn_rep_ptr = calloc(lattice->nspin * lattice->nbnds * lattice->nbnds,
                              sizeof(ELPH_cmplx));
        CHECK_ALLOC(Dkmn_rep_ptr);
    }

    // for computation of Dmats, we use all the nodes
    // ("nsym", "nk", "nspin", "nbndb", "nbnda")
    ND_int dmat_shift;
    ND_int ndmats =
        distribute_to_grps(nph_sym * nk_totalBZ, Comm->nqpools * Comm->nkpools,
                           Comm->commW_rank / Comm->commK_size, &dmat_shift);

    // start the progress bar for dmats
    struct progress_bar pbar[1];
    start_progressbar(pbar, Comm->commW_rank, ndmats);

    for (ND_int idmat = 0; idmat < ndmats; ++idmat)
    {
        ND_int isym = (idmat + dmat_shift) / nk_totalBZ;
        ND_int ikBZ = (idmat + dmat_shift) % nk_totalBZ;

        startp[0] = (hsize_t)isym;
        startp[1] = (hsize_t)ikBZ;

        // compute the dmats
        electronic_reps(wfcs, lattice, sym_data[isym].Rmat, sym_data[isym].tau,
                        sym_data[isym].time_rev, ikBZ, Dkmn_rep_ptr, Comm);

        if (Comm->commK_rank == 0)
        {
            // write data to file
            elph_h5_write_vara(dset_id, ELPH_H5_IO_FLOAT, startp, countp,
                               Dkmn_rep_ptr);
        }

        // update the progress bar
        print_progressbar(pbar);
    }
    if (Comm->commK_rank == 0)
    {
        free(Dkmn_rep_ptr);
        elph_h5_close_var(dset_id);
        elph_h5_close_file(file_id);
    }
}

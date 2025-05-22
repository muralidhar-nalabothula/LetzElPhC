#include <netcdf.h>
#include <netcdf_par.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "../common/dtypes.h"
#include "../common/error.h"
#include "../common/parallel.h"
#include "../common/progess_bar.h"
#include "../elphC.h"
#include "../interpolation/elph_lr_kernels.h"
#include "../io/io.h"
#include "../symmetries/symmetries.h"
#include "elph.h"

/*
 * This function contain the function to compute and
 * write long range part of el-ph matrix elements)
 */

void compute_and_write_elph_lr(int ncid, const struct WFC* wfcs,
                               const struct Lattice* lattice,
                               const struct Pseudo* pseudo,
                               struct Phonon* phonon, bool is_bare,
                               const struct ELPH_MPI_Comms* Comm)
{
    // NM: Note that ncid must be open only by the rank 0 of each k-pool
    ND_int nk_totalBZ = lattice->nkpts_BZ;
    ELPH_cmplx* elph_lr_ptr = NULL;

    ND_int nkpt_shift;
    ND_int nkpts_loc =
        distribute_to_grps(nk_totalBZ, Comm->nqpools * Comm->nkpools,
                           Comm->commW_rank / Comm->commK_size, &nkpt_shift);

    int varid, nc_err;
    size_t startp[4] = {nkpt_shift, 0, 0, 0};
    size_t countp[4] = {nkpts_loc, lattice->natom, 3, 2};

    if (Comm->commK_rank == 0)
    {
        def_ncVar(ncid, &varid, 4, ELPH_NC4_IO_FLOAT,
                  (ND_int[]){nk_totalBZ, lattice->natom, 3, 2}, "elph_lr",
                  (char*[]){"nk", "natoms", "pol", "re_im"}, NULL);

        // Make the access COLLECTIVE as  all can call the put_var function
        // simultaneously
        if ((nc_err = nc_var_par_access(ncid, varid, NC_COLLECTIVE)))
        {
            ERR(nc_err);
        }

        elph_lr_ptr =
            calloc(nkpts_loc * lattice->natom * 3, sizeof(ELPH_cmplx));
        CHECK_ALLOC(elph_lr_ptr);
    }

    ELPH_cmplx* tmp_buf = calloc(3 * lattice->natom, sizeof(*tmp_buf));
    CHECK_ALLOC(tmp_buf);
    //
    ELPH_float* Zvals = malloc(lattice->natom * sizeof(*Zvals));
    CHECK_ALLOC(Zvals);

    for (int ia = 0; ia < lattice->natom; ++ia)
    {
        int itype = lattice->atom_type[ia];
        Zvals[ia] = pseudo->loc_pseudo[itype].Zval;
    }

    // start the progress bar
    struct progress_bar pbar[1];
    start_progressbar(pbar, Comm->commW_rank, nkpts_loc);

    for (ND_int ik = 0; ik < nkpts_loc; ++ik)
    {
        ND_int ikBZ = ik + nkpt_shift;
        int ikibz = lattice->kmap[2 * ikBZ];
        //
        const ELPH_float* gvec = wfcs[ikibz].gvec;
        const ND_int npw_loc = wfcs[ikibz].npw_loc;
        const ELPH_float* qpt = lattice->kpt_fullBZ + 3 * ikBZ;
        //
        if (is_bare)
        {
            bare_lr_vertex(qpt, gvec, npw_loc, Zvals, lattice->natom,
                           lattice->atomic_pos, lattice->dimension,
                           lattice->volume, lattice->alat_vec[8], tmp_buf);
        }
        else
        {
            frohlich_lr_vertex(qpt, gvec, npw_loc, phonon->epsilon,
                               phonon->Zborn, phonon->Qpole, lattice->natom,
                               lattice->atomic_pos, lattice->dimension,
                               lattice->volume, lattice->alat_vec[8], tmp_buf);
        }
        // We nee to reduce over pws
        //
        ELPH_cmplx* recv_buf =
            elph_lr_ptr ? (elph_lr_ptr + ik * 3 * lattice->natom) : NULL;
        int mpi_error = MPI_Reduce(tmp_buf, recv_buf, 3 * lattice->natom,
                                   ELPH_MPI_cmplx, MPI_SUM, 0, Comm->commK);
        MPI_error_msg(mpi_error);
        //
        // update the progress bar
        print_progressbar(pbar);
    }
    //
    free(Zvals);
    free(tmp_buf);
    //
    if (Comm->commK_rank == 0)
    {
        if ((nc_err = nc_put_vara(ncid, varid, startp, countp, elph_lr_ptr)))
        {
            ERR(nc_err);
        }

        free(elph_lr_ptr);
    }
}

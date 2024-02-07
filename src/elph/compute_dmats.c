#include "elph.h"

/*
 * This function contain the wrapper functions to compute and 
 * write or read dmat functions to file)
*/


void compute_and_write_dmats(const char * file_name, const struct WFC * wfcs, \
        const struct Lattice * lattice, const ND_int nph_sym, \
        const ELPH_float * Rsym_mat,  const ELPH_float * tauR, \
        const bool *tim_revR, const struct ELPH_MPI_Comms * Comm)
{   
    ND_int nk_totalBZ = lattice->kmap->dims[0];
    ELPH_cmplx * Dkmn_rep_ptr = NULL;
    
    int ncid, varid, nc_err, dimids[6];

    size_t startp[6]={0,0,0,0,0,0};
    size_t countp[6]={1,1,lattice->nspin,lattice->nbnds,lattice->nbnds,2};
    
    if (Comm->commK_rank == 0)
    {
        if ((nc_err = nc_create_par(file_name, NC_NETCDF4, Comm->commR, MPI_INFO_NULL, &ncid))) ERR(nc_err);

        ND_function(def_ncVar, Nd_cmplxS) (ncid, &varid, 5, nd_idx{nph_sym,nk_totalBZ, \
                lattice->nspin,lattice->nbnds,lattice->nbnds}, "Dmats", \
                (char *[]){"nsym_ph","nkpts","nspin","nbands_b","nbands_a"}, \
                (size_t[]){1,1,lattice->nspin,lattice->nbnds,lattice->nbnds,2});
        
        // Make the access INDEPENDENT as not all can call the put_var function simultaneously
        if((nc_err = nc_var_par_access(ncid, varid, NC_INDEPENDENT))) ERR(nc_err);

        Dkmn_rep_ptr = calloc(lattice->nspin*lattice->nbnds*lattice->nbnds,sizeof(ELPH_cmplx));
    }

    // for computation of Dmats, we use all the nodes
    // ("nsym", "nk", "nspin", "nbndb", "nbnda")
    ND_int dmat_shift;
    ND_int ndmats = distribute_to_grps(nph_sym*nk_totalBZ, Comm->nqpools * Comm->nkpools, \
                    Comm->commW_rank/Comm->commK_size, &dmat_shift);
    
    for (ND_int idmat= 0; idmat < ndmats; ++idmat)
    {
        ND_int isym = (idmat + dmat_shift)/nk_totalBZ;
        ND_int ikBZ = (idmat + dmat_shift)%nk_totalBZ;

        startp[0] = isym; startp[1] = ikBZ;
        
        // compute the dmats
        electronic_reps(wfcs, lattice, Rsym_mat + 9*isym,  tauR + 3*isym, \
                                tim_revR[isym], ikBZ, Dkmn_rep_ptr, Comm);
        
        if (Comm->commK_rank == 0) 
        {   
            // write data to file
            if ((nc_err = nc_put_vara(ncid, varid, startp, countp, Dkmn_rep_ptr))) ERR(nc_err);
        }
    }
    if (Comm->commK_rank == 0)
    {
        free(Dkmn_rep_ptr);

        if ((nc_err = nc_close(ncid))) ERR(nc_err);
    }
} 




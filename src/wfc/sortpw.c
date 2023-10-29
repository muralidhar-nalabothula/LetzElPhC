/*
This file contains function which sorts gvectors
*/
#include "wfc.h"



int sort_cmp(const void *a, const void *b);
static void Sorted_gvecs_idxs(const ND_int npw, ELPH_float * gvecs, ND_int * indices);



void Sort_pw(const ND_int npw_tot, const ND_int npw_loc, const ND_int * fft_dims,\
                const ELPH_float * restrict gvec_in, const ELPH_cmplx * wfc_in, \
                const ND_int nsets, ND_int * npw_loc_out, ND_int * nG_xy,
                int ** gvec_out, ELPH_cmplx ** wfc_out, MPI_Comm mpi_comm)
{
    /*
    This function sorts pws and wfcs
    */
    /*
    First gather wfc and gvecs

    return local number of pws after sorting
    */

    int my_rank, Comm_size, mpi_error;
    mpi_error = MPI_Comm_size(mpi_comm, &Comm_size);
    mpi_error = MPI_Comm_rank(mpi_comm, &my_rank);

    ELPH_float * all_Gvecs         = NULL;
    ELPH_cmplx * wfc_root_buf      = NULL;
    ELPH_cmplx * wfc_root_buf_sort = NULL;

    int * gvec_buf  = NULL;
    int * counts_recv = NULL;
    int * displacements = NULL;
    int * counts_send = NULL;
    int * disp_send = NULL;

    ND_int * start_xy     = NULL;
    ND_int * indices      = NULL;
    ND_int * gxy_loc_comm = NULL; // local number of gxy in each cpu 
    
    if (npw_loc != get_mpi_local_size_idx(npw_tot, NULL,mpi_comm)) error_msg("npws mismatch while sorting");

    if (my_rank == 0)
    {   
        all_Gvecs    = malloc(sizeof(ELPH_float)*npw_tot*3);
        gvec_buf     = malloc(sizeof(int)*npw_tot*3);
        counts_recv  = malloc(4*sizeof(int)*Comm_size);
        indices      = malloc(sizeof(ND_int)*npw_tot);
        start_xy     = malloc(sizeof(ND_int)*fft_dims[0]*fft_dims[1]);
        gxy_loc_comm = malloc(sizeof(ND_int)*Comm_size);

        if (counts_recv == NULL) error_msg("Failed to allocate comm array");
        if (all_Gvecs == NULL)   error_msg("Failed to allocate gvec array");
        if (gvec_buf == NULL)    error_msg("Failed to allocate gvec int array");
        if (indices == NULL)     error_msg("Failed to allocate indices array");
        if (start_xy == NULL)    error_msg("Failed to allocate xy array");
        if (gxy_loc_comm == NULL)    error_msg("Failed to allocate xycomm array");

        displacements = counts_recv+Comm_size ;
        counts_send = counts_recv+2*Comm_size ;
        disp_send   = counts_recv+3*Comm_size ;
    }

    int pw_loc_int = npw_loc;
    mpi_error = MPI_Gather(&pw_loc_int, 1, MPI_INT, counts_recv, 1, MPI_INT, 0, mpi_comm);


    if (my_rank == 0)
    {
        int disp_rectemp = 0;
        for (int i = 0 ; i<Comm_size; ++i)
        {
            counts_recv[i] *= 3;
            displacements[i] = disp_rectemp;
            disp_rectemp += counts_recv[i];
        }
    }

    // gather to the root node
    mpi_error = MPI_Gatherv(gvec_in, 3*npw_loc, ELPH_MPI_float,
                all_Gvecs, counts_recv, displacements, 
                ELPH_MPI_float, 0, mpi_comm);

    Sorted_gvecs_idxs(npw_tot, all_Gvecs, indices);

    ND_int  xycount;
    if (my_rank == 0)
    {   
        xycount = 0;
        // initialize with gvecs which are out of the box
        int Gx_prev = fft_dims[0]+10;
        int Gy_prev = fft_dims[1]+10;

        for (ND_int ig = 0; ig < npw_tot; ++ig)
        {
            ND_int idx_temp     = indices[ig];
            gvec_buf[3*ig]   = rint(all_Gvecs[3*idx_temp]);
            gvec_buf[3*ig+1] = rint(all_Gvecs[3*idx_temp+1]);
            gvec_buf[3*ig+2] = rint(all_Gvecs[3*idx_temp+2]);
            if (gvec_buf[3*ig+1] != Gy_prev || gvec_buf[3*ig] != Gx_prev)
            {   
                Gx_prev = gvec_buf[3*ig];
                Gy_prev = gvec_buf[3*ig+1];
                start_xy[xycount] = ig;
                ++xycount;
            }
        }
    }   

    mpi_error = MPI_Bcast(&xycount, 1, ELPH_MPI_ND_INT, 0, mpi_comm);

    *nG_xy = get_mpi_local_size_idx(xycount, NULL,mpi_comm);

    mpi_error = MPI_Gather(nG_xy, 1, ELPH_MPI_ND_INT, gxy_loc_comm, 1, ELPH_MPI_ND_INT, 0, mpi_comm);

    if (my_rank == 0)
    {
        /* calc the npw for each cpu */
        ND_int xy_disp = 0;
        for (int i = 0 ; i<Comm_size; ++i)
        {
            // this is number of pw per code
            ND_int xy_in_this_cpu = gxy_loc_comm[i];
            ND_int xy_prev = xy_disp;
            xy_disp += xy_in_this_cpu;
            if (xy_disp < xycount)
            {
                counts_send[i] = start_xy[xy_disp]-start_xy[xy_prev];
            }
            else counts_send[i] = npw_tot-start_xy[xy_prev];
            counts_send[i] *= 3;
            disp_send[i] = 3*start_xy[xy_prev]; 
        }
    }

    int npw_this_cpu;
    mpi_error = MPI_Scatter(counts_send, 1, MPI_INT, &npw_this_cpu, 1, MPI_INT, 0,mpi_comm);
    // note we send 3*npws, so we divide by 3 now
    npw_this_cpu = npw_this_cpu/3;
    *npw_loc_out = npw_this_cpu;

    *gvec_out = malloc(sizeof(int)*3*npw_this_cpu);
    // scatter gvecs
    mpi_error = MPI_Scatterv(gvec_buf, counts_send, disp_send, \
                MPI_INT, *gvec_out, 3*npw_this_cpu, MPI_INT, 0, mpi_comm);

    free(gxy_loc_comm);
    free(start_xy);
    free(all_Gvecs);
    free(gvec_buf);
    
    *wfc_out  = malloc(sizeof(ELPH_cmplx)*npw_this_cpu*nsets);

    if (my_rank ==0)
    {   
        wfc_root_buf = malloc(sizeof(ELPH_cmplx)*npw_tot);
        wfc_root_buf_sort = malloc(sizeof(ELPH_cmplx)*npw_tot);
        for (int i = 0 ; i<Comm_size; ++i)
        {
            counts_send[i]   /= 3;
            disp_send[i]     /= 3;
            displacements[i] /= 3;
            counts_recv[i]   /= 3;
        }
    }

    for (ND_int iset = 0; iset < nsets; ++iset)
    {
        const ELPH_cmplx * wfin  = wfc_in   + iset*npw_loc;
        ELPH_cmplx * wfout = *wfc_out + iset*npw_this_cpu;

        mpi_error = MPI_Gatherv(wfin, npw_loc, ELPH_MPI_cmplx,
                wfc_root_buf, counts_recv, displacements, 
                ELPH_MPI_cmplx, 0, mpi_comm);
        
        /* sort */
        if (my_rank ==0)
        {
            for (ND_int ig = 0; ig < npw_tot; ++ig)
            {
                ND_int idx_temp     = indices[ig];
                wfc_root_buf_sort[ig] = wfc_root_buf[idx_temp];
            }
        }

        mpi_error = MPI_Scatterv(wfc_root_buf_sort, counts_send, disp_send, \
                ELPH_MPI_cmplx, wfout, npw_this_cpu, ELPH_MPI_cmplx, 0, mpi_comm);
    }

    free(indices);
    free(wfc_root_buf);
    free(wfc_root_buf_sort);
    free(counts_recv);

}






static void Sorted_gvecs_idxs(const ND_int npw, ELPH_float * gvecs, ND_int * indices )
{
    /*
    Sort gvecs and group all gvecs with same (x,y).

    // gvecs in crystal coordinates
    */
    if (gvecs == NULL) return ; 

    ELPH_float ** gvec_ptrs = malloc(sizeof(ELPH_float *)*npw);

    /* fill the struct */
    for (ND_int i = 0 ; i< npw; ++i)
    {
        gvec_ptrs[i] = gvecs + 3*i; 
    }

    qsort(gvec_ptrs, npw, sizeof(ELPH_float *), sort_cmp);

    // store the sorted indices
    for (ND_int i = 0 ; i< npw; ++i) indices[i] = (gvec_ptrs[i]-gvecs)/3;
    
    free(gvec_ptrs);

}

int sort_cmp(const void *a, const void *b)
{   
    /*
    First sorts along y, then x and finally z
    */
    ELPH_float * v1 = *(ELPH_float **)a;
    ELPH_float * v2 = *(ELPH_float **)b;

    ND_int Gax, Gay, Gaz;
    ND_int Gbx, Gby, Gbz;
    Gax = rint(v1[0]);
    Gay = rint(v1[1]);
    Gaz = rint(v1[2]);

    Gbx = rint(v2[0]);
    Gby = rint(v2[1]);
    Gbz = rint(v2[2]);

    if (Gay == Gby)
    {
        if (Gax == Gbx) return (Gaz > Gbz)-(Gaz < Gbz);
        else            return (Gax > Gbx)-(Gax < Gbx);
    }
    else return (Gay > Gby)-(Gay < Gby);
}



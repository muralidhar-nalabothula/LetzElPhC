/*
This is a routine to perform 3D FFT in parallel
*/
#include "fft.h"


void fft3D(struct ELPH_fft_plan * plan, const ND_int nsets, \
            ELPH_cmplx * wfcr, ELPH_cmplx * wfcG, const bool conjugate)
{
    /*
    a) FFT along xy and  (Nx,Ny,Nz_loc)-> (Gx, Gy, Nz_loc)
    
    b) scatter Gx,Gy to different cpus to perform fft along x and y.
        *Note that we use scatter instead of alltoallv so that some FFTs 
        can be computed while being communicating
    
    c) Perform FFT along Z
    
    d) box to sphere
    */

    // first some basic stuff
    int myrank, ncpus, mpi_error;
    mpi_error = MPI_Comm_rank(plan->comm, &myrank);
    mpi_error = MPI_Comm_size(plan->comm, &ncpus);

    ND_int Nx = plan->fft_dims[0];
    ND_int Ny = plan->fft_dims[1];
    ND_int Nz = plan->fft_dims[2];

    ELPH_float norm = Nx*Ny*Nz ;
    norm = 1.0/norm ;

    ND_int fft_buf_size = Nx*Ny*plan->nzloc; 
    
    for (ND_int iset = 0 ; iset < nsets; ++iset)
    {   
        ELPH_cmplx * wfcr_tmp = wfcr + iset*fft_buf_size;
        ELPH_cmplx * wfcG_tmp = wfcG + iset*plan->ngvecs_loc ;
        // copy to the buffer
        memcpy(plan->fft_data,wfcr_tmp,sizeof(ELPH_cmplx)*fft_buf_size);
        // a) perform FFTs along xy slab
        ND_function(fft_execute_plan, Nd_cmplxS) (plan->fplan_xy);
        
        // b) (i) pack the data for transpose 
        for (ND_int ixy = 0 ; ixy < plan->nGxy; ++ixy)
        {   
            ND_int Gx = plan->Gxy_total[2*ixy];
            ND_int Gy = plan->Gxy_total[2*ixy+1];
            // [Gx,Gy,0] =  Gx*Ny*Nzloc + Gy*Nzloc 
            ELPH_cmplx * xy_buf = plan->fft_data + plan->nzloc*(Gy + Gx*Ny);
            // (Nxy,Zloc)
            memcpy(plan->nz_buf + ixy*plan->nzloc, xy_buf, sizeof(ELPH_cmplx)*plan->nzloc);
        }
        // b) (ii)  transpose the data
        fwd_transpose(plan); // nz_buf has (NGxy_loc,Nz)
        
        // c) perform fft along z
        ND_function(fft_execute_plan, Nd_cmplxS) (plan->fplan_z);

        // d) box to sphere
        ND_int igvec = 0 ;
        for (ND_int ixy = 0 ; ixy < plan->nGxyloc; ++ixy)
        {
            ELPH_cmplx * zfft_ptr = plan->nz_buf + ixy*Nz ;
            for (ND_int ig = 0; ig < plan->ngxy_z[ixy]; ++ig)
            {
                ND_int Gz = plan->gvecs[3*igvec+2];
                if (Gz < 0 ) Gz += Nz ;
                wfcG_tmp[igvec]= zfft_ptr[Gz]*norm;
                if (conjugate) wfcG_tmp[igvec] = conj(wfcG_tmp[igvec]);
                ++igvec;
            }        
        }
        if (igvec != plan->ngvecs_loc) error_msg("Gvec mismatch in fwd execute. ");
    }
}






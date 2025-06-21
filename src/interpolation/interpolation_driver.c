#include <mpi.h>
#include <stdlib.h>

#include "common/ELPH_timers.h"
#include "common/dtypes.h"
#include "common/error.h"
#include "common/init_dtypes.h"
#include "common/parallel.h"
#include "common/print_info.h"
#include "dvloc/dvloc.h"
#include "elphC.h"
#include "io/io.h"
#include "io/qe/qe_io.h"

void interpolation_driver(const char* ph_save, const char* ph_save_interpolated,
                          enum ELPH_dft_code dft_code, int nqpools,
                          const ND_int* qgrid_new, MPI_Comm comm_world)
{
    //
    init_ELPH_clocks();

    struct ELPH_MPI_Comms* mpi_comms = malloc(sizeof(struct ELPH_MPI_Comms));
    CHECK_ALLOC(mpi_comms);

    // since kpoints donot makes sence, we must always set it to 1
    create_parallel_comms(nqpools, 1, comm_world, mpi_comms);

    print_ELPH_logo(mpi_comms->commW_rank, stdout);
    print_info_msg(mpi_comms->commW_rank,
                   "********** Interpolation Program started **********");

    struct Lattice* lattice = malloc(sizeof(struct Lattice));
    CHECK_ALLOC(lattice);
    init_lattice_type(lattice);

    struct Phonon* phonon = malloc(sizeof(struct Phonon));
    CHECK_ALLOC(phonon);
    init_phonon_type(phonon);

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

    int World_rank_tmp = mpi_comms->commW_rank;

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

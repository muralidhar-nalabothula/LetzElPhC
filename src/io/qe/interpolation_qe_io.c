#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "common/cwalk/cwalk.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/parallel.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"
#include "io/io.h"
#include "qe_io.h"

void get_interpolation_data_from_qe(struct Lattice* lattice,
                                    struct Phonon* phonon,
                                    const char* ph_save_dir, ELPH_float** Zvals,
                                    ELPH_float* alat,
                                    const struct ELPH_MPI_Comms* Comm)
{
    /*
       This functions gets basic ground state data to run interpolation. We
       donot have yambo SAVE, so we need to get few more quantities than in
       get_data_from_qe.c

       allocated memory must be freed outside this function
       */
    int mpi_error;
    char** pseudo_pots = NULL;  // must initialize to NULL else U.B
                                // only root will have pseudo_pots.

    get_data_from_qe(lattice, phonon, ph_save_dir, &pseudo_pots, Comm);
    // we need to get few more
    // 1) lattice vectors, -> volume, reciprocal lattice vectors
    // 2) natom -> nmodes
    // 3) Zvals
    // 4) nspinor
    // 5) set nnfftz_loc and nfftz_loc_shift
    // 6) atomic positions (in cart units)
    // 7) alat lengths

    // check if there is any assume_isolated tag in xml file
    //
    char* tmp_buffer = NULL;
    size_t tmp_buffer_len = 1;

    if (Comm->commW_rank == 0)
    {
        // sanity check
        //
        if (!pseudo_pots)
        {
            error_msg("Pseudo pots not set in the master CPU.");
        }
        //
        tmp_buffer_len = 1024 + strlen(ph_save_dir);
        tmp_buffer = calloc(tmp_buffer_len, 1);
        CHECK_ALLOC(tmp_buffer);
        //
        cwk_path_join(ph_save_dir, "data-file-schema.xml", tmp_buffer,
                      tmp_buffer_len);
        FILE* fp = fopen(tmp_buffer, "r");
        if (!fp)
        {
            error_msg("Error opening data-file-schema.xml file");
        }

        ezxml_t qexml = ezxml_parse_fp(fp);
        if (qexml == NULL)
        {
            error_msg("Error parsing data-file-schema.xml file");
        }

        // read all the remaining required stuff

        ezxml_t xml_tmp;
        const char* tmp_xml_str;
        // get lattice vectors
        for (int ilat = 0; ilat < 3; ++ilat)
        {
            ELPH_float a_tmp_read[3];
            char latname[4];
            snprintf(latname, sizeof(latname), "a%d", (int)(ilat + 1));
            //
            xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", 0,
                                "cell", 0, latname, -1);
            if (!xml_tmp)
            {
                error_msg(
                    "Parsing lattice vector from data-file-schema.xml file");
            }
            tmp_xml_str = xml_tmp->txt;
            if (!tmp_xml_str ||
                parser_doubles_from_string(tmp_xml_str, a_tmp_read, 3) != 3)
            {
                error_msg("Error parsing vector vec from data-file-schema.xml");
            }
            for (int ix = 0; ix < 3; ++ix)
            {
                lattice->alat_vec[3 * ix + ilat] = a_tmp_read[ix];  // a[:,i]
            }
        }

        // get number of atoms
        xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", -1);
        if (!xml_tmp)
        {
            error_msg(
                "Parsing atomic_structure from data-file-schema.xml file");
        }

        tmp_xml_str = ezxml_attr(xml_tmp, "nat");
        if (!tmp_xml_str)
        {
            error_msg("error nat attribute from data-file-schema.xml file");
        }
        lattice->natom = atoi(tmp_xml_str);
        int natom = lattice->natom;

        tmp_xml_str = ezxml_attr(xml_tmp, "alat");
        if (!tmp_xml_str)
        {
            error_msg("error alat attribute from data-file-schema.xml file");
        }
        alat[0] = atof(tmp_xml_str);
        alat[1] = alat[0];
        alat[2] = alat[0];

        ezxml_t atom_specs =
            ezxml_get(qexml, "output", 0, "atomic_species", -1);
        if (atom_specs == NULL)
        {
            error_msg(
                "error reading atomic spices from data-file-schema.xml file");
        }
        tmp_xml_str = ezxml_attr(atom_specs, "ntyp");
        if (!tmp_xml_str)
        {
            error_msg("error ntyp attribute from data-file-schema.xml file");
        }
        int ntype = atoi(tmp_xml_str);

        *Zvals = malloc(sizeof(ELPH_float) * natom);
        CHECK_ALLOC(*Zvals);

        ELPH_float* Zval_type = malloc(ntype * sizeof(*Zval_type));
        CHECK_ALLOC(Zval_type);

        char* atom_type_symbol = calloc(ntype * 16, sizeof(*atom_type_symbol));
        CHECK_ALLOC(atom_type_symbol);

        // read atom pos, Zvals
        for (int it = 0; it < ntype; ++it)
        {
            xml_tmp = ezxml_get(atom_specs, "species", it, "");
            if (xml_tmp == NULL)
            {
                error_msg(
                    "error reading atomic spices from data-file-schema.xml "
                    "file");
            }
            tmp_xml_str = ezxml_attr(xml_tmp, "name");
            if (!tmp_xml_str)
            {
                error_msg(
                    "error name attribute from data-file-schema.xml file");
            }
            if (strlen(tmp_xml_str) > 15)
            {
                error_msg("Atomic Name is very long.");
            }
            while (isspace((unsigned char)(*tmp_xml_str)))
            {
                ++tmp_xml_str;
            }
            strncpy_custom(atom_type_symbol + 16 * it, tmp_xml_str, 16);
            //
            cwk_path_join(ph_save_dir, pseudo_pots[it], tmp_buffer,
                          tmp_buffer_len);
            char atomic_sym[4];
            get_upf_element(tmp_buffer, atomic_sym, Zval_type + it);
            free(pseudo_pots[it]);
        }
        free(pseudo_pots);
        pseudo_pots = NULL;

        // read atomic pos
        lattice->atomic_pos = malloc(3 * natom * sizeof(*lattice->atomic_pos));
        CHECK_ALLOC(lattice->atomic_pos);

        for (int ia = 0; ia < natom; ++ia)
        {
            xml_tmp = ezxml_get(qexml, "output", 0, "atomic_structure", 0,
                                "atomic_positions", 0, "atom", ia, "");
            if (!xml_tmp)
            {
                error_msg(
                    "Parsing atomic_positions from data-file-schema.xml file");
            }
            tmp_xml_str = ezxml_attr(xml_tmp, "name");
            if (!tmp_xml_str)
            {
                error_msg(
                    "error name attribute from data-file-schema.xml file");
            }
            while (isspace((unsigned char)(*tmp_xml_str)))
            {
                ++tmp_xml_str;
            }

            int itype = -1;
            for (int it = 0; it < ntype; ++it)
            {
                if (strstr(atom_type_symbol + 16 * it, tmp_xml_str))
                {
                    itype = it;
                    break;
                }
            }
            if (itype < 0)
            {
                error_msg("Cannot find atomic type.");
            }

            *(*Zvals + ia) = Zval_type[itype];

            tmp_xml_str = xml_tmp->txt;
            if (parser_doubles_from_string(
                    tmp_xml_str, lattice->atomic_pos + 3 * ia, 3) != 3)
            {
                error_msg("Error parsing atomic positions");
            }
        }

        lattice->nspinor = 1;
        // read nspinor
        xml_tmp =
            ezxml_get(qexml, "output", 0, "magnetization", 0, "noncolin", -1);
        if (!xml_tmp)
        {
            error_msg("Parsing noncolin from data-file-schema.xml file");
        }
        tmp_xml_str = xml_tmp->txt;
        strncpy_custom(tmp_buffer, tmp_xml_str, tmp_buffer_len);
        lowercase_str(tmp_buffer);

        if (strstr(tmp_buffer, "true"))
        {
            lattice->nspinor = 2;
        }

        free(tmp_buffer);
        free(Zval_type);
        free(atom_type_symbol);
        ezxml_free(qexml);
        fclose(fp);
    }

    mpi_error = MPI_Bcast(lattice->alat_vec, 9, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&lattice->natom, 1, MPI_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&lattice->nspinor, 1, MPI_INT, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(alat, 3, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    if (Comm->commW_rank)
    {
        lattice->atomic_pos =
            malloc(3 * lattice->natom * sizeof(*lattice->atomic_pos));
        CHECK_ALLOC(lattice->atomic_pos);

        *Zvals = malloc(sizeof(ELPH_float) * lattice->natom);
        CHECK_ALLOC(*Zvals);
    }

    mpi_error =
        MPI_Bcast(*Zvals, lattice->natom, ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(lattice->atomic_pos, 3 * lattice->natom,
                          ELPH_MPI_float, 0, Comm->commW);
    MPI_error_msg(mpi_error);

    lattice->nfftz_loc = get_mpi_local_size_idx(
        lattice->fft_dims[2], &(lattice->nfftz_loc_shift), Comm->commK);

    if (lattice->nfftz_loc < 1)
    {
        error_msg(
            "Some cpus do not contain plane waves. Over parallelization !.");
    }

    lattice->nmodes = 3 * lattice->natom;
    lattice->volume = fabs(det3x3(lattice->alat_vec));
    reciprocal_vecs(lattice->alat_vec, lattice->blat_vec);
}

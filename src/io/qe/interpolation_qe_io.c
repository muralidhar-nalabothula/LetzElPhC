#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#include "common/cwalk/cwalk.h"
#include "common/numerical_func.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"
#include "qe_io.h"

void get_interpolation_data_from_qe(struct Lattice* lattice,
                                    struct Phonon* phonon,
                                    const char* ph_save_dir, ELPH_float** Zvals,
                                    const struct ELPH_MPI_Comms* Comm)
{
    /*
       This functions gets basic ground state data to run interpolation. We
       donot have yambo SAVE, so we need to get few more quantities than in
       get_data_from_qe.c
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

    // check if there is any assume_isolated tag in xml file
    //
    char* tmp_buffer = NULL;
    size_t temp_str_len = 1;

    if (Comm->commW_rank == 0)
    {
        temp_str_len = 1024 + strlen(ph_save_dir);
        tmp_buffer = calloc(temp_str_len, 1);
        CHECK_ALLOC(tmp_buffer);
        //
        cwk_path_join(ph_save_dir, "data-file-schema.xml", tmp_buffer,
                      temp_str_len);
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
        const char* tmp_str;
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
            tmp_str = xml_tmp->txt;
            if (!tmp_str ||
                parser_doubles_from_string(tmp_str, a_tmp_read) != 3)
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

        tmp_str = ezxml_attr(xml_tmp, "nat");
        if (!tmp_str)
        {
            error_msg("error nat attribute from data-file-schema.xml file");
        }
        lattice->natom = atoi(tmp_str);

        ezxml_free(qexml);
        fclose(fp);
    }
    // Bcast
    // 1) alat_vec
    // 2 natom
    //

    lattice->nmodes = 3 * lattice->natom;
    lattice->volume = fabs(det3x3(lattice->alat_vec));
    reciprocal_vecs(lattice->alat_vec, lattice->blat_vec);
}

/*
This file contains function that parses data-file-schema.xml file
*/

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/constants.h"
#include "../../common/error.h"
#include "../../common/string_func.h"
#include "../../elphC.h"
#include "../ezxml/ezxml.h"
#include "qe_io.h"

#define ELPH_XML_READ_LINE_SIZE 1000

void parse_qexml(const char* xml_file, ND_int* natoms, ELPH_float* lat_vec,
                 ELPH_float* alat, char* dim, bool* is_soc_present,
                 ND_int* nmag, ND_int* fft_dims, ND_int* nph_sym,
                 ELPH_float** ph_sym_mats, ELPH_float** ph_sym_tau,
                 bool** ph_trevs, bool* ph_mag_symm, bool* ph_tim_rev,
                 char** pseudo_dir, char*** pseudo_pots)
{
    /*
    get pseudo potential information, fft information from xml file

    need to free pseudo_dir and pseudo_pots pointers
    */

    // set defaults
    *is_soc_present = false;
    *dim = '3';
    *nmag = 1;

    // check if there is any assume_isolated tag in xml file
    FILE* fp = fopen(xml_file, "r");
    if (!fp)
    {
        error_msg("Error opening data-file-schema.xml file");
    }

    bool assume_isolated_found = false;

    char* tmp_read = malloc(ELPH_XML_READ_LINE_SIZE);
    CHECK_ALLOC(tmp_read);

    while (fgets(tmp_read, ELPH_XML_READ_LINE_SIZE, fp))
    {
        if (strstr(tmp_read, "assume_isolated"))
        {
            assume_isolated_found =
                true;  // if (strstr(tmp_read, "2D")) *dim = '2';
            break;
        }
    }

    rewind(fp);

    char* tmp_str;  // tmp variable for ezxml char pointer

    ezxml_t qexml = ezxml_parse_fp(fp);
    if (qexml == NULL)
    {
        error_msg("Error parsing data-file-schema.xml file");
    }

    if (assume_isolated_found)
    {
        char* assume_iso = ezxml_get(qexml, "output", 0, "boundary_conditions",
                                     0, "assume_isolated", -1)
                               ->txt;
        if (string_start_with(assume_iso, "2D", true))
        {
            *dim = '2';
        }
    }
    // next get the pseudo pot directory
    tmp_str =
        ezxml_get(qexml, "input", 0, "control_variables", 0, "pseudo_dir", -1)
            ->txt;

    *pseudo_dir = malloc(strlen(tmp_str) +
                         1);  // we need to free this outside of this function
    CHECK_ALLOC(*pseudo_dir);

    strcpy(*pseudo_dir, tmp_str);
    // printf("pseudo dir : %s",*pseudo_dir);
    //  get ntypes
    ezxml_t atom_specs = ezxml_get(qexml, "output", 0, "atomic_species", -1);
    if (atom_specs == NULL)
    {
        error_msg("error reading atomic spices from data-file-schema.xml file");
    }

    ND_int ntype = atoll(ezxml_attr(atom_specs, "ntyp"));
    *pseudo_pots = malloc(sizeof(char*) * ntype);
    CHECK_ALLOC(*pseudo_pots);

    char** pot_tmp = *pseudo_pots;

    for (ND_int itype = 0; itype < ntype; ++itype)
    {
        tmp_str =
            ezxml_get(atom_specs, "species", itype, "pseudo_file", -1)->txt;
        // printf("%d : %s \n",(int)itype, tmp_str);
        pot_tmp[itype] = malloc(1 + strlen(tmp_str));
        CHECK_ALLOC(pot_tmp[itype]);

        strcpy(pot_tmp[itype], tmp_str);
    }

    // get number of atoms
    *natoms = atoll(ezxml_attr(
        ezxml_get(qexml, "output", 0, "atomic_structure", -1), "nat"));
    //
    // get alat
    alat[0] = atof(ezxml_attr(
        ezxml_get(qexml, "output", 0, "atomic_structure", -1), "alat"));
    alat[1] = alat[0];
    alat[2] = alat[0];
    // get fft dims
    fft_dims[0] = atoll(ezxml_attr(
        ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1), "nr1"));
    fft_dims[1] = atoll(ezxml_attr(
        ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1), "nr2"));
    fft_dims[2] = atoll(ezxml_attr(
        ezxml_get(qexml, "output", 0, "basis_set", 0, "fft_grid", -1), "nr3"));

    // get lattice vectors
    ELPH_float a_tmp_read[3];

    tmp_str = ezxml_get(qexml, "output", 0, "atomic_structure", 0, "cell", 0,
                        "a1", -1)
                  ->txt;
    if (parser_doubles_from_string(tmp_str, a_tmp_read) != 3)
    {
        error_msg("Error parsing a1 vec from data-file-schema.xml");
    }
    for (int ix = 0; ix < 3; ++ix)
    {
        lat_vec[3 * ix + 0] = a_tmp_read[ix];  // a[:,i]
    }

    tmp_str = ezxml_get(qexml, "output", 0, "atomic_structure", 0, "cell", 0,
                        "a2", -1)
                  ->txt;
    if (parser_doubles_from_string(tmp_str, a_tmp_read) != 3)
    {
        error_msg("Error parsing a2 vec from data-file-schema.xml");
    }
    for (int ix = 0; ix < 3; ++ix)
    {
        lat_vec[3 * ix + 1] = a_tmp_read[ix];
    }

    tmp_str = ezxml_get(qexml, "output", 0, "atomic_structure", 0, "cell", 0,
                        "a3", -1)
                  ->txt;
    if (parser_doubles_from_string(tmp_str, a_tmp_read) != 3)
    {
        error_msg("Error parsing a3 vec from data-file-schema.xml");
    }
    for (int ix = 0; ix < 3; ++ix)
    {
        lat_vec[3 * ix + 2] = a_tmp_read[ix];
    }

    // check if soc is present
    tmp_str =
        ezxml_get(qexml, "output", 0, "magnetization", 0, "spinorbit", -1)->txt;
    strcpy(tmp_read, tmp_str);
    lowercase_str(tmp_read);

    if (strstr(tmp_read, "true"))
    {
        *is_soc_present = true;
    }

    // nmag
    // first check if this is a lsda calc
    bool lsda = false;
    tmp_str =
        ezxml_get(qexml, "output", 0, "magnetization", 0, "lsda", -1)->txt;
    strcpy(tmp_read, tmp_str);
    lowercase_str(tmp_read);

    if (strstr(tmp_read, "true"))
    {
        lsda = true;
    }

    *ph_tim_rev = true;

    // check if there is noinv flag
    bool no_inv = false;
    tmp_str =
        ezxml_get(qexml, "input", 0, "symmetry_flags", 0, "noinv", -1)->txt;
    strcpy(tmp_read, tmp_str);
    lowercase_str(tmp_read);

    if (strstr(tmp_read, "true"))
    {
        no_inv = true;
    }

    if (no_inv)
    {
        *ph_tim_rev = false;
    }

    if (lsda)
    {
        *nmag = 2;
    }
    else
    {
        bool is_non_collinear = false;
        bool mag_system = false;

        tmp_str =
            ezxml_get(qexml, "output", 0, "magnetization", 0, "noncolin", -1)
                ->txt;
        strcpy(tmp_read, tmp_str);
        lowercase_str(tmp_read);

        if (strstr(tmp_read, "true"))
        {
            is_non_collinear = true;
        }

        if (is_non_collinear)
        {
            tmp_str = ezxml_get(qexml, "output", 0, "magnetization", 0,
                                "do_magnetization", -1)
                          ->txt;
            strcpy(tmp_read, tmp_str);
            lowercase_str(tmp_read);

            if (strstr(tmp_read, "true"))
            {
                mag_system = true;
            }
        }
        //
        if (mag_system)
        {
            *nmag = 4;
        }
        else
        {
            *nmag = 1;
        }

        if (mag_system)
        {
            // NM : This is can be changed later in this function
            *ph_tim_rev = false;
        }
    }

    // Finally read the phonon symmetries from the xml file
    // first get number of symmetries
    tmp_str = ezxml_get(qexml, "output", 0, "symmetries", 0, "nsym", -1)->txt;

    *nph_sym = atoi(tmp_str);

    // create a temp buffers for reading
    // we create 2*nph_sym sets of symmetries (factor 2 to store the time rev
    // case) The second half is used as temporary buffer
    // and also to store time_rev symmetries in nmag == 1 case
    // Note: the second half([nph_sym:]) are only used when nmag == 1
    // and these are not computed in this function. but we allocate
    // the storage
    //
    // in case of nmag == 2 or nmag == 4, we donot have a pure time reversal
    // symmetry but we can have a rotation + time rev
    //
    *ph_sym_mats = malloc(sizeof(ELPH_float) * 3 * 3 * 2 * (*nph_sym));
    CHECK_ALLOC(*ph_sym_mats);

    *ph_sym_tau = malloc(sizeof(ELPH_float) * 3 * 2 * (*nph_sym));
    CHECK_ALLOC(*ph_sym_tau);

    *ph_trevs = calloc(2 * (*nph_sym), sizeof(bool));
    CHECK_ALLOC(*ph_trevs);
    bool* trev_ptr = *ph_trevs;

    for (ND_int isym = 0; isym < 2 * (*nph_sym); ++isym)
    {
        // a  smart compiller will remove this loop.
        // but lets stick to standard and leave these to compilers
        trev_ptr[isym] = false;
    }

    bool inversion_sym = false;
    bool mag_sym_found = false;

    ELPH_float I3x3[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for (ND_int isym = 0; isym < *nph_sym; ++isym)
    {
        ezxml_t sym_xml_tmp = ezxml_get(qexml, "output", 0, "symmetries", -1);
        if (sym_xml_tmp == NULL)
        {
            error_msg("Error reading in phonon symmetries");
        }
        //
        // rotation matrix (stored in transposed order) and frac .translation
        // are store in crystals coordinates parse rotation (in crystal units)
        // S_cart = (alat@S_crys^T@blat^T)
        tmp_str = ezxml_get(sym_xml_tmp, "symmetry", isym, "rotation", -1)->txt;
        if (parser_doubles_from_string(tmp_str, (*ph_sym_mats) + 9 * isym) != 9)
        {
            error_msg(
                "Error parsing symmetry matrices from data-file-schema.xml");
        }

        // parse translation (in crystal units) tau_cart = alat@tau_crys
        // It should be noted that we use Sx + v convention, but qe uses S(x+v)
        // so our v = S*tau_qe
        tmp_str = ezxml_get(sym_xml_tmp, "symmetry", isym,
                            "fractional_translation", -1)
                      ->txt;
        if (parser_doubles_from_string(tmp_str, (*ph_sym_tau) + 3 * isym) != 3)
        {
            error_msg(
                "Error parsing frac. trans. vecs from data-file-schema.xml");
        }
        // check if this symmetry corresponds to time_reversal symmetry for nmag
        // == 4
        const char* trev_tmp_str =
            ezxml_attr(ezxml_get(sym_xml_tmp, "symmetry", isym, "info", -1),
                       "time_reversal");
        if (trev_tmp_str)
        {
            strcpy(tmp_read, trev_tmp_str);
            lowercase_str(tmp_read);

            if (strstr(tmp_read, "true"))
            {
                mag_sym_found = true;
                *ph_tim_rev = true;
                trev_ptr[isym] = true;
            }
            else
            {
                trev_ptr[isym] = false;
            }
        }
        //

        // check if inversion is present
        ELPH_float sum = 0;
        ELPH_float* sym_mat_tmpp = (*ph_sym_mats) + 9 * isym;
        for (int ix = 0; ix < 9; ++ix)
        {
            sum += fabs((I3x3[ix] + sym_mat_tmpp[ix]) *
                        (I3x3[ix] + sym_mat_tmpp[ix]));
        }
        sum = sqrt(sum);
        if (sum < ELPH_EPS)
        {
            inversion_sym = true;
        }
    }

    // incase magnetic symmetries found, we need to rearrange the symmetries and
    // push back
    if (mag_sym_found)
    {
        // send all timerev symmetries to end
        // first make a temporary copy
        ELPH_float* sym_src_ptr = *ph_sym_mats + (*nph_sym) * 9;
        ELPH_float* sym_des_ptr = *ph_sym_mats;
        memcpy(sym_src_ptr, sym_des_ptr, sizeof(ELPH_float) * 9 * (*nph_sym));

        ELPH_float* tau_src_ptr = *ph_sym_tau + (*nph_sym) * 3;
        ELPH_float* tau_des_ptr = *ph_sym_tau;
        memcpy(tau_src_ptr, tau_des_ptr, sizeof(ELPH_float) * 3 * (*nph_sym));

        bool* trev_src_ptr = *ph_trevs + (*nph_sym);
        bool* trev_des_ptr = *ph_trevs;
        memcpy(trev_src_ptr, trev_des_ptr, sizeof(bool) * (*nph_sym));

        ND_int t_sym_copiled = 0;
        ND_int n_sym_copiled = 0;

        if (*nph_sym % 2)
        {
            error_msg("Wrong number of symmetries");
        }
        for (ND_int isym = 0; isym < *nph_sym; ++isym)
        {
            ND_int dest_shift = 0;
            if (trev_src_ptr[isym])
            {
                dest_shift = (*nph_sym) / 2 + t_sym_copiled;
                ++t_sym_copiled;
            }
            else
            {
                dest_shift = n_sym_copiled;
                ++n_sym_copiled;
            }
            ELPH_float* sym_src_ptr_tmp = sym_src_ptr + 9 * isym;
            ELPH_float* sym_des_ptr_tmp = sym_des_ptr + 9 * dest_shift;
            memcpy(sym_des_ptr_tmp, sym_src_ptr_tmp, sizeof(ELPH_float) * 9);

            ELPH_float* tau_src_ptr_tmp = tau_src_ptr + 3 * isym;
            ELPH_float* tau_des_ptr_tmp = tau_des_ptr + 3 * dest_shift;
            memcpy(tau_des_ptr_tmp, tau_src_ptr_tmp, sizeof(ELPH_float) * 3);

            trev_des_ptr[dest_shift] = trev_src_ptr[isym];
        }
        if (t_sym_copiled != n_sym_copiled)
        {
            error_msg("Wrong number of symmetries");
        }
    }

    *ph_mag_symm = mag_sym_found;
    // set time reversal to false if inversion is present (no longer needed) for
    // non mangetic materials.
    if (inversion_sym && !mag_sym_found)
    {
        *ph_tim_rev = false;
    }

    free(tmp_read);
    ezxml_free(qexml);
    fclose(fp);
}

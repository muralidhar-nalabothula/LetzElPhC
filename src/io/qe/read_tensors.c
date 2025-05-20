// reads tensors born charges and dielectric tensor
// from tensors.xml file
//

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../common/dtypes.h"
#include "../../common/error.h"
#include "../../common/numerical_func.h"
#include "../../common/string_func.h"
#include "../../elphC.h"
#include "../ezxml/ezxml.h"

void read_ph_tensors_qe(const char* tensor_xml_file, const ND_int natom,
                        struct Phonon* phonon)
{
    phonon->Zborn = NULL;
    phonon->epsilon = NULL;
    phonon->Qpole = NULL;
    // for now no Quadruple as q.e does not support it.

    FILE* fp = fopen(tensor_xml_file, "r");
    if (NULL == fp)
    {
        // file not found or not readable.
        // so we skip without reading them.
        return;
    }

    ezxml_t tensor_xml = ezxml_parse_fp(fp);
    if (NULL == tensor_xml)
    {
        error_msg("parsing tensor.xml file failed\n");
    }

    char* tmp_str;
    // first check if tensors exist !
    //
    bool epsilon_exists = false;
    bool Zeu_exists = false;

    tmp_str =
        ezxml_get(tensor_xml, "EF_TENSORS", 0, "DONE_EFFECTIVE_CHARGE_EU", -1)
            ->txt;
    if ('t' == tolower(*tmp_str))
    {
        Zeu_exists = true;
    }

    tmp_str =
        ezxml_get(tensor_xml, "EF_TENSORS", 0, "DONE_ELECTRIC_FIELD", -1)->txt;
    if ('t' == tolower(*tmp_str))
    {
        epsilon_exists = true;
    }

    if (epsilon_exists && Zeu_exists)
    {
        phonon->Zborn = malloc(9 * natom * sizeof(*phonon->Zborn));
        CHECK_ALLOC(phonon->Zborn);

        phonon->epsilon = malloc(9 * sizeof(*phonon->epsilon));
        CHECK_ALLOC(phonon->epsilon);
        // read dielectric tensor
        tmp_str =
            ezxml_get(tensor_xml, "EF_TENSORS", 0, "DIELECTRIC_CONSTANT", -1)
                ->txt;

        if (parser_doubles_from_string(tmp_str, phonon->epsilon) != 9)
        {
            error_msg("Parsing epsilon failed");
        }
        // read effective charges
        tmp_str =
            ezxml_get(tensor_xml, "EF_TENSORS", 0, "EFFECTIVE_CHARGES_EU", -1)
                ->txt;
        if (parser_doubles_from_string(tmp_str, phonon->Zborn) != 9 * natom)
        {
            error_msg("Parsing Born charges failed");
        }

        // transpose epsilon and born charges as they are stored in transpose
        // order
        transpose3x3f_inplace(phonon->epsilon);

        for (ND_int ia = 0; ia < natom; ++ia)
        {
            transpose3x3f_inplace(phonon->Zborn + 9 * ia);
        }
    }
    else
    {
        fprintf(stderr,
                "Warning : Tensor file exists but Born charges or dielectric "
                "tensor missing !\n");
    }

    ezxml_free(tensor_xml);
    fclose(fp);
}

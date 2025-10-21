#include <ctype.h>
#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io.h"
#include "io/inih/ini.h"

#define READ_STR_LEN 600

static void Bcast_interpolation_input_data(
    struct interpolation_usr_input* input, int root, MPI_Comm comm);
static int interpolation_input_handler(void* user, const char* section,
                                       const char* name, const char* value);

// function to alloc, initiate interpolation_usr_input
void init_interpolation_usr_input(struct interpolation_usr_input** input)
{
    // this function also sets defaults for the user input file
    *input = malloc(sizeof(struct interpolation_usr_input));
    CHECK_ALLOC(*input);

    struct interpolation_usr_input* inp = *input;

    inp->ph_save_dir = calloc(READ_STR_LEN * 3, 1);
    CHECK_ALLOC(inp->ph_save_dir);

    inp->ph_save_interpolation_dir = inp->ph_save_dir + READ_STR_LEN;
    inp->qlist_file = inp->ph_save_dir + 2 * READ_STR_LEN;

    strncpy_custom(inp->ph_save_dir, "ph_save", READ_STR_LEN);
    strncpy_custom(inp->ph_save_interpolation_dir, "ph_save_interpolation",
                   READ_STR_LEN);
    strcpy(inp->qlist_file, "");
    // qlist_file will be set to NULL later if no list is given.

    // defaults
    inp->interpolate_dvscf = true;
    //
    inp->asr = true;
    strncpy_custom(inp->asr_kind, "simple", sizeof(inp->asr_kind));
    //
    inp->loto = false;
    inp->loto_dir[0] = 0.0;
    inp->loto_dir[1] = 0.0;
    inp->loto_dir[2] = 0.0;
    //
    inp->qgrid_fine[0] = 1;
    inp->qgrid_fine[1] = 1;
    inp->qgrid_fine[2] = 1;
}

// function to free interpolation_usr_input struct data
void free_interpolation_usr_input(struct interpolation_usr_input* input)
{
    free(input->ph_save_dir);
    free(input);
}

static void Bcast_interpolation_input_data(
    struct interpolation_usr_input* input, int root, MPI_Comm comm)
{
    int mpi_error;

    // all char * will be bcasted in one single call
    mpi_error =
        MPI_Bcast(input->ph_save_dir, READ_STR_LEN * 3, MPI_CHAR, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->interpolate_dvscf, 1, MPI_C_BOOL, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->asr, 1, MPI_C_BOOL, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(input->asr_kind, sizeof(input->asr_kind), MPI_CHAR,
                          root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(&input->loto, 1, MPI_C_BOOL, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(input->loto_dir, ARRAY_LEN(input->loto_dir),
                          ELPH_MPI_float, root, comm);
    MPI_error_msg(mpi_error);

    mpi_error = MPI_Bcast(input->qgrid_fine, ARRAY_LEN(input->qgrid_fine),
                          ELPH_MPI_ND_INT, root, comm);
    MPI_error_msg(mpi_error);
}

static int interpolation_input_handler(void* user, const char* section,
                                       const char* name, const char* value)
{
    UNUSED_VAR(section);
    // All the new inputs are added here.
    // Note : Set the defaults in init_interpolation_usr_input function.
    struct interpolation_usr_input* inp = user;

    // check if value is just an empty string
    size_t nospace_len = 0;
    size_t val_len = strlen(value);
    for (size_t i = 0; i < val_len; ++i)
    {
        if (value[i] != ' ')
        {
            ++nospace_len;
        }
    }

    if (nospace_len == 0)
    {
        printf("Invalid input for %s ", name);
        error_msg("Invalid input");
    }

    // add inputs from here. use else if
    if (strcmp(name, "asr") == 0)
    {
        inp->asr = parse_bool_input(value);
    }
    else if (strcmp(name, "asr_kind") == 0)
    {
        strncpy_custom(inp->asr_kind, value, sizeof(inp->asr_kind));
        lowercase_str(inp->asr_kind);
    }
    else if (strcmp(name, "interpolate_dvscf") == 0)
    {
        inp->interpolate_dvscf = parse_bool_input(value);
    }
    else if (strcmp(name, "loto") == 0)
    {
        inp->loto = parse_bool_input(value);
    }
    else if (strcmp(name, "loto_dir") == 0)
    {
        ND_int lt_parsed = parse_floats_from_string(value, inp->loto_dir,
                                                    ARRAY_LEN(inp->loto_dir));
        if (lt_parsed != 3)
        {
            error_msg(
                "Wrong number of floats given in loto_dir. Expected 4 floats");
        }
    }
    else if (strcmp(name, "ph_save_dir") == 0)
    {
        strncpy_custom(inp->ph_save_dir, value, READ_STR_LEN);
    }
    else if (strcmp(name, "ph_save_interpolation_dir") == 0)
    {
        strncpy_custom(inp->ph_save_interpolation_dir, value, READ_STR_LEN);
    }
    else if (strcmp(name, "nq1") == 0)
    {
        inp->qgrid_fine[0] = atoll(value);
    }
    else if (strcmp(name, "nq2") == 0)
    {
        inp->qgrid_fine[1] = atoll(value);
    }
    else if (strcmp(name, "nq3") == 0)
    {
        inp->qgrid_fine[2] = atoll(value);
    }
    else if (strcmp(name, "qlist_file") == 0)
    {
        strncpy_custom(inp->qlist_file, value, READ_STR_LEN);
    }
    else
    {
        error_msg("Invalid variable in interpolation input file.");
    }
    return 1;
}

void read_interpolation_input_file(const char* input_file,
                                   struct interpolation_usr_input** input_data,
                                   MPI_Comm MPI_world_comm)
{
    // input_data must be free outside

    init_interpolation_usr_input(input_data);

    int mpi_world_rank, mpi_error;

    mpi_error = MPI_Comm_rank(MPI_world_comm, &mpi_world_rank);
    MPI_error_msg(mpi_error);

    if (mpi_world_rank == 0)
    {
        if (ini_parse(input_file, interpolation_input_handler, *input_data) < 0)
        {
            error_msg("Cannot open interpolation input file.");
        }
    }
    // broad cast;
    Bcast_interpolation_input_data(*input_data, 0, MPI_world_comm);
}

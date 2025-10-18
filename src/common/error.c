/**
 * @file
 * @brief Error handling and reporting functions
 *
 * Provides error reporting functionality with MPI support, including
 * formatted error messages with file location and automatic program
 * termination.
 */

#include "error.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Reports an error message and aborts MPI execution
 *
 * Prints formatted error information to stderr including the source location
 * (file, function, line number) and error message, then calls MPI_Abort to
 * terminate all MPI processes.
 *
 * @param error_msg Descriptive error message
 * @param file Source file where error occurred (typically __FILE__)
 * @param line Line number where error occurred (typically __LINE__)
 * @param func_name Function name where error occurred (typically __func__)
 *
 * @note This function does not return - it terminates the program via MPI_Abort
 * @note All MPI ranks will print the error, but abort is called on COMM_WORLD
 */
void elph_error_msg(const char* error_msg, const char* file,
                    const long long int line, const char* func_name)
{
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    fprintf(stderr, "*************************************\n");
    fprintf(stderr,
            "# [ Error !!!] :  File : %s, in function : %s at line : %lld \n"
            "Error msg : %s \n",
            file, func_name, line, error_msg);
    fprintf(stderr, "*************************************\n");
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
}

/**
 * @brief Reports an MPI error and aborts execution
 *
 * Converts an MPI error code to a human-readable string using MPI_Error_string,
 * then calls elph_error_msg to report and terminate.
 *
 * @param err_code MPI error code returned by MPI functions
 * @param file Source file where error occurred (typically __FILE__)
 * @param line Line number where error occurred (typically __LINE__)
 * @param func_name Function name where error occurred (typically __func__)
 *
 * @note This function does not return - it terminates the program
 */
void ELPH_MPI_error_msg(int err_code, const char* file,
                        const long long int line, const char* func_name)
{
    char message[MPI_MAX_ERROR_STRING + 1];
    int resultlen;
    MPI_Error_string(err_code, message, &resultlen);
    message[resultlen] = '\0';
    elph_error_msg(message, file, line, func_name);
}

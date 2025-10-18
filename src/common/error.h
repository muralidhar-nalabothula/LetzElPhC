/**
 * @file
 * @brief Error handling macros and function declarations
 *
 * Provides error codes, overflow checking, allocation verification,
 * and MPI/NetCDF error handling macros with automatic source location tracking.
 */

#pragma once

/**
 * @def ERR_DIR_DOES_NOT_EXIST
 * @brief Error code: Directory does not exist
 */
#define ERR_DIR_DOES_NOT_EXIST 1

/**
 * @def ERR_NOT_A_DIRECTORY
 * @brief Error code: Path exists but is not a directory
 */
#define ERR_NOT_A_DIRECTORY -1

/**
 * @def ERR_FILE_OPEN_READ
 * @brief Error code: Cannot open file for reading
 */
#define ERR_FILE_OPEN_READ 1

/**
 * @def ERR_FILE_OPEN_WRITE
 * @brief Error code: Cannot open file for writing
 */
#define ERR_FILE_OPEN_WRITE 2

/**
 * @def ERR_FILE_COPY_FAIL
 * @brief Error code: File copy operation failed
 */
#define ERR_FILE_COPY_FAIL 3

/**
 * @def CHECK_OVERFLOW_ERROR(var, max_var_limit)
 * @brief Checks if a variable exceeds maximum limit and aborts if true
 *
 * Verifies that var < max_var_limit. If the condition fails, reports an
 * overflow error with variable names and terminates execution.
 *
 * @param var Variable to check
 * @param max_var_limit Maximum allowed value (exclusive)
 *
 * Example:
 * @code
 * CHECK_OVERFLOW_ERROR(array_size, INT_MAX);
 * @endcode
 */
#define CHECK_OVERFLOW_ERROR(var, max_var_limit)                        \
    {                                                                   \
        if ((var) >= (max_var_limit))                                   \
        {                                                               \
            elph_error_msg("Over flow of variable " #var                \
                           " over the provided limit, " #max_var_limit, \
                           __FILE__, __LINE__, __func__);               \
        }                                                               \
    }

/**
 * @def MPI_error_msg(err_code)
 * @brief Checks MPI return code and reports error if not MPI_SUCCESS
 *
 * Verifies that an MPI function returned MPI_SUCCESS. If not, converts
 * the error code to a message and terminates execution.
 *
 * @param err_code MPI error code to check
 *
 * Example:
 * @code
 * int err = MPI_Send(...);
 * MPI_error_msg(err);
 * @endcode
 */
#define MPI_error_msg(err_code)                                           \
    {                                                                     \
        if ((err_code) != MPI_SUCCESS)                                    \
        {                                                                 \
            ELPH_MPI_error_msg((err_code), __FILE__, __LINE__, __func__); \
        }                                                                 \
    }

/**
 * @def error_msg(print_str)
 * @brief Reports a custom error message and terminates execution
 *
 * Prints the provided error message with source location and aborts via
 * MPI_Abort.
 *
 * @param print_str Error message string
 *
 * Example:
 * @code
 * if (invalid_input)
 *     error_msg("Invalid input parameter detected");
 * @endcode
 */
#define error_msg(print_str) \
    elph_error_msg((print_str), __FILE__, __LINE__, __func__);

/**
 * @def CHECK_ALLOC(ptr)
 * @brief Verifies memory allocation succeeded and aborts if NULL
 *
 * Checks if a pointer is non-NULL after allocation. If NULL, reports
 * allocation failure and terminates.
 *
 * @param ptr Pointer to check
 *
 * Example:
 * @code
 * double* array = malloc(n * sizeof(double));
 * CHECK_ALLOC(array);
 * @endcode
 */
#define CHECK_ALLOC(ptr)                                     \
    {                                                        \
        if (!(ptr))                                          \
        {                                                    \
            error_msg("Failed to allocate " #ptr " buffer"); \
        }                                                    \
    }

/**
 * @def ERR(e)
 * @brief Reports NetCDF error and terminates execution
 *
 * Converts NetCDF error code to string using nc_strerror and reports the error.
 *
 * @param e NetCDF error code
 *
 * Example:
 * @code
 * int status = nc_open(filename, NC_NOWRITE, &ncid);
 * if (status != NC_NOERR) ERR(status);
 * @endcode
 */
#define ERR(e)                                            \
    {                                                     \
        fprintf(stderr, "Error: %s\n", nc_strerror((e))); \
        error_msg("netcdf_error");                        \
    }

/**
 * @brief Reports an error message and aborts MPI execution
 *
 * @param error_msg Descriptive error message
 * @param file Source file where error occurred
 * @param line Line number where error occurred
 * @param func_name Function name where error occurred
 */
void elph_error_msg(const char* error_msg, const char* file,
                    const long long int line, const char* func_name);

/**
 * @brief Reports an MPI error and aborts execution
 *
 * @param err_code MPI error code
 * @param file Source file where error occurred
 * @param line Line number where error occurred
 * @param func_name Function name where error occurred
 */
void ELPH_MPI_error_msg(int err_code, const char* file,
                        const long long int line, const char* func_name);

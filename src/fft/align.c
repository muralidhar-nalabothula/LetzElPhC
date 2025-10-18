/**
 * @file
 * @brief SIMD alignment length determination for FFTW
 *
 * Determines the proper alignment length for SIMD operations with FFTW
 * by testing allocation alignment properties.
 */

#include <complex.h>
#include <fftw3.h>

#include "common/error.h"
#include "elphC.h"
#include "fft.h"

/**
 * @def ELPH_FFTW_SIMD_LEN
 * @brief Maximum SIMD alignment length to test
 */
#define ELPH_FFTW_SIMD_LEN 64

/**
 * @brief Determines SIMD alignment length for FFTW allocations
 *
 * Finds the alignment length in units of sizeof(ELPH_cmplx) by allocating
 * an FFTW-aligned buffer and testing pointer alignment at different offsets.
 * The alignment length is the smallest offset where the alignment matches
 * the base buffer alignment.
 *
 * Uses cached result on subsequent calls for efficiency.
 *
 * Algorithm:
 * 1. Allocate aligned buffer of ELPH_FFTW_SIMD_LEN complex elements
 * 2. Get alignment of base pointer
 * 3. Test alignment at each offset i = 1, 2, 3, ...
 * 4. First offset where alignment matches base is the SIMD length
 *
 * @return Alignment length in units of sizeof(ELPH_cmplx)
 *
 * @note Result is cached after first call (thread-safe for read after init)
 * @note Calls error_msg() if alignment cannot be determined
 */
ND_int alignment_len(void)
{
    /* get the aligment len in units of sizeof(ELPH_cmplx) bytes*/
    static ND_int align_len = -1;
    if (align_len > 0)
    {
        return align_len;
    }
    //
    align_len = ELPH_FFTW_SIMD_LEN + 10;
    ND_int temp_len = sizeof(ELPH_cmplx) * ELPH_FFTW_SIMD_LEN;
    ELPH_cmplx* temp = fftw_fun(malloc)(temp_len);
    CHECK_ALLOC(temp);
    // first get the alignment value of temp;
    int buf_alignment = fftw_fun(alignment_of)((void*)temp);
    for (int i = 1; i < ELPH_FFTW_SIMD_LEN; ++i)
    {
        ELPH_cmplx* ptr = temp + i;
        int ptr_align = fftw_fun(alignment_of)((void*)ptr);
        if (ptr_align == buf_alignment)
        {
            align_len = i;
            break;
        }
    }
    fftw_fun(free)(temp);
    if (align_len >= ELPH_FFTW_SIMD_LEN)
    {
        error_msg("Error determining simd len");
    }
    return align_len;
}

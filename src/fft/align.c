/* determines the simd aligment length */

#include <complex.h>
#include <fftw3.h>

#include "../common/error.h"
#include "../elphC.h"
#include "fft.h"

#define ELPH_FFTW_SIMD_LEN 64

ND_int alignment_len(void)
{
    /* get the aligment len in units of sizeof(ELPH_cmplx) bytes*/
    ND_int align_len = ELPH_FFTW_SIMD_LEN + 10;
    ND_int temp_len = 2 * sizeof(ELPH_float) * ELPH_FFTW_SIMD_LEN;
    // Complex numbers have same representation as two doubles
    // We are doing this to avoid strict aliasing violation when calling
    // fftw_alignment_of function

    ELPH_float* temp = fftw_fun(malloc)(temp_len);
    CHECK_ALLOC(temp);

    // first get the alignment value of temp;
    int buf_alignment = fftw_fun(alignment_of)(temp);

    for (int i = 1; i < ELPH_FFTW_SIMD_LEN; ++i)
    {
        ELPH_float* ptr = temp + 2 * i;
        int ptr_align = fftw_fun(alignment_of)(ptr);
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

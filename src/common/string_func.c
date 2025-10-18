/**
 * @file
 * @brief String manipulation utility functions
 *
 * Provides safe string operations including copying, case conversion,
 * numeric parsing, pattern matching, and character replacement.
 */

#include "string_func.h"

#include <ctype.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"

/**
 * @brief Safer version of strncpy with guaranteed null termination
 *
 * Copies up to (count-1) characters from src to dest and ensures
 * null termination. Unlike standard strncpy, always null-terminates
 * and doesn't pad with zeros.
 *
 * @param dest Destination buffer
 * @param src Source string
 * @param count Size of destination buffer (including space for null)
 * @return Pointer to dest
 *
 * @note If count is 1, only null terminator is written
 * @note If count is 0, no operation is performed
 */
char* strncpy_custom(char* dest, const char* src, size_t count)
{
    // this does strncpy(dest,src,n-1) and sets nth element as \0
    if (1 == count)
    {
        dest[0] = '\0';
    }
    else if (count > 1)
    {
        strncpy(dest, src, count - 1);
        dest[count - 1] = '\0';
    }
    return dest;
}

/**
 * @brief Converts all characters in string to lowercase
 *
 * Performs in-place lowercase conversion using locale-independent
 * character classification.
 *
 * @param str String to convert (modified in-place)
 */
void lowercase_str(char* str)
{
    // lower case all the chars in a string
    for (char* p = str; *p; ++p)
    {
        *p = tolower((unsigned char)(*p));
    }
}

/**
 * @brief Extracts all numeric values from a string
 *
 * Parses string for floating-point numbers (including signed values).
 * Can be called twice: first with out=NULL to count numbers, then
 * with allocated array to extract values.
 *
 * Recognizes numbers starting with digits or +/- followed by digit.
 *
 * @param str Input string to parse
 * @param out Output array for extracted values (can be NULL to count only)
 * @return Number of floating-point values found
 *
 * @note Uses strtod for parsing, which handles various formats (e.g., 1.5e-3)
 */
ND_int parser_doubles_from_string(const char* str, ELPH_float* out)
{
    const char* p = str;
    char* q;
    double temp_val;
    ND_int count = 0;
    while (*p)
    {
        if (isdigit((unsigned char)(*p)) ||
            ((*p == '-' || *p == '+') && isdigit((unsigned char)(*(p + 1)))))
        {
            temp_val = strtod(p, &q);
            if (p == q)
            {
                break;
            }
            else
            {
                p = q;
            }
            if (out != NULL)
            {
                out[count] = temp_val;
            }
            ++count;
        }
        else
        {
            ++p;
        }
    }
    return count;
}

/**
 * @brief Checks if string starts with given substring
 *
 * @param str String to check
 * @param compare_str Substring to match at beginning
 * @param trim If true, ignore leading whitespace in both strings
 * @return true if str starts with compare_str, false otherwise
 */
bool string_start_with(char* str, char* compare_str, bool trim)
{
    char* a;
    char* b;
    a = str;
    b = compare_str;
    if (trim)
    {
        while (isspace((unsigned char)(*a)))
        {
            ++a;
        }
        while (isspace((unsigned char)(*b)))
        {
            ++b;
        }
    }
    if (b[0] != a[0])
    {
        return false;
    }
    return !strncmp(a, b, strlen(b));
}

/**
 * @brief Checks if string ends with given substring
 *
 * Uses temporary buffer to reverse both strings, then checks if
 * reversed str starts with reversed compare_str.
 *
 * @param str String to check
 * @param compare_str Substring to match at end
 * @param trim If true, ignore trailing whitespace in both strings
 * @return true if str ends with compare_str, false otherwise
 *
 * @note Allocates temporary buffer internally (freed before return)
 */
bool string_end_with(char* str, char* compare_str, bool trim)
{
    char* temp_str =
        malloc(sizeof(char) * (strlen(str) + strlen(compare_str) + 2));
    CHECK_ALLOC(temp_str);
    char* a = temp_str;
    char* b = temp_str + strlen(str) + 1;
    strcpy(a, str);
    strcpy(b, compare_str);
    str_reverse_in_place(a);
    str_reverse_in_place(b);
    if (trim)
    {
        while (isspace((unsigned char)(*a)))
        {
            ++a;
        }
        while (isspace((unsigned char)(*b)))
        {
            ++b;
        }
    }
    bool ret_value = !strncmp(a, b, strlen(b));
    if (b[0] != a[0])
    {
        ret_value = false;
    }
    free(temp_str);
    return ret_value;
}

/**
 * @brief Reverses a string in-place
 *
 * Swaps characters from both ends moving toward center.
 *
 * @param str String to reverse (modified in-place)
 * @return Pointer to str
 */
char* str_reverse_in_place(char* str)
{
    ND_int len = strlen(str);
    if (len == 0)
    {
        return str;
    }
    char* p1 = str;
    char* p2 = str + len - 1;
    while (p1 < p2)
    {
        char tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
    return str;
}

/**
 * @brief Replaces delimiter characters with replacement characters
 *
 * Scans string for any character in delimiters and replaces with
 * corresponding character in replace_chars.
 *
 * @param str_in String to modify (modified in-place)
 * @param delimters Characters to search for
 * @param replace_chars Replacement characters (same length as delimters)
 *
 * @warning If strlen(delimters) != strlen(replace_chars), buffer overflow may
 * occur
 */
void str_replace_chars(char* str_in, const char* delimters,
                       const char* replace_chars)
{
    ND_int ndelimters = strlen(delimters);
    // if  ndelimters != strlen(replace_chars) buffer overflow
    ND_int str_in_len = strlen(str_in);
    for (ND_int i = 0; i < str_in_len; ++i)
    {
        for (ND_int j = 0; j < ndelimters; ++j)
        {
            if (str_in[i] == delimters[j])
            {
                str_in[i] = replace_chars[j];
                break;
            }
        }
    }
}

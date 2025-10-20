#include "string_func.h"

#include <ctype.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "error.h"
/*
This file contains some useful string functions
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

void lowercase_str(char* str)
{
    // lower case all the chars in a string
    for (char* p = str; *p; ++p)
    {
        *p = tolower((unsigned char)(*p));
    }
}

ND_int parse_floats_from_string(const char* str, ELPH_float* out,
                                ND_int out_size)
{
    /*
    Extract atmost out_size float values from given string

    out_size is size of the out buffer. In case more floats
    found than the buffer size, only the first out_size floats
    are written.

    returns  number of floats found (this can be different than out_size)
    when number of floats in the string are different than the buffer size.
    if out == NULL, it returns number of floats found in the string.
    */
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

            if (out != NULL && count < out_size)
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

bool string_start_with(char* str, char* compare_str, bool trim)
{
    /*
    Check if given string starts with a substring
    */
    if (str == NULL || compare_str == NULL)
    {
        return false;
    }
    //
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

bool string_end_with(char* str, char* compare_str, bool trim)
{
    /*
    Check if given string ends with a substring
    */
    if (str == NULL || compare_str == NULL)
    {
        return false;
    }
    //
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

char* str_reverse_in_place(char* str)
{
    /*
    Do a inplace reversing of string
    */
    if (str == NULL)
    {
        return NULL;
    }
    //
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

void str_replace_chars(char* str_in, const char* delimters,
                       const char* replace_chars)
{
    ND_int ndelimters = strlen(delimters);
    // if  ndelimters != strlen(replace_chars) buffer overflow
    if (ndelimters != strlen(replace_chars))
    {
        error_msg("Number of delimiters not equal to replace_chars.");
        return;
    }

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

/**
 * @file
 * @brief String manipulation utility function declarations
 */

#pragma once
#include <stdbool.h>
#include <stddef.h>

#include "elphC.h"

/**
 * @brief Safer strncpy with guaranteed null termination
 *
 * @param dest Destination buffer
 * @param src Source string
 * @param count Buffer size (including null terminator)
 * @return Pointer to dest
 */
char* strncpy_custom(char* dest, const char* src, size_t count);

/**
 * @brief Converts string to lowercase in-place
 *
 * @param str String to convert
 */
void lowercase_str(char* str);

/**
 * @brief Extracts floating-point numbers from string
 *
 * @param str Input string
 * @param out Output array (NULL to count only)
 * @return Number of values found
 */
ND_int parser_doubles_from_string(const char* str, ELPH_float* out);

/**
 * @brief Checks if string starts with substring
 *
 * @param str String to check
 * @param compare_str Substring to match
 * @param trim Ignore leading whitespace
 * @return true if match found
 */
bool string_start_with(char* str, char* compare_str, bool trim);

/**
 * @brief Checks if string ends with substring
 *
 * @param str String to check
 * @param compare_str Substring to match
 * @param trim Ignore trailing whitespace
 * @return true if match found
 */
bool string_end_with(char* str, char* compare_str, bool trim);

/**
 * @brief Reverses string in-place
 *
 * @param str String to reverse
 * @return Pointer to str
 */
char* str_reverse_in_place(char* str);

/**
 * @brief Replaces delimiter characters with replacements
 *
 * @param str_in String to modify
 * @param delimters Characters to find
 * @param replace_chars Replacement characters
 */
void str_replace_chars(char* str_in, const char* delimters,
                       const char* replace_chars);

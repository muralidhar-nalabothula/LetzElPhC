/**
 * @file
 * @brief Lightweight XML parsing library
 *
 * ezxml is a C library for parsing XML documents. It provides a simple
 * DOM-like interface for reading, manipulating, and writing XML data.
 * The parser modifies the input string in-place for efficiency.
 *
 * Copyright 2004-2006 Aaron Voisine <aaron@voisine.org>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _EZXML_H
#define _EZXML_H

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @def EZXML_BUFSIZE
 * @brief Size of internal memory buffers for reallocation
 */
#define EZXML_BUFSIZE 1024

/**
 * @def EZXML_NAMEM
 * @brief Flag: tag name is malloc'ed and must be freed
 */
#define EZXML_NAMEM 0x80

/**
 * @def EZXML_TXTM
 * @brief Flag: text content is malloc'ed and must be freed
 */
#define EZXML_TXTM 0x40

/**
 * @def EZXML_DUP
 * @brief Flag: attribute name and value are strdup'ed
 */
#define EZXML_DUP 0x20

    /**
     * @typedef ezxml_t
     * @brief Pointer to an ezxml structure
     */
    typedef struct ezxml* ezxml_t;

    /**
     * @struct ezxml
     * @brief XML element node structure
     *
     * Represents a single XML element with its attributes, text content,
     * and relationships to other elements in the document tree.
     */
    struct ezxml
    {
        char* name;  /**< Tag name */
        char** attr; /**< Tag attributes as array: name, value, name, value,
                        ..., NULL */
        char* txt;   /**< Tag character content (empty string if none) */
        size_t off;  /**< Offset from start of parent tag's character content */
        ezxml_t next;    /**< Next tag with same name at same depth */
        ezxml_t sibling; /**< Next tag with different name at same depth */
        ezxml_t ordered; /**< Next tag at same depth in document order */
        ezxml_t child;  /**< First child tag (one level deeper), NULL if none */
        ezxml_t parent; /**< Parent tag, NULL if this is root */
        short flags;    /**< Memory management flags (EZXML_NAMEM, EZXML_TXTM,
                           EZXML_DUP) */
    };

    /**
     * @brief Parses an XML string and creates an ezxml structure
     *
     * For efficiency, modifies the input string by adding null terminators
     * and decoding ampersand sequences. If you don't want the input modified,
     * pass a copy.
     *
     * @param s XML string to parse (will be modified in-place)
     * @param len Length of XML string
     * @return Root element of parsed XML tree, or NULL on failure
     *
     * @note The input string is modified during parsing
     * @note Use ezxml_error() to retrieve error messages on failure
     */
    ezxml_t ezxml_parse_str(char* s, size_t len);

    /**
     * @brief Parses XML from a file stream
     *
     * Reads entire stream into memory and parses it. Wrapper for
     * ezxml_parse_str().
     *
     * @param fp File stream to read from
     * @return Root element of parsed XML tree, or NULL on failure
     *
     * @note For XML files, consider using ezxml_parse_file() or
     * ezxml_parse_fd()
     */
    ezxml_t ezxml_parse_fp(FILE* fp);

    /**
     * @brief Returns first child tag with given name
     *
     * @param xml Parent element
     * @param name Tag name to search for
     * @return First child with matching name, or NULL if not found
     */
    ezxml_t ezxml_child(ezxml_t xml, const char* name);

/**
 * @def ezxml_next(xml)
 * @brief Returns next tag with same name at same depth
 *
 * @param xml Current element
 * @return Next sibling with same tag name, or NULL if none
 */
#define ezxml_next(xml) ((xml) ? xml->next : NULL)

    /**
     * @brief Returns Nth tag with same name at same depth
     *
     * @param xml Starting element
     * @param idx Index (0 returns xml itself, 1 returns next, etc.)
     * @return Nth element with same name, or NULL if not found
     */
    ezxml_t ezxml_idx(ezxml_t xml, int idx);

/**
 * @def ezxml_name(xml)
 * @brief Returns tag name of given element
 *
 * @param xml Element
 * @return Tag name, or NULL if xml is NULL
 */
#define ezxml_name(xml) ((xml) ? xml->name : NULL)

/**
 * @def ezxml_txt(xml)
 * @brief Returns text content of given element
 *
 * @param xml Element
 * @return Text content, or empty string if none
 */
#define ezxml_txt(xml) ((xml) ? xml->txt : "")

    /**
     * @brief Returns value of requested tag attribute
     *
     * @param xml Element to query
     * @param attr Attribute name
     * @return Attribute value, or NULL if not found
     */
    const char* ezxml_attr(ezxml_t xml, const char* attr);

    /**
     * @brief Traverses XML tree to retrieve specific subtag
     *
     * Takes variable-length list of tag names and indices. Argument list
     * must be terminated by either index -1 or empty string tag name.
     *
     * Example:
     * @code
     * title = ezxml_get(library, "shelf", 0, "book", 2, "title", -1);
     * @endcode
     * This retrieves title of 3rd book on 1st shelf of library.
     *
     * @param xml Starting element
     * @param ... Variable arguments: tag name (const char*), index (int), ...
     * @return Found element, or NULL if not found
     */
    ezxml_t ezxml_get(ezxml_t xml, ...);

    /**
     * @brief Converts ezxml structure back to XML string
     *
     * @param xml Root element to convert
     * @return XML string that must be freed by caller
     */
    char* ezxml_toxml(ezxml_t xml);

    /**
     * @brief Returns processing instructions for given target
     *
     * @param xml Element (typically root)
     * @param target PI target name
     * @return NULL-terminated array of PI strings
     */
    const char** ezxml_pi(ezxml_t xml, const char* target);

    /**
     * @brief Frees memory allocated for ezxml structure
     *
     * Recursively frees element and all its children.
     *
     * @param xml Element to free
     */
    void ezxml_free(ezxml_t xml);

    /**
     * @brief Returns parser error message
     *
     * @param xml Any element from parsed tree
     * @return Error message, or empty string if no error
     */
    const char* ezxml_error(ezxml_t xml);

    /**
     * @brief Creates new empty ezxml structure with given root tag name
     *
     * @param name Root tag name (not copied - must remain valid)
     * @return New root element
     */
    ezxml_t ezxml_new(const char* name);

/**
 * @def ezxml_new_d(name)
 * @brief Wrapper for ezxml_new() that strdup's name
 *
 * @param name Root tag name (will be duplicated)
 * @return New root element with duplicated name
 */
#define ezxml_new_d(name) ezxml_set_flag(ezxml_new(strdup(name)), EZXML_NAMEM)

    /**
     * @brief Adds child tag to element
     *
     * @param xml Parent element
     * @param name Child tag name (not copied - must remain valid)
     * @param off Offset from start of parent's character content
     * @return Newly created child element
     */
    ezxml_t ezxml_add_child(ezxml_t xml, const char* name, size_t off);

/**
 * @def ezxml_add_child_d(xml, name, off)
 * @brief Wrapper for ezxml_add_child() that strdup's name
 *
 * @param xml Parent element
 * @param name Child tag name (will be duplicated)
 * @param off Offset from start of parent's character content
 * @return Newly created child element with duplicated name
 */
#define ezxml_add_child_d(xml, name, off) \
    ezxml_set_flag(ezxml_add_child(xml, strdup(name), off), EZXML_NAMEM)

    /**
     * @brief Sets character content for given tag
     *
     * @param xml Element
     * @param txt Text content (not copied - must remain valid)
     * @return The element (xml parameter)
     */
    ezxml_t ezxml_set_txt(ezxml_t xml, const char* txt);

/**
 * @def ezxml_set_txt_d(xml, txt)
 * @brief Wrapper for ezxml_set_txt() that strdup's txt
 *
 * @param xml Element
 * @param txt Text content (will be duplicated)
 * @return The element with duplicated text
 */
#define ezxml_set_txt_d(xml, txt) \
    ezxml_set_flag(ezxml_set_txt(xml, strdup(txt)), EZXML_TXTM)

    /**
     * @brief Sets or adds tag attribute
     *
     * If attribute exists, updates its value. If not found, adds new attribute.
     * A value of NULL removes the attribute.
     *
     * @param xml Element
     * @param name Attribute name (not copied - must remain valid)
     * @param value Attribute value (not copied), or NULL to remove attribute
     * @return The element (xml parameter)
     */
    ezxml_t ezxml_set_attr(ezxml_t xml, const char* name, const char* value);

/**
 * @def ezxml_set_attr_d(xml, name, value)
 * @brief Wrapper for ezxml_set_attr() that strdup's name and value
 *
 * @param xml Element
 * @param name Attribute name (will be duplicated)
 * @param value Attribute value (will be duplicated, cannot be NULL)
 * @return The element with duplicated attribute
 */
#define ezxml_set_attr_d(xml, name, value) \
    ezxml_set_attr(ezxml_set_flag(xml, EZXML_DUP), strdup(name), strdup(value))

    /**
     * @brief Sets flag for given tag
     *
     * @param xml Element
     * @param flag Flag value (EZXML_NAMEM, EZXML_TXTM, or EZXML_DUP)
     * @return The element (xml parameter)
     */
    ezxml_t ezxml_set_flag(ezxml_t xml, short flag);

    /**
     * @brief Removes tag and its subtags without freeing memory
     *
     * Detaches element from tree structure but does not free it.
     *
     * @param xml Element to remove
     * @return The removed element (xml parameter)
     */
    ezxml_t ezxml_cut(ezxml_t xml);

    /**
     * @brief Inserts existing tag into ezxml structure
     *
     * @param xml Element to insert
     * @param dest Destination parent element
     * @param off Offset from start of dest's character content
     * @return The inserted element (xml parameter)
     */
    ezxml_t ezxml_insert(ezxml_t xml, ezxml_t dest, size_t off);

/**
 * @def ezxml_move(xml, dest, off)
 * @brief Moves existing tag to become subtag of dest
 *
 * Combines ezxml_cut() and ezxml_insert().
 *
 * @param xml Element to move
 * @param dest Destination parent element
 * @param off Offset from start of dest's character content
 * @return The moved element
 */
#define ezxml_move(xml, dest, off) ezxml_insert(ezxml_cut(xml), dest, off)

/**
 * @def ezxml_remove(xml)
 * @brief Removes tag and all its subtags, freeing memory
 *
 * Combines ezxml_cut() and ezxml_free().
 *
 * @param xml Element to remove and free
 */
#define ezxml_remove(xml) ezxml_free(ezxml_cut(xml))

#ifdef __cplusplus
}
#endif

#endif  // _EZXML_H

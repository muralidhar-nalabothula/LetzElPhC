/**
 * @file ELPH_hmap.h
 * @brief Generic hash map implementation with string keys
 * 
 * This header provides a type-safe hash map implementation for C using macros.
 * Supports string keys and any value type. Uses separate chaining for collision
 * resolution and automatic resizing.
 * 
 * Copyright (c) 2014 rxi
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license. See LICENSE for details.
 */

#ifndef MAP_H
#define MAP_H

#include <string.h>

struct map_node_t;
typedef struct map_node_t map_node_t;

/**
 * @struct map_base_t
 * @brief Base structure for hash map implementation
 * 
 * Contains the bucket array and metadata about map size.
 */
typedef struct
{
    map_node_t **buckets; /**< Array of pointers to bucket chain heads */
    unsigned nbuckets;    /**< Number of buckets in the hash table */
    unsigned nnodes;      /**< Total number of key-value pairs stored */
} map_base_t;

/**
 * @struct map_iter_t
 * @brief Iterator structure for traversing hash map entries
 */
typedef struct
{
    unsigned bucketidx;   /**< Current bucket index being iterated */
    map_node_t *node;     /**< Current node in the iteration */
} map_iter_t;

/**
 * @def map_t(T)
 * @brief Defines a type-safe map structure for value type T
 * 
 * Creates a map structure that can hold values of type T with string keys.
 * 
 * @param T The value type to store in the map
 */
#define map_t(T)         \
    struct               \
    {                    \
        map_base_t base; \
        T *ref;          \
        T tmp;           \
    }

/**
 * @def map_init(m)
 * @brief Initializes a map structure to empty state
 * 
 * Must be called before using any map operations.
 * 
 * @param m Pointer to map structure to initialize
 */
#define map_init(m) memset(m, 0, sizeof(*(m)))

/**
 * @def map_deinit(m)
 * @brief Deinitializes a map and frees all allocated memory
 * 
 * Should be called when done with a map to prevent memory leaks.
 * 
 * @param m Pointer to map structure to deinitialize
 */
#define map_deinit(m) map_deinit_(&(m)->base)

/**
 * @def map_get(m, key)
 * @brief Retrieves a value from the map by key
 * 
 * @param m Pointer to map structure
 * @param key String key to look up
 * @return Pointer to value if found, NULL otherwise. Result stored in m->ref
 */
#define map_get(m, key) ((m)->ref = map_get_(&(m)->base, key))

/**
 * @def map_set(m, key, value)
 * @brief Sets or updates a key-value pair in the map
 * 
 * If key exists, updates the value. Otherwise creates a new entry.
 * 
 * @param m Pointer to map structure
 * @param key String key for the entry
 * @param value Value to store (copied into map)
 * @return 0 on success, -1 on allocation failure
 */
#define map_set(m, key, value) \
    ((m)->tmp = (value), map_set_(&(m)->base, key, &(m)->tmp, sizeof((m)->tmp)))

/**
 * @def map_remove(m, key)
 * @brief Removes a key-value pair from the map
 * 
 * @param m Pointer to map structure
 * @param key String key of entry to remove
 */
#define map_remove(m, key) map_remove_(&(m)->base, key)

/**
 * @def map_iter(m)
 * @brief Creates an iterator for traversing the map
 * 
 * @param m Pointer to map structure (unused in initialization)
 * @return map_iter_t Initialized iterator
 */
#define map_iter(m) map_iter_()

/**
 * @def map_next(m, iter)
 * @brief Advances iterator to next entry
 * 
 * @param m Pointer to map structure
 * @param iter Pointer to iterator
 * @return Key string of next entry, or NULL when iteration complete
 */
#define map_next(m, iter) map_next_(&(m)->base, iter)

/**
 * @brief Internal function to deinitialize map base structure
 * 
 * @param m Pointer to map base structure
 */
void map_deinit_(map_base_t *m);

/**
 * @brief Internal function to retrieve value by key
 * 
 * @param m Pointer to map base structure
 * @param key Key string to look up
 * @return void* Pointer to value, or NULL if not found
 */
void *map_get_(map_base_t *m, const char *key);

/**
 * @brief Internal function to set or update a key-value pair
 * 
 * @param m Pointer to map base structure
 * @param key Key string
 * @param value Pointer to value data
 * @param vsize Size of value in bytes
 * @return int 0 on success, -1 on failure
 */
int map_set_(map_base_t *m, const char *key, void *value, int vsize);

/**
 * @brief Internal function to remove a key-value pair
 * 
 * @param m Pointer to map base structure
 * @param key Key string to remove
 */
void map_remove_(map_base_t *m, const char *key);

/**
 * @brief Internal function to create a new iterator
 * 
 * @return map_iter_t Initialized iterator
 */
map_iter_t map_iter_(void);

/**
 * @brief Internal function to advance iterator to next entry
 * 
 * @param m Pointer to map base structure
 * @param iter Pointer to iterator
 * @return const char* Key string of next entry, or NULL when complete
 */
const char *map_next_(map_base_t *m, map_iter_t *iter);

/**
 * @typedef map_void_t
 * @brief Map type for void pointer values
 */
typedef map_t(void *) map_void_t;

/**
 * @typedef map_str_t
 * @brief Map type for string (char*) values
 */
typedef map_t(char *) map_str_t;

/**
 * @typedef map_int_t
 * @brief Map type for integer values
 */
typedef map_t(int) map_int_t;

/**
 * @typedef map_char_t
 * @brief Map type for character values
 */
typedef map_t(char) map_char_t;

/**
 * @typedef map_float_t
 * @brief Map type for single-precision floating-point values
 */
typedef map_t(float) map_float_t;

/**
 * @typedef map_double_t
 * @brief Map type for double-precision floating-point values
 */
typedef map_t(double) map_double_t;

#endif

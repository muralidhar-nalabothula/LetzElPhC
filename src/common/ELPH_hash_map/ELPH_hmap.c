/**
 * @file ELPH_hmap.c
 * @brief Hash map implementation for string keys and generic values
 * 
 * This file implements a hash map data structure with string keys and void pointer
 * values. It uses separate chaining for collision resolution and dynamic resizing.
 * 
 * Copyright (c) 2014 rxi
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license. See LICENSE for details.
 */

#include "ELPH_hmap.h"

#include <stdlib.h>
#include <string.h>

/**
 * @struct map_node_t
 * @brief Internal node structure for hash map entries
 * 
 * Stores hash value, pointer to value data, next node pointer for chaining,
 * followed by variable-length key and value data in memory.
 */
struct map_node_t
{
    unsigned hash;      /**< Cached hash value of the key */
    void *value;        /**< Pointer to the value data */
    map_node_t *next;   /**< Next node in the collision chain */
    /* char key[]; */   /**< Flexible array member: key string follows this structure */
    /* char value[]; */ /**< Flexible array member: value data follows key string */
};

/**
 * @brief Computes hash value for a string key using DJB2 algorithm
 * 
 * Uses the DJB2 hash function (hash * 33 + char) which provides good
 * distribution for string keys.
 * 
 * @param str Null-terminated string to hash
 * @return unsigned Hash value for the string
 */
static unsigned map_hash(const char *str)
{
    unsigned hash = 5381;
    while (*str)
    {
        hash = ((hash << 5) + hash) ^ *str++;
    }
    return hash;
}

/**
 * @brief Creates a new map node with given key-value pair
 * 
 * Allocates memory for node structure plus key string and value data.
 * Memory layout: [node_struct][key_string][padding][value_data]
 * 
 * @param key Null-terminated string key
 * @param value Pointer to value data to copy
 * @param vsize Size of value data in bytes
 * @return map_node_t* Pointer to newly created node, or NULL on allocation failure
 */
static map_node_t *map_newnode(const char *key, void *value, int vsize)
{
    map_node_t *node;
    int ksize = strlen(key) + 1;
    int voffset = ksize + ((sizeof(void *) - ksize) % sizeof(void *));
    node = malloc(sizeof(*node) + voffset + vsize);
    if (!node)
    {
        return NULL;
    }
    memcpy(node + 1, key, ksize);
    node->hash = map_hash(key);
    node->value = ((char *)(node + 1)) + voffset;
    memcpy(node->value, value, vsize);
    return node;
}

/**
 * @brief Computes bucket index for a given hash value
 * 
 * Uses bitwise AND for fast modulo operation. Assumes bucket count
 * is a power of 2.
 * 
 * @param m Pointer to map base structure
 * @param hash Hash value to map to bucket
 * @return int Bucket index (0 to nbuckets-1)
 * 
 * @note If implementation changes to allow non-power-of-2 bucket counts,
 *       this should use modulo operator instead of AND
 */
static int map_bucketidx(map_base_t *m, unsigned hash)
{
    return hash & (m->nbuckets - 1);
}

/**
 * @brief Adds a node to the appropriate bucket in the hash map
 * 
 * Inserts node at the head of the collision chain for its bucket.
 * 
 * @param m Pointer to map base structure
 * @param node Node to add to the map
 */
static void map_addnode(map_base_t *m, map_node_t *node)
{
    int n = map_bucketidx(m, node->hash);
    node->next = m->buckets[n];
    m->buckets[n] = node;
}

/**
 * @brief Resizes the hash map to a new bucket count
 * 
 * Chains all existing nodes together, reallocates bucket array,
 * then redistributes nodes across new buckets.
 * 
 * @param m Pointer to map base structure
 * @param nbuckets New number of buckets (should be power of 2)
 * @return int 0 on success, -1 on allocation failure
 */
static int map_resize(map_base_t *m, int nbuckets)
{
    map_node_t *nodes, *node, *next;
    map_node_t **buckets;
    int i;
    
    /* Chain all nodes together */
    nodes = NULL;
    i = m->nbuckets;
    while (i--)
    {
        node = (m->buckets)[i];
        while (node)
        {
            next = node->next;
            node->next = nodes;
            nodes = node;
            node = next;
        }
    }
    
    /* Reset buckets */
    buckets = realloc(m->buckets, sizeof(*m->buckets) * nbuckets);
    if (buckets != NULL)
    {
        m->buckets = buckets;
        m->nbuckets = nbuckets;
    }
    if (m->buckets)
    {
        memset(m->buckets, 0, sizeof(*m->buckets) * m->nbuckets);
        /* Re-add nodes to buckets */
        node = nodes;
        while (node)
        {
            next = node->next;
            map_addnode(m, node);
            node = next;
        }
    }
    
    /* Return error code if realloc() failed */
    return (buckets == NULL) ? -1 : 0;
}

/**
 * @brief Gets reference to the node pointer for a given key
 * 
 * Searches for a node with matching key and returns pointer to the
 * pointer that references it (either in bucket array or previous node's next).
 * 
 * @param m Pointer to map base structure
 * @param key Key string to search for
 * @return map_node_t** Pointer to node pointer, or NULL if key not found
 */
static map_node_t **map_getref(map_base_t *m, const char *key)
{
    unsigned hash = map_hash(key);
    map_node_t **next;
    if (m->nbuckets > 0)
    {
        next = &m->buckets[map_bucketidx(m, hash)];
        while (*next)
        {
            if ((*next)->hash == hash && !strcmp((char *)(*next + 1), key))
            {
                return next;
            }
            next = &(*next)->next;
        }
    }
    return NULL;
}

/**
 * @brief Deinitializes and frees all memory used by the hash map
 * 
 * Frees all nodes and their associated data, then frees the bucket array.
 * 
 * @param m Pointer to map base structure to deinitialize
 */
void map_deinit_(map_base_t *m)
{
    map_node_t *next, *node;
    int i;
    i = m->nbuckets;
    while (i--)
    {
        node = m->buckets[i];
        while (node)
        {
            next = node->next;
            free(node);
            node = next;
        }
    }
    free(m->buckets);
}

/**
 * @brief Retrieves the value associated with a key
 * 
 * @param m Pointer to map base structure
 * @param key Key string to look up
 * @return void* Pointer to value data, or NULL if key not found
 */
void *map_get_(map_base_t *m, const char *key)
{
    map_node_t **next = map_getref(m, key);
    return next ? (*next)->value : NULL;
}

/**
 * @brief Sets or updates a key-value pair in the hash map
 * 
 * If key exists, updates its value. Otherwise, creates a new node.
 * Automatically resizes map when load factor exceeds 1.0.
 * 
 * @param m Pointer to map base structure
 * @param key Key string for the entry
 * @param value Pointer to value data to store
 * @param vsize Size of value data in bytes
 * @return int 0 on success, -1 on allocation failure
 */
int map_set_(map_base_t *m, const char *key, void *value, int vsize)
{
    int n, err;
    map_node_t **next, *node;
    
    /* Find & replace existing node */
    next = map_getref(m, key);
    if (next)
    {
        memcpy((*next)->value, value, vsize);
        return 0;
    }
    
    /* Add new node */
    node = map_newnode(key, value, vsize);
    if (node == NULL)
    {
        goto fail;
    }
    if (m->nnodes >= m->nbuckets)
    {
        n = (m->nbuckets > 0) ? (m->nbuckets << 1) : 1;
        err = map_resize(m, n);
        if (err)
        {
            goto fail;
        }
    }
    map_addnode(m, node);
    m->nnodes++;
    return 0;
    
fail:
    if (node)
    {
        free(node);
    }
    return -1;
}

/**
 * @brief Removes a key-value pair from the hash map
 * 
 * If key exists, removes its node and frees associated memory.
 * Decrements node count on successful removal.
 * 
 * @param m Pointer to map base structure
 * @param key Key string of entry to remove
 */
void map_remove_(map_base_t *m, const char *key)
{
    map_node_t *node;
    map_node_t **next = map_getref(m, key);
    if (next)
    {
        node = *next;
        *next = (*next)->next;
        free(node);
        m->nnodes--;
    }
}

/**
 * @brief Initializes a new map iterator
 * 
 * @return map_iter_t Initialized iterator positioned before first element
 */
map_iter_t map_iter_(void)
{
    map_iter_t iter;
    iter.bucketidx = -1;
    iter.node = NULL;
    return iter;
}

/**
 * @brief Advances iterator to next entry and returns its key
 * 
 * Iterates through all entries in the map in arbitrary order.
 * Call repeatedly until it returns NULL to iterate all entries.
 * 
 * @param m Pointer to map base structure
 * @param iter Pointer to iterator structure
 * @return const char* Key string of next entry, or NULL when iteration complete
 */
const char *map_next_(map_base_t *m, map_iter_t *iter)
{
    if (iter->node)
    {
        iter->node = iter->node->next;
        if (iter->node == NULL)
        {
            goto nextBucket;
        }
    }
    else
    {
    nextBucket:
        do
        {
            if (++iter->bucketidx >= m->nbuckets)
            {
                return NULL;
            }
            iter->node = m->buckets[iter->bucketidx];
        } while (iter->node == NULL);
    }
    return (char *)(iter->node + 1);
}

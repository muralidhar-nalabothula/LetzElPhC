/**
 * @file
 * @brief K-d tree implementation for efficient spatial searching
 *
 * Implements a k-dimensional tree data structure for fast nearest-neighbor
 * searches in k-dimensional space. Supports single nearest neighbor, k-nearest
 * neighbors, and range queries.
 *
 * This file is part of ``kdtree'', a library for working with kd-trees.
 * Copyright (C) 2007-2011 John Tsiombikas <nuclear@member.fsf.org>
 *
 * Single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "kdtree.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * @def SQ(x)
 * @brief Computes square of a value
 *
 * @param x Value to square
 * @return x * x
 */
#define SQ(x) ((x) * (x))

static void clear_rec(struct kdnode *node, void (*destr)(void *));
static int insert_rec(struct kdnode **node, const double *pos, void *data,
                      int dir, int dim);
static int rlist_insert(struct res_node *list, struct kdnode *item,
                        double dist_sq);
static struct res_node *rlist_pop_back(struct res_node *list);
static void clear_results(struct kdres *set);

static struct kdhyperrect *hyperrect_create(int dim, const double *min,
                                            const double *max);
static void hyperrect_free(struct kdhyperrect *rect);
static struct kdhyperrect *hyperrect_duplicate(const struct kdhyperrect *rect);
static void hyperrect_extend(struct kdhyperrect *rect, const double *pos);
static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos);

static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node *);
static void free_resnode_buffer();

/**
 * @brief Creates a new k-d tree
 *
 * @param k Dimensionality of the space
 * @return Pointer to newly created kdtree, or NULL on allocation failure
 */
struct kdtree *kd_create(int k)
{
    struct kdtree *tree;

    if (!(tree = malloc(sizeof *tree)))
    {
        return 0;
    }

    tree->dim = k;
    tree->root = 0;
    tree->destr = 0;
    tree->rect = 0;

    return tree;
}

/**
 * @brief Frees a k-d tree and all its resources
 *
 * Clears all nodes and deallocates the tree structure. Also frees
 * the result node buffer pool.
 *
 * @param tree Pointer to kdtree to free
 */
void kd_free(struct kdtree *tree)
{
    if (tree)
    {
        kd_clear(tree);
        free(tree);
    }

    free_resnode_buffer();
}

/**
 * @brief Recursively clears all nodes in subtree
 *
 * @param node Root of subtree to clear
 * @param destr Optional destructor function for data pointers
 */
static void clear_rec(struct kdnode *node, void (*destr)(void *))
{
    if (!node)
    {
        return;
    }

    clear_rec(node->left, destr);
    clear_rec(node->right, destr);

    if (destr)
    {
        destr(node->data);
    }
    free(node->pos);
    free(node);
}

/**
 * @brief Removes all elements from the tree
 *
 * @param tree Pointer to kdtree to clear
 */
void kd_clear(struct kdtree *tree)
{
    clear_rec(tree->root, tree->destr);
    tree->root = 0;

    if (tree->rect)
    {
        hyperrect_free(tree->rect);
        tree->rect = 0;
    }
}

/**
 * @brief Sets destructor function for data pointers
 *
 * The provided function will be called on data pointers when nodes
 * are removed from the tree.
 *
 * @param tree Pointer to kdtree
 * @param destr Destructor function (can be NULL)
 */
void kd_data_destructor(struct kdtree *tree, void (*destr)(void *))
{
    tree->destr = destr;
}

/**
 * @brief Recursively inserts node into k-d tree
 *
 * Uses standard k-d tree insertion: cycles through dimensions, splitting
 * on the current dimension and recursing into appropriate subtree.
 *
 * @param nptr Pointer to node pointer (for recursive insertion)
 * @param pos Position coordinates (dim elements)
 * @param data User data pointer
 * @param dir Current splitting dimension
 * @param dim Total dimensionality
 * @return 0 on success, -1 on allocation failure
 */
static int insert_rec(struct kdnode **nptr, const double *pos, void *data,
                      int dir, int dim)
{
    int new_dir;
    struct kdnode *node;

    if (!*nptr)
    {
        if (!(node = malloc(sizeof *node)))
        {
            return -1;
        }
        if (!(node->pos = malloc(dim * sizeof *node->pos)))
        {
            free(node);
            return -1;
        }
        memcpy(node->pos, pos, dim * sizeof *node->pos);
        node->data = data;
        node->dir = dir;
        node->left = node->right = 0;
        *nptr = node;
        return 0;
    }

    node = *nptr;
    new_dir = (node->dir + 1) % dim;
    if (pos[node->dir] < node->pos[node->dir])
    {
        return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
    }
    return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

/**
 * @brief Inserts a node into the k-d tree
 *
 * @param tree Pointer to kdtree
 * @param pos Position coordinates (tree->dim elements)
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insert(struct kdtree *tree, const double *pos, void *data)
{
    if (insert_rec(&tree->root, pos, data, 0, tree->dim))
    {
        return -1;
    }

    if (tree->rect == 0)
    {
        tree->rect = hyperrect_create(tree->dim, pos, pos);
    }
    else
    {
        hyperrect_extend(tree->rect, pos);
    }

    return 0;
}

/**
 * @brief Inserts a node with float coordinates
 *
 * Converts float array to double internally. Uses stack buffer for
 * dimensions <= 16, heap allocation for larger dimensions.
 *
 * @param tree Pointer to kdtree
 * @param pos Position coordinates as floats (tree->dim elements)
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insertf(struct kdtree *tree, const float *pos, void *data)
{
    static double sbuf[16];
    double *bptr, *buf = 0;
    int res, dim = tree->dim;

    if (dim > 16)
    {
        if (!(bptr = buf = malloc(dim * sizeof *bptr)))
        {
            return -1;
        }
    }
    else
    {
        bptr = buf = sbuf;
    }

    while (dim-- > 0)
    {
        *bptr++ = *pos++;
    }

    res = kd_insert(tree, buf, data);

    if (tree->dim > 16)
    {
        free(buf);
    }
    return res;
}

/**
 * @brief Inserts a 3D point (convenience function)
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insert3(struct kdtree *tree, double x, double y, double z, void *data)
{
    double buf[3];
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
    return kd_insert(tree, buf, data);
}

/**
 * @brief Inserts a 3D point with float coordinates (convenience function)
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insert3f(struct kdtree *tree, float x, float y, float z, void *data)
{
    double buf[3];
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
    return kd_insert(tree, buf, data);
}

/**
 * @brief Recursively finds all nodes within range of query point
 *
 * Traverses tree, adding nodes within range to result list.
 * Uses pruning: only explores subtrees that could contain points within range.
 *
 * Distance metric: Euclidean distance squared
 * d^2 = sum_i (node_i - pos_i)^2
 *
 * @param node Current node in traversal
 * @param pos Query position (dim elements)
 * @param range Search radius
 * @param list Result list to append to
 * @param ordered If true, maintain sorted order; if false, use -1.0 for dist_sq
 * @param dim Dimensionality
 * @return Number of nodes added, or -1 on error
 */
static int find_nearest(struct kdnode *node, const double *pos, double range,
                        struct res_node *list, int ordered, int dim)
{
    double dist_sq, dx;
    int i, ret, added_res = 0;

    if (!node)
    {
        return 0;
    }

    dist_sq = 0;
    for (i = 0; i < dim; i++)
    {
        dist_sq += SQ(node->pos[i] - pos[i]);
    }
    if (dist_sq <= SQ(range))
    {
        if (rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1)
        {
            return -1;
        }
        added_res = 1;
    }

    dx = pos[node->dir] - node->pos[node->dir];

    ret = find_nearest(dx <= 0.0 ? node->left : node->right, pos, range, list,
                       ordered, dim);
    if (ret >= 0 && fabs(dx) < range)
    {
        added_res += ret;
        ret = find_nearest(dx <= 0.0 ? node->right : node->left, pos, range,
                           list, ordered, dim);
    }
    if (ret == -1)
    {
        return -1;
    }
    added_res += ret;

    return added_res;
}

/**
 * @brief Recursively finds N nearest nodes to query point
 *
 * Maintains a bounded priority queue of size N. Updates dist_max as
 * queue fills to enable more aggressive pruning.
 *
 * @param node Current node in traversal
 * @param pos Query position (dim elements)
 * @param num Maximum number of neighbors to find
 * @param size Current number of results found
 * @param dist_max Maximum distance squared in current result set
 * @param list Result list (sorted by distance)
 * @param dim Dimensionality
 * @return 0 on success, -1 on error
 */
static int find_nearest_n(struct kdnode *node, const double *pos, int num,
                          int *size, double *dist_max, struct res_node *list,
                          int dim)
{
    double dist_sq, dx;
    int i, ret;

    if (!node)
    {
        return 0;
    }

    dist_sq = 0;
    for (i = 0; i < dim; i++)
    {
        dist_sq += SQ(node->pos[i] - pos[i]);
    }

    if (dist_sq < *dist_max)
    {
        if (*size < num)
        {
            ++(*size);
        }
        else
        {
            struct res_node *back = rlist_pop_back(list);
            if (back == 0)
            {
                return -1;
            }
            *dist_max = back->dist_sq > dist_sq ? back->dist_sq : dist_sq;
        }
        if (rlist_insert(list, node, dist_sq) == -1)
        {
            return -1;
        }
    }

    /* find signed distance from the splitting plane */
    dx = pos[node->dir] - node->pos[node->dir];

    ret = find_nearest_n(dx <= 0.0 ? node->left : node->right, pos, num, size,
                         dist_max, list, dim);
    if (ret >= 0 && SQ(dx) < *dist_max)
    {
        ret = find_nearest_n(dx <= 0.0 ? node->right : node->left, pos, num,
                             size, dist_max, list, dim);
    }
    return ret;
}

/**
 * @brief Internal recursive single nearest neighbor search
 *
 * Uses hyperrectangle pruning for efficient search. At each node:
 * 1. Recursively search nearer subtree
 * 2. Check current node
 * 3. If splitting plane distance < current best, search farther subtree
 *
 * @param node Current node in traversal
 * @param pos Query position (rect->dim elements)
 * @param result Output: pointer to nearest node found
 * @param result_dist_sq Output: squared distance to nearest node
 * @param rect Bounding hyperrectangle for current subtree
 */
static void kd_nearest_i(struct kdnode *node, const double *pos,
                         struct kdnode **result, double *result_dist_sq,
                         struct kdhyperrect *rect)
{
    int dir = node->dir;
    int i;
    double dummy, dist_sq;
    struct kdnode *nearer_subtree, *farther_subtree;
    double *nearer_hyperrect_coord, *farther_hyperrect_coord;

    /* Decide whether to go left or right in the tree */
    dummy = pos[dir] - node->pos[dir];
    if (dummy <= 0)
    {
        nearer_subtree = node->left;
        farther_subtree = node->right;
        nearer_hyperrect_coord = rect->max + dir;
        farther_hyperrect_coord = rect->min + dir;
    }
    else
    {
        nearer_subtree = node->right;
        farther_subtree = node->left;
        nearer_hyperrect_coord = rect->min + dir;
        farther_hyperrect_coord = rect->max + dir;
    }

    if (nearer_subtree)
    {
        /* Slice the hyperrect to get the hyperrect of the nearer subtree */
        dummy = *nearer_hyperrect_coord;
        *nearer_hyperrect_coord = node->pos[dir];
        /* Recurse down into nearer subtree */
        kd_nearest_i(nearer_subtree, pos, result, result_dist_sq, rect);
        /* Undo the slice */
        *nearer_hyperrect_coord = dummy;
    }

    /* Check the distance of the point at the current node, compare it
     * with our best so far */
    dist_sq = 0;
    for (i = 0; i < rect->dim; i++)
    {
        dist_sq += SQ(node->pos[i] - pos[i]);
    }
    if (dist_sq < *result_dist_sq)
    {
        *result = node;
        *result_dist_sq = dist_sq;
    }

    if (farther_subtree)
    {
        /* Get the hyperrect of the farther subtree */
        dummy = *farther_hyperrect_coord;
        *farther_hyperrect_coord = node->pos[dir];
        /* Check if we have to recurse down by calculating the closest
         * point of the hyperrect and see if it's closer than our
         * minimum distance in result_dist_sq. */
        if (hyperrect_dist_sq(rect, pos) < *result_dist_sq)
        {
            /* Recurse down into farther subtree */
            kd_nearest_i(farther_subtree, pos, result, result_dist_sq, rect);
        }
        /* Undo the slice on the hyperrect */
        *farther_hyperrect_coord = dummy;
    }
}

/**
 * @brief Finds single nearest neighbor to query point
 *
 * @param kd Pointer to kdtree
 * @param pos Query position (kd->dim elements)
 * @return Result set with at most one element, or NULL on error
 *
 * @note Caller must free result set with kd_res_free()
 */
struct kdres *kd_nearest(struct kdtree *kd, const double *pos)
{
    struct kdhyperrect *rect;
    struct kdnode *result;
    struct kdres *rset;
    double dist_sq;
    int i;

    if (!kd)
    {
        return 0;
    }
    if (!kd->rect)
    {
        return 0;
    }

    /* Allocate result set */
    if (!(rset = malloc(sizeof *rset)))
    {
        return 0;
    }
    if (!(rset->rlist = alloc_resnode()))
    {
        free(rset);
        return 0;
    }
    rset->rlist->next = 0;
    rset->tree = kd;

    /* Duplicate the bounding hyperrectangle, we will work on the copy */
    if (!(rect = hyperrect_duplicate(kd->rect)))
    {
        kd_res_free(rset);
        return 0;
    }

    /* Our first guesstimate is the root node */
    result = kd->root;
    dist_sq = 0;
    for (i = 0; i < kd->dim; i++)
    {
        dist_sq += SQ(result->pos[i] - pos[i]);
    }

    /* Search for the nearest neighbour recursively */
    kd_nearest_i(kd->root, pos, &result, &dist_sq, rect);

    /* Free the copy of the hyperrect */
    hyperrect_free(rect);

    /* Store the result */
    if (result)
    {
        if (rlist_insert(rset->rlist, result, -1.0) == -1)
        {
            kd_res_free(rset);
            return 0;
        }
        rset->size = 1;
        kd_res_rewind(rset);
        return rset;
    }
    else
    {
        kd_res_free(rset);
        return 0;
    }
}

/**
 * @brief Finds single nearest neighbor with float coordinates
 *
 * @param tree Pointer to kdtree
 * @param pos Query position as floats (tree->dim elements)
 * @return Result set with at most one element, or NULL on error
 */
struct kdres *kd_nearestf(struct kdtree *tree, const float *pos)
{
    static double sbuf[16];
    double *bptr, *buf = 0;
    int dim = tree->dim;
    struct kdres *res;

    if (dim > 16)
    {
        if (!(bptr = buf = malloc(dim * sizeof *bptr)))
        {
            return 0;
        }
    }
    else
    {
        bptr = buf = sbuf;
    }

    while (dim-- > 0)
    {
        *bptr++ = *pos++;
    }

    res = kd_nearest(tree, buf);

    if (tree->dim > 16)
    {
        free(buf);
    }
    return res;
}

/**
 * @brief Finds single nearest neighbor in 3D (convenience function)
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Result set with at most one element, or NULL on error
 */
struct kdres *kd_nearest3(struct kdtree *tree, double x, double y, double z)
{
    double pos[3];
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    return kd_nearest(tree, pos);
}

/**
 * @brief Finds single nearest neighbor in 3D with float coordinates
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Result set with at most one element, or NULL on error
 */
struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z)
{
    double pos[3];
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    return kd_nearest(tree, pos);
}

/**
 * @brief Finds N nearest neighbors to query point
 *
 * @param kd Pointer to kdtree
 * @param pos Query position (kd->dim elements)
 * @param num Maximum number of neighbors to find
 * @return Result set with at most num elements, or NULL on error
 *
 * @note Caller must free result set with kd_res_free()
 */
struct kdres *kd_nearest_n(struct kdtree *kd, const double *pos, int num)
{
    int ret, size = 0;
    struct kdres *rset;
    double dist_max = DBL_MAX;

    if (!(rset = malloc(sizeof *rset)))
    {
        return 0;
    }
    if (!(rset->rlist = alloc_resnode()))
    {
        free(rset);
        return 0;
    }
    rset->rlist->next = 0;
    rset->tree = kd;

    ret = find_nearest_n(kd->root, pos, num, &size, &dist_max, rset->rlist,
                         kd->dim);
    if (ret == -1)
    {
        kd_res_free(rset);
        return 0;
    }
    rset->size = size;
    kd_res_rewind(rset);
    return rset;
}

/**
 * @brief Finds N nearest neighbors with float coordinates
 *
 * @param tree Pointer to kdtree
 * @param pos Query position as floats (tree->dim elements)
 * @param num Maximum number of neighbors to find
 * @return Result set with at most num elements, or NULL on error
 */
struct kdres *kd_nearest_nf(struct kdtree *tree, const float *pos, int num)
{
    static double sbuf[16];
    double *bptr, *buf = 0;
    int dim = tree->dim;
    struct kdres *res;

    if (dim > 16)
    {
        if (!(bptr = buf = malloc(dim * sizeof *bptr)))
        {
            return 0;
        }
    }
    else
    {
        bptr = buf = sbuf;
    }

    while (dim-- > 0)
    {
        *bptr++ = *pos++;
    }

    res = kd_nearest_n(tree, buf, num);

    if (tree->dim > 16)
    {
        free(buf);
    }
    return res;
}

/**
 * @brief Finds N nearest neighbors in 3D (convenience function)
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param num Maximum number of neighbors to find
 * @return Result set with at most num elements, or NULL on error
 */
struct kdres *kd_nearest_n3(struct kdtree *tree, double x, double y, double z,
                            int num)
{
    double pos[3];
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    return kd_nearest_n(tree, pos, num);
}

/**
 * @brief Finds N nearest neighbors in 3D with float coordinates
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param num Maximum number of neighbors to find
 * @return Result set with at most num elements, or NULL on error
 */
struct kdres *kd_nearest_n3f(struct kdtree *tree, float x, float y, float z,
                             int num)
{
    double pos[3];
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    return kd_nearest_n(tree, pos, num);
}

/**
 * @brief Finds all nodes within range of query point
 *
 * @param kd Pointer to kdtree
 * @param pos Query position (kd->dim elements)
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 *
 * @note Caller must free result set with kd_res_free()
 */
struct kdres *kd_nearest_range(struct kdtree *kd, const double *pos,
                               double range)
{
    int ret;
    struct kdres *rset;

    if (!(rset = malloc(sizeof *rset)))
    {
        return 0;
    }
    if (!(rset->rlist = alloc_resnode()))
    {
        free(rset);
        return 0;
    }
    rset->rlist->next = 0;
    rset->tree = kd;

    if ((ret = find_nearest(kd->root, pos, range, rset->rlist, 0, kd->dim)) ==
        -1)
    {
        kd_res_free(rset);
        return 0;
    }
    rset->size = ret;
    kd_res_rewind(rset);
    return rset;
}

/**
 * @brief Finds all nodes within range with float coordinates
 *
 * @param kd Pointer to kdtree
 * @param pos Query position as floats (kd->dim elements)
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 */
struct kdres *kd_nearest_rangef(struct kdtree *kd, const float *pos,
                                float range)
{
    static double sbuf[16];
    double *bptr, *buf = 0;
    int dim = kd->dim;
    struct kdres *res;

    if (dim > 16)
    {
        if (!(bptr = buf = malloc(dim * sizeof *bptr)))
        {
            return 0;
        }
    }
    else
    {
        bptr = buf = sbuf;
    }

    while (dim-- > 0)
    {
        *bptr++ = *pos++;
    }

    res = kd_nearest_range(kd, buf, range);

    if (kd->dim > 16)
    {
        free(buf);
    }
    return res;
}

/**
 * @brief Finds all nodes within range in 3D (convenience function)
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 */
struct kdres *kd_nearest_range3(struct kdtree *tree, double x, double y,
                                double z, double range)
{
    double buf[3];
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
    return kd_nearest_range(tree, buf, range);
}

/**
 * @brief Finds all nodes within range in 3D with float coordinates
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 */
struct kdres *kd_nearest_range3f(struct kdtree *tree, float x, float y, float z,
                                 float range)
{
    double buf[3];
    buf[0] = x;
    buf[1] = y;
    buf[2] = z;
    return kd_nearest_range(tree, buf, range);
}

/**
 * @brief Frees a result set
 *
 * @param rset Result set to free
 */
void kd_res_free(struct kdres *rset)
{
    clear_results(rset);
    free_resnode(rset->rlist);
    free(rset);
}

/**
 * @brief Returns the number of elements in result set
 *
 * @param set Result set
 * @return Number of elements
 */
int kd_res_size(struct kdres *set) { return (set->size); }

/**
 * @brief Rewinds result set iterator to beginning
 *
 * @param rset Result set to rewind
 */
void kd_res_rewind(struct kdres *rset) { rset->riter = rset->rlist->next; }

/**
 * @brief Checks if result set iterator is at end
 *
 * @param rset Result set to check
 * @return Non-zero if at end, zero otherwise
 */
int kd_res_end(struct kdres *rset) { return rset->riter == 0; }

/**
 * @brief Advances result set iterator to next element
 *
 * @param rset Result set
 * @return Non-zero if advanced successfully, zero if no more elements
 */
int kd_res_next(struct kdres *rset)
{
    rset->riter = rset->riter->next;
    return rset->riter != 0;
}

/**
 * @brief Returns data pointer and position of current result element
 *
 * @param rset Result set
 * @param pos Output array for position (rset->tree->dim elements, can be NULL)
 * @return Data pointer associated with current element, or NULL if at end
 */
void *kd_res_item(struct kdres *rset, double *pos)
{
    if (rset->riter)
    {
        if (pos)
        {
            memcpy(pos, rset->riter->item->pos, rset->tree->dim * sizeof *pos);
        }
        return rset->riter->item->data;
    }
    return 0;
}

/**
 * @brief Returns data pointer and position with float output
 *
 * @param rset Result set
 * @param pos Output array for position as floats (can be NULL)
 * @return Data pointer associated with current element, or NULL if at end
 */
void *kd_res_itemf(struct kdres *rset, float *pos)
{
    if (rset->riter)
    {
        if (pos)
        {
            int i;
            for (i = 0; i < rset->tree->dim; i++)
            {
                pos[i] = rset->riter->item->pos[i];
            }
        }
        return rset->riter->item->data;
    }
    return 0;
}

/**
 * @brief Returns data pointer and 3D position (convenience function)
 *
 * @param rset Result set (must be 3-dimensional)
 * @param x Output pointer for X coordinate (can be NULL)
 * @param y Output pointer for Y coordinate (can be NULL)
 * @param z Output pointer for Z coordinate (can be NULL)
 * @return Data pointer associated with current element, or NULL if at end
 */
void *kd_res_item3(struct kdres *rset, double *x, double *y, double *z)
{
    if (rset->riter)
    {
        if (x)
        {
            *x = rset->riter->item->pos[0];
        }
        if (y)
        {
            *y = rset->riter->item->pos[1];
        }
        if (z)
        {
            *z = rset->riter->item->pos[2];
        }
        return rset->riter->item->data;
    }
    return 0;
}

/**
 * @brief Returns data pointer and 3D position with float output
 *
 * @param rset Result set (must be 3-dimensional)
 * @param x Output pointer for X coordinate (can be NULL)
 * @param y Output pointer for Y coordinate (can be NULL)
 * @param z Output pointer for Z coordinate (can be NULL)
 * @return Data pointer associated with current element, or NULL if at end
 */
void *kd_res_item3f(struct kdres *rset, float *x, float *y, float *z)
{
    if (rset->riter)
    {
        if (x)
        {
            *x = rset->riter->item->pos[0];
        }
        if (y)
        {
            *y = rset->riter->item->pos[1];
        }
        if (z)
        {
            *z = rset->riter->item->pos[2];
        }
        return rset->riter->item->data;
    }
    return 0;
}

/**
 * @brief Returns data pointer of current result element
 *
 * Equivalent to kd_res_item(set, NULL)
 *
 * @param set Result set
 * @return Data pointer, or NULL if at end
 */
void *kd_res_item_data(struct kdres *set) { return kd_res_item(set, 0); }

/**
 * @brief Returns distance from query point to current result element
 *
 * @param set Result set
 * @return Euclidean distance (not squared)
 */
double kd_res_dist(struct kdres *set) { return sqrt(set->riter->dist_sq); }

/* ---- hyperrectangle helpers ---- */

/**
 * @brief Creates a hyperrectangle with given bounds
 *
 * @param dim Dimensionality
 * @param min Minimum coordinates (dim elements)
 * @param max Maximum coordinates (dim elements)
 * @return Pointer to new hyperrectangle, or NULL on allocation failure
 */
static struct kdhyperrect *hyperrect_create(int dim, const double *min,
                                            const double *max)
{
    size_t size = dim * sizeof(double);
    struct kdhyperrect *rect = 0;

    if (!(rect = malloc(sizeof(struct kdhyperrect))))
    {
        return 0;
    }

    rect->dim = dim;
    if (!(rect->min = malloc(size)))
    {
        free(rect);
        return 0;
    }
    if (!(rect->max = malloc(size)))
    {
        free(rect->min);
        free(rect);
        return 0;
    }
    memcpy(rect->min, min, size);
    memcpy(rect->max, max, size);

    return rect;
}

/**
 * @brief Frees a hyperrectangle
 *
 * @param rect Hyperrectangle to free
 */
static void hyperrect_free(struct kdhyperrect *rect)
{
    free(rect->min);
    free(rect->max);
    free(rect);
}

/**
 * @brief Duplicates a hyperrectangle
 *
 * @param rect Hyperrectangle to duplicate
 * @return Pointer to new copy, or NULL on allocation failure
 */
static struct kdhyperrect *hyperrect_duplicate(const struct kdhyperrect *rect)
{
    return hyperrect_create(rect->dim, rect->min, rect->max);
}

/**
 * @brief Extends hyperrectangle to include given point
 *
 * Updates min/max bounds to ensure point is contained.
 *
 * @param rect Hyperrectangle to extend
 * @param pos Point coordinates (rect->dim elements)
 */
static void hyperrect_extend(struct kdhyperrect *rect, const double *pos)
{
    int i;

    for (i = 0; i < rect->dim; i++)
    {
        if (pos[i] < rect->min[i])
        {
            rect->min[i] = pos[i];
        }
        if (pos[i] > rect->max[i])
        {
            rect->max[i] = pos[i];
        }
    }
}

/**
 * @brief Computes squared distance from point to hyperrectangle
 *
 * Returns 0 if point is inside or on boundary of hyperrectangle.
 * Otherwise returns sum of squared distances along each axis where
 * point is outside bounds.
 *
 * Distance formula:
 * d^2 = sum_i [max(0, min_i - pos_i)]^2 + [max(0, pos_i - max_i)]^2
 *
 * @param rect Hyperrectangle
 * @param pos Point coordinates (rect->dim elements)
 * @return Squared distance to nearest point on hyperrectangle
 */
static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos)
{
    int i;
    double result = 0;

    for (i = 0; i < rect->dim; i++)
    {
        if (pos[i] < rect->min[i])
        {
            result += SQ(rect->min[i] - pos[i]);
        }
        else if (pos[i] > rect->max[i])
        {
            result += SQ(rect->max[i] - pos[i]);
        }
    }

    return result;
}

/* ---- static helpers ---- */

/**
 * @brief Static pool of free result nodes for reuse
 */
static struct res_node *free_nodes;

/**
 * @brief Allocates a result node from pool or heap
 *
 * Uses object pool pattern to reduce allocation overhead.
 *
 * @return Pointer to result node, or NULL on allocation failure
 */
static struct res_node *alloc_resnode(void)
{
    struct res_node *node;

    if (!free_nodes)
    {
        node = malloc(sizeof *node);
    }
    else
    {
        node = free_nodes;
        free_nodes = free_nodes->next;
        node->next = 0;
    }

    return node;
}

/**
 * @brief Returns result node to pool for reuse
 *
 * @param node Result node to free
 */
static void free_resnode(struct res_node *node)
{
    node->next = free_nodes;
    free_nodes = node;
}

/**
 * @brief Frees all nodes in the result node pool
 *
 * Should be called at program termination to avoid memory leaks.
 */
static void free_resnode_buffer()
{
    if (free_nodes)
    {
        struct res_node *ptr = free_nodes;
        while (ptr)
        {
            ptr = ptr->next;
            free(free_nodes);
            free_nodes = ptr;
        }
        free_nodes = 0;
    }
}

/**
 * @brief Inserts item into sorted result list
 *
 * If dist_sq >= 0, maintains sorted order by distance.
 * If dist_sq < 0, performs unsorted insertion at head.
 *
 * @param list Result list (head sentinel node)
 * @param item k-d tree node to insert
 * @param dist_sq Squared distance (-1.0 for unsorted insertion)
 * @return 0 on success, -1 on allocation failure
 *
 * @note TODO: Use heap sort for better performance with large result sets
 */
static int rlist_insert(struct res_node *list, struct kdnode *item,
                        double dist_sq)
{
    struct res_node *rnode;

    if (!(rnode = alloc_resnode()))
    {
        return -1;
    }
    rnode->item = item;
    rnode->dist_sq = dist_sq;

    if (dist_sq >= 0.0)
    {
        while (list->next && list->next->dist_sq < dist_sq)
        {
            list = list->next;
        }
    }
    rnode->next = list->next;
    list->next = rnode;
    return 0;
}

/**
 * @brief Removes and returns last element from result list
 *
 * @param list Result list (head sentinel node)
 * @return Previous-to-last node (new tail), or NULL if list empty
 *
 * @note The removed node is returned to the pool via free_resnode
 */
static struct res_node *rlist_pop_back(struct res_node *list)
{
    struct res_node *previous = 0;
    while (list->next)
    {
        previous = list;
        list = list->next;
    }
    if (previous)
    {
        previous->next = 0;
    }
    free_resnode(list);
    return previous;
}

/**
 * @brief Clears all elements from result list
 *
 * Returns all nodes to pool except the head sentinel.
 *
 * @param rset Result set containing list to clear
 */
static void clear_results(struct kdres *rset)
{
    struct res_node *tmp, *node = rset->rlist->next;

    while (node)
    {
        tmp = node;
        node = node->next;
        free_resnode(tmp);
    }

    rset->rlist->next = 0;
}

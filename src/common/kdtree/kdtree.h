/**
 * @file
 * @brief K-d tree data structure declarations
 *
 * This file is part of ``kdtree'', a library for working with kd-trees.
 * Copyright (C) 2007-2011 John Tsiombikas <nuclear@member.fsf.org>
 *
 * [License text same as .c file - omitted for brevity]
 */

#pragma once

/**
 * @struct kdhyperrect
 * @brief Axis-aligned hyperrectangle (bounding box)
 */
struct kdhyperrect
{
    int dim;           /**< Dimensionality */
    double *min, *max; /**< Minimum and maximum coordinates per dimension */
};

/**
 * @struct kdnode
 * @brief Node in k-d tree
 */
struct kdnode
{
    double *pos; /**< Position coordinates (dim elements) */
    int dir;     /**< Splitting dimension for this node */
    void *data;  /**< User data pointer */

    struct kdnode *left,
        *right; /**< Child nodes (negative/positive side of split) */
};

/**
 * @struct res_node
 * @brief Node in result list (linked list)
 */
struct res_node
{
    struct kdnode *item;   /**< Pointer to k-d tree node */
    double dist_sq;        /**< Squared distance to query point */
    struct res_node *next; /**< Next node in list */
};

/**
 * @struct kdtree
 * @brief K-d tree structure
 */
struct kdtree
{
    int dim;                  /**< Dimensionality of space */
    struct kdnode *root;      /**< Root node of tree */
    struct kdhyperrect *rect; /**< Bounding hyperrectangle of all points */
    void (*destr)(void *);    /**< Optional destructor for data pointers */
};

/**
 * @struct kdres
 * @brief Result set from nearest neighbor queries
 */
struct kdres
{
    struct kdtree *tree; /**< Pointer to source tree */
    struct res_node *rlist,
        *riter; /**< Result list head and current iterator */
    int size;   /**< Number of elements in result set */
};

/**
 * @brief Creates a k-d tree for k-dimensional data
 *
 * @param k Dimensionality
 * @return Pointer to new kdtree, or NULL on failure
 */
struct kdtree *kd_create(int k);

/**
 * @brief Frees the kdtree structure
 *
 * @param tree Pointer to kdtree to free
 */
void kd_free(struct kdtree *tree);

/**
 * @brief Removes all elements from the tree
 *
 * @param tree Pointer to kdtree to clear
 */
void kd_clear(struct kdtree *tree);

/**
 * @brief Sets destructor function for data pointers
 *
 * The provided function will be called on data pointers when nodes
 * are removed from the tree.
 *
 * @param tree Pointer to kdtree
 * @param destr Destructor function (can be NULL)
 */
void kd_data_destructor(struct kdtree *tree, void (*destr)(void *));

/**
 * @brief Inserts a node with specified position and optional data
 *
 * @param tree Pointer to kdtree
 * @param pos Position coordinates (tree->dim elements)
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insert(struct kdtree *tree, const double *pos, void *data);

/**
 * @brief Inserts a node with float coordinates
 *
 * @param tree Pointer to kdtree
 * @param pos Position coordinates as floats
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insertf(struct kdtree *tree, const float *pos, void *data);

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
int kd_insert3(struct kdtree *tree, double x, double y, double z, void *data);

/**
 * @brief Inserts a 3D point with float coordinates
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param data Optional user data pointer
 * @return 0 on success, -1 on failure
 */
int kd_insert3f(struct kdtree *tree, float x, float y, float z, void *data);

/**
 * @brief Finds the nearest node from a given point
 *
 * @param tree Pointer to kdtree
 * @param pos Query position
 * @return Result set with at most one element, or NULL on error
 *
 * @note Caller must free result with kd_res_free()
 */
struct kdres *kd_nearest(struct kdtree *tree, const double *pos);

/**
 * @brief Finds nearest node with float query position
 *
 * @param tree Pointer to kdtree
 * @param pos Query position as floats
 * @return Result set with at most one element, or NULL on error
 */
struct kdres *kd_nearestf(struct kdtree *tree, const float *pos);

/**
 * @brief Finds nearest node in 3D
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Result set with at most one element, or NULL on error
 */
struct kdres *kd_nearest3(struct kdtree *tree, double x, double y, double z);

/**
 * @brief Finds nearest node in 3D with float coordinates
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @return Result set with at most one element, or NULL on error
 */
struct kdres *kd_nearest3f(struct kdtree *tree, float x, float y, float z);

/**
 * @brief Finds the N nearest nodes from a given point
 *
 * @param tree Pointer to kdtree
 * @param pos Query position
 * @param num Maximum number of neighbors to find
 * @return Result set with at most N elements, or NULL on error
 *
 * @note Result set must be freed with kd_res_free()
 */
struct kdres *kd_nearest_n(struct kdtree *tree, const double *pos, int num);

/**
 * @brief Finds N nearest nodes with float query position
 *
 * @param tree Pointer to kdtree
 * @param pos Query position as floats
 * @param num Maximum number of neighbors to find
 * @return Result set with at most N elements, or NULL on error
 */
struct kdres *kd_nearest_nf(struct kdtree *tree, const float *pos, int num);

/**
 * @brief Finds N nearest nodes in 3D
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param num Maximum number of neighbors to find
 * @return Result set with at most N elements, or NULL on error
 */
struct kdres *kd_nearest_n3(struct kdtree *tree, double x, double y, double z,
                            int num);

/**
 * @brief Finds N nearest nodes in 3D with float coordinates
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param num Maximum number of neighbors to find
 * @return Result set with at most N elements, or NULL on error
 */
struct kdres *kd_nearest_n3f(struct kdtree *tree, float x, float y, float z,
                             int num);

/**
 * @brief Finds all nodes within a range from a given point
 *
 * @param tree Pointer to kdtree
 * @param pos Query position
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 *
 * @note Result set must be freed with kd_res_free()
 */
struct kdres *kd_nearest_range(struct kdtree *tree, const double *pos,
                               double range);

/**
 * @brief Finds all nodes within range with float query position
 *
 * @param tree Pointer to kdtree
 * @param pos Query position as floats
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 */
struct kdres *kd_nearest_rangef(struct kdtree *tree, const float *pos,
                                float range);

/**
 * @brief Finds all nodes within range in 3D
 *
 * @param tree Pointer to kdtree (must be 3-dimensional)
 * @param x X coordinate
 * @param y Y coordinate
 * @param z Z coordinate
 * @param range Search radius
 * @return Result set with all nodes within range, or NULL on error
 */
struct kdres *kd_nearest_range3(struct kdtree *tree, double x, double y,
                                double z, double range);

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
                                 float range);

/**
 * @brief Frees a result set
 *
 * @param set Result set to free
 */
void kd_res_free(struct kdres *set);

/**
 * @brief Returns the size of the result set in elements
 *
 * @param set Result set
 * @return Number of elements
 */
int kd_res_size(struct kdres *set);

/**
 * @brief Rewinds the result set iterator to beginning
 *
 * @param set Result set to rewind
 */
void kd_res_rewind(struct kdres *set);

/**
 * @brief Checks if iterator reached the end
 *
 * @param set Result set
 * @return Non-zero if at end, zero otherwise
 */
int kd_res_end(struct kdres *set);

/**
 * @brief Advances the result set iterator
 *
 * @param set Result set
 * @return Non-zero on success, zero if no more elements
 */
int kd_res_next(struct kdres *set);

/**
 * @brief Returns data pointer and position of current result item
 *
 * @param set Result set
 * @param pos Output array for position (can be NULL)
 * @return Data pointer, or NULL if at end
 */
void *kd_res_item(struct kdres *set, double *pos);

/**
 * @brief Returns data pointer and position with float output
 *
 * @param set Result set
 * @param pos Output array for position as floats (can be NULL)
 * @return Data pointer, or NULL if at end
 */
void *kd_res_itemf(struct kdres *set, float *pos);

/**
 * @brief Returns data pointer and 3D position
 *
 * @param set Result set (must be 3-dimensional)
 * @param x Output pointer for X (can be NULL)
 * @param y Output pointer for Y (can be NULL)
 * @param z Output pointer for Z (can be NULL)
 * @return Data pointer, or NULL if at end
 */
void *kd_res_item3(struct kdres *set, double *x, double *y, double *z);

/**
 * @brief Returns data pointer and 3D position with float output
 *
 * @param set Result set (must be 3-dimensional)
 * @param x Output pointer for X (can be NULL)
 * @param y Output pointer for Y (can be NULL)
 * @param z Output pointer for Z (can be NULL)
 * @return Data pointer, or NULL if at end
 */
void *kd_res_item3f(struct kdres *set, float *x, float *y, float *z);

/**
 * @brief Returns data pointer of current result item
 *
 * Equivalent to kd_res_item(set, NULL)
 *
 * @param set Result set
 * @return Data pointer, or NULL if at end
 */
void *kd_res_item_data(struct kdres *set);

/**
 * @brief Returns distance from query point to current result item
 *
 * @param set Result set
 * @return Euclidean distance (not squared)
 */
double kd_res_dist(struct kdres *set);

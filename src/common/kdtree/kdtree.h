/*
 * MIT License

Copyright (c) 2017 Leo Ma

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 * Copyright (C) 2017, Leo Ma <begeekmyfriend@gmail.com>
 */

#pragma once
#ifndef _KD_TREE_H
#define _KD_TREE_H

#include <stddef.h>

#define KDTREE_MAX_LEVEL 64
#define KDTREE_LEFT_INDEX 0
#define KDTREE_RIGHT_INDEX 1

typedef struct knn_list
{
    struct kdnode *node;
    double distance;
    struct knn_list *prev;
    struct knn_list *next;
} knn_list_t;

typedef struct kdnode
{
    long coord_index;
    double *coord;
    struct kdnode *left;
    struct kdnode *right;
    int r;
} kdnode_t;

typedef struct kdtree
{
    struct kdnode *root;
    size_t count;
    size_t capacity;
    double *coords;
    double **coord_table;
    long *coord_indexes;
    unsigned char *coord_deleted;
    unsigned char *coord_passed;
    struct knn_list knn_list_head;
    int dim;
    int knn_num;
    double eps;
    double eps_nearest_dist;
} kdtree_t;

struct kdtree *kdtree_init(int dim);
void kdtree_insert(struct kdtree *tree, double *coord);
void kdtree_rebuild(struct kdtree *tree);
void kdtree_knn_search(struct kdtree *tree, double *coord, int k);
void kdtree_eps_nearest_search(struct kdtree *tree, double *target, double EPS);
void kdtree_destroy(struct kdtree *tree);
void kdtree_dump(struct kdtree *tree);

#endif /* _KD_TREE_H */

/* k-D tree */
/*
  Copyright (C) 2015 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

#include "kdtree.h"

struct kd_Node {
    float *x;
    int i;
    struct kd_Node *prev, *next;
};

#ifndef _kdtree_h

typedef struct kd_Node *kd_node;
/*^*/

#endif
 
static float dist(const float *a, const float *b, int dim)
{
    float t, d;
    int j;

    d = 0.0f;
    for (j=0; j < dim; j++) {
        t = a[j] - b[j];
        d += t * t;
    }
    return d;
}

static void swap(kd_node x, kd_node y, int size) 
{
    int i;
    float tmp[SF_MAX_DIM];
    memcpy(tmp,  x->x, size);
    memcpy(x->x, y->x, size);
    memcpy(y->x, tmp,  size);
    i = x->i;
    x->i = y->i;
    y->i = i;
}

 
/* median by quickselect method */
static kd_node median(kd_node start, kd_node end, int idx, int size)
{
    kd_node p, store, md;
    float pivot;

    if (end <= start) return NULL;
    if (end == start + 1) return start;
 
    md = start + (end - start) / 2;
    
    while (1) {
        pivot = md->x[idx];
 
        swap(md, end - 1, size);
        for (store = p = start; p < end; p++) {
            if (p->x[idx] < pivot) {
                if (p != store)
                    swap(p, store, size);
                store++;
            }
        }
        swap(store, end - 1, size);
 
        /* median has duplicate values */
        if (store->x[idx] == md->x[idx])
            return md;
 
        if (store > md) end = store;
        else        start = store;
    }
}

static kd_node make_tree(kd_node t, int len, int i, int dim)
{
    kd_node n;

    if (0==len) return NULL;
 
    n = median(t, t + len, i, dim*sizeof(float));
 
    if (NULL != n) {
        i = (i + 1) % dim;
        n->prev  = make_tree(t, n - t, i, dim);
        n->next = make_tree(n + 1, t + len - (n + 1), i, dim);
    }

    return n;
}

float *kd_coord(kd_node n)
/*< node coordinates >*/
{
    return n->x;
}

kd_node kd_tree(float **data /* [len][dim] */, int len, int dim)
/*< make a k-D tree >*/
{
    kd_node n;
    int i, j;

    n = (kd_node) sf_alloc(len,sizeof(*n));
    
    for (i=0; i < len; i++) {
	n[i].i = i;
	n[i].x = sf_floatalloc(dim);
	for (j=0; j < dim; j++) {
	    n[i].x[j] = data[i][j];
	}
    }
    n = make_tree(n, len, 0, dim);

    return n;
}

void free_tree(kd_node n, int len)
/*< free allocated storage >*/
{
    int i;

    for (i=0; i < len; i++) {
	free(n[i].x);
    }
    free(n);
}

void kd_nearest(kd_node root, float* x, int i, int dim,
		kd_node *best, float *best_dist)
/*< find the nearest point in k-D tree >*/
{
    float d, dx, dx2;
 
    if (NULL == root) return;
    
    d = dist(root->x, x, dim);
    dx = root->x[i] - x[i];
    dx2 = dx * dx;
  
    if (d < *best_dist) {
        *best_dist = d;
        *best = root;
    }
    
    if (0.0f == *best_dist) return;
 
    i = (i + 1) % dim;
 
    kd_nearest(dx > 0 ? root->prev : root->next, x, i, dim, best, best_dist);
    if (dx2 >= *best_dist) return;
    kd_nearest(dx > 0 ? root->next : root->prev, x, i, dim, best, best_dist);
}

int kd_order(kd_node node, int* arr, int i)
/*< put node indeces in an array >*/
{
    if (NULL != node) {
	arr[i] = node->i;
	i++;
	
	i = kd_order(node->prev,arr,i);
	i = kd_order(node->next,arr,i);
    }

    return i;
}



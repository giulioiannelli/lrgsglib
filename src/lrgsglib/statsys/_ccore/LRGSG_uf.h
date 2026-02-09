/**
 * @file LRGSG_uf.h
 * @brief Union-Find (Disjoint Set Union) data structure.
 *
 * Used by Swendsen-Wang cluster algorithm for connected component
 * detection.  Supports path compression and union by rank.
 */
#ifndef LRGSG_UF_H
#define LRGSG_UF_H

#include <stddef.h>

typedef struct {
    size_t *parent;
    size_t *rank;
    size_t N;
} UnionFind;

/** Allocate and initialise a union-find structure for N elements. */
UnionFind *uf_create(size_t N);

/** Find the root of element x with path compression. */
size_t uf_find(UnionFind *uf, size_t x);

/** Union the sets containing x and y (union by rank). */
void uf_union(UnionFind *uf, size_t x, size_t y);

/** Reset all elements to singletons (reuse without realloc). */
void uf_reset(UnionFind *uf);

/** Free all memory. */
void uf_destroy(UnionFind *uf);

#endif /* LRGSG_UF_H */

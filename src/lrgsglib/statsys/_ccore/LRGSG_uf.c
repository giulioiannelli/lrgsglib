/**
 * @file LRGSG_uf.c
 * @brief Union-Find with path compression and union by rank.
 */
#include "LRGSG_uf.h"
#include <stdlib.h>
#include <string.h>

UnionFind *uf_create(size_t N) {
    UnionFind *uf = (UnionFind *)malloc(sizeof(UnionFind));
    if (!uf) return NULL;
    uf->N = N;
    uf->parent = (size_t *)malloc(N * sizeof(size_t));
    uf->rank   = (size_t *)calloc(N, sizeof(size_t));
    if (!uf->parent || !uf->rank) {
        free(uf->parent);
        free(uf->rank);
        free(uf);
        return NULL;
    }
    for (size_t i = 0; i < N; ++i)
        uf->parent[i] = i;
    return uf;
}

size_t uf_find(UnionFind *uf, size_t x) {
    while (uf->parent[x] != x) {
        uf->parent[x] = uf->parent[uf->parent[x]];  /* path halving */
        x = uf->parent[x];
    }
    return x;
}

void uf_union(UnionFind *uf, size_t x, size_t y) {
    size_t rx = uf_find(uf, x);
    size_t ry = uf_find(uf, y);
    if (rx == ry) return;
    if (uf->rank[rx] < uf->rank[ry]) {
        uf->parent[rx] = ry;
    } else if (uf->rank[rx] > uf->rank[ry]) {
        uf->parent[ry] = rx;
    } else {
        uf->parent[ry] = rx;
        uf->rank[rx]++;
    }
}

void uf_reset(UnionFind *uf) {
    for (size_t i = 0; i < uf->N; ++i)
        uf->parent[i] = i;
    memset(uf->rank, 0, uf->N * sizeof(size_t));
}

void uf_destroy(UnionFind *uf) {
    if (!uf) return;
    free(uf->parent);
    free(uf->rank);
    free(uf);
}

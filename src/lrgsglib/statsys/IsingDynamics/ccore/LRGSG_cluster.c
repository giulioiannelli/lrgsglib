/**
 * @file LRGSG_cluster.c
 * @brief Wolff and Swendsen-Wang cluster algorithms for Ising model.
 */
#include "LRGSG_cluster.h"
#include "LRGSG_uf.h"
#include "sfmtrng.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Wolff single-cluster step                                           */
/* ------------------------------------------------------------------ */

size_t wolff_step(size_t N, double T, spin_tp s, size_tp nlen, NodesEdges ne) {
    double inv_T = 1.0 / T;

    /* Pick random seed */
    size_t seed = (size_t)(RNG_dbl() * (double)N);
    if (seed >= N) seed = N - 1;

    /* BFS queue (worst case: entire lattice) */
    size_t *queue = (size_t *)malloc(N * sizeof(size_t));
    uint8_t *in_cluster = (uint8_t *)calloc(N, sizeof(uint8_t));
    if (!queue || !in_cluster) {
        free(queue);
        free(in_cluster);
        return 0;
    }

    queue[0] = seed;
    in_cluster[seed] = 1;
    size_t head = 0, tail = 1;

    while (head < tail) {
        size_t node = queue[head++];
        size_t deg = nlen[node];

        for (size_t j = 0; j < deg; ++j) {
            size_t nb = ne[node].neighbors[j];
            if (in_cluster[nb]) continue;

            double w = ne[node].weights[j];
            /* Bond activation: only when aligned relative to coupling */
            /* s[node]*s[nb]*sign(w) > 0  ⟺  s[node]*s[nb]*w > 0 */
            if ((double)s[node] * (double)s[nb] * w <= 0.0) continue;

            double p_add = 1.0 - exp(-2.0 * fabs(w) * inv_T);
            if (RNG_dbl() < p_add) {
                in_cluster[nb] = 1;
                queue[tail++] = nb;
            }
        }
    }

    /* Flip entire cluster */
    size_t cluster_size = tail;
    for (size_t i = 0; i < cluster_size; ++i) {
        s[queue[i]] = (int8_t)(-s[queue[i]]);
    }

    free(queue);
    free(in_cluster);
    return cluster_size;
}

/* ------------------------------------------------------------------ */
/* Wolff sweep: repeat until ~N spins flipped                          */
/* ------------------------------------------------------------------ */

size_t wolff_sweep(size_t N, double T, spin_tp s, size_tp nlen, NodesEdges ne) {
    size_t total_flipped = 0;
    while (total_flipped < N) {
        total_flipped += wolff_step(N, T, s, nlen, ne);
    }
    return total_flipped;
}

/* ------------------------------------------------------------------ */
/* Swendsen-Wang step                                                  */
/* ------------------------------------------------------------------ */

size_t sw_step(size_t N, double T, spin_tp s, size_tp nlen, NodesEdges ne) {
    double inv_T = 1.0 / T;

    UnionFind *uf = uf_create(N);
    if (!uf) return 0;

    /* Phase 1: Bond activation */
    for (size_t i = 0; i < N; ++i) {
        size_t deg = nlen[i];
        for (size_t j = 0; j < deg; ++j) {
            size_t nb = ne[i].neighbors[j];
            if (nb <= i) continue;  /* each edge once */

            double w = ne[i].weights[j];
            /* Only activate when aligned relative to coupling */
            if ((double)s[i] * (double)s[nb] * w <= 0.0) continue;

            double p_bond = 1.0 - exp(-2.0 * fabs(w) * inv_T);
            if (RNG_dbl() < p_bond) {
                uf_union(uf, i, nb);
            }
        }
    }

    /* Phase 2: Assign flip decisions per cluster root */
    /* Allocate a random bit for each potential root (lazy: one per node) */
    uint8_t *flip_decision = (uint8_t *)malloc(N * sizeof(uint8_t));
    if (!flip_decision) {
        uf_destroy(uf);
        return 0;
    }
    for (size_t i = 0; i < N; ++i) {
        flip_decision[i] = (RNG_dbl() < 0.5) ? 1 : 0;
    }

    /* Phase 3: Flip clusters and count them */
    size_t n_clusters = 0;
    uint8_t *seen_root = (uint8_t *)calloc(N, sizeof(uint8_t));
    if (!seen_root) {
        free(flip_decision);
        uf_destroy(uf);
        return 0;
    }

    for (size_t i = 0; i < N; ++i) {
        size_t root = uf_find(uf, i);
        if (!seen_root[root]) {
            seen_root[root] = 1;
            n_clusters++;
        }
        if (flip_decision[root]) {
            s[i] = (int8_t)(-s[i]);
        }
    }

    free(seen_root);
    free(flip_decision);
    uf_destroy(uf);
    return n_clusters;
}

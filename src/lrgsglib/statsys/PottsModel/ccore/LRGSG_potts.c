#include "LRGSG_potts.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

void potts_metropolis_sweep(size_t N, int32_t *s, int q, double T,
                            size_t *neigh_len, NodeEdges *node_edges) {
    for (size_t i = 0; i < N; ++i) {
        /* Pick random node */
        size_t nd = (size_t)(RNG_dbl() * N);
        if (nd >= N) nd = N - 1;

        int32_t current = s[nd];
        /* Propose new state != current */
        int32_t proposal = (int32_t)(RNG_dbl() * (q - 1));
        if (proposal >= current) proposal++;

        /* Energy change */
        double dE = 0.0;
        size_t n_nn = neigh_len[nd];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = node_edges[nd].neighbors[k];
            double w = node_edges[nd].weights[k];
            dE += w * (kronecker_delta_i32(current, s[j])
                     - kronecker_delta_i32(proposal, s[j]));
        }

        /* Accept/reject */
        if (dE <= 0.0 || (T > 0.0 && RNG_dbl() < exp(-dE / T))) {
            s[nd] = proposal;
        }
    }
}

double calc_potts_energy(size_t N, const int32_t *s, size_t *neigh_len,
                         NodeEdges *node_edges) {
    double E = 0.0;
    for (size_t i = 0; i < N; ++i) {
        size_t n_nn = neigh_len[i];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = node_edges[i].neighbors[k];
            if (j > i) {
                double w = node_edges[i].weights[k];
                E -= w * kronecker_delta_i32(s[i], s[j]);
            }
        }
    }
    return E;
}

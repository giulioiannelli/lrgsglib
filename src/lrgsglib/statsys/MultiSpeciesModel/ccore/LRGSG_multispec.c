#include "LRGSG_multispec.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

void multispec_metropolis_sweep(size_t N, int32_t *s, int k, int q, double T,
                                size_t *neigh_len, NodeEdges *node_edges) {
    /* Each "sweep" does N single-site updates */
    for (size_t step = 0; step < N; ++step) {
        size_t nd = (size_t)(RNG_dbl() * N);
        if (nd >= N) nd = N - 1;

        /* Choose random component to flip */
        int comp = (int)(RNG_dbl() * k);
        if (comp >= k) comp = k - 1;

        int32_t current = s[nd * k + comp];
        /* Propose new value != current */
        int32_t proposal = (int32_t)(RNG_dbl() * (q - 1));
        if (proposal >= current) proposal++;

        /* Energy change (identity interaction matrix for C backend) */
        double dE = 0.0;
        size_t n_nn = neigh_len[nd];
        for (size_t kk = 0; kk < n_nn; ++kk) {
            size_t j = node_edges[nd].neighbors[kk];
            double w = node_edges[nd].weights[kk];
            /* Only same-component interaction (identity matrix) */
            dE += w * (kronecker_delta_i32(current, s[j * k + comp])
                     - kronecker_delta_i32(proposal, s[j * k + comp]));
        }

        if (dE <= 0.0 || (T > 0.0 && RNG_dbl() < exp(-dE / T))) {
            s[nd * k + comp] = proposal;
        }
    }
}

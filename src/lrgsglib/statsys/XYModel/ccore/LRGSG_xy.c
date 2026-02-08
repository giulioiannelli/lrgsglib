#include "LRGSG_xy.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

void xy_metropolis_sweep(size_t N, double *theta, double T, double delta,
                         size_t *neigh_len, NodeEdges *node_edges) {
    for (size_t step = 0; step < N; ++step) {
        size_t nd = (size_t)(genrand_real2() * N);
        if (nd >= N) nd = N - 1;

        double current = theta[nd];
        double proposal = fmod(current + delta * (2.0 * genrand_real2() - 1.0), 2.0 * M_PI);
        if (proposal < 0.0) proposal += 2.0 * M_PI;

        double dE = 0.0;
        size_t n_nn = neigh_len[nd];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = node_edges[nd].neighbors[k];
            double w = node_edges[nd].weights[k];
            dE += w * (cos(current - theta[j]) - cos(proposal - theta[j]));
        }

        if (dE <= 0.0 || (T > 0.0 && genrand_real2() < exp(-dE / T))) {
            theta[nd] = proposal;
        }
    }
}

double calc_xy_energy(size_t N, const double *theta, size_t *neigh_len,
                      NodeEdges *node_edges) {
    double E = 0.0;
    for (size_t i = 0; i < N; ++i) {
        size_t n_nn = neigh_len[i];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = node_edges[i].neighbors[k];
            if (j > i) {
                double w = node_edges[i].weights[k];
                E -= w * cos(theta[i] - theta[j]);
            }
        }
    }
    return E;
}

double calc_xy_magnetisation(size_t N, const double *theta) {
    return calc_kuramoto_order_param(N, theta);
}

#include "LRGSG_heisenberg.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
#include <math.h>

void heisenberg_metropolis_sweep(size_t N, double *spin, double T, double delta,
                                 size_t *neigh_len, NodeEdges *node_edges) {
    for (size_t step = 0; step < N; ++step) {
        size_t nd = (size_t)(RNG_dbl() * N);
        if (nd >= N) nd = N - 1;

        double *current = &spin[nd * 3];
        double proposal[3];
        random_rotation_3d(current, proposal, delta);

        /* Energy change */
        double dE = 0.0;
        size_t n_nn = neigh_len[nd];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = node_edges[nd].neighbors[k];
            double w = node_edges[nd].weights[k];
            double *sj = &spin[j * 3];
            dE += w * (dot3(current, sj) - dot3(proposal, sj));
        }

        if (dE <= 0.0 || (T > 0.0 && RNG_dbl() < exp(-dE / T))) {
            current[0] = proposal[0];
            current[1] = proposal[1];
            current[2] = proposal[2];
        }
    }
}

double calc_heisenberg_energy(size_t N, const double *spin, size_t *neigh_len,
                              NodeEdges *node_edges) {
    double E = 0.0;
    for (size_t i = 0; i < N; ++i) {
        const double *si = &spin[i * 3];
        size_t n_nn = neigh_len[i];
        for (size_t k = 0; k < n_nn; ++k) {
            size_t j = node_edges[i].neighbors[k];
            if (j > i) {
                double w = node_edges[i].weights[k];
                E -= w * dot3(si, &spin[j * 3]);
            }
        }
    }
    return E;
}

double calc_heisenberg_magnetisation(size_t N, const double *spin) {
    double mx = 0.0, my = 0.0, mz = 0.0;
    for (size_t i = 0; i < N; ++i) {
        mx += spin[i * 3];
        my += spin[i * 3 + 1];
        mz += spin[i * 3 + 2];
    }
    mx /= (double)N;
    my /= (double)N;
    mz /= (double)N;
    return sqrt(mx * mx + my * my + mz * mz);
}

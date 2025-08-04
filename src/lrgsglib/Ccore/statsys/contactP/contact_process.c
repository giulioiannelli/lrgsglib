#include "LRGSG_utils.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double infection_rate(size_t j, spin_tp s, size_tp nlen, NodeEdges ne) {
    double r = 0.0;
    for (size_t k = 0; k < nlen; k++) {
        size_t i = ne.neighbors[k];
        double w = ne.weights[k];
        if (w > 0.0 && s[i])
            r += w;
    }
    return r;
}

static double recovery_rate(size_t j, double mu, spin_tp s, size_tp nlen,
                            NodeEdges ne) {
    double r = mu;
    for (size_t k = 0; k < nlen; k++) {
        size_t i = ne.neighbors[k];
        double w = ne.weights[k];
        if (w < 0.0 && s[i])
            r += -w;
    }
    return r;
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s N mu steps edge_file\n", argv[0]);
        return EXIT_FAILURE;
    }

    size_t N = strtozu(argv[1]);
    double mu = strtod(argv[2], NULL);
    size_t steps = strtozu(argv[3]);
    const char *edge_file = argv[4];

    srand((unsigned)time(NULL));

    spin_tp state = __chMalloc(N * sizeof(*state));
    for (size_t i = 0; i < N; i++)
        state[i] = rand() & 1U;

    Edges edges;
    NodesEdges node_edges;
    size_tp neigh_len;
    process_edges(edge_file, N, &edges, &node_edges, &neigh_len);

    for (size_t t = 0; t < steps; t++) {
        for (size_t j = 0; j < N; j++) {
            double prob;
            if (state[j]) {
                double r = recovery_rate(j, mu, state, neigh_len[j],
                                         node_edges[j]);
                prob = 1.0 - exp(-r);
                if (((double)rand() / RAND_MAX) < prob)
                    state[j] = 0;
            } else {
                double r = infection_rate(j, state, neigh_len[j],
                                          node_edges[j]);
                prob = 1.0 - exp(-r);
                if (((double)rand() / RAND_MAX) < prob)
                    state[j] = 1;
            }
        }
    }

    fwrite(state, sizeof(*state), N, stdout);

    free(neigh_len);
    free(edges);
    for (size_t i = 0; i < N; i++) {
        free(node_edges[i].neighbors);
        free(node_edges[i].weights);
    }
    free(node_edges);
    free(state);

    return 0;
}


#include "LRGSG_cp.h"

double cp_infection_rate(size_t node, spin_tp state, size_tp degree, NodeEdges edges) {
    (void)node;
    double rate = 0.0;
    for (size_t idx = 0; idx < degree; ++idx) {
        size_t neighbour = edges.neighbors[idx];
        double weight = edges.weights[idx];
        if (weight > 0.0 && state[neighbour]) {
            rate += weight;
        }
    }
    return rate;
}

double cp_recovery_rate(size_t node, double mu, spin_tp state, size_tp degree, NodeEdges edges) {
    (void)node;
    double rate = mu;
    for (size_t idx = 0; idx < degree; ++idx) {
        size_t neighbour = edges.neighbors[idx];
        double weight = edges.weights[idx];
        if (weight < 0.0 && state[neighbour]) {
            rate += -weight;
        }
    }
    return rate;
}

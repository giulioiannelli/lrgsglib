#include "LRGSG_cp.h"

#include <math.h>
#include <strings.h>

double cp_infection_rate(size_t node, spin_tp state, size_t degree, NodeEdges edges) {
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

cp_activation_t cp_activation_from_string(const char *name) {
    if (!name) {
        return CP_ACTIVATION_TANH;
    }
    if (strcasecmp(name, "relu") == 0) {
        return CP_ACTIVATION_RELU;
    }
    return CP_ACTIVATION_TANH;
}

double cp_linear_input(double gamma, spin_tp state, size_t degree, NodeEdges edges) {
    double total = 0.0;
    for (size_t idx = 0; idx < degree; ++idx) {
        size_t neighbour = edges.neighbors[idx];
        double weight = edges.weights[idx];
        total += weight * (double)state[neighbour];
    }
    return gamma * total;
}

double cp_activation_probability(double lambda, cp_activation_t activation) {
    double prob;
    switch (activation) {
        case CP_ACTIVATION_RELU:
            prob = lambda > 0.0 ? lambda : 0.0;
            break;
        case CP_ACTIVATION_TANH:
        default:
            prob = 0.5 * (tanh(lambda) + 1.0);
            break;
    }
    if (prob < 0.0) {
        prob = 0.0;
    } else if (prob > 1.0) {
        prob = 1.0;
    }
    return prob;
}

double cp_recovery_rate(size_t node, double mu, spin_tp state, size_t degree, NodeEdges edges) {
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

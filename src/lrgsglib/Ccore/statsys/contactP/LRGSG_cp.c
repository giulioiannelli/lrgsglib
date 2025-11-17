#include "LRGSG_cp.h"
#include "LRGSG_utils.h"

#include <math.h>
#include <stdlib.h>
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

/* Activation function: RELU - clipped to [0,1] */
static double cp_activation_relu(double lambda) {
    if (lambda <= 0.0) return 0.0;
    if (lambda >= 1.0) return 1.0;
    return lambda;
}

/* Activation function: TANH - naturally in [0,1] */
static double cp_activation_tanh(double lambda) {
    return 0.5 * (tanh(lambda) + 1.0);
}

/* Get activation function pointer from enum */
cp_activation_func_t cp_get_activation_function(cp_activation_t activation) {
    switch (activation) {
        case CP_ACTIVATION_RELU:
            return cp_activation_relu;
        case CP_ACTIVATION_TANH:
        default:
            return cp_activation_tanh;
    }
}

/* DEPRECATED: Use cp_get_activation_function() instead */
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

/* Active border management functions */

void cp_init_active_border(ActiveBorder *ab, size_t N) {
    ab->capacity = N;
    ab->count = 0;
    ab->nodes = __chMalloc(N * sizeof(*ab->nodes));
    ab->in_border = __chCalloc(N, sizeof(*ab->in_border));
    ab->position = __chMalloc(N * sizeof(*ab->position));
}

void cp_free_active_border(ActiveBorder *ab) {
    free(ab->nodes);
    free(ab->in_border);
    free(ab->position);
}

void cp_add_to_border(ActiveBorder *ab, size_t node) {
    if (!ab->in_border[node]) {
        ab->position[node] = ab->count;  // Store position
        ab->nodes[ab->count++] = node;
        ab->in_border[node] = 1;
    }
}

void cp_remove_from_border(ActiveBorder *ab, size_t idx) {
    size_t node = ab->nodes[idx];
    ab->in_border[node] = 0;
    ab->count--;
    if (idx < ab->count) {
        // Move last element to this position
        size_t last_node = ab->nodes[ab->count];
        ab->nodes[idx] = last_node;
        ab->position[last_node] = idx;  // Update position of moved node
    }
}

int cp_has_active_neighbor(const spin_tp state, size_t degree, const NodeEdges edges_node) {
    for (size_t j = 0; j < degree; ++j) {
        if (state[edges_node.neighbors[j]] == 1) {
            return 1;
        }
    }
    return 0;
}

void cp_build_initial_border(ActiveBorder *ab, const spin_tp state, size_t N,
                             const NodesEdges node_edges, const size_tp neigh_len) {
    // First pass: add all active nodes
    for (size_t i = 0; i < N; ++i) {
        if (state[i] == 1) {
            cp_add_to_border(ab, i);
        }
    }
    
    // Second pass: add all neighbors of active nodes
    for (size_t i = 0; i < N; ++i) {
        if (state[i] == 1) {
            size_t degree = neigh_len[i];
            NodeEdges edges_node = node_edges[i];
            for (size_t j = 0; j < degree; ++j) {
                size_t neighbor = edges_node.neighbors[j];
                cp_add_to_border(ab, neighbor);
            }
        }
    }
}

void cp_update_border_after_flip(ActiveBorder *ab, size_t node, int8_t new_state,
                                 const spin_tp state __attribute__((unused)), size_t N __attribute__((unused)),
                                 const NodesEdges node_edges, const size_tp neigh_len) {
    if (new_state == 1) {
        // Node became active: add all its neighbors to border
        size_t degree = neigh_len[node];
        NodeEdges edges_node = node_edges[node];
        for (size_t j = 0; j < degree; ++j) {
            size_t neighbor = edges_node.neighbors[j];
            cp_add_to_border(ab, neighbor);
        }
    }
    // When node becomes inactive, we DON'T remove anything from border
    // This trades some wasted samples for much faster border maintenance
    // The border will naturally shrink as activity moves away
}

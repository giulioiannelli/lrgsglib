#include "LRGSG_cp.h"
#include "LRGSG_utils.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>

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

void cp_lambda_init(double gamma, const spin_tp state, size_t N,
                    const NodesEdges node_edges, const size_tp neigh_len,
                    double *lambda) {
    for (size_t i = 0; i < N; ++i) {
        double total = 0.0;
        size_t degree = neigh_len[i];
        NodeEdges edges_node = node_edges[i];
        for (size_t j = 0; j < degree; ++j) {
            size_t neighbor = edges_node.neighbors[j];
            total += edges_node.weights[j] * (double)state[neighbor];
        }
        lambda[i] = gamma * total;
    }
}

void cp_lambda_update_neighbors(double gamma, size_t node, int delta_state,
                                const NodesEdges node_edges, const size_tp neigh_len,
                                double *lambda) {
    if (delta_state == 0) {
        return;
    }
    size_t degree = neigh_len[node];
    NodeEdges edges_node = node_edges[node];
    double delta = gamma * (double)delta_state;
    for (size_t j = 0; j < degree; ++j) {
        size_t neighbor = edges_node.neighbors[j];
        lambda[neighbor] += delta * edges_node.weights[j];
    }
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

int cp_reached_absorbing_state(size_t sum, size_t N, size_t t, size_t steps) {
    if (sum == 0 || sum == N) {
        fprintf(stderr, "Absorbing state reached (sum=%zu) at step %zu/%zu\n", sum, t, steps);
        return 1;
    }
    return 0;
}

/* Frontier tracking implementation (three-list system) */

static inline void cp_frontier_add_active_node(cp_frontier_t *frontier, size_t node) {
    frontier->pos_active[node] = frontier->active_count;
    frontier->active_list[frontier->active_count++] = node;
    frontier->node_status[node] = NODE_ACTIVE;
}

static inline void cp_frontier_remove_active_node(cp_frontier_t *frontier, size_t node) {
    size_t idx = frontier->pos_active[node];
    size_t last_idx = --frontier->active_count;
    if (idx < last_idx) {
        size_t swapped = frontier->active_list[last_idx];
        frontier->active_list[idx] = swapped;
        frontier->pos_active[swapped] = idx;
    }
}

static inline void cp_frontier_add_boundary_node(cp_frontier_t *frontier, size_t node) {
    frontier->pos_boundary[node] = frontier->boundary_count;
    frontier->boundary_list[frontier->boundary_count++] = node;
    frontier->node_status[node] = NODE_BOUNDARY;
}

static inline void cp_frontier_remove_boundary_node(cp_frontier_t *frontier, size_t node) {
    size_t idx = frontier->pos_boundary[node];
    size_t last_idx = --frontier->boundary_count;
    if (idx < last_idx) {
        size_t swapped = frontier->boundary_list[last_idx];
        frontier->boundary_list[idx] = swapped;
        frontier->pos_boundary[swapped] = idx;
    }
    frontier->node_status[node] = NODE_INACTIVE;
}

void cp_frontier_init(cp_frontier_t *frontier, size_t N) {
    frontier->active_list = __chMalloc(N * sizeof(size_t));
    frontier->boundary_list = __chMalloc(N * sizeof(size_t));
    frontier->pos_active = __chMalloc(N * sizeof(size_t));
    frontier->pos_boundary = __chMalloc(N * sizeof(size_t));
    frontier->node_status = __chCalloc(N, sizeof(cp_node_status_t));
    frontier->active_neighbor_count = __chCalloc(N, sizeof(size_t));
    frontier->N = N;
    frontier->active_count = 0;
    frontier->boundary_count = 0;
    frontier->use_lambda_boundary = 0;
}

void cp_frontier_free(cp_frontier_t *frontier) {
    free(frontier->active_list);
    free(frontier->boundary_list);
    free(frontier->pos_active);
    free(frontier->pos_boundary);
    free(frontier->node_status);
    free(frontier->active_neighbor_count);
}

void cp_frontier_build(cp_frontier_t *frontier, const spin_tp state, const double *lambda,
                       size_t N, const NodesEdges node_edges, const size_tp neigh_len) {
    frontier->N = N;
    frontier->active_count = 0;
    frontier->boundary_count = 0;
    memset(frontier->node_status, 0, frontier->N * sizeof(cp_node_status_t));
    memset(frontier->active_neighbor_count, 0, frontier->N * sizeof(size_t));

    /* First pass: build active list and count active neighbors */
    for (size_t i = 0; i < N; ++i) {
        if (state[i]) {
            cp_frontier_add_active_node(frontier, i);

            size_t degree = neigh_len[i];
            NodeEdges edges_node = node_edges[i];
            for (size_t j = 0; j < degree; ++j) {
                size_t neighbor = edges_node.neighbors[j];
                frontier->active_neighbor_count[neighbor]++;
            }
        }
    }

    /* Second pass: boundary = inactive nodes with active neighbors (and positive lambda if enabled) */
    for (size_t i = 0; i < N; ++i) {
        if (!state[i] && frontier->active_neighbor_count[i] > 0) {
            if (!frontier->use_lambda_boundary || lambda[i] > 0.0) {
                cp_frontier_add_boundary_node(frontier, i);
            }
        }
    }
}

void cp_frontier_node_activate(cp_frontier_t *frontier, size_t node,
                               const NodesEdges node_edges, const size_tp neigh_len,
                               const spin_tp state, const double *lambda) {
    /* Remove from boundary list if present */
    if (frontier->node_status[node] == NODE_BOUNDARY) {
        cp_frontier_remove_boundary_node(frontier, node);
    }

    /* Add to active list */
    cp_frontier_add_active_node(frontier, node);

    /* Update neighbors: increment count and add inactive ones to boundary */
    size_t degree = neigh_len[node];
    NodeEdges edges_node = node_edges[node];
    for (size_t j = 0; j < degree; ++j) {
        size_t neighbor = edges_node.neighbors[j];
        frontier->active_neighbor_count[neighbor]++;

        if (!state[neighbor] && frontier->node_status[neighbor] == NODE_INACTIVE) {
            if (!frontier->use_lambda_boundary || lambda[neighbor] > 0.0) {
                cp_frontier_add_boundary_node(frontier, neighbor);
            }
        }
    }
}

void cp_frontier_node_deactivate(cp_frontier_t *frontier, size_t node,
                                 const NodesEdges node_edges, const size_tp neigh_len,
                                 const spin_tp state, const double *lambda) {
    /* Remove from active list */
    if (frontier->node_status[node] == NODE_ACTIVE) {
        cp_frontier_remove_active_node(frontier, node);
    }

    /* Update neighbors: decrement count and remove from boundary if needed */
    size_t degree = neigh_len[node];
    NodeEdges edges_node = node_edges[node];
    for (size_t j = 0; j < degree; ++j) {
        size_t neighbor = edges_node.neighbors[j];
        frontier->active_neighbor_count[neighbor]--;

        if (!state[neighbor] && frontier->node_status[neighbor] == NODE_BOUNDARY) {
            if (frontier->active_neighbor_count[neighbor] == 0 ||
                (frontier->use_lambda_boundary && lambda[neighbor] <= 0.0)) {
                cp_frontier_remove_boundary_node(frontier, neighbor);
            }
        }
    }

    /* Update this node's status */
    if (frontier->active_neighbor_count[node] > 0 &&
        (!frontier->use_lambda_boundary || lambda[node] > 0.0)) {
        cp_frontier_add_boundary_node(frontier, node);
    } else {
        frontier->node_status[node] = NODE_INACTIVE;
    }
}

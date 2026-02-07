#ifndef HOLME_KIM_HH_
#define HOLME_KIM_HH_

#include "graph_tool.hh"
#include <random>
#include <vector>
#include <algorithm>
#include <set>

using namespace boost;
using namespace graph_tool;

/**
 * Create a Holme-Kim powerlaw cluster graph.
 *
 * Algorithm:
 * 1. Start with m+1 fully connected nodes
 * 2. For each new node:
 *    a. Connect to first target via preferential attachment
 *    b. For remaining m-1 connections:
 *       - With probability p: try to connect to a neighbor of last target (triangle)
 *       - Otherwise: connect via preferential attachment
 *
 * This produces scale-free networks with tunable clustering.
 *
 * @param g The graph to populate (should be empty)
 * @param n Total number of nodes
 * @param m Number of edges per new node
 * @param p Probability of triad formation (0=pure BA, 1=max clustering)
 * @param seed Random seed (0 = use random device)
 */
template <class Graph>
void create_holme_kim(Graph& g, int n, int m, double p, unsigned long seed) {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (n <= 0)
        throw std::invalid_argument("n must be positive");
    if (m < 1)
        throw std::invalid_argument("m must be >= 1");
    if (m >= n)
        throw std::invalid_argument("m must be < n");
    if (p < 0.0 || p > 1.0)
        throw std::invalid_argument("p must be in [0, 1]");

    // Initialize RNG
    std::mt19937_64 rng;
    if (seed == 0) {
        std::random_device rd;
        rng.seed(rd());
    } else {
        rng.seed(seed);
    }

    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // Create initial m+1 fully connected nodes
    int initial_nodes = m + 1;
    std::vector<vertex_t> vertices;
    vertices.reserve(n);

    for (int i = 0; i < initial_nodes; ++i) {
        vertices.push_back(add_vertex(g));
    }

    // Connect initial nodes (complete graph)
    for (int i = 0; i < initial_nodes; ++i) {
        for (int j = i + 1; j < initial_nodes; ++j) {
            add_edge(vertices[i], vertices[j], g);
        }
    }

    // Degree list for preferential attachment (each node appears degree times)
    std::vector<int> degree_list;
    degree_list.reserve(2 * m * n);

    // Initial nodes each have degree = initial_nodes - 1
    for (int i = 0; i < initial_nodes; ++i) {
        for (int d = 0; d < initial_nodes - 1; ++d) {
            degree_list.push_back(i);
        }
    }

    // Track neighbors for triangle formation
    std::vector<std::set<int>> neighbors(n);
    for (int i = 0; i < initial_nodes; ++i) {
        for (int j = 0; j < initial_nodes; ++j) {
            if (i != j) {
                neighbors[i].insert(j);
            }
        }
    }

    // Add remaining nodes
    for (int new_idx = initial_nodes; new_idx < n; ++new_idx) {
        vertex_t v_new = add_vertex(g);
        vertices.push_back(v_new);

        std::set<int> targets;  // Set to avoid duplicates
        int last_target = -1;

        // We need m connections
        int needed = m;

        while (static_cast<int>(targets.size()) < needed) {
            bool use_triangle = false;

            // First connection always via PA
            // For subsequent: with probability p, try triangle formation
            if (!targets.empty() && last_target >= 0 && uniform(rng) < p) {
                use_triangle = true;
            }

            if (use_triangle && !neighbors[last_target].empty()) {
                // Try to find a neighbor of last_target not already connected
                std::vector<int> available;
                for (int nbr : neighbors[last_target]) {
                    if (targets.find(nbr) == targets.end()) {
                        available.push_back(nbr);
                    }
                }

                if (!available.empty()) {
                    std::uniform_int_distribution<size_t> dist(0, available.size() - 1);
                    int target = available[dist(rng)];
                    targets.insert(target);
                    last_target = target;
                    continue;
                }
                // If no available neighbors, fall through to PA
            }

            // Preferential attachment
            bool found = false;
            int max_attempts = 1000;
            for (int attempt = 0; attempt < max_attempts && !found; ++attempt) {
                std::uniform_int_distribution<size_t> dist(0, degree_list.size() - 1);
                int target = degree_list[dist(rng)];
                if (targets.find(target) == targets.end()) {
                    targets.insert(target);
                    last_target = target;
                    found = true;
                }
            }
            if (!found) {
                // Fallback: pick any node not yet connected
                for (int i = 0; i < new_idx && static_cast<int>(targets.size()) < needed; ++i) {
                    if (targets.find(i) == targets.end()) {
                        targets.insert(i);
                        last_target = i;
                        break;
                    }
                }
            }
        }

        // Add edges
        for (int target : targets) {
            add_edge(v_new, vertices[target], g);
            // Update degree list
            degree_list.push_back(new_idx);
            degree_list.push_back(target);
            // Update neighbor sets
            neighbors[new_idx].insert(target);
            neighbors[target].insert(new_idx);
        }
    }
}

#endif // HOLME_KIM_HH_

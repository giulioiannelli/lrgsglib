#ifndef EXTENDED_BARABASI_ALBERT_HH_
#define EXTENDED_BARABASI_ALBERT_HH_

#include "graph_tool.hh"
#include <random>
#include <vector>

using namespace boost;
using namespace graph_tool;

/**
 * Create an Extended Barabasi-Albert graph with initial attractiveness.
 *
 * Attachment probability: P(i) ~ (k_i + a)
 * When a=0, this reduces to standard BA.
 * Higher a leads to more uniform degree distribution.
 *
 * @param g The graph to populate (should be empty)
 * @param n Total number of nodes
 * @param m Number of edges per new node
 * @param a Initial attractiveness parameter (>= 0)
 * @param seed Random seed (0 = use random device)
 */
template <class Graph>
void create_extended_barabasi_albert(Graph& g, int n, int m, double a, unsigned long seed) {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (n <= 0)
        throw std::invalid_argument("n must be positive");
    if (m < 1)
        throw std::invalid_argument("m must be >= 1");
    if (m >= n)
        throw std::invalid_argument("m must be < n");
    if (a < 0)
        throw std::invalid_argument("a must be >= 0");

    // Initialize RNG
    std::mt19937_64 rng;
    if (seed == 0) {
        std::random_device rd;
        rng.seed(rd());
    } else {
        rng.seed(seed);
    }

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

    // Track degrees for weighted sampling
    std::vector<double> degrees(n, 0.0);
    for (int i = 0; i < initial_nodes; ++i) {
        degrees[i] = static_cast<double>(initial_nodes - 1);
    }

    // Add remaining nodes with preferential attachment + attractiveness
    for (int new_idx = initial_nodes; new_idx < n; ++new_idx) {
        vertex_t v_new = add_vertex(g);
        vertices.push_back(v_new);

        // Compute weights: degree[i] + a for all existing nodes
        std::vector<double> weights(new_idx);
        double total_weight = 0.0;
        for (int i = 0; i < new_idx; ++i) {
            weights[i] = degrees[i] + a;
            total_weight += weights[i];
        }

        // Select m unique targets via weighted sampling
        std::vector<int> targets;
        targets.reserve(m);
        std::vector<bool> selected(new_idx, false);

        while (static_cast<int>(targets.size()) < m) {
            // Weighted random selection
            std::uniform_real_distribution<double> dist(0.0, total_weight);
            double r = dist(rng);

            double cumsum = 0.0;
            for (int i = 0; i < new_idx; ++i) {
                if (selected[i]) continue;
                cumsum += weights[i];
                if (r <= cumsum) {
                    targets.push_back(i);
                    selected[i] = true;
                    total_weight -= weights[i];
                    break;
                }
            }
        }

        // Add edges to selected targets
        for (int target : targets) {
            add_edge(v_new, vertices[target], g);
            degrees[target] += 1.0;
        }
        degrees[new_idx] = static_cast<double>(m);
    }
}

#endif // EXTENDED_BARABASI_ALBERT_HH_

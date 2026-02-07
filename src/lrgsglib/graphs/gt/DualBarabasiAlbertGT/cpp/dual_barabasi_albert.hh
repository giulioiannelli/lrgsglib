#ifndef DUAL_BARABASI_ALBERT_HH_
#define DUAL_BARABASI_ALBERT_HH_

#include "graph_tool.hh"
#include <random>
#include <vector>
#include <algorithm>

using namespace boost;
using namespace graph_tool;

/**
 * Create a Dual Barabasi-Albert graph with two attachment modes.
 *
 * With probability p, new nodes use m1 edges.
 * With probability 1-p, new nodes use m2 edges.
 *
 * @param g The graph to populate (should be empty)
 * @param n Total number of nodes
 * @param m1 Number of edges in mode 1
 * @param m2 Number of edges in mode 2
 * @param p Probability of using m1 (vs m2)
 * @param seed Random seed (0 = use random device)
 */
template <class Graph>
void create_dual_barabasi_albert(Graph& g, int n, int m1, int m2, double p, unsigned long seed) {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (n <= 0)
        throw std::invalid_argument("n must be positive");
    if (m1 < 1 || m2 < 1)
        throw std::invalid_argument("m1 and m2 must be >= 1");
    if (std::max(m1, m2) >= n)
        throw std::invalid_argument("max(m1, m2) must be < n");
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

    // Create initial max(m1, m2)+1 fully connected nodes
    int max_m = std::max(m1, m2);
    int initial_nodes = max_m + 1;
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

    // Degree list for weighted sampling
    std::vector<int> degree_list;
    degree_list.reserve(2 * max_m * n);

    for (int i = 0; i < initial_nodes; ++i) {
        for (int d = 0; d < initial_nodes - 1; ++d) {
            degree_list.push_back(i);
        }
    }

    std::uniform_real_distribution<double> p_dist(0.0, 1.0);

    // Add remaining nodes
    for (int new_idx = initial_nodes; new_idx < n; ++new_idx) {
        vertex_t v_new = add_vertex(g);
        vertices.push_back(v_new);

        // Choose m based on probability p
        int m = (p_dist(rng) < p) ? m1 : m2;

        // Select m unique targets via preferential attachment
        std::vector<int> targets;
        targets.reserve(m);

        while (static_cast<int>(targets.size()) < m) {
            std::uniform_int_distribution<size_t> dist(0, degree_list.size() - 1);
            int target = degree_list[dist(rng)];

            if (std::find(targets.begin(), targets.end(), target) == targets.end()) {
                targets.push_back(target);
            }
        }

        // Add edges to selected targets
        for (int target : targets) {
            add_edge(v_new, vertices[target], g);
            degree_list.push_back(new_idx);
            degree_list.push_back(target);
        }
    }
}

#endif // DUAL_BARABASI_ALBERT_HH_

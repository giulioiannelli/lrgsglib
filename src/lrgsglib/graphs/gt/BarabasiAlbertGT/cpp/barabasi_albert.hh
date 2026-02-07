#ifndef BARABASI_ALBERT_HH_
#define BARABASI_ALBERT_HH_

#include "graph_tool.hh"
#include <random>
#include <vector>
#include <algorithm>

using namespace boost;
using namespace graph_tool;

/**
 * Create a Barabasi-Albert preferential attachment graph.
 *
 * Algorithm:
 * 1. Start with m+1 fully connected nodes
 * 2. Add new nodes one at a time
 * 3. Each new node connects to m existing nodes
 * 4. Connection probability is proportional to degree
 *
 * Uses degree-list trick for O(1) preferential attachment sampling.
 *
 * @param g The graph to populate (should be empty)
 * @param n Total number of nodes
 * @param m Number of edges per new node
 * @param seed Random seed (0 = use random device)
 */
template <class Graph>
void create_barabasi_albert(Graph& g, int n, int m, unsigned long seed) {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (n <= 0)
        throw std::invalid_argument("n must be positive");
    if (m < 1)
        throw std::invalid_argument("m must be >= 1");
    if (m >= n)
        throw std::invalid_argument("m must be < n");

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

    // Degree list for weighted sampling (each node appears degree times)
    // This is the standard trick for O(1) preferential attachment sampling
    std::vector<int> degree_list;
    degree_list.reserve(2 * m * n);  // Pre-allocate

    // Initial nodes each have degree = initial_nodes - 1
    for (int i = 0; i < initial_nodes; ++i) {
        for (int d = 0; d < initial_nodes - 1; ++d) {
            degree_list.push_back(i);
        }
    }

    // Add remaining nodes with preferential attachment
    for (int new_idx = initial_nodes; new_idx < n; ++new_idx) {
        vertex_t v_new = add_vertex(g);
        vertices.push_back(v_new);

        // Select m unique targets via preferential attachment
        std::vector<int> targets;
        targets.reserve(m);

        while (static_cast<int>(targets.size()) < m) {
            // Sample from degree list (O(1) preferential attachment)
            std::uniform_int_distribution<size_t> dist(0, degree_list.size() - 1);
            int target = degree_list[dist(rng)];

            // Check for duplicates (usually fast since m is small)
            if (std::find(targets.begin(), targets.end(), target) == targets.end()) {
                targets.push_back(target);
            }
        }

        // Add edges to selected targets
        for (int target : targets) {
            add_edge(v_new, vertices[target], g);
            // Update degree list: both endpoints gain degree 1
            degree_list.push_back(new_idx);
            degree_list.push_back(target);
        }
    }
}

#endif // BARABASI_ALBERT_HH_

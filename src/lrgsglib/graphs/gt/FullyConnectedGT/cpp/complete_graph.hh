#ifndef COMPLETE_GRAPH_HH_
#define COMPLETE_GRAPH_HH_

#include "graph_tool.hh"
#include <vector>

using namespace boost;
using namespace graph_tool;

/**
 * Create a complete (fully connected) graph.
 *
 * Every node is connected to every other node.
 * Total edges: n*(n-1)/2
 *
 * @param g The graph to populate (should be empty)
 * @param n Number of nodes
 */
template <class Graph>
void create_complete_graph(Graph& g, int n) {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (n <= 0)
        throw std::invalid_argument("n must be positive");

    // Add all vertices
    std::vector<vertex_t> vertices;
    vertices.reserve(n);

    for (int i = 0; i < n; ++i) {
        vertices.push_back(add_vertex(g));
    }

    // Add all edges (complete graph)
    // For large n, batch edge addition is faster
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            add_edge(vertices[i], vertices[j], g);
        }
    }
}

#endif // COMPLETE_GRAPH_HH_

#ifndef MULTIPLICATIVE_CASCADE_HH_
#define MULTIPLICATIVE_CASCADE_HH_

#include "graph_tool.hh"
#include <vector>
#include <unordered_map>

using namespace boost;
using namespace graph_tool;

/**
 * Build grid edges between selected nodes on a 2D lattice.
 *
 * For each selected node, checks the right and down neighbors.
 * If the neighbor is also selected, adds an edge. This finds
 * each edge exactly once (no duplicates).
 *
 * @param g         The graph (must already have n_selected vertices)
 * @param rows      Row coordinates of selected nodes
 * @param cols      Column coordinates of selected nodes
 * @param n_selected Number of selected nodes
 * @param size      Grid side length
 * @param periodic  Whether to use periodic boundary conditions
 */
template <class Graph>
void build_cascade_edges(
    Graph& g,
    const std::vector<int>& rows,
    const std::vector<int>& cols,
    int n_selected,
    int size,
    bool periodic
) {
    // Build hash map: (row * size + col) -> vertex index
    std::unordered_map<long long, int> coord_to_idx;
    coord_to_idx.reserve(n_selected * 2);

    for (int i = 0; i < n_selected; ++i) {
        long long key = static_cast<long long>(rows[i]) * size + cols[i];
        coord_to_idx[key] = i;
    }

    // Check right and down neighbors for each selected node
    for (int i = 0; i < n_selected; ++i) {
        int r = rows[i];
        int c = cols[i];

        // Right neighbor (c + 1)
        int nc_right = c + 1;
        if (periodic) {
            nc_right = nc_right % size;
        }
        if (nc_right < size) {
            long long key_right = static_cast<long long>(r) * size + nc_right;
            auto it = coord_to_idx.find(key_right);
            if (it != coord_to_idx.end()) {
                add_edge(vertex(i, g), vertex(it->second, g), g);
            }
        }

        // Down neighbor (r + 1)
        int nr_down = r + 1;
        if (periodic) {
            nr_down = nr_down % size;
        }
        if (nr_down < size) {
            long long key_down = static_cast<long long>(nr_down) * size + c;
            auto it = coord_to_idx.find(key_down);
            if (it != coord_to_idx.end()) {
                add_edge(vertex(i, g), vertex(it->second, g), g);
            }
        }
    }
}

#endif // MULTIPLICATIVE_CASCADE_HH_

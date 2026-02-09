#ifndef TRIANGULAR_LATTICE_HH_
#define TRIANGULAR_LATTICE_HH_

#include "graph_tool.hh"

using namespace boost;
using namespace graph_tool;

// Template function to create a triangular lattice graph
template <class Graph>
void create_triangular_lattice(Graph& g, int width, int height) {
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

    if (width <= 0 || height <= 0)
        throw std::invalid_argument("Width and height must be positive integers");

    std::vector<vertex_t> last_row(width);

    // Create vertices and edges for the triangular lattice
    for (int i = 0; i < height; ++i) {
        std::vector<vertex_t> this_row(width);
        for (int j = 0; j < width; ++j) {
            this_row[j] = add_vertex(g);
            if (j > 0) {
                add_edge(this_row[j - 1], this_row[j], g);  // horizontal edge
                if (i > 0) {
                    add_edge(last_row[j - 1], this_row[j], g);  // diagonal edge
                }
            }
            if (i > 0) {
                add_edge(last_row[j], this_row[j], g);  // vertical edge
            }
        }
        last_row = this_row;
    }
}

#endif // TRIANGULAR_LATTICE_HH_

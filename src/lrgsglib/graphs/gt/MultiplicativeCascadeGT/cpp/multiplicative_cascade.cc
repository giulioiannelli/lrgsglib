#include "multiplicative_cascade.hh"
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>

// Binding function: extracts Python lists -> C++ vectors, dispatches
void build_cascade_edges_bind(
    GraphInterface& gi,
    boost::python::list py_rows,
    boost::python::list py_cols,
    int n_selected,
    int size,
    bool periodic
) {
    // Convert Python lists to C++ vectors
    std::vector<int> rows(n_selected);
    std::vector<int> cols(n_selected);
    for (int i = 0; i < n_selected; ++i) {
        rows[i] = boost::python::extract<int>(py_rows[i]);
        cols[i] = boost::python::extract<int>(py_cols[i]);
    }

    gt_dispatch<>()(
        [&](auto& g) {
            build_cascade_edges(g, rows, cols, n_selected, size, periodic);
        },
        all_graph_views()
    )(gi.get_graph_view());
}

BOOST_PYTHON_MODULE(multiplicative_cascade)
{
    using namespace boost::python;

    def("build_cascade_edges", build_cascade_edges_bind,
        "Build grid edges between selected nodes on a 2D lattice.\n\n"
        "Args:\n"
        "    gi: GraphInterface (internal graph object)\n"
        "    rows: list of row coordinates\n"
        "    cols: list of column coordinates\n"
        "    n_selected: number of selected nodes\n"
        "    size: grid side length\n"
        "    periodic: whether to use periodic boundary conditions\n");
}

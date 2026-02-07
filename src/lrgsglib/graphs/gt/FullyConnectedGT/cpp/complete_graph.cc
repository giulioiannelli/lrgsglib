#include "complete_graph.hh"
#include <boost/python.hpp>

// Binding function for complete graph
void complete_graph_bind(GraphInterface& gi, int n)
{
    gt_dispatch<>()(
        [&](auto& g) { create_complete_graph(g, n); },
        all_graph_views()
    )(gi.get_graph_view());
}

// Define the Python module
BOOST_PYTHON_MODULE(complete_graph)
{
    using namespace boost::python;

    def("create_complete_graph", complete_graph_bind,
        "Create a complete (fully connected) graph.\n\n"
        "Every node is connected to every other node.\n"
        "Total edges: n*(n-1)/2\n\n"
        "Args:\n"
        "    gi: GraphInterface (internal graph object)\n"
        "    n: Number of nodes\n");
}

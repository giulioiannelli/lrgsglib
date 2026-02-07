#include "barabasi_albert.hh"
#include <boost/python.hpp>

// Binding function for standard Barabasi-Albert
void barabasi_albert_bind(GraphInterface& gi, int n, int m, unsigned long seed)
{
    gt_dispatch<>()(
        [&](auto& g) { create_barabasi_albert(g, n, m, seed); },
        all_graph_views()
    )(gi.get_graph_view());
}

// Define the Python module
BOOST_PYTHON_MODULE(barabasi_albert)
{
    using namespace boost::python;

    def("create_barabasi_albert", barabasi_albert_bind,
        "Create a Barabasi-Albert preferential attachment graph.\n\n"
        "Args:\n"
        "    gi: GraphInterface (internal graph object)\n"
        "    n: Total number of nodes\n"
        "    m: Number of edges per new node\n"
        "    seed: Random seed (0 = use random device)\n");
}

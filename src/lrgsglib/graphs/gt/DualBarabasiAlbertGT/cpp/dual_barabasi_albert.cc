#include "dual_barabasi_albert.hh"
#include <boost/python.hpp>

// Binding function for Dual Barabasi-Albert
void dual_barabasi_albert_bind(GraphInterface& gi, int n, int m1, int m2, double p, unsigned long seed)
{
    gt_dispatch<>()(
        [&](auto& g) { create_dual_barabasi_albert(g, n, m1, m2, p, seed); },
        all_graph_views()
    )(gi.get_graph_view());
}

// Define the Python module
BOOST_PYTHON_MODULE(dual_barabasi_albert)
{
    using namespace boost::python;

    def("create_dual_barabasi_albert", dual_barabasi_albert_bind,
        "Create a Dual Barabasi-Albert graph with two attachment modes.\n\n"
        "With probability p, new nodes use m1 edges; otherwise m2.\n\n"
        "Args:\n"
        "    gi: GraphInterface (internal graph object)\n"
        "    n: Total number of nodes\n"
        "    m1: Number of edges in mode 1\n"
        "    m2: Number of edges in mode 2\n"
        "    p: Probability of using m1 (vs m2)\n"
        "    seed: Random seed (0 = use random device)\n");
}

#include "extended_barabasi_albert.hh"
#include <boost/python.hpp>

// Binding function for Extended Barabasi-Albert (with attractiveness)
void extended_barabasi_albert_bind(GraphInterface& gi, int n, int m, double a, unsigned long seed)
{
    gt_dispatch<>()(
        [&](auto& g) { create_extended_barabasi_albert(g, n, m, a, seed); },
        all_graph_views()
    )(gi.get_graph_view());
}

// Define the Python module
BOOST_PYTHON_MODULE(extended_barabasi_albert)
{
    using namespace boost::python;

    def("create_extended_barabasi_albert", extended_barabasi_albert_bind,
        "Create an Extended Barabasi-Albert graph with initial attractiveness.\n\n"
        "Attachment probability: P(i) ~ (k_i + a)\n"
        "When a=0, reduces to standard BA. Higher a = more uniform degrees.\n\n"
        "Args:\n"
        "    gi: GraphInterface (internal graph object)\n"
        "    n: Total number of nodes\n"
        "    m: Number of edges per new node\n"
        "    a: Initial attractiveness parameter (>= 0)\n"
        "    seed: Random seed (0 = use random device)\n");
}

#include "holme_kim.hh"
#include <boost/python.hpp>

// Binding function for Holme-Kim graph
void holme_kim_bind(GraphInterface& gi, int n, int m, double p, unsigned long seed)
{
    gt_dispatch<>()(
        [&](auto& g) { create_holme_kim(g, n, m, p, seed); },
        all_graph_views()
    )(gi.get_graph_view());
}

// Define the Python module
BOOST_PYTHON_MODULE(holme_kim)
{
    using namespace boost::python;

    def("create_holme_kim", holme_kim_bind,
        "Create a Holme-Kim powerlaw cluster graph.\n\n"
        "Args:\n"
        "    gi: GraphInterface (internal graph object)\n"
        "    n: Total number of nodes\n"
        "    m: Number of edges per new node\n"
        "    p: Probability of triad formation (0=pure BA, 1=max clustering)\n"
        "    seed: Random seed (0 = use random device)\n");
}

#include "triangular_lattice.hh"
#include <boost/python.hpp>

// Function to be called from Python
void triangular_lattice_bind(GraphInterface& gi, int width, int height)
{
    // We need to invoke the appropriate instance of the template function
    // at runtime using the gt_dispatch<>() function, which handles the type
    // resolution. Here, we're considering all possible graph views and
    // the most basic vertex property maps (which we don't actually use here).
    gt_dispatch<>()(
        [&](auto& g) { create_triangular_lattice(g, width, height); },
        all_graph_views()
    )(gi.get_graph_view());
}

// Define the Python module that exposes the triangular lattice creation function
BOOST_PYTHON_MODULE(triangular_lattice)
{
    using namespace boost::python;
    def("create_triangular_lattice", triangular_lattice_bind);
}

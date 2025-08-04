import graph_tool.all as gt

# Import the C++ module for creating a triangular lattice (called triangular_lattice.so)
import triangular_lattice

def create_triangular_lattice(width, height):
    """
    Create a triangular lattice graph with given width and height using the
    graph_tool library and a C++ extension for efficient implementation.

    Args:
        width (int): The width of the lattice.
        height (int): The height of the lattice.

    Returns:
        Graph: A graph_tool Graph object representing the triangular lattice.
    """
    # Create a new empty Graph using graph_tool
    g = gt.Graph(directed=False)

    # For graph objects, we need to pass the internal GraphInterface which is
    # accessed via g._Graph__graph.

    # The C++ function needs no property maps, so we only pass the graph and dimensions
    triangular_lattice.create_triangular_lattice(g._Graph__graph, width, height)

    return g

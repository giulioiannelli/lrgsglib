"""
Lattice2DGT: graph-tool implementation of 2D lattice signed graphs.

Mirrors the API of nx_patches.Lattice2D for easy switching between backends.
"""
from __future__ import annotations

from typing import Literal, Optional, Tuple, Union
import numpy as np

import graph_tool.all as gt
import graph_tool.generation as gen

from ..SignedGraphGT import SignedGraphGT


# Try to import C++ triangular lattice extension
try:
    from .cpp.lattice_generators import create_triangular_lattice as _cpp_triangular
    _CPP_TRIANGULAR_AVAILABLE = True
except ImportError:
    _CPP_TRIANGULAR_AVAILABLE = False


class Lattice2DGT(SignedGraphGT):
    """
    2D lattice graph using graph-tool backend.

    Mirrors nx_patches.Lattice2D API for easy backend switching.

    Parameters
    ----------
    side1 : int
        Size in first dimension.
    side2 : int, optional
        Size in second dimension. If None, uses side1 (square lattice).
    geo : str
        Lattice geometry: 'sqr' (square), 'tri' (triangular), 'hex' (hexagonal).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    periodic : bool
        Whether to use periodic boundary conditions.
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> lat = Lattice2DGT(side1=50, geo='sqr', pflip=0.3)
    >>> lat.flip_random_fract_edges()
    >>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")

    >>> # Triangular lattice (fast C++ backend)
    >>> lat_tri = Lattice2DGT(side1=100, geo='tri')
    """

    # Supported geometries
    GEOMETRIES = {'sqr', 'tri', 'hex'}

    def __init__(
        self,
        side1: int,
        side2: Optional[int] = None,
        geo: Literal['sqr', 'tri', 'hex'] = 'sqr',
        pflip: float = 0.0,
        periodic: bool = False,
        seed: Optional[int] = None,
    ):
        # Validate inputs
        if geo not in self.GEOMETRIES:
            raise ValueError(f"geo must be one of {self.GEOMETRIES}, got '{geo}'")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        self.side1 = side1
        self.side2 = side2 if side2 is not None else side1
        self.geo = geo
        self.periodic = periodic

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Create the graph based on geometry
        if geo == 'sqr':
            G, pos = self._create_square_lattice()
        elif geo == 'tri':
            G, pos = self._create_triangular_lattice()
        elif geo == 'hex':
            G, pos = self._create_hexagonal_lattice()

        # Store position property
        self._pos = pos

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

    def _create_square_lattice(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Create square lattice using GT's native function."""
        G = gen.lattice([self.side1, self.side2], periodic=self.periodic)

        # Create position property map
        pos = G.new_vertex_property("vector<double>")
        for v in G.vertices():
            idx = int(v)
            x = idx % self.side1
            y = idx // self.side1
            pos[v] = [float(x), float(y)]

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, pos

    def _create_triangular_lattice(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Create triangular lattice using C++ extension or fallback."""
        if _CPP_TRIANGULAR_AVAILABLE:
            # Fast C++ path
            G = _cpp_triangular(self.side1, self.side2)
        else:
            # Fallback: create points and triangulate
            G = self._triangular_lattice_fallback()

        # Create position property map (triangular grid)
        pos = G.new_vertex_property("vector<double>")
        for v in G.vertices():
            idx = int(v)
            x = idx % self.side1
            y = idx // self.side1
            # Offset odd rows for triangular geometry
            x_offset = 0.5 if y % 2 == 1 else 0.0
            pos[v] = [float(x) + x_offset, float(y) * np.sqrt(3) / 2]

        # Add sign property
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, pos

    def _triangular_lattice_fallback(self) -> gt.Graph:
        """Create triangular lattice without C++ extension."""
        # Create regular grid points
        points = []
        for y in range(self.side2):
            for x in range(self.side1):
                x_offset = 0.5 if y % 2 == 1 else 0.0
                points.append([x + x_offset, y * np.sqrt(3) / 2])

        points = np.array(points)
        G, _ = gen.triangulation(points, type='delaunay')
        return G

    def _create_hexagonal_lattice(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Create hexagonal (honeycomb) lattice."""
        # Hexagonal lattice: create triangular and remove every third vertex
        # Or build directly from hex grid coordinates

        G = gt.Graph(directed=False)
        pos = G.new_vertex_property("vector<double>")

        # Create vertices in honeycomb pattern
        # Each unit cell has 2 vertices (A and B sublattice)
        vertex_map = {}  # (i, j, sublattice) -> vertex

        for j in range(self.side2):
            for i in range(self.side1):
                # A sublattice
                v_a = G.add_vertex()
                x_a = i * 1.5
                y_a = j * np.sqrt(3) + (0 if i % 2 == 0 else np.sqrt(3) / 2)
                pos[v_a] = [x_a, y_a]
                vertex_map[(i, j, 'A')] = v_a

                # B sublattice (offset)
                v_b = G.add_vertex()
                x_b = x_a + 0.5
                y_b = y_a + np.sqrt(3) / 2 if i % 2 == 0 else y_a - np.sqrt(3) / 2
                pos[v_b] = [x_b, y_b]
                vertex_map[(i, j, 'B')] = v_b

                # Connect A-B within cell
                G.add_edge(v_a, v_b)

        # Connect between cells
        for j in range(self.side2):
            for i in range(self.side1):
                v_a = vertex_map[(i, j, 'A')]
                v_b = vertex_map[(i, j, 'B')]

                # Horizontal connections
                if i < self.side1 - 1:
                    v_a_next = vertex_map[(i + 1, j, 'A')]
                    G.add_edge(v_b, v_a_next)
                elif self.periodic:
                    v_a_next = vertex_map[(0, j, 'A')]
                    G.add_edge(v_b, v_a_next)

                # Vertical connections
                if j < self.side2 - 1:
                    if i % 2 == 0:
                        v_b_up = vertex_map[(i, j + 1, 'B')]
                        G.add_edge(v_a, v_b_up)
                    else:
                        v_a_up = vertex_map[(i, j + 1, 'A')]
                        G.add_edge(v_b, v_a_up)

        # Add sign property
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, pos

    @property
    def pos(self) -> gt.VertexPropertyMap:
        """Vertex position property map for visualization."""
        return self._pos

    def draw(self, ax=None, **kwargs):
        """
        Draw the lattice with edge colors based on sign.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on. If None, creates new figure.
        **kwargs : dict
            Additional arguments passed to gt.graph_draw.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))

        # Color edges by sign
        edge_color = self.G.new_edge_property("vector<double>")
        sign = self.G.edge_properties.get("sign")

        for e in self.G.edges():
            if sign is not None and sign[e] == -1:
                edge_color[e] = [0.8, 0.2, 0.2, 0.8]  # red for negative
            else:
                edge_color[e] = [0.2, 0.6, 0.2, 0.8]  # green for positive

        # Default drawing options
        draw_opts = {
            'pos': self._pos,
            'edge_color': edge_color,
            'vertex_size': 4,
            'edge_pen_width': 1.0,
            'mplfig': ax,
        }
        draw_opts.update(kwargs)

        gt.graph_draw(self.G, **draw_opts)
        ax.set_aspect('equal')

        return ax

__all__ = ["Lattice2DGT"]

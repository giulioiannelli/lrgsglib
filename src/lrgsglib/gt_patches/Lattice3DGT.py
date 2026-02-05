"""
Lattice3DGT: graph-tool implementation of 3D lattice signed graphs.

Mirrors the API of nx_patches.Lattice3D for easy switching between backends.
"""
from __future__ import annotations

from typing import Literal, Optional, Tuple, Union

import numpy as np

import graph_tool.all as gt
import graph_tool.generation as gen

from .signed_graph_gt import SignedGraphGT


class Lattice3DGT(SignedGraphGT):
    """
    3D lattice graph using graph-tool backend.

    Mirrors nx_patches.Lattice3D API for easy backend switching.

    Parameters
    ----------
    dim : int or tuple of 3 ints
        Lattice dimensions. If int, used for all axes (cubic lattice).
    geo : str
        Lattice geometry: 'sc' (simple cubic), 'bcc' (body-centered cubic),
        'fcc' (face-centered cubic).
    pflip : float
        Fraction of edges to flip to negative (0.0 to 1.0).
    periodic : bool
        Whether to use periodic boundary conditions.
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> lat = Lattice3DGT(dim=10, geo='sc', pflip=0.3)
    >>> lat.flip_random_fract_edges()
    >>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")

    >>> # BCC lattice (body-centered cubic)
    >>> lat_bcc = Lattice3DGT(dim=5, geo='bcc', periodic=True)
    """

    # Supported geometries
    GEOMETRIES = {"sc", "bcc", "fcc"}

    # Node multiplier per unit cell
    NODE_MULTIPLIERS = {"sc": 1, "bcc": 2, "fcc": 4}

    def __init__(
        self,
        dim: Union[int, Tuple[int, int, int]] = 10,
        geo: Literal["sc", "bcc", "fcc"] = "sc",
        pflip: float = 0.0,
        periodic: bool = False,
        seed: Optional[int] = None,
    ):
        # Validate inputs
        if geo not in self.GEOMETRIES:
            raise ValueError(f"geo must be one of {self.GEOMETRIES}, got '{geo}'")
        if not 0.0 <= pflip <= 1.0:
            raise ValueError(f"pflip must be in [0, 1], got {pflip}")

        # Parse dimensions
        if isinstance(dim, int):
            self.dim = (dim, dim, dim)
        elif isinstance(dim, tuple) and len(dim) == 3:
            self.dim = tuple(sorted(dim, reverse=True))
        else:
            raise ValueError("dim must be int or tuple of 3 ints")

        self.geo = geo
        self.periodic = periodic
        self.node_multiplier = self.NODE_MULTIPLIERS[geo]

        # Set seed
        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Create the graph based on geometry
        if geo == "sc":
            G, pos = self._create_simple_cubic()
        elif geo == "bcc":
            G, pos = self._create_bcc()
        elif geo == "fcc":
            G, pos = self._create_fcc()

        # Store position property
        self._pos = pos

        # Initialize parent class
        super().__init__(G=G, pflip=pflip, seed=seed)

    def _create_simple_cubic(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Create simple cubic lattice using GT's native function."""
        G = gen.lattice(list(self.dim), periodic=self.periodic)

        # Create position property map (3D coordinates)
        pos = G.new_vertex_property("vector<double>")
        nx, ny, nz = self.dim
        for v in G.vertices():
            idx = int(v)
            x = idx % nx
            y = (idx // nx) % ny
            z = idx // (nx * ny)
            pos[v] = [float(x), float(y), float(z)]

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, pos

    def _create_bcc(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Create body-centered cubic lattice."""
        G = gt.Graph(directed=False)
        nx, ny, nz = self.dim

        # BCC has two atoms per unit cell:
        # - Corner atoms at (x, y, z)
        # - Center atoms at (x+0.5, y+0.5, z+0.5)
        offsets = [(0, 0, 0), (0.5, 0.5, 0.5)]

        # Create position property
        pos = G.new_vertex_property("vector<double>")

        # Create all nodes
        node_map = {}  # (x, y, z) -> vertex
        vertex_idx = 0

        range_max = [n if self.periodic else n - 1 for n in self.dim]

        for oz, oy, ox in offsets:
            for z in range(nz if self.periodic else nz - (1 if oz > 0 else 0)):
                for y in range(ny if self.periodic else ny - (1 if oy > 0 else 0)):
                    for x in range(nx if self.periodic else nx - (1 if ox > 0 else 0)):
                        coord = (x + ox, y + oy, z + oz)
                        v = G.add_vertex()
                        pos[v] = [float(coord[0]), float(coord[1]), float(coord[2])]
                        node_map[coord] = v

        # BCC connectivity: center connected to 8 corners
        corner_to_center = [
            (0.5, 0.5, 0.5),
            (0.5, 0.5, -0.5),
            (0.5, -0.5, 0.5),
            (-0.5, 0.5, 0.5),
            (0.5, -0.5, -0.5),
            (-0.5, -0.5, -0.5),
            (-0.5, 0.5, -0.5),
            (-0.5, -0.5, 0.5),
        ]

        def wrap_coord(coord, dims):
            """Wrap coordinate for periodic boundary conditions."""
            return (coord[0] % dims[0], coord[1] % dims[1], coord[2] % dims[2])

        # Add edges
        for coord, v in node_map.items():
            # Connect to neighbors based on offset type
            if all(c == int(c) for c in coord):
                # Corner atom - connect to nearest center atoms
                for dx, dy, dz in corner_to_center:
                    neighbor = (coord[0] + dx, coord[1] + dy, coord[2] + dz)
                    if self.periodic:
                        neighbor = wrap_coord(neighbor, self.dim)
                    if neighbor in node_map:
                        G.add_edge(v, node_map[neighbor])

        # Add sign property
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, pos

    def _create_fcc(self) -> Tuple[gt.Graph, gt.VertexPropertyMap]:
        """Create face-centered cubic lattice."""
        G = gt.Graph(directed=False)
        nx, ny, nz = self.dim

        # FCC has four atoms per unit cell:
        # - Corner at (0, 0, 0)
        # - Face centers at (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)
        offsets = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]

        # Create position property
        pos = G.new_vertex_property("vector<double>")

        # Create all nodes
        node_map = {}

        for oz, oy, ox in offsets:
            for z in range(nz if self.periodic else nz - (1 if oz > 0 else 0)):
                for y in range(ny if self.periodic else ny - (1 if oy > 0 else 0)):
                    for x in range(nx if self.periodic else nx - (1 if ox > 0 else 0)):
                        coord = (x + ox, y + oy, z + oz)
                        v = G.add_vertex()
                        pos[v] = [float(coord[0]), float(coord[1]), float(coord[2])]
                        node_map[coord] = v

        def wrap_coord(coord, dims):
            """Wrap coordinate for periodic boundary conditions."""
            return (coord[0] % dims[0], coord[1] % dims[1], coord[2] % dims[2])

        # FCC connectivity: each atom connects to 12 nearest neighbors
        # (atoms at distance sqrt(2)/2 in the same lattice)
        # Neighbors differ by one of: (+/- 0.5, +/- 0.5, 0) and permutations
        neighbor_deltas = [
            (0.5, 0.5, 0),
            (0.5, -0.5, 0),
            (-0.5, 0.5, 0),
            (-0.5, -0.5, 0),
            (0.5, 0, 0.5),
            (0.5, 0, -0.5),
            (-0.5, 0, 0.5),
            (-0.5, 0, -0.5),
            (0, 0.5, 0.5),
            (0, 0.5, -0.5),
            (0, -0.5, 0.5),
            (0, -0.5, -0.5),
        ]

        # Add edges
        added_edges = set()
        for coord, v in node_map.items():
            for dx, dy, dz in neighbor_deltas:
                neighbor = (coord[0] + dx, coord[1] + dy, coord[2] + dz)
                if self.periodic:
                    neighbor = wrap_coord(neighbor, self.dim)
                if neighbor in node_map:
                    # Avoid duplicate edges
                    edge_key = tuple(sorted([coord, neighbor]))
                    if edge_key not in added_edges:
                        G.add_edge(v, node_map[neighbor])
                        added_edges.add(edge_key)

        # Add sign property
        sign = G.new_edge_property("int")
        for e in G.edges():
            sign[e] = 1
        G.edge_properties["sign"] = sign

        return G, pos

    @property
    def side1(self) -> int:
        """First dimension size (for compatibility)."""
        return self.dim[0]

    @property
    def side2(self) -> int:
        """Second dimension size (for compatibility)."""
        return self.dim[1]

    @property
    def side3(self) -> int:
        """Third dimension size."""
        return self.dim[2]

    def get_positions(self) -> np.ndarray:
        """
        Get 3D positions of all nodes.

        Returns
        -------
        ndarray
            N x 3 array of node positions.
        """
        positions = np.zeros((self.N, 3))
        for v in self.G.vertices():
            positions[int(v)] = self._pos[v]
        return positions

    def draw(
        self,
        ax=None,
        projection: Literal["xy", "xz", "yz", "3d"] = "xy",
        **kwargs,
    ):
        """
        Draw the lattice.

        Parameters
        ----------
        ax : matplotlib axes, optional
            Axes to draw on.
        projection : str
            2D projection plane or '3d' for 3D plot.
        **kwargs
            Passed to graph-tool's graph_draw.
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))

        pos = self.get_positions()

        # Project to 2D
        if projection == "xy":
            pos_2d = pos[:, :2]
        elif projection == "xz":
            pos_2d = pos[:, [0, 2]]
        elif projection == "yz":
            pos_2d = pos[:, 1:]
        else:
            raise ValueError(f"projection must be 'xy', 'xz', 'yz', got {projection}")

        # Create 2D position property for drawing
        pos_prop = self.G.new_vertex_property("vector<double>")
        for v in self.G.vertices():
            pos_prop[v] = pos_2d[int(v)].tolist()

        # Color edges by sign
        sign_prop = self.G.edge_properties["sign"]
        edge_color = self.G.new_edge_property("vector<double>")
        for e in self.G.edges():
            if sign_prop[e] == 1:
                edge_color[e] = [0.2, 0.2, 0.8, 0.7]  # Blue for positive
            else:
                edge_color[e] = [0.8, 0.2, 0.2, 0.7]  # Red for negative

        gt.graph_draw(
            self.G,
            pos=pos_prop,
            edge_color=edge_color,
            mplfig=ax.figure,
            **kwargs,
        )

        return ax

    def __repr__(self) -> str:
        neg = self.count_negative_edges()
        return (
            f"Lattice3DGT(dim={self.dim}, geo='{self.geo}', "
            f"N={self.N}, edges={self.num_edges}, negative={neg})"
        )

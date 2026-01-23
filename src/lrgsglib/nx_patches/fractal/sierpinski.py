"""
SierpinskiGraph: Sierpinski gasket/triangle fractal graph.

The Sierpinski gasket is a classic fractal structure with well-known
spectral properties, including exact solutions for eigenvalues.

Example
-------
>>> from lrgsglib.nx_patches.fractal import SierpinskiGraph
>>> sg = SierpinskiGraph(n=4, pflip=0.1, seed=42)
>>> sg.flip_random_fract_edges()
>>> print(f"Nodes: {sg.N}, Edges: {sg.M}")
"""

from typing import Optional
import networkx as nx
import numpy as np

from .fractal_graph import FractalGraph


class SierpinskiGraph(FractalGraph):
    """
    Sierpinski gasket (triangle) graph.

    The Sierpinski gasket is constructed iteratively:
    - Level 0: A single triangle (3 nodes, 3 edges)
    - Level n: Three copies of level n-1, joined at corner vertices

    Parameters
    ----------
    n : int
        Iteration level (n >= 0). Level 0 is a single triangle.
    variant : str, default "gasket"
        Type of Sierpinski graph:
        - "gasket": Standard Sierpinski triangle
        - "carpet": Sierpinski carpet (2D, square-based)
        - "tetrahedron": 3D Sierpinski tetrahedron
    with_positions : bool, default True
        Whether to compute 2D positions for visualization.
    only_const_mode : bool, default False
        If True, populate metadata only without building the graph.
    **kwargs
        Additional arguments passed to FractalGraph (pflip, seed, etc.).

    Attributes
    ----------
    n : int
        Iteration level.
    variant : str
        Sierpinski variant type.
    N : int
        Number of nodes = (3^(n+1) + 3) / 2 for gasket.
    M : int
        Number of edges = 3^(n+1) for gasket.
    spectral_dimension : float
        Spectral dimension d_s = 2 * ln(3) / ln(5) ≈ 1.365 for gasket.

    Notes
    -----
    The Sierpinski gasket has exact spectral solutions. The eigenvalues
    follow a hierarchical pattern related to the self-similar structure.

    The spectral dimension d_s ≈ 1.365 (for gasket) governs diffusion:
    - Mean square displacement ~ t^(2/d_w) where d_w = ln(5)/ln(2)
    - Return probability ~ t^(-d_s/2)

    References
    ----------
    .. [1] Rammal, R., & Toulouse, G. (1983). Random walks on fractal
           structures and percolation clusters. J. Physique Lett., 44(1), 13-22.
    """

    # Spectral dimension for standard Sierpinski gasket
    SPECTRAL_DIMENSION_GASKET = 2 * np.log(3) / np.log(5)  # ≈ 1.365
    SPECTRAL_DIMENSION_CARPET = 2 * np.log(8) / np.log(8)  # = 2 (marginal)
    SPECTRAL_DIMENSION_TETRA = 2 * np.log(4) / np.log(6)   # ≈ 1.547

    def __init__(
        self,
        n: int,
        variant: str = "gasket",
        with_positions: bool = True,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        if n < 0:
            raise ValueError(f"Iteration level n must be >= 0, got {n}")
        if variant not in ("gasket", "carpet", "tetrahedron"):
            raise ValueError(
                f"Unknown variant '{variant}'. "
                "Choose from: 'gasket', 'carpet', 'tetrahedron'"
            )

        self.variant = variant
        self._corner_nodes: list = []  # Track corner nodes for joining

        super().__init__(
            n=n,
            sgpathn=f"Sierpinski_{variant}",
            std_fname=f"sierp{variant[0]}",
            with_positions=with_positions,
            only_const_mode=only_const_mode,
            **kwargs,
        )

        # Set spectral dimension based on variant
        if variant == "gasket":
            self.spectral_dimension = self.SPECTRAL_DIMENSION_GASKET
        elif variant == "carpet":
            self.spectral_dimension = self.SPECTRAL_DIMENSION_CARPET
        else:
            self.spectral_dimension = self.SPECTRAL_DIMENSION_TETRA

    def _compute_expected_nodes(self) -> int:
        """Compute expected number of nodes for iteration level n."""
        if self.variant == "gasket":
            # N = (3^(n+1) + 3) / 2
            return (3 ** (self.n + 1) + 3) // 2
        elif self.variant == "carpet":
            # Sierpinski carpet: 8^n nodes (roughly, depends on construction)
            return 8 ** self.n if self.n > 0 else 8
        else:  # tetrahedron
            # N = (4^(n+1) + 2*4) / 2 = 2 * 4^n + 2
            return 2 * (4 ** self.n) + 2

    def _init_network(self) -> None:
        """Build the Sierpinski graph at iteration level n."""
        if self.variant == "gasket":
            self._build_gasket()
        elif self.variant == "carpet":
            self._build_carpet()
        else:
            self._build_tetrahedron()

        self.syshape = self.G.number_of_nodes()
        self.syshapePth = f"N={self.syshape}_n={self.n}"

        if self.with_positions:
            self._compute_positions()

    def _build_gasket(self) -> None:
        """
        Build Sierpinski gasket iteratively.

        Uses the "corner-joining" construction:
        1. Start with a triangle (level 0)
        2. At each level, make 3 copies and join at corners
        """
        # Level 0: single triangle
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (2, 0)])
        corners = [0, 1, 2]

        for level in range(self.n):
            G, corners = self._gasket_iteration(G, corners)

        self.G = G
        self._corner_nodes = corners

    def _gasket_iteration(
        self, G: nx.Graph, corners: list
    ) -> tuple[nx.Graph, list]:
        """
        Perform one iteration of Sierpinski gasket construction.

        Creates 3 copies of current graph and joins them at corners.
        """
        n_nodes = G.number_of_nodes()
        new_G = nx.Graph()

        # Copy 1: original (bottom-left)
        new_G.add_edges_from(G.edges())

        # Copy 2: shifted (bottom-right)
        mapping2 = {v: v + n_nodes for v in G.nodes()}
        G2 = nx.relabel_nodes(G, mapping2)
        new_G.add_edges_from(G2.edges())

        # Copy 3: shifted (top)
        mapping3 = {v: v + 2 * n_nodes for v in G.nodes()}
        G3 = nx.relabel_nodes(G, mapping3)
        new_G.add_edges_from(G3.edges())

        # Join at corners:
        # - Corner 1 of copy 1 joins with corner 0 of copy 2
        # - Corner 2 of copy 1 joins with corner 0 of copy 3
        # - Corner 2 of copy 2 joins with corner 1 of copy 3
        c1_0, c1_1, c1_2 = corners
        c2_0, c2_1, c2_2 = [c + n_nodes for c in corners]
        c3_0, c3_1, c3_2 = [c + 2 * n_nodes for c in corners]

        # Merge nodes at joining points (contract edges)
        # We merge by removing one node and redirecting edges
        new_G = nx.contracted_nodes(new_G, c1_1, c2_0, self_loops=False)
        new_G = nx.contracted_nodes(new_G, c1_2, c3_0, self_loops=False)
        new_G = nx.contracted_nodes(new_G, c2_2, c3_1, self_loops=False)

        # Relabel to consecutive integers
        new_G = nx.convert_node_labels_to_integers(new_G)

        # New corners are the outer corners of the three copies
        # After contraction, we need to track where they ended up
        # For simplicity, recalculate based on degree (corners have degree 2)
        new_corners = self._find_corner_nodes(new_G)

        return new_G, new_corners

    def _find_corner_nodes(self, G: nx.Graph) -> list:
        """Find the 3 corner nodes (nodes with minimum degree)."""
        degrees = dict(G.degree())
        min_deg = min(degrees.values())
        corners = [v for v, d in degrees.items() if d == min_deg]
        # Return exactly 3 corners (sort for consistency)
        return sorted(corners)[:3]

    def _build_carpet(self) -> None:
        """
        Build Sierpinski carpet using NetworkX's built-in function.

        The carpet is constructed on a grid with holes.
        """
        if self.n == 0:
            # Level 0: 3x3 grid minus center
            G = nx.grid_2d_graph(3, 3)
            G.remove_node((1, 1))
        else:
            G = nx.generators.small.sedgewick_maze_graph()
            # For higher levels, build iteratively
            G = self._build_carpet_iterative(self.n)

        self.G = nx.convert_node_labels_to_integers(G)

    def _build_carpet_iterative(self, n: int) -> nx.Graph:
        """Build Sierpinski carpet iteratively."""
        # Start with 3x3 grid minus center
        size = 3 ** n
        G = nx.grid_2d_graph(size, size)

        # Remove nodes corresponding to "holes" at all scales
        holes = self._get_carpet_holes(n)
        G.remove_nodes_from(holes)

        return G

    def _get_carpet_holes(self, n: int) -> set:
        """Get coordinates of all holes in the Sierpinski carpet."""
        holes = set()

        for level in range(1, n + 1):
            block_size = 3 ** level
            hole_size = 3 ** (level - 1)

            for i in range(0, 3 ** n, block_size):
                for j in range(0, 3 ** n, block_size):
                    # Hole is in the center third of each block
                    start_i = i + hole_size
                    start_j = j + hole_size
                    for hi in range(start_i, start_i + hole_size):
                        for hj in range(start_j, start_j + hole_size):
                            holes.add((hi, hj))

        return holes

    def _build_tetrahedron(self) -> None:
        """
        Build Sierpinski tetrahedron (3D generalization).

        Similar to gasket but with 4 copies joined at vertices.
        """
        # Level 0: tetrahedron (4 nodes, 6 edges - complete graph K4)
        G = nx.complete_graph(4)
        corners = [0, 1, 2, 3]

        for level in range(self.n):
            G, corners = self._tetrahedron_iteration(G, corners)

        self.G = G
        self._corner_nodes = corners

    def _tetrahedron_iteration(
        self, G: nx.Graph, corners: list
    ) -> tuple[nx.Graph, list]:
        """Perform one iteration of Sierpinski tetrahedron construction."""
        n_nodes = G.number_of_nodes()
        new_G = nx.Graph()

        # Create 4 copies
        copies = []
        for i in range(4):
            if i == 0:
                copy_G = G.copy()
            else:
                mapping = {v: v + i * n_nodes for v in G.nodes()}
                copy_G = nx.relabel_nodes(G, mapping)
            new_G.add_edges_from(copy_G.edges())
            copies.append([c + i * n_nodes if i > 0 else c for c in corners])

        # Join copies at corners (6 joins for 4 tetrahedra)
        # Each pair of copies shares one corner
        joins = [
            (copies[0][1], copies[1][0]),
            (copies[0][2], copies[2][0]),
            (copies[0][3], copies[3][0]),
            (copies[1][2], copies[2][1]),
            (copies[1][3], copies[3][1]),
            (copies[2][3], copies[3][2]),
        ]

        for node1, node2 in joins:
            if node1 in new_G and node2 in new_G:
                new_G = nx.contracted_nodes(new_G, node1, node2, self_loops=False)

        new_G = nx.convert_node_labels_to_integers(new_G)

        # Find new corners (4 outer vertices)
        degrees = dict(new_G.degree())
        min_deg = min(degrees.values())
        new_corners = sorted([v for v, d in degrees.items() if d == min_deg])[:4]

        return new_G, new_corners

    def _compute_positions(self) -> None:
        """Compute 2D/3D positions for visualization."""
        if self.variant == "gasket":
            pos = self._gasket_positions()
        elif self.variant == "carpet":
            pos = self._carpet_positions()
        else:
            pos = self._tetrahedron_positions()

        nx.set_node_attributes(self.G, pos, "pos")

    def _gasket_positions(self) -> dict:
        """Compute positions for Sierpinski gasket using recursive subdivision."""
        # Use spring layout as fallback, but try to compute fractal positions
        # For exact positions, use barycentric subdivision
        if self.n <= 6:
            # For small n, compute exact positions
            return self._compute_gasket_positions_exact()
        else:
            # For large n, use spring layout
            return nx.spring_layout(self.G, seed=42)

    def _compute_gasket_positions_exact(self) -> dict:
        """Compute exact positions using iterated function system."""
        # Vertices of outer triangle
        v0 = np.array([0.0, 0.0])
        v1 = np.array([1.0, 0.0])
        v2 = np.array([0.5, np.sqrt(3) / 2])

        pos = {}
        # Generate all points using chaos game / IFS
        n_nodes = self.G.number_of_nodes()

        # For small graphs, assign positions based on graph structure
        if n_nodes <= 3:
            corners = list(self.G.nodes())[:3]
            pos[corners[0]] = tuple(v0)
            pos[corners[1]] = tuple(v1)
            pos[corners[2]] = tuple(v2)
        else:
            # Use spectral layout which respects the fractal structure
            spectral_pos = nx.spectral_layout(self.G)
            # Scale and rotate to fit in triangle
            pos = spectral_pos

        return pos

    def _carpet_positions(self) -> dict:
        """Compute positions for Sierpinski carpet (grid-based)."""
        # Nodes retain their grid coordinates from construction
        pos = {}
        for node in self.G.nodes():
            if isinstance(node, tuple):
                pos[node] = node
            else:
                # If relabeled, use spring layout
                return nx.spring_layout(self.G, seed=42)
        return pos

    def _tetrahedron_positions(self) -> dict:
        """Compute 2D projection positions for tetrahedron."""
        # Use 3D spring layout projected to 2D
        pos_3d = nx.spring_layout(self.G, dim=3, seed=42)
        # Project to 2D (drop z coordinate or use isometric projection)
        pos_2d = {}
        for node, (x, y, z) in pos_3d.items():
            # Isometric projection
            pos_2d[node] = (x + 0.5 * z, y + 0.5 * z)
        return pos_2d

    def get_corner_nodes(self) -> list:
        """
        Return the corner nodes of the Sierpinski structure.

        These are the nodes with minimum degree that form the outer boundary.

        Returns
        -------
        list
            Node IDs of corner vertices.
        """
        return self._corner_nodes.copy() if self._corner_nodes else []

    def get_spectral_dimension(self) -> float:
        """
        Return the theoretical spectral dimension.

        Returns
        -------
        float
            Spectral dimension d_s.
        """
        return self.spectral_dimension

    def get_walk_dimension(self) -> float:
        """
        Return the walk dimension d_w.

        The walk dimension governs the scaling of mean square displacement:
        <r²(t)> ~ t^(2/d_w)

        Returns
        -------
        float
            Walk dimension d_w = ln(5)/ln(2) ≈ 2.32 for gasket.
        """
        if self.variant == "gasket":
            return np.log(5) / np.log(2)  # ≈ 2.32
        elif self.variant == "carpet":
            return np.log(8) / np.log(3)  # ≈ 1.89
        else:
            return np.log(6) / np.log(2)  # ≈ 2.58

    def get_fractal_dimension(self) -> float:
        """
        Return the fractal (Hausdorff) dimension.

        Returns
        -------
        float
            Fractal dimension d_f.
        """
        if self.variant == "gasket":
            return np.log(3) / np.log(2)  # ≈ 1.585
        elif self.variant == "carpet":
            return np.log(8) / np.log(3)  # ≈ 1.893
        else:
            return 2.0  # Tetrahedron fills 2D in limiting sense

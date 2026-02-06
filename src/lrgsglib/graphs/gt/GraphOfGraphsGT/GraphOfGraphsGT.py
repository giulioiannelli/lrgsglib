"""
GraphOfGraphsGT - Generalized hierarchical graph structure using graph-tool backend.

This module provides a flexible "graph of graphs" class that accepts arbitrary
SignedGraph types as base and fiber components, generalizing the Dirac lattice
pattern (DiracComb, DiracBrush).

Example
-------
>>> from lrgsglib.graphs.gt.GraphOfGraphsGT import GraphOfGraphsGT
>>> gog = GraphOfGraphsGT(
...     base_graph_type='Lattice2D',
...     base_params={'side1': 5, 'geo': 'sqr'},
...     fiber_graph_type='ErdosRenyi',
...     fiber_params={'n': 10, 'p': 0.2},
...     anchor_policy='first',
...     seed=42
... )
>>> print(f"N={gog.N}, N_base={gog.N_base}, N_fiber={gog.N_fiber}")
"""

from __future__ import annotations

from typing import Any, Callable, Dict, Literal, Optional, Union

import numpy as np

try:
    import graph_tool as gt
    from graph_tool import Graph

    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = object

from ....gt_patches.signed_graph_gt import SignedGraphGT
from ._generators import (
    extract_edges_from_gt_graph,
    generate_anchor_edges,
    generate_fiber_edges,
)
from ._policies import AnchorPolicy, get_anchor_index, resolve_anchor_indices


__all__ = ["GraphOfGraphsGT", "GraphOfGraphs"]


# Graph type registry for name -> class resolution
_GRAPH_TYPE_REGISTRY: Dict[str, Callable[..., Any]] = {}


def _register_graph_types():
    """Lazily populate the graph type registry."""
    if _GRAPH_TYPE_REGISTRY:
        return

    # Import graph types on demand to avoid circular imports
    from ..Lattice2DGT import Lattice2DGT
    from ..Lattice3DGT import Lattice3DGT
    from ..ErdosRenyiGT import ErdosRenyiGT
    from ..BarabasiAlbertGT import BarabasiAlbertGT
    from ..WattsStrogatzGT import WattsStrogatzGT
    from ..StochasticBlockModelGT import StochasticBlockModelGT

    _GRAPH_TYPE_REGISTRY.update({
        # Lattice types
        "Lattice2D": Lattice2DGT,
        "Lattice2DGT": Lattice2DGT,
        "Lattice3D": Lattice3DGT,
        "Lattice3DGT": Lattice3DGT,
        # Random types
        "ErdosRenyi": ErdosRenyiGT,
        "ErdosRenyiGT": ErdosRenyiGT,
        "BarabasiAlbert": BarabasiAlbertGT,
        "BarabasiAlbertGT": BarabasiAlbertGT,
        "WattsStrogatz": WattsStrogatzGT,
        "WattsStrogatzGT": WattsStrogatzGT,
        "StochasticBlockModel": StochasticBlockModelGT,
        "StochasticBlockModelGT": StochasticBlockModelGT,
    })


def _resolve_graph_type(type_name: str) -> Callable[..., Any]:
    """Resolve a graph type name to its class."""
    _register_graph_types()
    if type_name not in _GRAPH_TYPE_REGISTRY:
        available = ", ".join(sorted(_GRAPH_TYPE_REGISTRY.keys()))
        raise ValueError(
            f"Unknown graph type: {type_name}. Available: {available}"
        )
    return _GRAPH_TYPE_REGISTRY[type_name]


class GraphOfGraphsGT(SignedGraphGT):
    """Generalized graph of graphs structure using graph-tool backend.

    A GraphOfGraphs consists of a base graph where each node has a fiber graph
    attached to it. The fiber is connected to the base at an anchor node
    determined by the anchor policy.

    Parameters
    ----------
    base_graph_type : str
        Type name of the base graph (e.g., 'Lattice2D', 'ErdosRenyi')
    base_params : dict
        Parameters to pass to the base graph constructor
    fiber_graph_type : str
        Type name of the fiber graph
    fiber_params : dict or callable
        Parameters for fiber graph construction. If callable, should be
        a function(base_idx) -> dict for heterogeneous fibers.
    anchor_policy : str, AnchorPolicy, or callable, default 'first'
        How to select anchor nodes:
        - 'first': Use fiber node 0
        - 'center': Use middle node
        - 'last': Use last node
        - 'random': Random node (seeded)
        - callable: function(base_idx, n_fiber) -> anchor_idx
    pflip : float, default 0.0
        Fraction of edges to flip to negative
    seed : int, optional
        Random seed for reproducibility

    Attributes
    ----------
    gog_structure : dict
        Metadata about the graph-of-graphs structure
    N_base : int
        Number of nodes in the base graph
    N_fiber : int
        Number of nodes in each fiber graph

    Example
    -------
    >>> gog = GraphOfGraphsGT(
    ...     base_graph_type='Lattice2D',
    ...     base_params={'side1': 4, 'geo': 'sqr'},
    ...     fiber_graph_type='ErdosRenyi',
    ...     fiber_params={'n': 8, 'p': 0.3},
    ...     seed=42
    ... )
    >>> print(f"Total nodes: {gog.N}")  # 16 + 16*8 = 144
    """

    def __init__(
        self,
        base_graph_type: str,
        base_params: dict,
        fiber_graph_type: str,
        fiber_params: Union[dict, Callable[[int], dict]],
        anchor_policy: Union[str, AnchorPolicy, Callable[[int, int], int]] = "first",
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        # Store configuration
        self._base_graph_type = base_graph_type
        self._base_params = base_params
        self._fiber_graph_type = fiber_graph_type
        self._fiber_params = fiber_params
        self._anchor_policy = anchor_policy
        self._seed = seed

        # Initialize RNG before generation
        self._rng = np.random.default_rng(seed)

        # Determine if fibers are heterogeneous
        self._heterogeneous = callable(fiber_params)

        # Generate base and fiber graphs to get dimensions
        self._base_graph_instance = self._instantiate_base_graph()
        self._fiber_graph_instance = self._instantiate_fiber_graph(0)

        self._n_base = self._base_graph_instance.N
        self._n_fiber = self._fiber_graph_instance.N

        # Resolve anchor indices for all fibers
        self._anchor_indices = resolve_anchor_indices(
            anchor_policy, self._n_base, self._n_fiber, self._rng
        )

        # Build metadata structure
        self._gog_structure = self._build_metadata()

        # Generate composite graph
        G = self._generate()

        # Initialize base class
        super().__init__(G=G, pflip=pflip, seed=seed)

        # Apply sign flips if requested
        if pflip > 0:
            self.flip_random_fract_edges()

    def _instantiate_base_graph(self) -> SignedGraphGT:
        """Create the base graph instance."""
        cls = _resolve_graph_type(self._base_graph_type)
        # Use a separate seed for base graph
        base_seed = None
        if self._seed is not None:
            base_seed = int(self._rng.integers(0, 2**31))
        return cls(**self._base_params, seed=base_seed)

    def _instantiate_fiber_graph(self, base_idx: int) -> SignedGraphGT:
        """Create a fiber graph instance for a given base node."""
        cls = _resolve_graph_type(self._fiber_graph_type)

        if self._heterogeneous:
            params = self._fiber_params(base_idx)
        else:
            params = self._fiber_params

        # Use a separate seed for each fiber
        fiber_seed = None
        if self._seed is not None:
            fiber_seed = int(self._rng.integers(0, 2**31))

        return cls(**params, seed=fiber_seed)

    def _generate(self) -> "Graph":
        """Generate the composite graph-of-graphs structure.

        Vertex layout:
        [--- Base nodes ---][--- Fiber 0 ---][--- Fiber 1 ---] ... [--- Fiber N_base-1 ---]
           0 to N_base-1         N_fiber          N_fiber                  N_fiber

        Returns
        -------
        Graph
            The generated graph-tool Graph with sign property
        """
        total_nodes = self._n_base + self._n_base * self._n_fiber

        # Create graph
        G = gt.Graph(directed=False)
        G.add_vertex(total_nodes)

        # Collect all edges
        edge_list = []
        edge_signs = []

        # 1. Base edges
        base_edges = extract_edges_from_gt_graph(self._base_graph_instance.G)
        edge_list.extend(base_edges)
        # Get base edge signs
        if "sign" in self._base_graph_instance.G.edge_properties:
            sign_prop = self._base_graph_instance.G.edge_properties["sign"]
            for e in self._base_graph_instance.G.edges():
                edge_signs.append(sign_prop[e])
        else:
            edge_signs.extend([1] * len(base_edges))

        # 2. Fiber edges (one copy per base node)
        for base_idx in range(self._n_base):
            if self._heterogeneous:
                # Create new fiber graph for heterogeneous case
                fiber_graph = self._instantiate_fiber_graph(base_idx)
            else:
                # Reuse prototype for homogeneous case
                fiber_graph = self._fiber_graph_instance

            fiber_edges = generate_fiber_edges(
                fiber_graph.G, base_idx, self._n_base, self._n_fiber
            )
            edge_list.extend(fiber_edges)

            # Get fiber edge signs
            if "sign" in fiber_graph.G.edge_properties:
                sign_prop = fiber_graph.G.edge_properties["sign"]
                for e in fiber_graph.G.edges():
                    edge_signs.append(sign_prop[e])
            else:
                edge_signs.extend([1] * len(fiber_edges))

        # 3. Anchor edges (connect fiber anchors to base nodes)
        anchor_edges = generate_anchor_edges(
            self._n_base, self._n_fiber, self._anchor_indices
        )
        edge_list.extend(anchor_edges)
        edge_signs.extend([1] * len(anchor_edges))  # Anchor edges are positive

        # Batch add edges
        if edge_list:
            G.add_edge_list(edge_list)

        # Add sign property
        sign_prop = G.new_edge_property("int")
        for i, e in enumerate(G.edges()):
            sign_prop[e] = edge_signs[i]
        G.edge_properties["sign"] = sign_prop

        return G

    def _build_metadata(self) -> dict:
        """Build the gog_structure metadata dictionary."""
        return {
            # Graph types
            "base_graph_class": self._base_graph_type,
            "base_params": self._base_params,
            "fiber_graph_class": self._fiber_graph_type,
            "fiber_params": (
                "<callable>" if self._heterogeneous else self._fiber_params
            ),
            # Dimensions
            "N_base": self._n_base,
            "N_fiber": self._n_fiber,
            "N_total": self._n_base + self._n_base * self._n_fiber,
            # Anchor config
            "anchor_policy": (
                self._anchor_policy
                if isinstance(self._anchor_policy, str)
                else str(self._anchor_policy)
            ),
            "anchor_indices": self._anchor_indices,
            # Flags
            "heterogeneous": self._heterogeneous,
            "structure_type": "graph_of_graphs",
            "can_use_separated_spectrum": self.can_use_separated_spectrum(),
        }

    @property
    def gog_structure(self) -> dict:
        """Metadata about the graph-of-graphs structure."""
        return self._gog_structure

    @property
    def dirac_structure(self) -> dict:
        """Backward compatibility with Dirac lattice spectral methods."""
        return {
            "base_nodes": self._n_base,
            "fiber_nodes": self._n_fiber,
            "total_nodes": self._n_base + self._n_base * self._n_fiber,
            "structure": "graph_of_graphs",
        }

    @property
    def N_base(self) -> int:
        """Number of nodes in the base graph."""
        return self._n_base

    @property
    def N_fiber(self) -> int:
        """Number of nodes in each fiber graph."""
        return self._n_fiber

    def base_vertex_indices(self) -> list[int]:
        """Return indices of all base vertices."""
        return list(range(self._n_base))

    def fiber_vertex_indices(self, base_idx: int) -> list[int]:
        """Return indices of fiber vertices attached to a base node.

        Parameters
        ----------
        base_idx : int
            Index of the base node

        Returns
        -------
        list[int]
            Vertex indices in the fiber
        """
        if base_idx < 0 or base_idx >= self._n_base:
            raise ValueError(f"base_idx must be in [0, {self._n_base})")
        fiber_start = self._n_base + base_idx * self._n_fiber
        return list(range(fiber_start, fiber_start + self._n_fiber))

    def anchor_vertex(self, base_idx: int) -> int:
        """Return the anchor vertex index for a fiber.

        Parameters
        ----------
        base_idx : int
            Index of the base node

        Returns
        -------
        int
            Global vertex index of the anchor node
        """
        if base_idx < 0 or base_idx >= self._n_base:
            raise ValueError(f"base_idx must be in [0, {self._n_base})")
        fiber_start = self._n_base + base_idx * self._n_fiber
        return fiber_start + self._anchor_indices[base_idx]

    def vertex_layer(self, v: int) -> Literal["base", "fiber"]:
        """Determine which layer a vertex belongs to.

        Parameters
        ----------
        v : int
            Vertex index

        Returns
        -------
        str
            'base' or 'fiber'
        """
        if v < 0 or v >= self.N:
            raise ValueError(f"Vertex index out of range: {v}")
        return "base" if v < self._n_base else "fiber"

    def can_use_separated_spectrum(self) -> bool:
        """Check if separated spectrum computation is valid.

        The Dirac spectral formula applies when:
        1. All fibers are identical (homogeneous)
        2. Anchor index is constant across all fibers

        Returns
        -------
        bool
            True if separated spectrum optimization is valid
        """
        if self._heterogeneous:
            return False
        # Check if all anchor indices are the same
        return len(set(self._anchor_indices)) == 1

    def get_base_graph(self) -> SignedGraphGT:
        """Return the base graph instance."""
        return self._base_graph_instance

    def get_fiber_graph(self) -> SignedGraphGT:
        """Return the prototype fiber graph instance.

        Note: For heterogeneous fibers, this returns the fiber at index 0.
        """
        return self._fiber_graph_instance

    def __repr__(self) -> str:
        return (
            f"GraphOfGraphsGT("
            f"base={self._base_graph_type}[{self._n_base}], "
            f"fiber={self._fiber_graph_type}[{self._n_fiber}], "
            f"N={self.N})"
        )


# Backward compatibility alias
GraphOfGraphs = GraphOfGraphsGT

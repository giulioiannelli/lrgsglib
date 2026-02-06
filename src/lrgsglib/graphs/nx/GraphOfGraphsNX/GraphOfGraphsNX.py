"""
GraphOfGraphsNX - Generalized hierarchical graph structure using NetworkX backend.

This module provides a flexible "graph of graphs" class that accepts arbitrary
SignedGraph types as base and fiber components, generalizing the Dirac lattice
pattern (DiracComb, DiracBrush).

Example
-------
>>> from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX
>>> gog = GraphOfGraphsNX(
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

import networkx as nx
import numpy as np

from ..MultispectralGraphNX import MultispectralGraphNX
from ._generators import (
    build_composite_nx_graph,
    extract_edge_signs_from_nx_graph,
    extract_edges_from_nx_graph,
    generate_anchor_edges,
    generate_fiber_edges,
)
from ._policies import AnchorPolicy, get_anchor_index, resolve_anchor_indices
from . import _spectral


__all__ = ["GraphOfGraphsNX", "GraphOfGraphs"]


# Graph type registry for name -> class resolution
_GRAPH_TYPE_REGISTRY: Dict[str, Callable[..., Any]] = {}


def _register_graph_types():
    """Lazily populate the graph type registry."""
    if _GRAPH_TYPE_REGISTRY:
        return

    # Import graph types on demand to avoid circular imports
    from ..Lattice2DNX import Lattice2DNX
    from ..Lattice3DNX import Lattice3DNX
    from ..ErdosRenyiNX import ErdosRenyiNX
    from ..BarabasiAlbertNX import BarabasiAlbertNX
    from ..WattsStrogatzNX import WattsStrogatzNX
    from ..StochasticBlockModelNX import StochasticBlockModelNX

    _GRAPH_TYPE_REGISTRY.update({
        # Lattice types
        "Lattice2D": Lattice2DNX,
        "Lattice2DNX": Lattice2DNX,
        "Lattice3D": Lattice3DNX,
        "Lattice3DNX": Lattice3DNX,
        # Random types
        "ErdosRenyi": ErdosRenyiNX,
        "ErdosRenyiNX": ErdosRenyiNX,
        "BarabasiAlbert": BarabasiAlbertNX,
        "BarabasiAlbertNX": BarabasiAlbertNX,
        "WattsStrogatz": WattsStrogatzNX,
        "WattsStrogatzNX": WattsStrogatzNX,
        "StochasticBlockModel": StochasticBlockModelNX,
        "StochasticBlockModelNX": StochasticBlockModelNX,
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


class GraphOfGraphsNX(MultispectralGraphNX):
    """Generalized graph of graphs structure using NetworkX backend.

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
    >>> gog = GraphOfGraphsNX(
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
        # Store configuration (before super().__init__ calls _generate)
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

        # These will be set by _generate()
        self._base_graph_instance: Any = None
        self._fiber_graph_instance: Any = None
        self._n_base: int = 0
        self._n_fiber: int = 0
        self._anchor_indices: list[int] = []
        self._gog_structure: dict = {}

        # Spectral caching
        self._base_eigenvalues = None
        self._fiber_laplacian = None

        # Initialize base class (which calls _generate)
        super().__init__(pflip=pflip, seed=seed, **kwargs)

    def _instantiate_base_graph(self) -> Any:
        """Create the base graph instance."""
        cls = _resolve_graph_type(self._base_graph_type)
        # Use a separate seed for base graph
        base_seed = None
        if self._seed is not None:
            base_seed = int(self._rng.integers(0, 2**31))
        return cls(**self._base_params, seed=base_seed)

    def _instantiate_fiber_graph(self, base_idx: int) -> Any:
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

    def _generate(self) -> nx.Graph:
        """Generate the composite graph-of-graphs structure.

        Vertex layout:
        [--- Base nodes ---][--- Fiber 0 ---][--- Fiber 1 ---] ... [--- Fiber N_base-1 ---]
           0 to N_base-1         N_fiber          N_fiber                  N_fiber

        Returns
        -------
        nx.Graph
            The generated NetworkX graph with 'sign' edge attribute
        """
        # Generate base and fiber graphs to get dimensions
        self._base_graph_instance = self._instantiate_base_graph()
        self._fiber_graph_instance = self._instantiate_fiber_graph(0)

        self._n_base = self._base_graph_instance.N
        self._n_fiber = self._fiber_graph_instance.N

        # Resolve anchor indices for all fibers
        self._anchor_indices = resolve_anchor_indices(
            self._anchor_policy, self._n_base, self._n_fiber, self._rng
        )

        total_nodes = self._n_base + self._n_base * self._n_fiber

        # Collect all edges
        base_edges = extract_edges_from_nx_graph(self._base_graph_instance.G)
        base_signs = extract_edge_signs_from_nx_graph(self._base_graph_instance.G)

        fiber_edges_list = []
        fiber_signs_list = []

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
            fiber_edges_list.append(fiber_edges)
            fiber_signs_list.append(extract_edge_signs_from_nx_graph(fiber_graph.G))

        anchor_edges = generate_anchor_edges(
            self._n_base, self._n_fiber, self._anchor_indices
        )

        # Build composite graph
        H = build_composite_nx_graph(
            n_total=total_nodes,
            base_edges=base_edges,
            fiber_edges_list=fiber_edges_list,
            anchor_edges=anchor_edges,
            base_signs=base_signs,
            fiber_signs_list=fiber_signs_list,
            n_base=self._n_base,
            n_fiber=self._n_fiber,
        )

        # Store for compatibility (MultispectralGraphNX expects H)
        self.H = H

        # Build metadata after generation
        self._gog_structure = self._build_metadata()

        return H

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
            "can_use_separated_spectrum": self._can_use_separated_spectrum(),
        }

    def _can_use_separated_spectrum(self) -> bool:
        """Internal check for separated spectrum validity."""
        if self._heterogeneous:
            return False
        return len(set(self._anchor_indices)) == 1

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

    @property
    def base_graph(self) -> nx.Graph:
        """The base graph instance (NX compatibility with DiracLattice)."""
        return self._base_graph_instance.G

    @property
    def fiber_graph(self) -> nx.Graph:
        """The fiber graph instance (NX compatibility with DiracLattice)."""
        return self._fiber_graph_instance.G

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
        return self._can_use_separated_spectrum()

    def is_dirac_structure(self) -> bool:
        """True if structure supports separated spectrum (like Dirac lattices)."""
        return self.can_use_separated_spectrum()

    def get_base_graph(self) -> Any:
        """Return the base graph instance."""
        return self._base_graph_instance

    def get_fiber_graph(self) -> Any:
        """Return the prototype fiber graph instance.

        Note: For heterogeneous fibers, this returns the fiber at index 0.
        """
        return self._fiber_graph_instance

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        return self._n_base + self._n_base * self._n_fiber

    def __repr__(self) -> str:
        return (
            f"GraphOfGraphsNX("
            f"base={self._base_graph_type}[{self._n_base}], "
            f"fiber={self._fiber_graph_type}[{self._n_fiber}], "
            f"N={self.N})"
        )


# Attach spectral methods to class
GraphOfGraphsNX.compute_base_spectrum = _spectral.compute_base_spectrum
GraphOfGraphsNX.compute_fiber_laplacian = _spectral.compute_fiber_laplacian
GraphOfGraphsNX.compute_separated_spectrum = _spectral.compute_separated_spectrum
GraphOfGraphsNX.get_base_fiber_dimensions = _spectral.get_base_fiber_dimensions


# Backward compatibility alias
GraphOfGraphs = GraphOfGraphsNX

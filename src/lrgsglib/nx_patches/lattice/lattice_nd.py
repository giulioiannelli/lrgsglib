"""
LatticeND: Abstract base class for N-dimensional lattices.

This module provides a unified interface for lattice graphs of arbitrary
dimensionality, with common functionality shared between 2D and 3D lattices.

Example
-------
>>> from lrgsglib.nx_patches.lattice import Lattice2D
>>> lat = Lattice2D(side1=10, geo='sqr', pflip=0.1)
>>> lat.get_expected_num_nodes()
100
"""

from abc import abstractmethod
from typing import Any, Sequence
import random
import warnings

import networkx as nx
import numpy as np

from ...config.const import SG_GRAPHINT_REPR
from ..SignedGraph.SignedGraph import SignedGraph


class LatticeND(SignedGraph):
    """
    Abstract base class for N-dimensional signed lattices.

    This class provides the common interface and functionality shared by
    Lattice2D and Lattice3D. It should not be instantiated directly.

    Parameters
    ----------
    dimensions : tuple[int, ...]
        Lattice dimensions (e.g., (10, 10) for 2D, (8, 8, 8) for 3D).
    geo : str
        Lattice geometry identifier.
    pbc : bool, default True
        Whether to use periodic boundary conditions.
    fbc_val : float, default 1.0
        Fixed boundary value for non-periodic boundaries.
    with_positions : bool, default False
        Whether to compute and store node positions.
    only_const_mode : bool, default False
        If True, populate metadata only without building the graph.
    **kwargs
        Additional arguments passed to SignedGraph.

    Attributes
    ----------
    ndim : int
        Number of dimensions.
    dimensions : tuple[int, ...]
        Lattice dimensions.
    syshape : tuple[int, ...]
        System shape (may differ from dimensions for certain geometries).
    syshapePth : str
        Path string representation of system shape.
    pbc : bool
        Periodic boundary conditions flag.
    geo : str
        Geometry identifier.
    z : int
        Coordination number (average node degree).
    node_multiplier : int
        Multiplier for node count (>1 for basis atoms like BCC, FCC).
    H : nx.Graph
        Coordinate-labeled graph representation.
    G : nx.Graph
        Integer-labeled graph representation.

    Notes
    -----
    Subclasses must implement:
    - `_get_geometry_registry()`: Return dict mapping geometry names to configs
    - `_init_lattice()`: Build the actual lattice graph
    - `get_central_edge()`: Return the central edge of the lattice
    """

    def __init__(
        self,
        dimensions: tuple[int, ...],
        geo: str,
        pbc: bool = True,
        fbc_val: float = 1.0,
        with_positions: bool = False,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.only_const_mode = only_const_mode
        self.ndim = len(dimensions)
        self.dimensions = self._validate_dimensions(dimensions)
        self.pbc = pbc
        self.fbc_val = fbc_val
        self.with_positions = with_positions

        # Initialize geometry (must be implemented by subclass)
        self._init_geometry(geo)

        # Build or stub the graph
        if not only_const_mode:
            self._init_lattice()
        else:
            self.G = nx.Graph()
            self._init_metadata_only()

        # Initialize parent SignedGraph
        super().__init__(self.G, **kwargs)

    def _validate_dimensions(self, dimensions: tuple[int, ...]) -> tuple[int, ...]:
        """Validate and normalize dimensions."""
        if not all(isinstance(d, int) and d > 0 for d in dimensions):
            raise ValueError("All dimensions must be positive integers.")
        return tuple(dimensions)

    @abstractmethod
    def _get_geometry_registry(self) -> dict[str, dict]:
        """
        Return geometry configuration registry.

        Returns
        -------
        dict
            Mapping from geometry name to configuration dict with keys:
            - 'z': coordination number
            - 'gen': generator function
            - 'node_mult': node multiplier (optional, default 1)

        Example
        -------
        >>> def _get_geometry_registry(self):
        ...     return {
        ...         'sqr': {'z': 4, 'gen': squared_lattice_graph},
        ...         'tri': {'z': 6, 'gen': triangular_lattice_graph},
        ...     }
        """
        raise NotImplementedError

    @abstractmethod
    def _init_geometry(self, geo: str) -> None:
        """
        Initialize geometry-specific attributes.

        Parameters
        ----------
        geo : str
            Geometry identifier (e.g., 'sqr', 'tri', 'sc', 'bcc').
        """
        raise NotImplementedError

    @abstractmethod
    def _init_lattice(self) -> None:
        """Build the lattice graph with coordinate (H) and integer (G) labels."""
        raise NotImplementedError

    def _init_metadata_only(self) -> None:
        """Initialize metadata without building the graph."""
        self.node_multiplier = getattr(self, 'node_multiplier', 1)
        total_nodes = self.node_multiplier * int(np.prod(self.dimensions))

        if all(d == self.dimensions[0] for d in self.dimensions):
            self.syshapePth = f"N={total_nodes}"
        else:
            dim_part = '_'.join([f"L{i}={d}" for i, d in enumerate(self.dimensions)])
            self.syshapePth = f"{dim_part}_N={total_nodes}"

        self.syshape = self.dimensions

    def get_expected_num_nodes(self) -> int:
        """
        Return the expected number of nodes for this lattice.

        Returns
        -------
        int
            Total number of nodes.
        """
        return int(self.node_multiplier * np.prod(self.dimensions))

    @abstractmethod
    def get_central_edge(self, on_g: str = SG_GRAPHINT_REPR) -> tuple:
        """
        Return the central edge of the lattice.

        Parameters
        ----------
        on_g : str, default 'G'
            Graph representation to use ('G' for integer, 'H' for coordinate).

        Returns
        -------
        tuple
            Edge as (node1, node2) tuple.
        """
        raise NotImplementedError


class LatticeNwContainerBase(dict):
    """
    Base class for lattice network container patterns.

    This provides common functionality for both 2D and 3D lattice
    nwContainer classes, including star patterns and random edge selection.

    Parameters
    ----------
    lattice : LatticeND
        The parent lattice object.
    iterable : list, optional
        Initial keys to populate.
    constant : Any, optional
        Constant value to assign to initial keys.
    **kwargs
        Additional dict initialization arguments.
    """

    def __init__(
        self,
        lattice: LatticeND,
        iterable: Sequence = (),
        constant: Any = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.update((key, constant) for key in iterable)
        self.l = lattice
        self.rd = self.l.graph_reprs
        self.rNodeFlip = {
            g: random.sample(list(self.l.gr[g].nodes()), self.l.nflip)
            for g in self.rd
        }
        self.centedge = {g: self.l.get_central_edge(g) for g in self.rd}

        # Initialize common patterns
        self['single'] = {g: [self.centedge[g]] for g in self.rd}
        self['singleXERR'] = {
            g: self.get_links_XERR(self.centedge[g][0], g)
            for g in self.rd
        }
        self['rand'] = {g: list(self.l.fleset[g]) for g in self.rd}
        self['randXERR'] = {
            g: self._get_rand_xerr_pattern(g)
            for g in self.rd
        }

    def get_links_XERR(self, node: Any, on_g: str) -> list[tuple]:
        """
        Get star pattern: all edges incident to a node.

        Parameters
        ----------
        node : Any
            The center node.
        on_g : str
            Graph representation ('G' or 'H').

        Returns
        -------
        list[tuple]
            List of edges as (node, neighbor) tuples.
        """
        return [(node, nn) for nn in self.l.get_graph_neighbors(node, on_g)]

    def _get_rand_xerr_pattern(self, on_g: str) -> list[tuple]:
        """Get random XERR pattern avoiding fully negative neighborhoods."""
        from ...config.const import COUNT_XERR_PATTERNS

        if COUNT_XERR_PATTERNS:
            return [
                k for i in self.rNodeFlip[on_g]
                for k in self.get_links_XERR(i, on_g)
            ]

        tmplst = list(self.rNodeFlip[on_g])
        grph = self.l.gr[on_g]
        idx = 0
        patternList = []

        while idx < len(tmplst):
            leval = [
                all(nnn['weight'] == -1 for nnn in grph[nn].values())
                for nn in grph.neighbors(tmplst[idx])
            ]
            if any(leval):
                tmplst.pop(idx)
            else:
                patternList.extend(self.get_links_XERR(tmplst[idx], on_g))
                idx += 1

        return list(set(patternList))

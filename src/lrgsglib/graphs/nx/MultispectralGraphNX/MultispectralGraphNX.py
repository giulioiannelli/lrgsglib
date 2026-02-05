"""Abstract base class for multispectral graph structures.

This module provides MultispectralGraphNX, an abstract base class for graph types
with multiple spectral dimensions. Subclasses implement specific generators:
  - MultiplicativeCascadeGraph
  - VicsekGraph
  - DiracLatticeGraph (with DiracCombGraph and DiracBrushGraph subclasses)
"""

from __future__ import annotations

import networkx as nx

from ....config.const import (
    MSG_STDFN,
    MSG_SGPATH,
    MSG_ONLY_CONST_MODE,
    MSG_PHTABB,
)
from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


__all__ = ["MultispectralGraphNX"]


class MultispectralGraphNX(SignedGraphNX):
    """Abstract base class for multispectral graph generators.

    Subclasses must implement:
    - _generate() -> nx.Graph
    - get_expected_num_nodes() -> int
    """

    def __init__(
        self,
        stdFnameSFFX: str = MSG_STDFN,
        sgpathn: str = MSG_SGPATH,
        only_const_mode: bool = MSG_ONLY_CONST_MODE,
        with_positions: bool = False,
        **kwargs,
    ) -> None:
        """Initialize MultispectralGraphNX base class.

        Parameters
        ----------
        stdFnameSFFX : str
            Standard filename suffix
        sgpathn : str
            Signed graph path name
        only_const_mode : bool
            If True, skip graph generation (for metadata-only usage)
        with_positions : bool
            If True, attach node positions from grid coordinates
        **kwargs
            Additional arguments passed to SignedGraphNX
        """
        self.only_const_mode = only_const_mode
        self.with_positions = with_positions
        self.positions: dict[int, tuple[float, float]] | None = None
        self.H: nx.Graph | None = None
        self.std_fname = f"{MSG_PHTABB}{stdFnameSFFX}" if stdFnameSFFX else MSG_PHTABB
        self.sgpathn = sgpathn or MSG_SGPATH

        self._verify_pflip(kwargs.get("pflip", 0.0))

        # Generate graph or create empty placeholder
        if not self.only_const_mode:
            self.G = self._generate()
            self.G = self._relabel_with_metadata(self.G, with_positions)
        else:
            self.G = nx.Graph()

        # Set system shape path (only if not already set by subclass)
        if not hasattr(self, 'syshapePth') or self.syshapePth is None:
            try:
                self.syshapePth = f"N={self.G.number_of_nodes()}"
            except Exception:
                self.syshapePth = ""

        super().__init__(self.G, **kwargs)

    def _generate(self) -> nx.Graph:
        """Generate the graph. Must be implemented by subclasses.

        Returns
        -------
        nx.Graph
            Generated graph with arbitrary node labels

        Raises
        ------
        NotImplementedError
            If not implemented by subclass
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _generate()"
        )

    def get_expected_num_nodes(self) -> int:
        """Return expected node count from generator parameters.

        Must be implemented by subclasses.

        Returns
        -------
        int
            Expected number of nodes

        Raises
        ------
        NotImplementedError
            If not implemented by subclass
        """
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement get_expected_num_nodes()"
        )

    def _relabel_with_metadata(self, H: nx.Graph, with_positions: bool) -> nx.Graph:
        """Relabel nodes to integers and optionally attach positions.

        Parameters
        ----------
        H : nx.Graph
            Graph with arbitrary node labels
        with_positions : bool
            If True, extract positions from tuple node labels

        Returns
        -------
        nx.Graph
            Graph with integer node labels and metadata
        """
        mapping = {node: idx for idx, node in enumerate(sorted(H.nodes()))}
        G_int = nx.relabel_nodes(H, mapping, copy=True)
        nx.set_node_attributes(
            G_int,
            {idx: node for node, idx in mapping.items()},
            "original_label"
        )

        if with_positions:
            pos_attr: dict[int, tuple[float, float]] = {}
            for node, idx in mapping.items():
                if isinstance(node, tuple) and len(node) == 2:
                    x, y = float(node[1]), float(-node[0])
                    pos_attr[idx] = (x, y)
            if pos_attr:
                nx.set_node_attributes(G_int, pos_attr, "pos")
                self.positions = pos_attr

        return G_int

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(N={self.G.number_of_nodes()})"

    def __str__(self) -> str:
        return self.__repr__()

"""
Unified position handling mixin for graph types.

Provides standardized methods for computing and storing node positions
across different graph families (lattices, random graphs, etc.).
"""
from typing import Any, Callable, Optional, TYPE_CHECKING

import networkx as nx

if TYPE_CHECKING:
    from ..SignedGraphNX.SignedGraphNX import SignedGraphNX


class PositionMixin:
    """
    Standardized position handling for graph types.

    This mixin provides common methods for computing, storing, and managing
    node positions for visualization. It supports multiple position computation
    strategies: from coordinate labels, custom layout functions, or spring layout.

    Attributes
    ----------
    with_positions : bool
        Whether positions are enabled for this graph.

    Methods
    -------
    _init_positions(with_positions, layout_func=None)
        Initialize position handling for the graph.
    _compute_positions(layout_func=None)
        Compute and set node positions.
    _positions_from_labels()
        Extract positions from tuple-labeled nodes (H graph).
    _sync_positions_to_G()
        Synchronize positions from H to G representation.

    Examples
    --------
    >>> class MyGraph(PositionMixin, SignedGraph):
    ...     def __init__(self, ..., with_positions=False, **kwargs):
    ...         # Build graph...
    ...         super().__init__(self.G, **kwargs)
    ...         self._init_positions(with_positions)
    """

    with_positions: bool

    def _init_positions(
        self,
        with_positions: bool,
        layout_func: Optional[Callable[[nx.Graph], dict]] = None,
    ) -> None:
        """
        Initialize position handling for the graph.

        Parameters
        ----------
        with_positions : bool
            Whether to compute and store node positions.
        layout_func : Callable[[nx.Graph], dict], optional
            Custom layout function that takes a graph and returns a dict
            mapping nodes to (x, y) positions. If None, uses automatic
            detection based on graph structure.

        Notes
        -----
        This method should be called after graph construction but can be
        called either before or after SignedGraph initialization depending
        on the subclass needs.
        """
        self.with_positions = with_positions
        if with_positions and not getattr(self, "only_const_mode", False):
            self._compute_positions(layout_func)

    def _compute_positions(
        self,
        layout_func: Optional[Callable[[nx.Graph], dict]] = None,
    ) -> None:
        """
        Compute and set node positions on the graph.

        Parameters
        ----------
        layout_func : Callable[[nx.Graph], dict], optional
            Custom layout function. If None:
            - Uses coordinate labels from H if available
            - Falls back to spring layout otherwise

        Notes
        -----
        Positions are stored as 'pos' node attribute on both G and H
        representations (if H exists).
        """
        if layout_func is not None:
            pos = layout_func(self.G)
        elif hasattr(self, "H") and self.H is not None:
            pos = self._positions_from_labels()
        else:
            pos = nx.spring_layout(self.G)

        nx.set_node_attributes(self.G, pos, "pos")

        # Sync to H if available
        if hasattr(self, "H") and self.H is not None:
            self._sync_positions_to_H(pos)

    def _positions_from_labels(self) -> dict[Any, tuple[float, float]]:
        """
        Extract positions from tuple-labeled nodes in H representation.

        Returns
        -------
        dict
            Dictionary mapping G nodes to (x, y) positions.

        Notes
        -----
        Assumes H has tuple-labeled nodes where the tuple represents
        coordinates. For 2D: (x, y), for 3D: projects to 2D.
        """
        # Get H nodes and their labels
        h_nodes = list(self.H.nodes())

        if not h_nodes:
            return {}

        first_node = h_nodes[0]

        # Handle different coordinate dimensions
        if isinstance(first_node, tuple):
            if len(first_node) == 2:
                # 2D coordinates - use directly
                pos_h = {node: (float(node[0]), float(node[1])) for node in h_nodes}
            elif len(first_node) == 3:
                # 3D coordinates - project to 2D
                from ....utils.basic.linalg import project_3d_to_2d

                theta = getattr(self, "theta", 0.5236)  # pi/6
                phi = getattr(self, "phi", 0.5236)
                pos_h = {
                    node: project_3d_to_2d(node[0], node[1], node[2], theta, phi)
                    for node in h_nodes
                }
            else:
                # Higher dimensions - use first 2 coords
                pos_h = {node: (float(node[0]), float(node[1])) for node in h_nodes}
        else:
            # Non-tuple nodes - fall back to spring layout
            return nx.spring_layout(self.G)

        # Map H positions to G nodes
        if hasattr(self, "map_node") and "G" in self.map_node:
            node_mapping = self.map_node["G"].get("H", {})
            pos_g = {node_mapping.get(h_node, h_node): position
                     for h_node, position in pos_h.items()}
        else:
            # No mapping available yet, use H positions directly
            pos_g = pos_h

        return pos_g

    def _sync_positions_to_H(self, pos_g: dict) -> None:
        """
        Synchronize positions from G to H representation.

        Parameters
        ----------
        pos_g : dict
            Dictionary mapping G nodes to positions.
        """
        if not hasattr(self, "H") or self.H is None:
            return

        if hasattr(self, "map_node") and "H" in self.map_node:
            node_mapping = self.map_node["H"].get("G", {})
            pos_h = {node_mapping.get(g_node, g_node): position
                     for g_node, position in pos_g.items()}
        else:
            pos_h = pos_g

        nx.set_node_attributes(self.H, pos_h, "pos")

    def _sync_positions_to_G(self) -> None:
        """
        Synchronize positions from H to G representation.

        Notes
        -----
        Useful when H has been modified and G needs to be updated.
        Called after SignedGraph initialization when mappings are available.
        """
        if not hasattr(self, "H") or self.H is None:
            return

        pos_h = nx.get_node_attributes(self.H, "pos")
        if not pos_h:
            return

        if hasattr(self, "map_node") and "G" in self.map_node:
            node_mapping = self.map_node["G"].get("H", {})
            pos_g = {node_mapping.get(h_node, h_node): position
                     for h_node, position in pos_h.items()}
        else:
            pos_g = pos_h

        nx.set_node_attributes(self.G, pos_g, "pos")

        # Update graph representations dict if available
        if hasattr(self, "gr"):
            self.gr["G"] = self.G
            if "H" in self.gr:
                self.gr["H"] = self.H

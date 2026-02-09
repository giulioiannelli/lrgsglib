"""
MultispectralGraphGT - Abstract base class for multispectral graph-tool graphs.

Provides shared infrastructure for graph types that embed nodes in a
probability/position space (e.g. MultiplicativeCascade, Vicsek).
"""

from __future__ import annotations

from abc import abstractmethod
from typing import Optional

import numpy as np
from numpy.typing import NDArray

try:
    import graph_tool as gt
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False

from ..SignedGraphGT import SignedGraphGT


__all__ = ["MultispectralGraphGT"]


class MultispectralGraphGT(SignedGraphGT):
    """Abstract base for multispectral graphs with position data.

    Subclasses must implement ``_generate()`` and
    ``get_expected_num_nodes()``.

    Attributes
    ----------
    probability_matrix : np.ndarray or None
        The computed probability matrix.
    sample_fraction : float
        Fraction of nodes sampled (1.0 for Vicsek).
    syshapePth : str
        Parameter string for file paths.
    stdFnameSFFX : str
        Standard filename suffix.
    sgpathn : str
        Signed-graph path name.
    """

    # Subclasses set these before calling super().__init__
    probability_matrix: Optional[NDArray[np.floating]] = None
    sample_fraction: float = 1.0
    syshapePth: str = ""
    stdFnameSFFX: str = ""
    sgpathn: str = ""

    # ------------------------------------------------------------------ #
    #  Position infrastructure                                            #
    # ------------------------------------------------------------------ #

    @property
    def pos(self) -> Optional["gt.VertexPropertyMap"]:
        """Vertex position property map (``vector<double>``)."""
        return self.G.vertex_properties.get("pos")

    def get_positions(self) -> NDArray[np.floating]:
        """Return positions as an (N, 2) numpy array.

        Returns
        -------
        np.ndarray
            Shape ``(N, 2)`` with ``(x, y)`` coordinates per vertex.

        Raises
        ------
        RuntimeError
            If no position property is stored on the graph.
        """
        pos = self.pos
        if pos is None:
            raise RuntimeError(
                "No 'pos' vertex property stored. "
                "Positions are set during graph generation."
            )
        return pos.get_2d_array([0, 1]).T  # shape (N, 2)

    # ------------------------------------------------------------------ #
    #  Drawing                                                            #
    # ------------------------------------------------------------------ #

    def draw(self, ax=None, layout: Optional[str] = None, **kwargs):
        """Draw the graph with edge colors based on sign.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.  Creates a new figure if *None*.
        layout : str, optional
            Force a layout algorithm (``"sfdp"``).  By default the stored
            grid positions are used; if none exist, ``sfdp_layout`` is
            used automatically.
        **kwargs
            Forwarded to ``gt.graph_draw``.

        Returns
        -------
        matplotlib.axes.Axes
        """
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 8))

        from graph_tool.draw import graph_draw, sfdp_layout

        # Resolve positions
        pos_prop = self.pos
        if layout == "sfdp" or pos_prop is None:
            pos_prop = sfdp_layout(self.G)

        # Edge colours by sign
        edge_color = self.G.new_edge_property("vector<double>")
        sign = self.G.edge_properties.get("sign")
        for e in self.G.edges():
            if sign is not None and sign[e] == -1:
                edge_color[e] = [0.8, 0.2, 0.2, 0.8]  # red
            else:
                edge_color[e] = [0.2, 0.6, 0.2, 0.8]  # green

        draw_opts = {
            "pos": pos_prop,
            "edge_color": edge_color,
            "vertex_size": 2,
            "edge_pen_width": 0.3,
            "mplfig": ax,
        }
        draw_opts.update(kwargs)

        graph_draw(self.G, **draw_opts)
        ax.set_aspect("equal")
        return ax

    # ------------------------------------------------------------------ #
    #  Abstract interface                                                 #
    # ------------------------------------------------------------------ #

    @abstractmethod
    def _generate(self) -> "gt.Graph":
        """Build and return the graph-tool Graph."""

    @abstractmethod
    def get_expected_num_nodes(self) -> int:
        """Return the expected number of nodes."""

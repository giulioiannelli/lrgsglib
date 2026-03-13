"""N-dimensional lattice using graph-tool, inheriting SignedGraphGT."""
from __future__ import annotations

from typing import Optional, Tuple, Union

import numpy as np

try:
    import graph_tool.all as gt
    import graph_tool.generation as gen
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False

from ..SignedGraphGT import SignedGraphGT


class LatticeNDGT(SignedGraphGT):
    """
    N-dimensional cubic lattice graph using graph-tool backend.

    Inherits from SignedGraphGT for full signed graph functionality
    (spectral analysis, entropy, clustering, dynamics integration).

    Parameters
    ----------
    shape : tuple of int
        Dimensions of the lattice (e.g., (10, 10, 10) for 3D cubic).
    periodic : bool, default True
        Whether to use periodic boundary conditions.
    pflip : float, default 0.0
        Fraction of edges to flip to negative (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.

    Examples
    --------
    >>> lat = LatticeNDGT(shape=(10, 10, 10), periodic=True)
    >>> lat.flip_random_fract_edges()
    >>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")
    """

    def __init__(
        self,
        shape: Union[tuple, list],
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ):
        if not GT_AVAILABLE:
            raise ImportError(
                "graph-tool is not installed. Install with: "
                "conda install -c conda-forge graph-tool"
            )

        self.shape = tuple(shape)
        self.periodic = periodic

        if seed is not None:
            np.random.seed(seed)
            gt.seed_rng(seed)

        # Create lattice
        G = gen.lattice(list(self.shape), periodic=periodic)

        # Add sign property (all positive initially)
        sign = G.new_edge_property("int", val=1)
        G.edge_properties["sign"] = sign

        # Set syshapePth for file naming
        shape_str = "x".join(map(str, self.shape))
        self._syshapePth = f"N={G.num_vertices()}"

        # Pop any NX-specific kwargs that might leak through facades
        for key in ("sgpathn", "stdFnameSFFX", "only_const_mode",
                     "path_data", "path_plot", "init_nw_dict"):
            kwargs.pop(key, None)

        # Initialize SignedGraphGT
        super().__init__(
            G=G,
            pflip=pflip,
            seed=seed,
            sgpathn=f"lnd_{shape_str}_gt",
        )

    @property
    def ndim(self) -> int:
        """Number of dimensions."""
        return len(self.shape)


__all__ = ["LatticeNDGT"]

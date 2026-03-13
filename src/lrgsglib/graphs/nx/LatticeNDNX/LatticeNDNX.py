"""
LatticeNDNX: Generic N-dimensional cubic lattice using NetworkX.

For specialized geometries (tri, hex, bcc, fcc, ...) use the
dimension-specific classes Lattice2DNX or Lattice3DNX directly.

Example
-------
>>> from lrgsglib.graphs.nx import LatticeNDNX
>>> lat = LatticeNDNX(shape=(5, 5, 5, 5), periodic=True, pflip=0.1)
"""

from __future__ import annotations

from typing import Optional, Tuple, Union

import networkx as nx

from ..SignedGraphNX.SignedGraphNX import SignedGraphNX
from ..funcs.lattice import LatticeND_graph_FastPatch


class LatticeNDNX(SignedGraphNX):
    """
    Generic N-dimensional cubic lattice using NetworkX backend.

    Creates a hypercubic lattice with nearest-neighbor connections.
    For specialized 2D/3D geometries (triangular, hexagonal, BCC, FCC, ...),
    use ``Lattice2DNX`` or ``Lattice3DNX`` instead.

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
    >>> lat = LatticeNDNX(shape=(8, 8, 8), periodic=True)
    >>> print(f"Nodes: {lat.N}, Edges: {lat.Ne}")
    Nodes: 512, Edges: 1536

    >>> lat4d = LatticeNDNX(shape=(4, 4, 4, 4), periodic=False)
    >>> print(f"4D lattice: {lat4d.N} nodes")
    4D lattice: 256 nodes
    """

    sgpathn = "lattice_nd"

    def __init__(
        self,
        shape: Union[tuple, list],
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        **kwargs,
    ) -> None:
        self.shape = tuple(shape)
        self.periodic = periodic

        # Build lattice using existing fast generator (tuple-labeled nodes)
        H = LatticeND_graph_FastPatch(self.shape, periodic=periodic)

        # Convert to integer-labeled graph
        G = nx.convert_node_labels_to_integers(H)

        # Store coordinate mapping
        self.H = H
        self.syshape = self.shape
        shape_str = "x".join(map(str, self.shape))
        self.syshapePth = f"N={G.number_of_nodes()}"

        # Initialize SignedGraphNX
        kwargs.setdefault("make_dir_tree", False)
        super().__init__(G, pflip=pflip, seed=seed, **kwargs)

    @property
    def ndim(self) -> int:
        """Number of dimensions."""
        return len(self.shape)


__all__ = ["LatticeNDNX"]

"""Loader methods for SignedGraphGT (NX-compatible binary I/O)."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ....config.const import BIN
from ....config.funcs import build_p_fname

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphGT import SignedGraphGT


def load_eigV_all(
    self: "SignedGraphGT",
    exName: str = "",
    ext: str = BIN,
    path_sgdata: Optional[Path] = None,
    binarize: bool = True,
    reshape: bool = True,
) -> NDArray:
    """Load eigenvectors from a binary file.

    Parameters
    ----------
    exName : str
        Filename suffix.
    ext : str
        File extension (default ``.bin``).
    path_sgdata : Path, optional
        Override directory (default ``self.path_graph``).
    binarize : bool
        If True, read as int8; else float64.
    reshape : bool
        If True, reshape flat array to ``(-1, N)``.

    Returns
    -------
    ndarray
        Loaded eigenvectors.
    """
    self._ensure_paths()
    fname = build_p_fname("eigV", self.pflip, out_suffix=exName, ext=ext)
    base = path_sgdata or self.path_graph
    fpath = base / fname

    if not fpath.exists():
        raise FileNotFoundError(f"Eigenvector file not found: {fpath}")

    dtype = np.int8 if binarize else np.float64
    data = np.fromfile(fpath, dtype=dtype)

    if reshape and data.size > 0:
        data = data.reshape(-1, self.N)

    self.eigVPname = fpath
    return data


def set_edgel_from_bin(
    self: "SignedGraphGT",
    file_path: str | Path,
    mode: str = "numpy",
) -> None:
    """Load an edge list from binary and add edges to the GT graph.

    Binary format: structured numpy array ``(uint64 i, uint64 j, float64 w_ij)``.

    Parameters
    ----------
    file_path : str or Path
        Path to the binary edge list file.
    mode : str
        ``'numpy'`` (structured array read).
    """
    fpath = Path(file_path)
    if not fpath.exists():
        raise FileNotFoundError(f"Edge list file not found: {fpath}")

    dtype = [("i", np.uint64), ("j", np.uint64), ("w_ij", np.float64)]
    edges = np.fromfile(fpath, dtype=dtype)

    # Determine number of nodes
    max_node = max(edges["i"].max(), edges["j"].max()) if len(edges) else 0
    n_nodes = int(max_node) + 1

    # Build a fresh graph
    from graph_tool import Graph
    G = Graph(directed=False)
    G.add_vertex(n_nodes)
    sign_prop = G.new_edge_property("int", val=1)
    G.edge_properties["sign"] = sign_prop

    for row in edges:
        e = G.add_edge(int(row["i"]), int(row["j"]))
        sign_prop[e] = int(np.sign(row["w_ij"]))

    self.G = G
    self.invalidate_cache()

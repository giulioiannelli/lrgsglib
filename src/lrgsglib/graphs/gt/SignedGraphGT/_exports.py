"""Export methods for SignedGraphGT (NX-compatible binary I/O).

File formats are identical to NX so C backends consume the same files.
"""

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


def export_eigV_all(
    self: "SignedGraphGT",
    exName: str = "",
    ext: str = BIN,
    binarize: bool = True,
    path_sgdata: Optional[Path] = None,
) -> None:
    """Export all eigenvectors to a binary file.

    Parameters
    ----------
    exName : str
        Suffix for the filename.
    ext : str
        File extension (default ``.bin``).
    binarize : bool
        If True, store sign(eigV) as int8; else store float64.
    path_sgdata : Path, optional
        Override output directory (default ``self.path_graph``).
    """
    if self.eigV is None:
        self.compute_laplacian_spectrum_weigV()

    self._ensure_paths()
    fname = build_p_fname("eigV", self.pflip, out_suffix=exName, ext=ext)
    base = path_sgdata or self.path_graph
    base.mkdir(parents=True, exist_ok=True)
    self.eigVPname = base / fname

    if binarize:
        # Column-major eigV: eigV[:, i] is i-th eigenvector
        out = np.array([
            np.sign(self.eigV[:, i]).astype(np.int8)
            for i in range(self.eigV.shape[1])
        ])
    else:
        # Row-major export: each row is an eigenvector
        out = self.eigV.T.astype(np.float64)

    with open(self.eigVPname, "wb") as f:
        out.tofile(f)


def export_adj_bin(
    self: "SignedGraphGT",
    exName: str = "",
    upper_triangle: bool = True,
) -> None:
    """Export signed adjacency matrix to binary file.

    Parameters
    ----------
    exName : str
        Suffix for the filename.
    upper_triangle : bool
        If True, export only upper triangle (saves ~50% disk).
    """
    self._ensure_paths()
    fname = build_p_fname("adj", self.pflip, out_suffix=exName, ext=BIN)
    base = self.path_graph
    base.mkdir(parents=True, exist_ok=True)
    self.path_exp_adj = base / fname

    A = self.get_signed_adjacency().astype(np.float64)
    if upper_triangle:
        A = np.triu(A)

    with open(self.path_exp_adj, "wb") as f:
        A.tofile(f)


def export_ising_clust(
    self: "SignedGraphGT",
    NoClust: int = 1,
    exName: str = "",
) -> None:
    """Export Ising cluster node indices to binary files.

    One file per cluster, named ``cl0_..., cl1_...``.

    Parameters
    ----------
    NoClust : int
        Number of clusters to export.
    exName : str
        Suffix for filenames.
    """
    self._ensure_paths()
    base = getattr(self, "path_ising", self.path_graph)
    base.mkdir(parents=True, exist_ok=True)

    clusters = getattr(self, "biggestClSet", [])
    if not clusters:
        logger.warning("No clusters to export.")
        return

    self.clPname = []
    n_export = min(NoClust, len(clusters))
    for i in range(n_export):
        fname = build_p_fname(f"cl{i}", self.pflip, out_suffix=exName, ext=BIN)
        fpath = base / fname
        node_arr = np.array(sorted(clusters[i]), dtype=np.uint64)
        with open(fpath, "wb") as f:
            node_arr.tofile(f)
        self.clPname.append(fpath)

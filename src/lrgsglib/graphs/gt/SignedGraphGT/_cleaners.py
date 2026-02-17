"""Cleanup methods for SignedGraphGT (NX-compatible file removal)."""

from __future__ import annotations

import logging
import os
from typing import TYPE_CHECKING

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphGT import SignedGraphGT


def remove_ising_clust_files(self: "SignedGraphGT") -> None:
    """Remove Ising cluster binary files."""
    for fpath in getattr(self, "clPname", []):
        try:
            if os.path.exists(fpath):
                os.remove(fpath)
        except OSError as e:
            logger.warning(f"Could not remove cluster file {fpath}: {e}")


def remove_edgl_file(self: "SignedGraphGT") -> None:
    """Remove exported edge list file."""
    fpath = getattr(self, "path_exp_edgl", None)
    if fpath and os.path.exists(fpath):
        try:
            os.remove(fpath)
        except OSError as e:
            logger.warning(f"Could not remove edge list file {fpath}: {e}")


def remove_eigV_file(self: "SignedGraphGT") -> None:
    """Remove exported eigenvector file."""
    fpath = getattr(self, "eigVPname", None)
    if fpath and os.path.exists(fpath):
        try:
            os.remove(fpath)
        except OSError as e:
            logger.warning(f"Could not remove eigV file {fpath}: {e}")


def remove_adj_file(self: "SignedGraphGT") -> None:
    """Remove exported adjacency file."""
    fpath = getattr(self, "path_exp_adj", None)
    if fpath and os.path.exists(fpath):
        try:
            os.remove(fpath)
        except OSError as e:
            logger.warning(f"Could not remove adj file {fpath}: {e}")


def remove_exported_files(self: "SignedGraphGT") -> None:
    """Remove all exported binary files."""
    if hasattr(self, "clPname"):
        remove_ising_clust_files(self)
    remove_edgl_file(self)
    remove_eigV_file(self)
    remove_adj_file(self)


def clean_gclutil(
    self: "SignedGraphGT",
    k: str = "",
    val: str = "",
    on_g: str = "graph",
) -> None:
    """Clear a graph clustering utility cache entry."""
    gyn = getattr(self, "_graph_yn", None)
    if gyn and (k, val) in gyn:
        gyn[(k, val)] = None

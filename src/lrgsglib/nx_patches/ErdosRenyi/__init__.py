"""ErdosRenyi module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.ErdosRenyi is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.ErdosRenyiNX import ErdosRenyiNX as ErdosRenyi

__all__ = ["ErdosRenyi"]

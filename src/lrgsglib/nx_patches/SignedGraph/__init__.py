"""Signed graph module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import sys
import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.SignedGraph is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.SignedGraphNX import SignedGraphNX as SignedGraph
from ...graphs.nx.SignedGraphNX import _backend

__all__ = ["SignedGraph", "_backend"]

# Fix submodule shadowing: when this module is imported as a submodule,
# it shadows the SignedGraph class in the parent nx_patches.__init__.py.
# We need to ensure that `from lrgsglib.nx_patches import SignedGraph` returns
# the class, not this module.
#
# Python binds submodules to parent attributes AFTER __init__.py runs, so we
# need to use an import hook or delay the fix. Instead, we'll remove ourselves
# from the parent's __dict__ so that __getattr__ can handle the lookup.
_parent = sys.modules.get("lrgsglib.nx_patches")
if _parent is not None and hasattr(_parent, "__dict__"):
    # Remove the module binding so __getattr__ can return the class
    _parent.__dict__.pop("SignedGraph", None)

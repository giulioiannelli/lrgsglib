"""
gt_patches.random — DEPRECATED.

All random graph GT implementations have moved to ``lrgsglib.graphs.gt``.
This shim re-exports for backward compatibility.
"""
import warnings as _warnings

_WARNED = set()


def __getattr__(name):
    _PUBLIC = {
        "ErdosRenyiGT", "BarabasiAlbertGT",
        "WattsStrogatzGT", "StochasticBlockModelGT",
    }
    if name in _PUBLIC:
        if name not in _WARNED:
            _warnings.warn(
                f"Importing {name!r} from 'lrgsglib.gt_patches.random' is "
                f"deprecated. Use 'lrgsglib.graphs.gt' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            _WARNED.add(name)

        if name == "ErdosRenyiGT":
            from lrgsglib.graphs.gt.ErdosRenyiGT import ErdosRenyiGT
            return ErdosRenyiGT
        elif name == "BarabasiAlbertGT":
            from lrgsglib.graphs.gt.BarabasiAlbertGT import BarabasiAlbertGT
            return BarabasiAlbertGT
        elif name == "WattsStrogatzGT":
            from lrgsglib.graphs.gt.WattsStrogatzGT import WattsStrogatzGT
            return WattsStrogatzGT
        elif name == "StochasticBlockModelGT":
            from lrgsglib.graphs.gt.StochasticBlockModelGT import (
                StochasticBlockModelGT,
            )
            return StochasticBlockModelGT

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__: list[str] = []

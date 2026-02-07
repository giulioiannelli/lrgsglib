"""
BipartiteFromDegreeSequenceGT: Backward-compatibility wrapper.

Use ``BipartiteGraphGT(top_degrees=..., bottom_degrees=...)`` instead.
"""

from __future__ import annotations

from typing import Optional, Sequence

from ..BipartiteGraphGT.BipartiteGraphGT import BipartiteGraphGT


__all__ = ["BipartiteFromDegreeSequenceGT"]


class BipartiteFromDegreeSequenceGT(BipartiteGraphGT):
    """Backward-compat wrapper.

    Use ``BipartiteGraphGT(top_degrees=..., bottom_degrees=...)``
    instead.

    Parameters
    ----------
    top_degrees : Sequence[int]
        Degree sequence for top nodes.
    bottom_degrees : Sequence[int]
        Degree sequence for bottom nodes.
    **kwargs
        Forwarded to ``BipartiteGraphGT``.
    """

    def __init__(
        self,
        top_degrees: Sequence[int],
        bottom_degrees: Sequence[int],
        **kwargs,
    ) -> None:
        super().__init__(
            top_degrees=top_degrees,
            bottom_degrees=bottom_degrees,
            **kwargs,
        )

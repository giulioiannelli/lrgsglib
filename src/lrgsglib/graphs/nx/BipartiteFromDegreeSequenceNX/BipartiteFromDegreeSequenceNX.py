"""
BipartiteFromDegreeSequenceNX: Backward-compatibility wrapper.

Use ``BipartiteGraphNX(top_degrees=..., bottom_degrees=...)`` instead.
"""

from typing import Sequence

from ..BipartiteGraphNX.BipartiteGraphNX import BipartiteGraphNX


__all__ = ["BipartiteFromDegreeSequenceNX"]


class BipartiteFromDegreeSequenceNX(BipartiteGraphNX):
    """Backward-compat wrapper.

    Use ``BipartiteGraphNX(top_degrees=..., bottom_degrees=...)``
    instead.

    Parameters
    ----------
    top_degrees : Sequence[int]
        Degree sequence for top nodes.
    bottom_degrees : Sequence[int]
        Degree sequence for bottom nodes.
    **kwargs
        Forwarded to ``BipartiteGraphNX``.
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

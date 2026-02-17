"""Order parameters — gap from shared, compute_pinf stays GT-local.

The gap implementation lives in ``graphs._shared._ordparams``.
``compute_pinf`` is defined in ``SignedGraphGT.py`` directly (uses
GT-specific cluster logic).
"""

from ..._shared._ordparams import (  # noqa: F401
    compute_gap,
    compute_gap_between,
    get_gap,
)

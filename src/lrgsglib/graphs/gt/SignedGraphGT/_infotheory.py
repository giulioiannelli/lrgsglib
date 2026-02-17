"""Information theory — re-exports from shared location.

The canonical implementation lives in ``graphs._shared._infotheory``.
"""

from ..._shared._infotheory import (  # noqa: F401
    compute_signed_laplacian_entropy,
    compute_renyi_entropy_profile,
    get_entropy,
    get_specific_heat,
    get_entropy_derivative,
    get_renyi_results,
)

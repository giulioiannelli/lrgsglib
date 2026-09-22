"""Information theory — re-exports from shared location.

The canonical implementation lives in ``graphs._shared._infotheory``.
"""

from ..._shared._infotheory import (  # noqa: F401
    compute_renyi_entropy_profile,
    compute_signed_laplacian_entropy,
    get_entropy,
    get_entropy_derivative,
    get_renyi_results,
    get_specific_heat,
)

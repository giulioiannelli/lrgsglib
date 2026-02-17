"""Order parameters — gap from shared, compute_pinf NX-specific."""

import numpy as np
from typing import Optional, TYPE_CHECKING
from ....config.const import SG_REPR

# Gap methods from shared engine-agnostic location
from ..._shared._ordparams import (  # noqa: F401
    compute_gap,
    compute_gap_between,
    get_gap,
)

if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX


def compute_pinf(
    self: "SignedGraphNX",
    which: int = 0,
    val: Optional[object] = 1,
    on_g: str = SG_REPR,
) -> None:
    """Compute the infinite-cluster probability (Pinf) for ``self``."""
    if val is None:
        val = 1

    clustd = np.array(self.get_eigV_cluster_sizes(which, True, val, on_g))
    self.Pinf = clustd[0] / self.N

    denominator = np.sum(clustd) - clustd[0]
    if denominator > 0:
        self.Pinf_var = np.sum(clustd @ clustd - clustd[0] ** 2) / denominator
    else:
        self.Pinf_var = 0.0

    if hasattr(self, "Pinf_dict"):
        self.Pinf_dict[which] = (self.Pinf, self.Pinf_var)
    else:
        self.Pinf_dict = {which: (self.Pinf, self.Pinf_var)}

import numpy as np
from typing import Optional, TYPE_CHECKING
from numpy.typing import NDArray

from ...utils.basic import dtype_numerical_precision
from ...utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

# reuse spectral building blocks
from ._spectral import (
    compute_laplacian_spectrum_weigV,
    make_eigV_transposed,
    make_eigV_column_major,
)


if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def compute_signed_laplacian_entropy(
    self,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    w_thresh: float = None,
    typf: type = np.float64,
    backend: str = 'numpy',
    transpose: Optional[bool] = None,
) -> None:
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    desired_transpose = (
        transpose if transpose is not None else getattr(self, "_eigV_is_transposed", True)
    )

    spectrum_missing = (
        not hasattr(self, "eigv") or self.eigv is None or len(self.eigv) != self.N
    )

    if spectrum_missing:
        compute_laplacian_spectrum_weigV(self, typf=typf, transpose=desired_transpose, backend=backend)
    else:
        current_transposed = getattr(self, "_eigV_is_transposed", False)
        if transpose is not None:
            if transpose and not current_transposed:
                make_eigV_transposed(self)
            elif not transpose and current_transposed:
                make_eigV_column_major(self)

    eigenvalues = np.asarray(self.eigv, dtype=typf)
    threshold = w_thresh if w_thresh is not None else dtype_numerical_precision(typf)

    normalized_entropy, entropy_derivative, variance_profile, time_grid = (
        compute_entropy_observables_from_eigenvalues(
            eigenvalues=eigenvalues,
            num_nodes=self.N,
            steps=steps,
            t1=t1,
            t2=t2,
            typf=typf,
            threshold=threshold,
            pad_last=False,
        )
    )

    self.entropy = normalized_entropy
    self.specific_heat = entropy_derivative
    self.variance_profile = variance_profile
    self.tauscale = time_grid
    self.entropy_params = {
        "steps": steps,
        "t1": t1,
        "t2": t2,
        "w_thresh": threshold,
        "typf": typf,
        "backend": backend,
        "transpose": transpose,
    }


def get_entropy(self: "SignedGraph") -> NDArray:
    if not hasattr(self, "entropy") or self.entropy is None:
        compute_signed_laplacian_entropy(self)
    return self.entropy


def get_specific_heat(self: "SignedGraph") -> NDArray:
    has_entropy = hasattr(self, "entropy_derivative") and self.entropy_derivative is not None
    if not has_entropy:
        compute_signed_laplacian_entropy(self)
    return self.entropy_derivative


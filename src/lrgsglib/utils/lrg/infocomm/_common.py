"""Common utilities shared across entropy computation modules."""

from numpy.typing import NDArray


def _normalize_entropy_profile(entropy_profile: NDArray, entropy_norm: str) -> NDArray:
    """Normalize entropy profile according to specified mode.

    Parameters
    ----------
    entropy_profile : NDArray
        Raw entropy profile.
    entropy_norm : str
        Normalization mode: "standard" or "complement".

    Returns
    -------
    NDArray
        Normalized entropy profile.

    Raises
    ------
    ValueError
        If entropy_norm is not "standard" or "complement".
    """
    norm = entropy_norm.lower()
    if norm == "standard":
        return entropy_profile
    if norm == "complement":
        return 1 - entropy_profile
    raise ValueError("entropy_norm must be 'standard' or 'complement'.")

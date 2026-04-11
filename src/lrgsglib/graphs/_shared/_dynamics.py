"""Engine-agnostic dynamics energy mixin functions.

Computes RBIM and spherical SK-model energies from eigenvectors.
Uses ``get_edges_with_weights()`` (protocol method on both NX and GT)
for edge access.  The ``on_g`` parameter is forwarded when the
concrete implementation accepts it (NX does, GT ignores it).
"""

import numpy as np
from ...utils.lrg import compute_ising_pairwise_energy


def _get_edges(self, on_g=None):
    """Get edges as list of (u, v, weight) using the protocol method."""
    if on_g is not None:
        try:
            return self.get_edges_with_weights(on_g=on_g)
        except TypeError:
            pass
    return self.get_edges_with_weights()


def compute_rbim_energy_eigV(
    self, which: int = 0, use_gpu: bool = False, on_g=None, **kwargs
):
    """Compute Random Bond Ising Model energy for the binarized eigenvector."""
    spins = self.get_eigV(which, binarize=True)
    edges = _get_edges(self, on_g)
    if not hasattr(self, "energy_eigV_RBIM"):
        self.energy_eigV_RBIM = {}
    self.energy_eigV_RBIM[which] = compute_ising_pairwise_energy(
        spins, edges, use_gpu=use_gpu
    )


def compute_sksph_energy_eigV(
    self, which: int = 0, use_gpu: bool = False, on_g=None, **kwargs
):
    """Compute spherical SK-model energy for the continuous eigenvector."""
    state = self.get_eigV(which, binarize=False)
    edges = _get_edges(self, on_g)
    if not hasattr(self, "energy_eigV_SKSPH"):
        self.energy_eigV_SKSPH = {}
    self.energy_eigV_SKSPH[which] = compute_ising_pairwise_energy(
        state, edges, use_gpu=use_gpu
    )


def compute_sksph_energy_eigV_all(self, on_g=None, **kwargs):
    """Compute SK-spherical energy for all eigenvectors."""
    if not hasattr(self, "eigV") or self.eigV is None:
        self.compute_laplacian_spectrum_weigV()
    if not hasattr(self, "energy_eigV_SKSPH"):
        self.energy_eigV_SKSPH = {}
    n_eigv = self.eigV.shape[1] if self.eigV.ndim == 2 else len(self.eigV)
    for which in range(n_eigv):
        if which not in self.energy_eigV_SKSPH:
            compute_sksph_energy_eigV(self, which, on_g=on_g)


def get_sksph_energy_eigV(self, which: int = 0):
    """Return cached SK-spherical energy for eigenvector ``which``."""
    has = hasattr(self, "energy_eigV_SKSPH") and which in getattr(
        self, "energy_eigV_SKSPH", {}
    )
    if not has:
        compute_sksph_energy_eigV(self, which)
    return self.energy_eigV_SKSPH[which]


def get_all_sksph_energy_eigV(self, as_dict: bool = False, on_g=None, **kwargs):
    """Return SK-spherical energies for all eigenvectors."""
    if not hasattr(self, "energy_eigV_SKSPH"):
        compute_sksph_energy_eigV_all(self, on_g=on_g, **kwargs)
    if as_dict:
        return self.energy_eigV_SKSPH
    return np.array(list(self.energy_eigV_SKSPH.values()))


def compute_rbim_energy_eigV_all(self, on_g=None, **kwargs):
    """Compute RBIM energy for every eigenvector currently stored in ``self.eigV``.

    Does NOT force a full-spectrum computation: it operates on whatever
    eigenvectors are already cached (populated lazily by
    ``compute_k_eigvV`` / ``get_eigV_check``). Callers that need more
    modes must request them first via ``compute_k_eigvV(k=...)``.
    """
    if not hasattr(self, "eigV") or self.eigV is None:
        return
    if not hasattr(self, "energy_eigV_RBIM"):
        self.energy_eigV_RBIM = {}
    if self.eigV.ndim == 2:
        transposed = getattr(self, "_eigV_is_transposed", False)
        n_eigv = self.eigV.shape[0] if transposed else self.eigV.shape[1]
    else:
        n_eigv = len(self.eigV)
    for which in range(n_eigv):
        if which not in self.energy_eigV_RBIM:
            compute_rbim_energy_eigV(self, which, on_g=on_g)


def get_rbim_energy_eigV(self, which: int = 0):
    """Return cached RBIM energy for eigenvector ``which``."""
    has = hasattr(self, "energy_eigV_RBIM") and which in getattr(
        self, "energy_eigV_RBIM", {}
    )
    if not has:
        compute_rbim_energy_eigV(self, which)
    return self.energy_eigV_RBIM[which]


def get_all_rbim_energy_eigV(self, as_dict: bool = False, on_g=None, **kwargs):
    """Return RBIM energies for all eigenvectors."""
    if not hasattr(self, "energy_eigV_RBIM"):
        compute_rbim_energy_eigV_all(self, on_g=on_g, **kwargs)
    if as_dict:
        return self.energy_eigV_RBIM
    return np.array(list(self.energy_eigV_RBIM.values()))

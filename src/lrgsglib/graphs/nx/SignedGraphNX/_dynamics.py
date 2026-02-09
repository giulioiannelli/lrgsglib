import numpy as np

from ....config.const import *
from ....utils.lrg import compute_ising_pairwise_energy
# Reuse spectral utilities
from ._spectral import (
    get_eigV_bin_check,
    get_eigV_check,
    compute_laplacian_spectrum_weigV,
)

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX


def compute_rbim_energy_eigV(self: "SignedGraphNX", which: int = 0, use_gpu: bool = False, on_g: str = SG_REPR):
    """Compute Random Bond Ising Model energy for the binarized eigenvector."""
    spins = get_eigV_bin_check(self, which)
    edges = list(self.gr[on_g].edges(data='weight'))
    if not hasattr(self, "energy_eigV_RBIM"):
        self.energy_eigV_RBIM = {}
    self.energy_eigV_RBIM[which] = compute_ising_pairwise_energy(spins, edges, use_gpu=use_gpu)


def compute_sksph_energy_eigV(self: "SignedGraphNX", which: int = 0, use_gpu: bool = False, on_g: str = SG_REPR):
    """Compute spherical SK-model energy for the continuous eigenvector."""
    state = get_eigV_check(self, which, binarize=False)
    edges = list(self.gr[on_g].edges(data='weight'))
    if not hasattr(self, "energy_eigV_SKSPH"):
        self.energy_eigV_SKSPH = {}
    self.energy_eigV_SKSPH[which] = compute_ising_pairwise_energy(state, edges, use_gpu=use_gpu)


def compute_sksph_energy_eigV_all(self: "SignedGraphNX", on_g: str = SG_REPR, **kw_laplspect):
    if not hasattr(self, "eigV"):
        compute_laplacian_spectrum_weigV(self, **kw_laplspect)
    if not hasattr(self, "energy_eigV_SKSPH"):
        self.energy_eigV_SKSPH = {}
    for which in range(len(self.eigV)):
        if which not in self.energy_eigV_SKSPH:
            compute_sksph_energy_eigV(self, which, on_g=on_g)


def get_sksph_energy_eigV(self: "SignedGraphNX", which: int = 0):
    """Return cached spherical SK energy for eigenvector `which`. Compute if missing."""
    has_energy = hasattr(self, "energy_eigV_SKSPH") and which in getattr(self, "energy_eigV_SKSPH", {})
    if not has_energy:
        compute_sksph_energy_eigV(self, which)
    return self.energy_eigV_SKSPH[which]


def get_all_sksph_energy_eigV(self: "SignedGraphNX", as_dict: bool = False, on_g: str = SG_REPR, **kw_laplspect):
    """Return spherical SK energies for all eigenvectors, as dict or array."""
    if not hasattr(self, "energy_eigV_SKSPH"):
        compute_sksph_energy_eigV_all(self, on_g, **kw_laplspect)
    if as_dict:
        return self.energy_eigV_SKSPH
    else:
        return np.array(list(self.energy_eigV_SKSPH.values()))


def compute_rbim_energy_eigV_all(self: "SignedGraphNX", on_g: str = SG_REPR, **kw_laplspect):
    if not hasattr(self, "eigV"):
        compute_laplacian_spectrum_weigV(self, **kw_laplspect)
    if not hasattr(self, "energy_eigV_RBIM"):
        self.energy_eigV_RBIM = {}
    for which in range(len(self.eigV)):
        if which not in self.energy_eigV_RBIM:
            compute_rbim_energy_eigV(self, which, on_g=on_g)


def get_rbim_energy_eigV(self: "SignedGraphNX", which: int = 0):
    has_rbim = hasattr(self, "energy_eigV_RBIM")
    if not has_rbim:
        compute_rbim_energy_eigV(self, which)
    elif which not in self.energy_eigV_RBIM.keys():
        compute_rbim_energy_eigV(self, which)
    return self.energy_eigV_RBIM[which]


def get_all_rbim_energy_eigV(self: "SignedGraphNX", as_dict: bool = False, on_g: str = SG_REPR, **kw_laplspect):
    if not hasattr(self, "energy_eigV_RBIM"):
        compute_rbim_energy_eigV_all(self, on_g, **kw_laplspect)
    if as_dict:
        return self.energy_eigV_RBIM
    else:
        return np.array(list(self.energy_eigV_RBIM.values()))

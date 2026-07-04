"""Spectral proposal / field construction shared by dynamics models.

Bodies extracted from the legacy Ising monolith's ``topo_*`` paths (Phase 2
of the Ising refactor), dtype-agnostic so the same helpers serve
``BinDynSys``- and ``VecDynSys``-based models.

Two capabilities live here (Ising plan §4):

- **eigenvector-subspace machinery** — the ``(M, N)`` basis of the *M*
  lowest signed-Laplacian eigenvectors, the RBIM energy of each binarized
  eigenvector, and the best-eigenvector seed. Consumed by the spectral
  Metropolis move (``move='spectral'``) and the CEM optimizer, which both
  search in the coefficient space ``s = sign(c · V)``.
- **spectral field construction** — the softmax-weighted eigenvector field
  ``h_i = strength · Σ_k c_k v_k[i]`` with ``c = softmax(−E_RBIM/τ)``
  (today's TFCA field, ``field='spectral'`` on the annealing scheme).

Everything consumes the graph through its spectral interface
(``compute_k_eigvV`` / ``get_eigV`` / ``compute_rbim_energy_eigV_all``); no
model internals leak in here.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ..graphs.protocols import DynamicsGraphProtocol

__all__ = [
    "SPECTRAL_SPARSE_N_MIN",
    "ensure_spectral_subspace",
    "eigvec_subspace_matrix",
    "rbim_subspace_energies",
    "best_eigenvector_seed",
    "softmax_mode_weights",
    "build_spectral_field",
    "greedy_polish_csr",
    "SpectralMetropolisEngine",
]

#: Above this many sites the subspace computation forces the sparse ``eigsh``
#: solver (when the graph backend exposes one) — dense diagonalisation is
#: O(N^3) and impractical there. Inherited from the legacy topo paths.
SPECTRAL_SPARSE_N_MIN: int = 500


def ensure_spectral_subspace(sg: "DynamicsGraphProtocol", n_modes: int) -> None:
    """Compute the partial Laplacian spectrum on *sg* if not already cached.

    Recomputes when the cached ``eigV`` holds fewer than *n_modes* rows —
    e.g. after an init populated only ``eigV[0]`` for a ground-state spin
    start. Forces the sparse solver (``backend='scipy'``) on large graphs
    where dense diagonalisation is impractical.
    """
    n_sites = int(sg.N)
    cached = getattr(sg, "eigV", None)
    if cached is not None:
        cached = np.asarray(cached)
        if cached.ndim == 2:
            transposed = getattr(sg, "_eigV_is_transposed", False)
            n_cached = cached.shape[0] if transposed else cached.shape[1]
        else:
            n_cached = len(cached)
        if n_cached >= n_modes:
            return
    kwargs: dict = {"k": n_modes}
    if n_sites > SPECTRAL_SPARSE_N_MIN and n_modes < n_sites // 2:
        import inspect

        sig = inspect.signature(sg.compute_k_eigvV)
        if "backend" in sig.parameters:
            kwargs["backend"] = "scipy"
    sg.compute_k_eigvV(**kwargs)


def eigvec_subspace_matrix(
    sg: "DynamicsGraphProtocol", n_modes: int
) -> NDArray[np.float64]:
    """Return the ``(n_modes, N)`` matrix of continuous eigenvectors.

    Row *k* is the *k*-th signed-Laplacian eigenvector (not binarized);
    ``sg.get_eigV`` handles the NX/GT storage layout. The subspace is
    computed on demand via :func:`ensure_spectral_subspace`.
    """
    ensure_spectral_subspace(sg, n_modes)
    vecs = [sg.get_eigV(which=k, binarize=False) for k in range(n_modes)]
    return np.vstack(vecs).astype(np.float64)


def rbim_subspace_energies(
    sg: "DynamicsGraphProtocol", n_modes: int
) -> NDArray[np.float64]:
    """Return the ``(n_modes,)`` RBIM energies of the binarized eigenvectors."""
    ensure_spectral_subspace(sg, n_modes)
    sg.compute_rbim_energy_eigV_all()
    # Lazily-created instance dict (not a static protocol member).
    rbim_by_mode = getattr(sg, "energy_eigV_RBIM")
    return np.array([rbim_by_mode[k] for k in range(n_modes)], dtype=np.float64)


def best_eigenvector_seed(
    sg: "DynamicsGraphProtocol", n_modes: int
) -> tuple[NDArray[np.int8], int]:
    """Return (binarized ±1 spins, mode index) of the lowest-RBIM-energy
    eigenvector — the standard seed of the spectral moves."""
    rbim_E = rbim_subspace_energies(sg, n_modes)
    best_k = int(np.argmin(rbim_E))
    spins = np.asarray(sg.get_eigV(which=best_k, binarize=True)).copy()
    # Ensure {-1, +1}: sign(0) → 0 must become +1.
    spins[spins == 0] = 1
    return spins.astype(np.int8), best_k


def softmax_mode_weights(
    rbim_energies: NDArray[np.float64], tau: float
) -> NDArray[np.float64]:
    """Softmax weights ``c_k ∝ exp(−E_k/(|E_min| τ))`` over the subspace.

    Energies are normalized by ``|E_min|`` first so the gaps are O(0.01–0.1)
    regardless of system size; lower RBIM energy → higher weight.
    """
    E_min_abs = float(np.abs(rbim_energies.min()))
    normed = rbim_energies / E_min_abs if E_min_abs > 0 else rbim_energies
    logits = -normed / tau
    logits -= logits.max()  # numerical stability
    weights = np.exp(logits)
    return weights / weights.sum()


def build_spectral_field(
    subspace_vecs: NDArray[np.float64],
    rbim_energies: NDArray[np.float64],
    tau: float,
    field_strength: float,
) -> NDArray[np.float64]:
    """Construct the spectral external field (the TFCA field).

    ``h_i = field_strength · Σ_k c_k v_k[i]`` with the softmax weights of
    :func:`softmax_mode_weights` — the low-RBIM-energy eigenvectors dominate,
    biasing the quench toward their sign structure.

    Parameters
    ----------
    subspace_vecs : ndarray, shape ``(M, N)``
        Row-major continuous eigenvectors.
    rbim_energies : ndarray, shape ``(M,)``
        RBIM energy per binarized eigenvector.
    tau : float
        Softmax temperature (lower → sharper weighting).
    field_strength : float
        Overall field magnitude scale.
    """
    weights = softmax_mode_weights(rbim_energies, tau)
    return field_strength * (weights @ subspace_vecs)  # (M,)@(M,N) → (N,)


def greedy_polish_csr(
    spins: NDArray,
    nbr_idx: NDArray,
    nbr_J: NDArray,
    nbr_ptr: NDArray,
    field: NDArray | None = None,
    max_sweeps: int = 50,
) -> NDArray[np.int8]:
    """Zero-temperature sequential quench: flip while any flip lowers energy.

    Sweeps sites in fixed order, flipping every site whose
    ``ΔE = 2 s_k (Σ_j J_kj s_j + h_k)`` is negative, until a sweep makes no
    flip (or *max_sweeps* is hit). Works on the model's symmetric-J CSR, so
    the polish minimizes the SAME Hamiltonian as ``delta_E``/``total_energy``
    (plan §3b consistency). Returns a new ±1 int8 array; the input is not
    modified.
    """
    s = np.asarray(spins, dtype=np.float64).copy()
    n_sites = len(s)
    h = None if field is None else np.asarray(field, dtype=np.float64)
    for _ in range(max_sweeps):
        any_flip = False
        for node in range(n_sites):
            lo, hi = int(nbr_ptr[node]), int(nbr_ptr[node + 1])
            h_eff = float(nbr_J[lo:hi] @ s[nbr_idx[lo:hi]])
            if h is not None:
                h_eff += float(h[node])
            if 2.0 * s[node] * h_eff < 0.0:  # flipping lowers energy
                s[node] = -s[node]
                any_flip = True
        if not any_flip:
            break
    return s.astype(np.int8)


class SpectralMetropolisEngine:
    """Metropolis in the eigenvector-coefficient space (``move='spectral'``).

    The legacy "topological Metropolis" as a configure-then-run engine with
    the same ``compile_step() -> sweep() -> ΔE`` interface as the site
    engines in ``statsys._thermal``: the walk lives on a continuous field
    ``f = Σ_k c_k v_k`` in the span of the *M* lowest signed-Laplacian
    eigenvectors; a proposal perturbs one softmax-selected mode's
    coefficient, the trial state is ``sign(f)`` (zeros keep their old
    sign), and the FULL §3b Hamiltonian of the trial state decides the M1
    acceptance. One sweep = N proposals. Per-mode proposal widths σ_m adapt
    toward a ~30% acceptance every *chunk_size* sweeps, optionally followed
    by a zero-temperature greedy polish that tracks the best-so-far state.

    The engine consumes the model through a narrow surface: ``s`` (±1 int8
    state, reassigned on accept), ``total_energy(s)``, and the CSR
    couplings + field passed explicitly for the polish.
    """

    def __init__(
        self,
        model: Any,
        *,
        T: float,
        subspace_vecs: NDArray[np.float64],
        rbim_energies: NDArray[np.float64],
        seed_spins: NDArray[np.int8],
        seed_mode: int,
        sigma_init: float,
        chunk_size: int,
        polish: bool,
        polish_sweeps: int,
        nbr_idx: NDArray,
        nbr_J: NDArray,
        nbr_ptr: NDArray,
        field: NDArray | None = None,
        rng: Any = None,
    ) -> None:
        if T < 0.0:
            raise ValueError(f"Temperature must be >= 0; got {T}.")
        self.model = model
        self.T = float(T)
        self.V = np.asarray(subspace_vecs, dtype=np.float64)
        n_modes = self.V.shape[0]
        self.coeffs = np.zeros(n_modes, dtype=np.float64)
        self.coeffs[seed_mode] = 1.0
        self._field_walk = self.V[seed_mode].copy()
        self.sigma = np.full(n_modes, float(sigma_init), dtype=np.float64)
        self.chunk_size = int(chunk_size)
        self.polish = bool(polish)
        self.polish_sweeps = int(polish_sweeps)
        self._csr = (nbr_idx, nbr_J, nbr_ptr)
        self._h = field
        self.rng = rng if rng is not None else np.random

        # Mode-selection weights ∝ |RBIM energy| (legacy semantics).
        w = np.abs(np.asarray(rbim_energies, dtype=np.float64)).copy()
        if w.sum() == 0.0:
            w[:] = 1.0
        self.mode_weights = w / w.sum()

        # Seed state: the engine OWNS the walk state and mirrors it into
        # the model on construction (the spectral move starts from the
        # best eigenvector, not from the model's ic).
        model.s = np.asarray(seed_spins, dtype=np.int8).copy()
        self._E = float(model.total_energy())
        self.best_spins = model.s.copy()
        self.best_energy = self._E
        self._accept = np.zeros(n_modes, dtype=np.int64)
        self._total = np.zeros(n_modes, dtype=np.int64)
        self._sweeps_done = 0

    def _maybe_adapt(self) -> None:
        """σ adaptation + optional greedy polish, every chunk_size sweeps."""
        if self._sweeps_done % self.chunk_size != 0:
            return
        mask = self._total > 0
        rate = np.zeros_like(self.sigma)
        rate[mask] = self._accept[mask] / self._total[mask]
        self.sigma[mask & (rate > 0.35)] *= 1.1
        self.sigma[mask & (rate < 0.25)] *= 0.9
        self._accept[:] = 0
        self._total[:] = 0
        if self.polish:
            idx, J, ptr = self._csr
            polished = greedy_polish_csr(
                self.model.s,
                idx,
                J,
                ptr,
                field=self._h,
                max_sweeps=self.polish_sweeps,
            )
            E_pol = float(self.model.total_energy(polished))
            if E_pol < self.best_energy:
                self.best_energy = E_pol
                self.best_spins = polished.copy()

    def compile_step(self):
        model = self.model
        rng = self.rng
        V = self.V
        n_modes = V.shape[0]
        mode_w = self.mode_weights
        T = self.T
        n_proposals = int(len(model.s))

        def sweep() -> float:
            E_start = self._E
            spins = model.s.astype(np.float64, copy=False)
            for _ in range(n_proposals):
                m = int(rng.choice(n_modes, p=mode_w))
                delta = float(rng.normal(0.0, self.sigma[m]))
                field_new = self._field_walk + delta * V[m]
                s_new = np.sign(field_new)
                zeros = s_new == 0.0
                s_new[zeros] = spins[zeros]
                E_new = float(model.total_energy(s_new))
                dE = E_new - self._E
                if dE <= 0.0 or (T > 0.0 and rng.random() < math.exp(-dE / T)):
                    self.coeffs[m] += delta
                    self._field_walk = field_new
                    model.s = s_new.astype(np.int8)
                    spins = s_new
                    self._E = E_new
                    self._accept[m] += 1
                self._total[m] += 1
            if self._E < self.best_energy:
                self.best_energy = self._E
                self.best_spins = model.s.copy()
            self._sweeps_done += 1
            self._maybe_adapt()
            return self._E - E_start

        return sweep

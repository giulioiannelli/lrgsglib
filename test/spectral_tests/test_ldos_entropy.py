"""
Unit tests for LDOS-based entropy framework.

Tests compute_ldos_entropy and compute_ldos_specific_heat from
lrgsglib.utils.lrg.quantum, verifying:
  - Output shapes and value ranges
  - Reduction to spectral entropy at tau -> 0
  - Reduction to classical LRG entropy for global S(tau)
  - Correct behavior on unsigned vs signed graphs
  - Specific heat shape and peak detection
"""

import numpy as np
import pytest
from lrgsglib.graphs.nx import Lattice2DNX
from lrgsglib.utils.lrg.quantum import (
    compute_ldos_entropy,
    compute_ldos_specific_heat,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_lattice(side: int, pflip: float = 0.0) -> Lattice2DNX:
    """Build square lattice, compute full spectrum."""
    np.random.seed(42)
    lat = Lattice2DNX(side1=side, geo="sqr", pflip=pflip)
    if pflip > 0:
        lat.flip_random_fract_edges()
    lat.compute_laplacian_spectrum_weigV()
    return lat


def _eigV_col(lat: Lattice2DNX) -> np.ndarray:
    """Column-major eigenvectors: shape (N, N), col k = k-th eigvec."""
    return lat.eigV.T


def _two_community_bridge_graph():
    """Create a two-clique graph joined by a single bridge edge.

    Returns (eigenvalues, eigvecs_col, bridge_nodes, non_bridge_nodes).
    bridge_nodes are the two endpoints of the bridge.
    """
    import networkx as nx

    G = nx.Graph()
    # Two cliques of size 6
    clique_a = list(range(6))
    clique_b = list(range(6, 12))
    for nodes in (clique_a, clique_b):
        for i in nodes:
            for j in nodes:
                if i < j:
                    G.add_edge(i, j, sign=1, weight=1.0)
    # Bridge: node 5 -- node 6
    G.add_edge(5, 6, sign=1, weight=1.0)

    N = G.number_of_nodes()
    L = nx.laplacian_matrix(G).toarray().astype(float)
    eigvals, eigvecs = np.linalg.eigh(L)  # eigvecs column-major

    bridge_nodes = np.array([5, 6])
    non_bridge = np.array([n for n in range(N) if n not in (5, 6)])
    return eigvals, eigvecs, bridge_nodes, non_bridge


# ---------------------------------------------------------------------------
# TestLdosEntropyBasic
# ---------------------------------------------------------------------------

class TestLdosEntropyBasic:
    """Basic shape and range checks for compute_ldos_entropy."""

    def test_shape(self):
        """S_local has shape (N,), S_global is scalar."""
        lat = _build_lattice(8)
        S_local, S_global = compute_ldos_entropy(
            lat.eigv, _eigV_col(lat), tau=1.0
        )
        assert S_local.shape == (lat.N,)
        assert np.isscalar(S_global) or S_global.ndim == 0

    def test_normalized_range(self):
        """0 <= S_local <= 1 and 0 <= S_global <= 1."""
        lat = _build_lattice(8)
        S_local, S_global = compute_ldos_entropy(
            lat.eigv, _eigV_col(lat), tau=1.0
        )
        assert np.all(S_local >= -1e-12), f"min S_local = {S_local.min()}"
        assert np.all(S_local <= 1 + 1e-12), f"max S_local = {S_local.max()}"
        assert S_global >= -1e-12
        assert S_global <= 1 + 1e-12

    def test_identity_eigenvectors(self):
        """For V=I, each node sees only one mode -> S_local should be low."""
        N = 10
        eigvals = np.arange(N, dtype=float)
        eigvecs = np.eye(N)  # already column-major
        S_local, S_global = compute_ldos_entropy(eigvals, eigvecs, tau=1.0)

        # With identity eigenvectors, each node has LDOS concentrated on
        # one mode, so S_local should be near 0
        assert np.all(S_local < 0.1), (
            f"Expected low S_local for identity eigvecs, got {S_local}"
        )

    def test_uniform_eigenvectors(self):
        """For uniform |v_k(i)|^2 = 1/N, S_local approx S_global."""
        N = 16
        eigvals = np.arange(N, dtype=float)
        # Orthogonal matrix with uniform absolute values: Hadamard-like
        # Use DFT matrix (unitary with |entry|^2 = 1/N)
        k = np.arange(N)
        V = np.exp(2j * np.pi * np.outer(k, k) / N) / np.sqrt(N)
        # Take real part of absolute values squared -- they are all 1/N
        # Actually use the real orthogonal: V[i,k] = sqrt(2/N)*cos(...)
        # Simplest: just use columns of |V|^2 all equal
        # Build: random orthogonal with |entry|^2 ~ 1/N (Haar random, large N)
        # For exact test: V = F/sqrt(N) DFT: |V[i,k]|^2 = 1/N exactly
        V_col = V.real  # won't be orthogonal, use the DFT directly
        # DFT V: |V[i,k]|^2 = 1/N for all i,k, but V is complex
        # compute_ldos_entropy squares: V2 = eigvecs**2, but for complex
        # the code uses ** 2 not |.|^2. Let's use a real construction.
        # Uniform real eigvecs: all entries +/-1/sqrt(N). Use Hadamard.
        from scipy.linalg import hadamard
        H = hadamard(N) / np.sqrt(N)  # |H[i,k]|^2 = 1/N
        S_local, S_global = compute_ldos_entropy(eigvals, H, tau=1.0)

        # With uniform |v_k(i)|^2, all nodes identical -> S_local = S_global
        assert np.allclose(S_local, S_global, atol=1e-10), (
            f"S_local should equal S_global for uniform eigvecs.\n"
            f"S_local range: [{S_local.min():.6f}, {S_local.max():.6f}], "
            f"S_global = {S_global:.6f}"
        )


# ---------------------------------------------------------------------------
# TestLdosEntropyReductions
# ---------------------------------------------------------------------------

class TestLdosEntropyReductions:
    """Test limiting-case reductions of the LDOS entropy."""

    def test_tau_zero_is_spectral_entropy(self):
        """At very small tau, S(tau,i) approx S_spec(i).

        S_spec(i) = -sum_k |v_k(i)|^2 log|v_k(i)|^2 / log(N)
        """
        lat = _build_lattice(10)
        eigV_col = _eigV_col(lat)
        N = lat.N

        S_local, _ = compute_ldos_entropy(lat.eigv, eigV_col, tau=0.001)

        # Compute spectral entropy directly
        V2 = eigV_col ** 2  # shape (N, N)
        with np.errstate(divide="ignore", invalid="ignore"):
            log_V2 = np.where(V2 > 1e-30, np.log(V2), 0.0)
        S_spec = -np.sum(V2 * log_V2, axis=1) / np.log(N)

        # Should match closely
        max_dev = np.max(np.abs(S_local - S_spec))
        assert max_dev < 0.05, (
            f"tau->0 should recover spectral entropy. Max deviation = {max_dev}"
        )

    def test_global_equals_classical(self):
        """S_global should equal classical -sum_k p_k log p_k / log(N).

        where p_k = exp(-tau|lambda_k|) / Z.
        """
        lat = _build_lattice(10)
        eigV_col = _eigV_col(lat)
        N = lat.N
        tau = 1.0

        _, S_global = compute_ldos_entropy(lat.eigv, eigV_col, tau=tau)

        # Classical computation
        lam = np.abs(lat.eigv)
        boltz = np.exp(-tau * lam)
        Z = np.sum(boltz)
        p = boltz / Z
        mask = p > 1e-30
        S_classical = -np.sum(p[mask] * np.log(p[mask])) / np.log(N)

        assert np.isclose(S_global, S_classical, atol=1e-12), (
            f"S_global = {S_global}, S_classical = {S_classical}"
        )

    def test_unsigned_reduction(self):
        """For a graph with all lambda >= 0, use_abs=True gives same as False."""
        lat = _build_lattice(8, pflip=0.0)
        eigV_col = _eigV_col(lat)

        # All eigenvalues should be >= 0 for unsigned graph
        assert np.all(lat.eigv >= -1e-10), "Unsigned graph should have non-negative eigvals"

        S_local_abs, S_global_abs = compute_ldos_entropy(
            lat.eigv, eigV_col, tau=1.0, use_abs=True
        )
        S_local_raw, S_global_raw = compute_ldos_entropy(
            lat.eigv, eigV_col, tau=1.0, use_abs=False
        )

        assert np.allclose(S_local_abs, S_local_raw, atol=1e-12)
        assert np.isclose(S_global_abs, S_global_raw, atol=1e-12)


# ---------------------------------------------------------------------------
# TestLdosEntropyOnLattice
# ---------------------------------------------------------------------------

class TestLdosEntropyOnLattice:
    """Physics-motivated tests on actual graph structures."""

    def test_unsigned_lattice_low_variance(self):
        """On pflip=0 lattice, std(S_local) should be small.

        Translation invariance implies near-uniform eigenvector amplitudes.
        """
        lat = _build_lattice(12, pflip=0.0)
        S_local, _ = compute_ldos_entropy(lat.eigv, _eigV_col(lat), tau=1.0)
        sigma = np.std(S_local)
        assert sigma < 0.05, (
            f"Unsigned lattice should have low S_local variance, got std = {sigma}"
        )

    def test_signed_lattice_higher_variance(self):
        """On pflip=0.15, std(S_local) > std for pflip=0."""
        lat_clean = _build_lattice(12, pflip=0.0)
        lat_dirty = _build_lattice(12, pflip=0.15)

        S_clean, _ = compute_ldos_entropy(
            lat_clean.eigv, _eigV_col(lat_clean), tau=1.0
        )
        S_dirty, _ = compute_ldos_entropy(
            lat_dirty.eigv, _eigV_col(lat_dirty), tau=1.0
        )

        sigma_clean = np.std(S_clean)
        sigma_dirty = np.std(S_dirty)

        assert sigma_dirty > sigma_clean, (
            f"Signed lattice should have higher variance. "
            f"Clean: {sigma_clean:.6f}, Dirty: {sigma_dirty:.6f}"
        )

    def test_bridge_nodes_lower_entropy(self):
        """On a two-community graph with bridge, bridge nodes have lower S_local.

        Bridge nodes participate in fewer large-scale modes, leading to
        more concentrated LDOS and lower entropy.
        """
        eigvals, eigvecs, bridge, non_bridge = _two_community_bridge_graph()

        # Use moderate tau to probe mesoscale structure
        S_local, _ = compute_ldos_entropy(eigvals, eigvecs, tau=1.0)

        mean_bridge = np.mean(S_local[bridge])
        mean_interior = np.mean(S_local[non_bridge])

        assert mean_bridge < mean_interior, (
            f"Bridge nodes should have lower S_local. "
            f"Bridge: {mean_bridge:.4f}, Interior: {mean_interior:.4f}"
        )


# ---------------------------------------------------------------------------
# TestLdosSpecificHeat
# ---------------------------------------------------------------------------

class TestLdosSpecificHeat:
    """Tests for compute_ldos_specific_heat."""

    def test_shape(self):
        """C_local has shape (n_tau, N), C_global has shape (n_tau,)."""
        lat = _build_lattice(8)
        tau_grid = np.logspace(-2, 2, 50)

        S_loc, S_glob, C_loc, C_glob = compute_ldos_specific_heat(
            lat.eigv, _eigV_col(lat), tau_grid
        )

        n_tau = len(tau_grid)
        assert S_loc.shape == (n_tau, lat.N)
        assert S_glob.shape == (n_tau,)
        assert C_loc.shape == (n_tau, lat.N)
        assert C_glob.shape == (n_tau,)

    def test_global_has_peaks(self):
        """On an unsigned lattice, C_global should have at least one peak."""
        lat = _build_lattice(12, pflip=0.0)
        tau_grid = np.logspace(-2, 3, 200)

        _, _, _, C_glob = compute_ldos_specific_heat(
            lat.eigv, _eigV_col(lat), tau_grid
        )

        # Find local maxima (value larger than both neighbors)
        peaks = []
        for i in range(1, len(C_glob) - 1):
            if C_glob[i] > C_glob[i - 1] and C_glob[i] > C_glob[i + 1]:
                peaks.append(i)

        assert len(peaks) >= 1, (
            f"Expected at least one peak in C_global, found {len(peaks)}. "
            f"C_global range: [{C_glob.min():.4f}, {C_glob.max():.4f}]"
        )

    def test_peaks_match_classical(self):
        """C_global peaks should be at the same tau as classical C(tau)."""
        from lrgsglib.utils.lrg.infocomm import (
            compute_entropy_observables_from_eigenvalues,
        )

        lat = _build_lattice(12, pflip=0.0)
        tau_grid = np.logspace(-2, 3, 200)

        # LDOS-based
        _, _, _, C_glob = compute_ldos_specific_heat(
            lat.eigv, _eigV_col(lat), tau_grid
        )

        # Classical
        _, C_classical, _, tau_classical = (
            compute_entropy_observables_from_eigenvalues(
                eigenvalues=lat.eigv,
                num_nodes=lat.N,
                steps=200,
                t1=-2,
                t2=3,
                specific_heat_scale="none",
            )
        )

        # Find peak indices
        def _find_peak(arr):
            return np.argmax(arr)

        peak_ldos = _find_peak(C_glob)
        peak_classical = _find_peak(C_classical)

        # Peak positions should be within a few grid points
        tau_peak_ldos = tau_grid[peak_ldos]
        tau_peak_classical = tau_classical[peak_classical]

        # Allow factor-of-2 tolerance in tau (they use the same grid)
        ratio = tau_peak_ldos / tau_peak_classical
        assert 0.5 < ratio < 2.0, (
            f"Peak mismatch: LDOS tau={tau_peak_ldos:.3f}, "
            f"classical tau={tau_peak_classical:.3f}, ratio={ratio:.3f}"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

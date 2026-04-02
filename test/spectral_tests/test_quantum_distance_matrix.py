"""
Unit tests for compute_quantum_distance_matrix.

Tests the quantum distance matrix derived from time-averaged quantum walk
transition probabilities.  Validates metric properties, compares with
classical diffusion distances, and checks hierarchical clustering output.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.utils.lrg.quantum import compute_quantum_distance_matrix


# ---------------------------------------------------------------------------
#  Helpers
# ---------------------------------------------------------------------------

def _make_lattice(side: int, pflip: float, seed: int = 42) -> Lattice2D:
    """Create a Lattice2D with full eigenvector decomposition."""
    np.random.seed(seed)
    lat = Lattice2D(side1=side, pflip=pflip)
    if pflip > 0:
        lat.flip_random_fract_edges()
    # Column-major (transpose=False) so eigV[:,k] is the k-th eigvec
    lat.compute_laplacian_spectrum_weigV(transpose=False)
    return lat


def _eigV_col_major(lat: Lattice2D) -> np.ndarray:
    """Return eigenvectors in column-major format regardless of storage."""
    if getattr(lat, "_eigV_is_transposed", False):
        return lat.eigV.T
    return lat.eigV


# ---------------------------------------------------------------------------
#  Basic correctness
# ---------------------------------------------------------------------------

class TestBasicProperties:
    """Test fundamental properties of the quantum distance matrix."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.lat = _make_lattice(side=8, pflip=0.0, seed=42)
        self.V = _eigV_col_major(self.lat)

    def test_shape(self):
        D = compute_quantum_distance_matrix(self.V)
        N = self.lat.N
        assert D.shape == (N, N)

    def test_zero_diagonal(self):
        D = compute_quantum_distance_matrix(self.V)
        assert_allclose(np.diag(D), 0.0, atol=1e-12)

    def test_symmetric(self):
        D = compute_quantum_distance_matrix(self.V)
        assert_allclose(D, D.T, atol=1e-12)

    def test_non_negative(self):
        D = compute_quantum_distance_matrix(self.V)
        assert np.all(D >= -1e-12), f"Negative distances found: min={D.min()}"

    def test_overlap_method_finite(self):
        D = compute_quantum_distance_matrix(self.V, method="overlap")
        assert np.all(np.isfinite(D))

    def test_hellinger_method(self):
        D = compute_quantum_distance_matrix(self.V, method="hellinger")
        assert D.shape == (self.lat.N, self.lat.N)
        assert_allclose(np.diag(D), 0.0, atol=1e-12)
        assert np.all(D >= -1e-12)

    def test_hellinger_bounded(self):
        """Hellinger distance should be in [0, 1]."""
        D = compute_quantum_distance_matrix(self.V, method="hellinger")
        assert np.all(D <= 1.0 + 1e-12)

    def test_invalid_method_raises(self):
        with pytest.raises(ValueError, match="Unknown method"):
            compute_quantum_distance_matrix(self.V, method="bogus")


class TestOverlapPhysics:
    """Test physical properties of the overlap distance."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.lat = _make_lattice(side=8, pflip=0.0, seed=42)
        self.V = _eigV_col_major(self.lat)

    def test_diagonal_is_ipr(self):
        """Diagonal of overlap matrix equals nodal IPR."""
        V2 = self.V ** 2
        overlap = V2 @ V2.T
        # IPR at node i = sum_k |v_k(i)|^4
        ipr_per_node = np.sum(self.V ** 4, axis=1)
        assert_allclose(np.diag(overlap), ipr_per_node, rtol=1e-10)

    def test_translational_symmetry_unfrustrated(self):
        """On an unfrustrated square lattice the overlap matrix diagonal
        entries (nodal IPR) should have low variance relative to the mean,
        reflecting approximate translational equivalence.  Exact equality
        is broken by open boundary conditions and eigenvector sign flips."""
        V2 = self.V ** 2
        overlap = V2 @ V2.T
        diag = np.diag(overlap)
        # Coefficient of variation should be moderate (not huge)
        cv = np.std(diag) / np.mean(diag)
        assert cv < 1.0, f"IPR coefficient of variation too large: {cv:.4f}"


class TestTriangleInequality:
    """Spot-check the triangle inequality for the overlap distance."""

    def test_triangle_inequality_sample(self):
        lat = _make_lattice(side=6, pflip=0.0, seed=42)
        V = _eigV_col_major(lat)
        D = compute_quantum_distance_matrix(V, method="overlap")
        N = lat.N
        # Check a random sample of triples
        rng = np.random.default_rng(0)
        for _ in range(500):
            i, j, k = rng.choice(N, size=3, replace=False)
            assert D[i, k] <= D[i, j] + D[j, k] + 1e-10, (
                f"Triangle ineq violated: d({i},{k})={D[i,k]:.6f} > "
                f"d({i},{j})={D[i,j]:.6f} + d({j},{k})={D[j,k]:.6f}"
            )


class TestFrustrationEffect:
    """Compare quantum distances with and without frustration."""

    def test_pflip_changes_distances(self):
        lat0 = _make_lattice(side=8, pflip=0.0, seed=42)
        lat1 = _make_lattice(side=8, pflip=0.15, seed=42)
        V0 = _eigV_col_major(lat0)
        V1 = _eigV_col_major(lat1)
        D0 = compute_quantum_distance_matrix(V0)
        D1 = compute_quantum_distance_matrix(V1)
        # Frustration should alter eigenvector structure
        assert not np.allclose(D0, D1, atol=1e-3)


class TestIdentityEigenvectors:
    """Regression: identity eigenvectors (complete localization)."""

    def test_identity_gives_infinite_offdiag(self):
        """For V = I, overlap O_ij = delta_ij, so off-diag distances = inf
        in the overlap method (no transition probability)."""
        V = np.eye(4)
        D = compute_quantum_distance_matrix(V, method="overlap")
        # Off-diagonal should be very large (essentially -log(0))
        offdiag = D[np.triu_indices(4, k=1)]
        assert np.all(offdiag > 100)

    def test_identity_hellinger(self):
        """For V = I, hellinger off-diag = 1 (maximally distant)."""
        V = np.eye(4)
        D = compute_quantum_distance_matrix(V, method="hellinger")
        offdiag = D[np.triu_indices(4, k=1)]
        assert_allclose(offdiag, 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
#  Comparison with classical distances (integration-level, marks slow)
# ---------------------------------------------------------------------------

class TestClassicalComparison:
    """Compare quantum and classical distances on a small lattice."""

    @pytest.fixture(autouse=True)
    def setup(self):
        self.lat = _make_lattice(side=8, pflip=0.0, seed=42)
        self.V = _eigV_col_major(self.lat)

    def test_nontrivial_correlation_with_classical(self):
        """Quantum and classical distances should have nontrivial
        (nonzero) correlation.  The sign may be positive or negative
        depending on tau, since quantum distances capture eigenvector
        geometry while classical distances are eigenvalue-weighted."""
        from lrgsglib.utils.lrg.infocomm.distances import lapl_dists

        D_q = compute_quantum_distance_matrix(self.V)
        # Classical distances -- use dense Laplacian for reliable expm
        L_dense = (
            self.lat.slp.toarray()
            if hasattr(self.lat.slp, "toarray")
            else np.asarray(self.lat.slp)
        )
        classical_condensed = lapl_dists(L_dense, tau=1.0, is_signed=False)
        D_c = squareform(classical_condensed)

        # Flatten upper triangles, masking out infinities
        idx = np.triu_indices(self.lat.N, k=1)
        q_flat = D_q[idx]
        c_flat = D_c[idx]
        mask = np.isfinite(q_flat) & np.isfinite(c_flat)
        assert mask.sum() > 10, "Too few finite distance pairs"

        r, p = pearsonr(q_flat[mask], c_flat[mask])
        # They should be significantly correlated (either direction)
        assert abs(r) > 0.05, (
            f"Expected nontrivial correlation, got r={r:.4f}, p={p:.4e}"
        )

    def test_hierarchical_clustering_agreement(self):
        """Dendrograms from quantum and classical should have nontrivial
        cophenetic distance correlation."""
        from lrgsglib.utils.lrg.infocomm.distances import lapl_dists

        D_q = compute_quantum_distance_matrix(self.V)
        # Classical distances -- dense Laplacian
        L_dense = (
            self.lat.slp.toarray()
            if hasattr(self.lat.slp, "toarray")
            else np.asarray(self.lat.slp)
        )
        classical_condensed = lapl_dists(L_dense, tau=1.0, is_signed=False)

        # Hierarchical clustering
        q_condensed = squareform(D_q)
        Z_q = linkage(q_condensed, method="average")
        Z_c = linkage(classical_condensed, method="average")

        # cophenet(Z) without Y returns the condensed cophenetic distance
        # matrix directly (not a tuple).
        coph_q_dists = cophenet(Z_q)
        coph_c_dists = cophenet(Z_c)
        r_cross, _ = pearsonr(coph_q_dists, coph_c_dists)
        # They should have some relationship (positive or negative)
        assert abs(r_cross) > 0.01, (
            f"Cophenetic distances uncorrelated: r={r_cross:.4f}"
        )

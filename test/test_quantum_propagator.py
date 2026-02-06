"""
Unit tests for quantum propagator functions.

Tests quantum propagator unitarity, probability conservation, entropy bounds,
and other quantum mechanical properties.
"""

import numpy as np
import pytest
from lrgsglib.nx_patches.Lattice2D import Lattice2D
from lrgsglib.nx_patches.ErdosRenyi import ErdosRenyi
from lrgsglib.utils.lrg.quantum import (
    compute_quantum_propagator_spectral,
    compute_quantum_propagator_matrix,
    quantum_density_matrix_evolution,
    quantum_probability_distribution,
    compute_quantum_coherence,
    von_neumann_entropy,
    quantum_classical_divergence,
    compute_interference_visibility,
    compute_quantum_observables_from_eigenvalues,
)


class TestQuantumPropagatorSpectral:
    """Test spectral method quantum propagator computation."""

    def test_unitarity(self):
        """Test U(t)^† U(t) = I (unitarity)."""
        l = Lattice2D(side1=10, pflip=0.1, seed=42)
        l.flip_random_fract_edges()
        l.compute_laplacian_spectrum_weigV()

        U_t = compute_quantum_propagator_spectral(
            l.eigv, l.eigV, t=1.0, full_matrix=True
        )

        # Check unitarity
        identity = U_t @ U_t.conj().T
        assert np.allclose(identity, np.eye(l.N), atol=1e-10)

    def test_return_types(self):
        """Test return types for full_matrix flag."""
        l = Lattice2D(side1=8, seed=42)
        l.compute_laplacian_spectrum_weigV()

        # Test with full_matrix=False
        phases, eigvecs = compute_quantum_propagator_spectral(
            l.eigv, l.eigV, t=1.0, full_matrix=False
        )
        assert phases.shape == (l.N,)
        assert phases.dtype == np.complex128
        assert eigvecs.shape == (l.N, l.N)

        # Test with full_matrix=True
        U_t = compute_quantum_propagator_spectral(
            l.eigv, l.eigV, t=1.0, full_matrix=True
        )
        assert U_t.shape == (l.N, l.N)
        assert U_t.dtype == np.complex128

    def test_time_reversibility(self):
        """Test U(-t) = U(t)^†."""
        l = Lattice2D(side1=8, seed=42)
        l.compute_laplacian_spectrum_weigV()

        U_t = compute_quantum_propagator_spectral(
            l.eigv, l.eigV, t=1.0, full_matrix=True
        )
        U_minus_t = compute_quantum_propagator_spectral(
            l.eigv, l.eigV, t=-1.0, full_matrix=True
        )

        assert np.allclose(U_minus_t, U_t.conj().T, atol=1e-10)


class TestQuantumProbabilityDistribution:
    """Test quantum walk probability distribution."""

    def test_probability_conservation(self):
        """Test ∑_j P_j(t) = 1."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        P_t = quantum_probability_distribution(
            l.eigv, l.eigV, t=1.0, init_node=50
        )

        assert np.isclose(np.sum(P_t), 1.0, atol=1e-10)

    def test_probability_positive(self):
        """Test P_j(t) >= 0 for all j."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        P_t = quantum_probability_distribution(
            l.eigv, l.eigV, t=1.0, init_node=50
        )

        assert np.all(P_t >= 0)

    def test_initial_condition(self):
        """Test P_j(0) = δ_ij."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        init_node = 50
        P_0 = quantum_probability_distribution(
            l.eigv, l.eigV, t=0.0, init_node=init_node
        )

        expected = np.zeros(l.N)
        expected[init_node] = 1.0

        assert np.allclose(P_0, expected, atol=1e-10)


class TestQuantumDensityMatrix:
    """Test quantum density matrix evolution."""

    def test_hermiticity(self):
        """Test ρ = ρ† (Hermiticity)."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        rho_t = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=1.0, init_type="localized", init_node=50
        )

        assert np.allclose(rho_t, rho_t.conj().T, atol=1e-10)

    def test_trace_preservation(self):
        """Test Tr(ρ(t)) = 1."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        rho_t = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=1.0, init_type="localized", init_node=50
        )

        assert np.isclose(np.trace(rho_t), 1.0, atol=1e-10)

    def test_positive_semidefinite(self):
        """Test ρ is positive semidefinite (eigenvalues >= 0)."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        rho_t = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=1.0, init_type="localized", init_node=50
        )

        eigenvals = np.linalg.eigvalsh(rho_t)
        assert np.all(eigenvals >= -1e-10)  # Allow small numerical errors

    def test_init_types(self):
        """Test different initial state types."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        # Localized
        rho_localized = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=0.0, init_type="localized", init_node=50
        )
        assert np.isclose(rho_localized[50, 50], 1.0, atol=1e-10)

        # Uniform
        rho_uniform = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=0.0, init_type="uniform"
        )
        assert np.allclose(rho_uniform, np.eye(l.N) / l.N, atol=1e-10)

        # Thermal
        rho_thermal = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=0.0, init_type="thermal"
        )
        assert np.isclose(np.trace(rho_thermal), 1.0, atol=1e-10)


class TestVonNeumannEntropy:
    """Test von Neumann entropy computation."""

    def test_pure_state_entropy(self):
        """Test S_vN = 0 for pure states."""
        rho_pure = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
        entropy = von_neumann_entropy(rho_pure)
        assert np.isclose(entropy, 0.0, atol=1e-10)

    def test_maximally_mixed_entropy(self):
        """Test S_vN = log(N) for maximally mixed state."""
        N = 10
        rho_mixed = np.eye(N, dtype=complex) / N
        entropy = von_neumann_entropy(rho_mixed)
        assert np.isclose(entropy, np.log(N), atol=1e-10)

    def test_entropy_bounds(self):
        """Test 0 <= S_vN <= log(N)."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        rho_t = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=1.0, init_type="localized", init_node=50
        )

        entropy = von_neumann_entropy(rho_t)
        assert 0 <= entropy <= np.log(l.N) + 1e-10


class TestQuantumCoherence:
    """Test quantum coherence computation."""

    def test_classical_state_zero_coherence(self):
        """Test C = 0 for diagonal (classical) states."""
        rho_classical = np.diag([0.5, 0.5])
        coherence = compute_quantum_coherence(rho_classical)
        assert np.isclose(coherence, 0.0, atol=1e-10)

    def test_superposition_coherence(self):
        """Test C > 0 for quantum superpositions."""
        rho_superposition = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        coherence = compute_quantum_coherence(rho_superposition)
        assert coherence > 0

    def test_coherence_positive(self):
        """Test C >= 0 always."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        rho_t = quantum_density_matrix_evolution(
            l.eigv, l.eigV, t=1.0, init_type="localized", init_node=50
        )

        coherence = compute_quantum_coherence(rho_t)
        assert coherence >= 0


class TestQuantumClassicalDivergence:
    """Test quantum-classical divergence computation."""

    def test_divergence_positive(self):
        """Test D >= 0."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        div = quantum_classical_divergence(
            l.eigv, l.eigV, t=1.0, init_node=50
        )

        assert div >= 0

    def test_different_metrics(self):
        """Test different distance metrics."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        div_fro = quantum_classical_divergence(
            l.eigv, l.eigV, t=1.0, init_node=50, metric="frobenius"
        )
        div_trace = quantum_classical_divergence(
            l.eigv, l.eigV, t=1.0, init_node=50, metric="trace"
        )
        div_max = quantum_classical_divergence(
            l.eigv, l.eigV, t=1.0, init_node=50, metric="max"
        )

        assert div_fro >= 0
        assert div_trace >= 0
        assert div_max >= 0


class TestInterferenceVisibility:
    """Test interference visibility computation."""

    def test_visibility_bounds(self):
        """Test 0 <= V <= 1."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        t_array = np.linspace(0, 10, 50)
        prob_t, visibility = compute_interference_visibility(
            l.eigv, l.eigV, t_array, node_i=50, node_j=60
        )

        assert 0 <= visibility <= 1

    def test_transition_probability_range(self):
        """Test transition probabilities are valid."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        t_array = np.linspace(0, 10, 50)
        prob_t, visibility = compute_interference_visibility(
            l.eigv, l.eigV, t_array, node_i=50, node_j=60
        )

        assert np.all(prob_t >= 0)
        assert np.all(prob_t <= 1)


class TestQuantumObservablesBatch:
    """Test batch computation of quantum observables."""

    def test_observable_keys(self):
        """Test that all expected keys are present."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        t_array = np.linspace(0, 5, 50)
        obs = compute_quantum_observables_from_eigenvalues(
            l.eigv, l.eigV, l.N, t_array, init_node=50
        )

        expected_keys = {
            "time",
            "von_neumann_entropy",
            "coherence",
            "purity",
            "participation_ratio",
            "probability_distribution",
        }
        assert set(obs.keys()) == expected_keys

    def test_observable_shapes(self):
        """Test shapes of returned arrays."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        n_times = 50
        t_array = np.linspace(0, 5, n_times)
        obs = compute_quantum_observables_from_eigenvalues(
            l.eigv, l.eigV, l.N, t_array, init_node=50
        )

        assert obs["time"].shape == (n_times,)
        assert obs["von_neumann_entropy"].shape == (n_times,)
        assert obs["coherence"].shape == (n_times,)
        assert obs["purity"].shape == (n_times,)
        assert obs["participation_ratio"].shape == (n_times,)
        assert obs["probability_distribution"].shape == (n_times, l.N)

    def test_purity_bounds(self):
        """Test 1/N <= Tr(ρ²) <= 1."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        t_array = np.linspace(0, 5, 50)
        obs = compute_quantum_observables_from_eigenvalues(
            l.eigv, l.eigV, l.N, t_array, init_node=50
        )

        assert np.all(obs["purity"] >= 1.0 / l.N - 1e-10)
        assert np.all(obs["purity"] <= 1.0 + 1e-10)


class TestSignedGraphQuantumMethods:
    """Test quantum methods integrated into SignedGraph class."""

    def test_compute_quantum_propagator_method(self):
        """Test SignedGraph.compute_quantum_propagator()."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        U_t = l.compute_quantum_propagator(t=1.0, full_matrix=True)

        # Check unitarity
        identity = U_t @ U_t.conj().T
        assert np.allclose(identity, np.eye(l.N), atol=1e-10)

    def test_quantum_walk_probabilities_method(self):
        """Test SignedGraph.quantum_walk_probabilities()."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        P_t = l.quantum_walk_probabilities(t=1.0, init_node=50)

        assert np.isclose(np.sum(P_t), 1.0, atol=1e-10)
        assert np.all(P_t >= 0)

    def test_quantum_observables_time_series_method(self):
        """Test SignedGraph.quantum_observables_time_series()."""
        l = Lattice2D(side1=10, seed=42)
        l.compute_laplacian_spectrum_weigV()

        t_array = np.linspace(0, 5, 50)
        obs = l.quantum_observables_time_series(t_array, init_node=50)

        assert "von_neumann_entropy" in obs
        assert "coherence" in obs
        assert obs["time"].shape == (50,)

    def test_no_eigendecomposition_error(self):
        """Test error when eigendecomposition not computed."""
        l = Lattice2D(side1=10, seed=42)

        with pytest.raises(ValueError, match="Eigendecomposition required"):
            l.compute_quantum_propagator(t=1.0)

        with pytest.raises(ValueError, match="Eigendecomposition required"):
            l.quantum_walk_probabilities(t=1.0, init_node=50)

        with pytest.raises(ValueError, match="Eigendecomposition required"):
            l.quantum_observables_time_series(np.linspace(0, 5, 50))


class TestSignedGraphs:
    """Test quantum propagators on signed graphs."""

    def test_signed_lattice_quantum_walk(self):
        """Test quantum walk on signed lattice."""
        l = Lattice2D(side1=10, pflip=0.2, seed=42)
        l.flip_random_fract_edges()
        l.compute_laplacian_spectrum_weigV()

        P_t = quantum_probability_distribution(
            l.eigv, l.eigV, t=1.0, init_node=50
        )

        # Probability conservation holds even with negative edges
        assert np.isclose(np.sum(P_t), 1.0, atol=1e-10)

    def test_signed_erdos_renyi_quantum_walk(self):
        """Test quantum walk on signed Erdos-Renyi graph."""
        er = ErdosRenyi(n=50, p=0.1, pflip=0.2, seed=42)
        er.compute_laplacian_spectrum_weigV()

        P_t = quantum_probability_distribution(
            er.eigv, er.eigV, t=1.0, init_node=25
        )

        assert np.isclose(np.sum(P_t), 1.0, atol=1e-10)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

"""Parity tests: verify NX and GT produce identical results for shared methods.

Strategy: create a graph in NX, convert to GT via ``from_networkx()``,
compute the same observable on both, assert numerical equivalence.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT


# ---------- fixtures ----------

@pytest.fixture
def paired_graphs():
    """Create matching NX and GT 8x8 square lattice with frustration."""
    nx_sg = Lattice2DNX(side1=8, geo="sqr", pflip=0.25, seed=42)
    nx_sg.flip_random_fract_edges()

    # NX stores signs as 'weight' attribute, converter expects 'sign' by default
    gt_sg = SignedGraphGT.from_networkx(
        nx_sg.gr["G"], pflip=0.25, seed=42
    )
    # Transfer signs: NX uses 'weight', GT uses 'sign' property
    # from_networkx uses nx_to_gt which looks for 'sign' by default.
    # NX Lattice2D stores edge signs in the 'weight' attribute.
    # Re-convert with correct attribute:
    from lrgsglib.graphs.gt._converters import nx_to_gt
    G_gt = nx_to_gt(nx_sg.gr["G"], sign_attr="weight")
    gt_sg = SignedGraphGT(G=G_gt, pflip=0.25, seed=42)

    # Verify same topology and signs
    assert nx_sg.N == gt_sg.N
    assert gt_sg.count_negative_edges() == sum(
        1 for _, _, d in nx_sg.gr["G"].edges(data=True) if d.get("weight", 1) < 0
    )
    return nx_sg, gt_sg


@pytest.fixture
def paired_with_spectrum(paired_graphs):
    """Paired graphs with full eigendecomposition computed."""
    nx_sg, gt_sg = paired_graphs
    nx_sg.compute_laplacian_spectrum_weigV()
    gt_sg.compute_laplacian_spectrum_weigV()
    return nx_sg, gt_sg


# ---------- Info theory parity ----------

class TestInfoTheoryParity:
    """Shannon and Renyi entropy must match between engines."""

    def test_shannon_entropy(self, paired_graphs):
        nx_sg, gt_sg = paired_graphs
        nx_sg.compute_signed_laplacian_entropy(steps=100, t1=-2, t2=5)
        gt_sg.compute_signed_laplacian_entropy(steps=100, t1=-2, t2=5)

        np.testing.assert_allclose(
            nx_sg.entropy, gt_sg.entropy, atol=1e-10,
            err_msg="Shannon entropy mismatch"
        )
        np.testing.assert_allclose(
            nx_sg.specific_heat, gt_sg.specific_heat, atol=1e-10,
            err_msg="Specific heat mismatch"
        )
        np.testing.assert_allclose(
            nx_sg.tauscale, gt_sg.tauscale, atol=1e-10,
            err_msg="Tau scale mismatch"
        )

    def test_renyi_entropy(self, paired_graphs):
        nx_sg, gt_sg = paired_graphs
        for q in [0.5, 2.0, 3.0]:
            nx_sg.compute_renyi_entropy_profile(q=q, steps=100)
            gt_sg.compute_renyi_entropy_profile(q=q, steps=100)

            nx_res = nx_sg.renyi_results[q]
            gt_res = gt_sg.renyi_results[q]
            np.testing.assert_allclose(
                nx_res["renyi_entropy"], gt_res["renyi_entropy"],
                atol=1e-10, err_msg=f"Renyi entropy mismatch for q={q}"
            )
            np.testing.assert_allclose(
                nx_res["specific_heat"], gt_res["specific_heat"],
                atol=1e-10, err_msg=f"Renyi specific heat mismatch for q={q}"
            )

    def test_get_entropy_cached(self, paired_graphs):
        """get_entropy() should compute lazily and match."""
        nx_sg, gt_sg = paired_graphs
        nx_e = nx_sg.get_entropy()
        gt_e = gt_sg.get_entropy()
        np.testing.assert_allclose(nx_e, gt_e, atol=1e-10)

    def test_get_specific_heat_cached(self, paired_graphs):
        nx_sg, gt_sg = paired_graphs
        nx_c = nx_sg.get_specific_heat()
        gt_c = gt_sg.get_specific_heat()
        np.testing.assert_allclose(nx_c, gt_c, atol=1e-10)


# ---------- Dynamics energy parity ----------

class TestDynamicsEnergyParity:
    """RBIM and SK-spherical energies from eigenvectors must match."""

    def test_rbim_energy_single(self, paired_with_spectrum):
        nx_sg, gt_sg = paired_with_spectrum
        nx_sg.compute_rbim_energy_eigV(which=0)
        gt_sg.compute_rbim_energy_eigV(which=0)

        assert nx_sg.energy_eigV_RBIM[0] == pytest.approx(
            gt_sg.energy_eigV_RBIM[0], abs=1e-10
        )

    def test_sksph_energy_single(self, paired_with_spectrum):
        nx_sg, gt_sg = paired_with_spectrum
        nx_sg.compute_sksph_energy_eigV(which=0)
        gt_sg.compute_sksph_energy_eigV(which=0)

        assert nx_sg.energy_eigV_SKSPH[0] == pytest.approx(
            gt_sg.energy_eigV_SKSPH[0], abs=1e-10
        )

    def test_rbim_energy_all(self, paired_with_spectrum):
        nx_sg, gt_sg = paired_with_spectrum
        nx_arr = nx_sg.get_all_rbim_energy_eigV()
        gt_arr = gt_sg.get_all_rbim_energy_eigV()
        np.testing.assert_allclose(nx_arr, gt_arr, atol=1e-10)

    def test_sksph_energy_all(self, paired_with_spectrum):
        nx_sg, gt_sg = paired_with_spectrum
        nx_arr = nx_sg.get_all_sksph_energy_eigV()
        gt_arr = gt_sg.get_all_sksph_energy_eigV()
        np.testing.assert_allclose(nx_arr, gt_arr, atol=1e-10)


# ---------- Spectral gap parity ----------

class TestSpectralGapParity:
    """Spectral gap must match between engines."""

    def test_gap(self, paired_graphs):
        nx_sg, gt_sg = paired_graphs
        nx_sg.compute_gap()
        gt_sg.compute_gap()
        assert nx_sg.gap == pytest.approx(gt_sg.gap, abs=1e-10)

    def test_gap_between(self, paired_graphs):
        nx_sg, gt_sg = paired_graphs
        nx_sg.compute_gap_between(low=0, high=2)
        gt_sg.compute_gap_between(low=0, high=2)
        assert nx_sg.gap == pytest.approx(gt_sg.gap, abs=1e-10)

    def test_get_gap_cached(self, paired_graphs):
        nx_sg, gt_sg = paired_graphs
        assert nx_sg.get_gap() == pytest.approx(gt_sg.get_gap(), abs=1e-10)


# ---------- Quantum parity ----------

class TestQuantumParity:
    """Quantum propagator must match between engines."""

    def test_quantum_walk_probs(self, paired_with_spectrum):
        nx_sg, gt_sg = paired_with_spectrum
        for t in [0.1, 1.0, 5.0]:
            nx_probs = nx_sg.quantum_walk_probabilities(t=t, init_node=0)
            gt_probs = gt_sg.quantum_walk_probabilities(t=t, init_node=0)
            np.testing.assert_allclose(
                nx_probs, gt_probs, atol=1e-10,
                err_msg=f"Quantum walk mismatch at t={t}"
            )

    def test_quantum_observables(self, paired_with_spectrum):
        nx_sg, gt_sg = paired_with_spectrum
        t_arr = np.linspace(0.1, 5.0, 20)
        nx_obs = nx_sg.quantum_observables_time_series(t_arr, init_node=0)
        gt_obs = gt_sg.quantum_observables_time_series(t_arr, init_node=0)
        for key in nx_obs:
            np.testing.assert_allclose(
                nx_obs[key], gt_obs[key], atol=1e-10,
                err_msg=f"Quantum observable '{key}' mismatch"
            )

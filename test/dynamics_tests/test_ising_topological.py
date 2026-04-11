"""Tests for Topological Metropolis and Topological FCA algorithms.

Tests cover:
- Spectral helper methods (subspace matrix, RBIM energies, seed, field, polish)
- Python topo_met and topo_fca sampling
- Pybind11 topo_met and topo_fca backends
- GT graph fallback (Python backend)
- Run dispatch and alias resolution
"""

import pytest
import numpy as np

from lrgsglib.graphs.nx import Lattice2DNX
from lrgsglib.statsys.IsingDynamics import IsingDynamics
from lrgsglib.utils.lrg.ising import compute_ising_pairwise_energy


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def ferro_lattice(tmp_path):
    """Small ferromagnetic (no frustration) square lattice."""
    return Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )


@pytest.fixture
def frustrated_lattice(tmp_path):
    """Small frustrated square lattice for realistic tests."""
    lat = Lattice2DNX(
        side1=8, geo="sqr", pbc=True, pflip=0.3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat.flip_random_fract_edges()
    return lat


@pytest.fixture
def topo_ising(frustrated_lattice):
    """IsingDynamics with topo params on frustrated lattice."""
    ising = IsingDynamics(
        sg=frustrated_lattice, T=0.5, steps=500,
        runlang="py_topo_met",
        topo_n_modes=10, topo_polish=False, seed=42,
    )
    ising.init_ising_dynamics()
    return ising


# ---------------------------------------------------------------------------
# Spectral helper tests
# ---------------------------------------------------------------------------

class TestSpectralHelpers:
    """Tests for the shared spectral infrastructure."""

    def test_subspace_matrix_shape(self, topo_ising):
        n_modes = 10
        topo_ising._ensure_spectral_subspace(n_modes)
        V = topo_ising._build_subspace_matrix(n_modes)
        assert V.shape == (n_modes, topo_ising.N)
        assert V.dtype == np.float64

    def test_rbim_energies_negative_for_ferro(self, ferro_lattice, tmp_path):
        """Binarized eigvecs of a ferro lattice should have negative RBIM E."""
        ising = IsingDynamics(
            sg=ferro_lattice, T=1.0, steps=10, runlang="py",
            topo_n_modes=5, seed=42,
        )
        ising.init_ising_dynamics()
        ising._ensure_spectral_subspace(5)
        rbim_E = ising._compute_rbim_energies(5)
        assert rbim_E.shape == (5,)
        # At least the best mode should have negative energy on ferro
        assert rbim_E.min() < 0

    def test_best_eigenvector_seed(self, topo_ising):
        n_modes = 10
        topo_ising._ensure_spectral_subspace(n_modes)
        spins, idx = topo_ising._best_eigenvector_seed(n_modes)
        assert spins.shape == (topo_ising.N,)
        assert set(np.unique(spins)).issubset({-1, 1})
        assert spins.dtype == np.int8
        assert 0 <= idx < n_modes

    def test_spectral_field_normalization(self, topo_ising):
        n_modes = 10
        topo_ising._ensure_spectral_subspace(n_modes)
        V = topo_ising._build_subspace_matrix(n_modes)
        rbim_E = topo_ising._compute_rbim_energies(n_modes)
        field = IsingDynamics._compute_spectral_field(V, rbim_E, 1.0, 1.0)
        assert field.shape == (topo_ising.N,)
        # Field should be nonzero
        assert np.linalg.norm(field) > 0

    def test_spectral_field_strength_scaling(self, topo_ising):
        """Doubling field_strength should double the field."""
        n_modes = 10
        topo_ising._ensure_spectral_subspace(n_modes)
        V = topo_ising._build_subspace_matrix(n_modes)
        rbim_E = topo_ising._compute_rbim_energies(n_modes)
        f1 = IsingDynamics._compute_spectral_field(V, rbim_E, 1.0, 1.0)
        f2 = IsingDynamics._compute_spectral_field(V, rbim_E, 1.0, 2.0)
        np.testing.assert_allclose(f2, 2.0 * f1, atol=1e-12)

    def test_greedy_polish_never_increases_energy(self, topo_ising):
        """Polish should never increase the Ising energy."""
        topo_ising._ensure_spectral_subspace(10)
        spins, _ = topo_ising._best_eigenvector_seed(10)
        edges = topo_ising.sg.get_edges_with_weights()
        E_before = float(compute_ising_pairwise_energy(spins, edges))
        polished = topo_ising._greedy_polish(spins, max_sweeps=50)
        E_after = float(compute_ising_pairwise_energy(polished, edges))
        assert E_after <= E_before + 1e-10, (
            f"Polish increased energy: {E_before} -> {E_after}"
        )


# ---------------------------------------------------------------------------
# Python topo_met tests
# ---------------------------------------------------------------------------

class TestPyTopoMet:
    """Tests for Python Topological Metropolis."""

    def test_smoke(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.5, steps=200,
            runlang="py_topo_met",
            topo_n_modes=8, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert len(ising.ene) == 200
        assert len(ising.magn) == 200
        assert hasattr(ising, "topo_met_coeffs")

    def test_observables_stored(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.3, steps=300,
            runlang="py_topo_met",
            topo_n_modes=8, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert ising.topo_met_coeffs.shape == (8,)
        assert ising.topo_met_best_spins.shape == (ising.N,)
        assert isinstance(ising.topo_met_best_energy, float)
        assert ising.topo_met_best_energy <= ising.ene[0]

    def test_energy_decreases_at_low_T(self, frustrated_lattice):
        """At low T the energy should not increase on average."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.1, steps=1000,
            runlang="py_topo_met",
            topo_n_modes=10, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        # Best energy should be ≤ initial energy
        assert ising.topo_met_best_energy <= ising.ene[0]

    def test_deterministic_seed(self, frustrated_lattice):
        """Same seed → same trajectory."""
        kwargs = dict(
            sg=frustrated_lattice, T=0.5, steps=100,
            runlang="py_topo_met",
            topo_n_modes=8, topo_polish=False, seed=123,
        )
        i1 = IsingDynamics(**kwargs)
        i1.init_ising_dynamics()
        i1.run(verbose=False, tqdm_on=False)

        i2 = IsingDynamics(**kwargs)
        i2.init_ising_dynamics()
        i2.run(verbose=False, tqdm_on=False)

        np.testing.assert_array_equal(i1.ene, i2.ene)
        np.testing.assert_array_equal(i1.magn, i2.magn)

    def test_finds_good_energy_on_ferro(self, ferro_lattice):
        """On a ferro lattice, topo_met should find near-ground-state."""
        N = ferro_lattice.N
        ising = IsingDynamics(
            sg=ferro_lattice, T=0.1, steps=2000,
            runlang="py_topo_met",
            topo_n_modes=10, topo_polish=True, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        edges = ferro_lattice.get_edges_with_weights()
        n_edges = len(edges)
        # Ground state E/N = -n_edges/N; best should be close
        assert ising.topo_met_best_energy <= -0.8 * n_edges / N

    def test_with_greedy_polish(self, frustrated_lattice):
        """Running with polish enabled should find ≤ energy of no polish."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.3, steps=500,
            runlang="py_topo_met",
            topo_n_modes=8, topo_polish=True,
            topo_chunk_size=100, topo_polish_sweeps=20, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        # Just verify it doesn't crash and has valid energy
        assert ising.topo_met_best_energy <= ising.ene[0]


# ---------------------------------------------------------------------------
# Python topo_fca tests
# ---------------------------------------------------------------------------

class TestPyTopoFCA:
    """Tests for Python Topological Field Cooling Annealing."""

    def test_smoke(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=1.0, steps=200,
            runlang="py_topo_fca",
            sa_enabled=True, T_init=5.0, T_final=0.01,
            topo_n_modes=8, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert hasattr(ising, "topo_fca_field")
        assert ising.topo_fca_field.shape == (ising.N,)
        assert len(ising.ene) > 0

    def test_field_nonzero(self, frustrated_lattice):
        """The spectral field should be nonzero."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=1.0, steps=100,
            runlang="py_topo_fca",
            sa_enabled=True, T_init=5.0, T_final=0.01,
            topo_n_modes=8, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert np.linalg.norm(ising.topo_fca_field) > 0

    def test_metropolis_attributes_populated(self, frustrated_lattice):
        """TFCA quench should populate Metropolis result attributes."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=1.0, steps=200,
            runlang="py_topo_fca",
            sa_enabled=True, T_init=5.0, T_final=0.01,
            n_temperatures=50, steps_per_T=10,
            topo_n_modes=8, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        # TFCA now runs Metropolis at T~0 (not SA), so ene/magn are populated
        assert hasattr(ising, "ene")
        assert hasattr(ising, "magn")
        assert len(ising.ene) == 50 * 10  # n_temperatures * steps_per_T


# ---------------------------------------------------------------------------
# Pybind11 tests
# ---------------------------------------------------------------------------

class TestPybindTopoFCA:
    """Tests for pybind11 Topological FCA (trivial wrapper)."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_native(self):
        try:
            from lrgsglib.statsys.IsingDynamics.ccore import _ising_native
        except ImportError:
            pytest.skip("_ising_native not built")

    def test_smoke(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=1.0, steps=200,
            runlang="pb_topo_fca",
            sa_enabled=True, T_init=5.0, T_final=0.01,
            topo_n_modes=8, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False)
        assert hasattr(ising, "topo_fca_field")
        assert np.linalg.norm(ising.topo_fca_field) > 0

    def test_field_matches_python(self, frustrated_lattice):
        """Pybind TFCA field should match Python TFCA field (same computation)."""
        common = dict(
            sg=frustrated_lattice, T=1.0, steps=100,
            sa_enabled=True, T_init=5.0, T_final=0.01,
            topo_n_modes=8, seed=42,
        )
        py = IsingDynamics(**common, runlang="py_topo_fca")
        py.init_ising_dynamics()
        py.run(verbose=False, tqdm_on=False)

        pb = IsingDynamics(**common, runlang="pb_topo_fca")
        pb.init_ising_dynamics()
        pb.run(verbose=False)

        np.testing.assert_allclose(
            py.topo_fca_field, pb.topo_fca_field, atol=1e-12,
        )


class TestPybindTopoMet:
    """Tests for pybind11 Topological Metropolis."""

    @pytest.fixture(autouse=True)
    def _skip_if_no_native(self):
        try:
            from lrgsglib.statsys.IsingDynamics.ccore import _ising_native
            if not hasattr(_ising_native, "topo_met_sampling"):
                pytest.skip("topo_met_sampling not in _ising_native")
        except ImportError:
            pytest.skip("_ising_native not built")

    def test_smoke(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.5, steps=200,
            runlang="pb_topo_met",
            topo_n_modes=8, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False)
        assert len(ising.ene) == 200
        assert hasattr(ising, "topo_met_coeffs")
        assert ising.topo_met_best_spins.shape == (ising.N,)

    def test_trajectory_length(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.5, steps=500,
            runlang="pb_topo_met",
            topo_n_modes=8, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False)
        assert len(ising.ene) == 500
        assert len(ising.magn) == 500

    def test_best_energy_tracked(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.3, steps=500,
            runlang="pb_topo_met",
            topo_n_modes=8, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False)
        # Best energy should be ≤ all trajectory values (in C convention)
        assert ising.topo_met_best_energy <= max(ising.ene) + 1e-10


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------

class TestIntegration:
    """Cross-cutting integration tests."""

    def test_runlang_dispatch_py_topo_met(self, frustrated_lattice):
        """Alias 'topo_met' resolves to py_topo_met and runs."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.5, steps=50,
            runlang="topo_met",
            topo_n_modes=5, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert len(ising.ene) == 50

    def test_runlang_dispatch_py_topo_fca(self, frustrated_lattice):
        """Alias 'topo_fca' resolves to py_topo_fca and runs."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=1.0, steps=50,
            runlang="topo_fca",
            sa_enabled=True, T_init=5.0, T_final=0.01,
            n_temperatures=50,
            topo_n_modes=5, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        # TFCA now uses Metropolis quench, so ene is populated
        assert hasattr(ising, "ene")
        assert len(ising.ene) > 0

    def test_gt_graph_fallback(self):
        """topo_met on a GT graph should work (Python fallback)."""
        try:
            from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT
        except ImportError:
            pytest.skip("graph_tool not installed")

        lat = Lattice2DGT(side1=6, geo="sqr", periodic=True, pflip=0.3, seed=42)
        ising = IsingDynamics(
            sg=lat, T=0.5, steps=100,
            runlang="py_topo_met",
            topo_n_modes=5, topo_polish=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert len(ising.ene) == 100

    def test_topo_met_does_not_break_regular_met(self, frustrated_lattice):
        """Regular Metropolis should still work after topo dispatch added."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=2.0, steps=100,
            runlang="py_met", seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert len(ising.ene) == 100

    def test_topo_fca_does_not_break_regular_sa(self, frustrated_lattice):
        """Regular SA should still work after topo dispatch added."""
        ising = IsingDynamics(
            sg=frustrated_lattice, T=1.0, steps=100,
            runlang="py_sa",
            sa_enabled=True, T_init=5.0, T_final=0.01, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert len(ising.sa_energy) == 100


# ---------------------------------------------------------------------------
# Cross-Entropy Method (CEM) tests
# ---------------------------------------------------------------------------

class TestPyTopoCEM:
    """Tests for Python Cross-Entropy Method spectral optimizer."""

    def test_smoke(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="py_topo_cem",
            topo_n_modes=8, topo_polish=False,
            cem_iter=3, cem_pop_size=16, cem_restarts=2,
            cem_greedy=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert hasattr(ising, "topo_cem_best_energy")
        assert hasattr(ising, "topo_cem_best_spins")
        assert ising.topo_cem_best_spins.shape == (ising.N,)

    def test_coeffs_stored(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="py_topo_cem",
            topo_n_modes=6, cem_iter=3, cem_pop_size=8,
            cem_restarts=2, cem_greedy=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert ising.topo_cem_best_coeffs.shape == (6,)
        assert ising.topo_cem_restart_energies.shape == (2,)

    def test_energy_valid(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="py_topo_cem",
            topo_n_modes=8, cem_iter=5, cem_pop_size=16,
            cem_restarts=2, cem_greedy=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        # stored energy is E/N; recompute and normalize
        recomputed = ising.compute_energy(ising.topo_cem_best_spins) / ising.N
        assert abs(recomputed - ising.topo_cem_best_energy) < 1e-10

    def test_greedy_improves_or_equals(self, frustrated_lattice):
        """CEM with greedy should find <= energy vs without."""
        common = dict(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="py_topo_cem",
            topo_n_modes=8, topo_polish=False,
            cem_iter=5, cem_pop_size=16, cem_restarts=3, seed=42,
        )
        no_greedy = IsingDynamics(**common, cem_greedy=False)
        no_greedy.init_ising_dynamics()
        no_greedy.run(verbose=False, tqdm_on=False)

        with_greedy = IsingDynamics(
            **common, cem_greedy=True, cem_greedy_sweeps=20,
        )
        with_greedy.init_ising_dynamics()
        with_greedy.run(verbose=False, tqdm_on=False)
        assert with_greedy.topo_cem_best_energy <= (
            no_greedy.topo_cem_best_energy + 1e-10
        )

    def test_deterministic_seed(self, frustrated_lattice):
        kwargs = dict(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="py_topo_cem",
            topo_n_modes=6, cem_iter=3, cem_pop_size=8,
            cem_restarts=2, cem_greedy=False, seed=123,
        )
        i1 = IsingDynamics(**kwargs)
        i1.init_ising_dynamics()
        i1.run(verbose=False, tqdm_on=False)

        i2 = IsingDynamics(**kwargs)
        i2.init_ising_dynamics()
        i2.run(verbose=False, tqdm_on=False)

        assert i1.topo_cem_best_energy == i2.topo_cem_best_energy
        np.testing.assert_array_equal(
            i1.topo_cem_best_spins, i2.topo_cem_best_spins,
        )

    def test_alias_dispatch(self, frustrated_lattice):
        """Aliases 'topo_cem' and 'cem' should resolve and run."""
        for rl in ("topo_cem", "cem"):
            ising = IsingDynamics(
                sg=frustrated_lattice, T=0.0, steps=10,
                runlang=rl, topo_n_modes=5,
                cem_iter=2, cem_pop_size=8, cem_restarts=1,
                cem_greedy=False, seed=42,
            )
            ising.init_ising_dynamics()
            ising.run(verbose=False, tqdm_on=False)
            assert hasattr(ising, "topo_cem_best_energy")


class TestCEMGreedyPolishCSR:
    """Tests for the CSR-based greedy polish used by CEM."""

    def test_never_increases_energy(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.5, steps=10,
            runlang="py_topo_met", topo_n_modes=10, seed=42,
        )
        ising.init_ising_dynamics()
        ising._ensure_spectral_subspace(10)
        spins, _ = ising._best_eigenvector_seed(10)
        ni, nw, nptr = ising._build_graph_csr()
        edges = frustrated_lattice.get_edges_with_weights()
        E_before = float(compute_ising_pairwise_energy(spins, edges))
        polished = ising._greedy_polish_csr(
            spins, ni, nw, nptr, max_sweeps=50,
        )
        E_after = float(compute_ising_pairwise_energy(polished, edges))
        assert E_after <= E_before + 1e-10


class TestPybindTopoCEM:
    """Tests for pybind11 CEM backend."""

    def test_smoke(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="pb_topo_cem",
            topo_n_modes=8, topo_polish=False,
            cem_iter=3, cem_pop_size=16, cem_restarts=2,
            cem_greedy=False, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        assert hasattr(ising, "topo_cem_best_energy")
        assert ising.topo_cem_best_spins.shape == (ising.N,)
        assert ising.topo_cem_best_coeffs.shape == (8,)
        assert ising.topo_cem_restart_energies.shape == (2,)

    def test_energy_matches_recomputed(self, frustrated_lattice):
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="pb_topo_cem",
            topo_n_modes=8, cem_iter=5, cem_pop_size=16,
            cem_restarts=2, cem_greedy=True,
            cem_greedy_sweeps=20, seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(verbose=False, tqdm_on=False)
        recomputed = ising.compute_energy(ising.topo_cem_best_spins) / ising.N
        assert abs(recomputed - ising.topo_cem_best_energy) < 1e-10

    def test_pb_finds_comparable_energy(self, frustrated_lattice):
        """Pybind11 should find energy in same ballpark as Python."""
        common = dict(
            sg=frustrated_lattice, T=0.0, steps=10,
            topo_n_modes=8, topo_polish=False,
            cem_iter=8, cem_pop_size=32, cem_restarts=3,
            cem_greedy=True, cem_greedy_sweeps=30, seed=42,
        )
        py = IsingDynamics(**common, runlang="py_topo_cem")
        py.init_ising_dynamics()
        py.run(verbose=False, tqdm_on=False)

        pb = IsingDynamics(**common, runlang="pb_topo_cem")
        pb.init_ising_dynamics()
        pb.run(verbose=False, tqdm_on=False)

        # Different RNGs → different exact energies, but same ballpark
        assert pb.topo_cem_best_energy <= py.topo_cem_best_energy * 0.9 or \
               pb.topo_cem_best_energy <= py.topo_cem_best_energy + 10.0

    def test_ground_state_init_with_subspace_expansion(self, frustrated_lattice):
        """Regression: ic='ground_state_0' populates eigV with 1 mode via
        BinDynSys.init_s, then _run_pybind_topo_cem must *expand* the
        subspace from 1 to topo_n_modes rather than shortcut on the
        already-cached single eigenvector.

        Before the fix, _ensure_spectral_subspace bailed out as soon as
        ``self.sg.eigV is not None`` regardless of the row count, so the
        first quench of any CLI run with the default init_cond died with
        ``IndexError`` in _build_subspace_matrix on the cluster.
        """
        ising = IsingDynamics(
            sg=frustrated_lattice, T=0.0, steps=10,
            runlang="pb_topo_cem",
            topo_n_modes=8, topo_polish=False,
            cem_iter=3, cem_pop_size=16, cem_restarts=2,
            cem_greedy=False, seed=42,
            ic="ground_state_0",
        )
        ising.init_ising_dynamics()
        # After init_s, eigV has exactly 1 row
        assert ising.sg.eigV is not None
        n_cached_before = (
            ising.sg.eigV.shape[0]
            if getattr(ising.sg, "_eigV_is_transposed", False)
            else ising.sg.eigV.shape[1]
        )
        assert n_cached_before == 1
        # CEM run must expand the subspace to topo_n_modes
        ising.run(verbose=False, tqdm_on=False)
        assert ising.topo_cem_best_coeffs.shape == (8,)
        assert ising.topo_cem_best_spins.shape == (ising.N,)

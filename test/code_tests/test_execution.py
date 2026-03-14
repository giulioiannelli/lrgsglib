"""
Test that core methods execute without error (not checking correctness).

Validates:
- Edge sign flipping
- Matrix updates
- Spectral computation with different backends
- Dynamics init + run (Python backend)
- Dynamics init + run (C backend)
"""

import pytest
import numpy as np


# ===================================================================
# Graph methods
# ===================================================================


@pytest.mark.code
def test_flip_random_fract_edges(tmp_path):
    """flip_random_fract_edges() runs on a lattice with pflip > 0."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat.flip_random_fract_edges()
    # Verify some edges are negative
    G = lat.gr["G"]
    signs = [G[u][v].get("weight", 1.0) for u, v in G.edges()]
    assert any(s < 0 for s in signs), "pflip=0.3 should produce negative edges"


@pytest.mark.code
def test_upd_graph_matrices(tmp_path):
    """upd_graph_matrices() produces adjacency and Laplacian."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    lat.upd_graph_matrices()
    assert hasattr(lat, "adj") and lat.adj is not None
    assert hasattr(lat, "lap") and lat.lap is not None
    assert lat.adj.shape == (lat.N, lat.N)
    assert lat.lap.shape == (lat.N, lat.N)


# ===================================================================
# Spectral methods
# ===================================================================


@pytest.mark.code
@pytest.mark.parametrize("backend", ["numpy", "scipy"])
def test_get_graph_lspectrum(backend, tmp_path):
    """get_graph_lspectrum() returns (L, eigenvalues) for each backend."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    L, eigvals = get_graph_lspectrum(lat.gr["G"], library=backend)
    assert L.shape == (lat.N, lat.N)
    assert len(eigvals) == lat.N


@pytest.mark.code
def test_compute_entropy_observables(tmp_path):
    """compute_entropy_observables_from_eigenvalues() returns 4-tuple."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
    from lrgsglib.utils.lrg.infocomm import (
        compute_entropy_observables_from_eigenvalues,
    )

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    _, eigvals = get_graph_lspectrum(lat.gr["G"], library="numpy")
    result = compute_entropy_observables_from_eigenvalues(
        eigenvalues=eigvals, num_nodes=lat.N,
    )
    assert len(result) == 4
    S, C, var, tau = result
    assert len(S) == len(tau)
    assert len(C) == len(tau)


# ===================================================================
# Dynamics: Python backend
# ===================================================================


@pytest.mark.code
def test_ising_python_runs(tmp_path):
    """Ising Python dynamics runs without error."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(sg=lat, T=2.0, steps=5, runlang="py", seed=42)
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False)
    assert len(ising.s) == lat.N


@pytest.mark.code
def test_contact_process_python_runs(tmp_path):
    """Contact Process EI Python dynamics runs without error."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    cp = ContactProcessEI(sg=lat, gamma=1.5, steps=5, runlang="py", seed=42)
    cp.init_contact_dynamics()
    cp.run(tqdm_on=False)
    assert len(cp.s) == lat.N


@pytest.mark.code
def test_voter_python_runs(tmp_path):
    """Voter Model Python dynamics runs without error."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter = VoterModel(sg=lat, steps=5, runlang="py", seed=42)
    voter.init_voter_dynamics()
    voter.run(tqdm_on=False)
    assert len(voter.s) == lat.N


# ===================================================================
# Dynamics: C backend
# ===================================================================


@pytest.mark.code
@pytest.mark.integration
def test_ising_c_backend_runs(tmp_path):
    """Ising C1b backend runs without error."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(
        sg=lat, T=2.0, steps=5, runlang="C0E",
        seed=42, save_magnetization=True,
    )
    ising.init_ising_dynamics()
    ising.run(tqdm_on=False, verbose=False)
    assert len(ising.s) == lat.N


@pytest.mark.code
@pytest.mark.integration
def test_contact_process_c_backend_runs(tmp_path):
    """Contact Process C1c backend runs without error."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    cp = ContactProcessEI(
        sg=lat, gamma=1.5, steps=5, runlang="C1D", seed=42,
    )
    cp.init_contact_dynamics()
    cp.run(tqdm_on=False, verbose=False)
    assert len(cp.s) == lat.N


@pytest.mark.code
@pytest.mark.integration
def test_voter_c_backend_runs(tmp_path):
    """Voter C backend runs without error (snapshot mode)."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    lat = Lattice2DNX(
        side1=6, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter = VoterModel(
        sg=lat, steps=5, runlang="C0S", seed=42,
    )
    voter.init_voter_dynamics()
    voter.run(tqdm_on=False, verbose=False)
    assert len(voter.s) == lat.N

"""
Test that all graph types and dynamics classes can be instantiated.

Uses small systems (side=4) to keep runtime low. Tests verify that
objects are created successfully with N > 0 and Ne > 0.
"""

import pytest
import numpy as np


# ===================================================================
# Graph instantiation tests
# ===================================================================


# --- Lattice 2D (parametrize geometry) ---

@pytest.mark.code
@pytest.mark.parametrize("geo", ["sqr", "tri", "hex"])
@pytest.mark.parametrize("pbc", [True, False])
def test_lattice2d(geo, pbc, tmp_path):
    from lrgsglib.graphs.nx import Lattice2DNX
    # tri+PBC needs m>2 and n>4, so use side=6 for safety
    side = 6 if (geo == "tri" and pbc) else 4
    lat = Lattice2DNX(
        side1=side, geo=geo, pbc=pbc, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert lat.N > 0
    assert lat.Ne > 0


@pytest.mark.code
@pytest.mark.parametrize("pflip", [0.0, 0.3, 0.5])
def test_lattice2d_pflip(pflip, tmp_path):
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=pflip,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    if pflip > 0:
        lat.flip_random_fract_edges()
    assert lat.N == 16


# --- Lattice 3D (parametrize geometry) ---

@pytest.mark.code
@pytest.mark.parametrize("geo", ["sc", "bcc", "fcc"])
def test_lattice3d(geo, tmp_path):
    from lrgsglib.graphs.nx import Lattice3DNX
    lat = Lattice3DNX(
        dim=4, geo=geo, pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert lat.N > 0
    assert lat.Ne > 0


# --- Random graph models ---

@pytest.mark.code
def test_erdos_renyi(tmp_path):
    from lrgsglib.graphs.nx import ErdosRenyiNX
    g = ErdosRenyiNX(
        n=20, p=0.3, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N == 20
    assert g.Ne > 0


@pytest.mark.code
def test_barabasi_albert(tmp_path):
    from lrgsglib.graphs.nx import BarabasiAlbertNX
    g = BarabasiAlbertNX(
        n=20, m=2, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N == 20
    assert g.Ne > 0


@pytest.mark.code
def test_watts_strogatz(tmp_path):
    from lrgsglib.graphs.nx import WattsStrogatzNX
    g = WattsStrogatzNX(
        n=20, k=4, p=0.3, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N == 20
    assert g.Ne > 0


@pytest.mark.code
def test_k_regular(tmp_path):
    from lrgsglib.graphs.nx import kRegularGraphNX
    g = kRegularGraphNX(
        n=20, k=4, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N == 20


@pytest.mark.code
def test_complete_graph(tmp_path):
    # CompleteGraphNX is abstract; FullyConnectedNX is the concrete impl
    from lrgsglib.graphs.nx import FullyConnectedNX
    g = FullyConnectedNX(
        N=10, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N == 10
    # Complete graph: N*(N-1)/2 edges
    assert g.Ne == 45


@pytest.mark.code
def test_fully_connected(tmp_path):
    from lrgsglib.graphs.nx import FullyConnectedNX
    g = FullyConnectedNX(
        N=10, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N == 10


# --- Multispectral / Hierarchical ---

@pytest.mark.code
def test_multiplicative_cascade(tmp_path):
    from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX
    g = MultiplicativeCascadeGraphNX(
        p1=0.8, p2=0.6, fraction=0.4, iterations=3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N > 0
    assert g.Ne > 0


@pytest.mark.code
def test_vicsek(tmp_path):
    from lrgsglib.graphs.nx import VicsekGraphNX
    g = VicsekGraphNX(
        N=20, k=2,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N > 0
    assert g.Ne > 0


# --- Fractal ---

@pytest.mark.code
def test_sierpinski(tmp_path):
    from lrgsglib.graphs.nx import SierpinskiNX
    g = SierpinskiNX(
        n=3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N > 0
    assert g.Ne > 0


# --- Dirac ---

@pytest.mark.code
def test_dirac_comb(tmp_path):
    from lrgsglib.graphs.nx import DiracCombGraphNX
    g = DiracCombGraphNX(
        base_nodes=16, fiber_nodes=3,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.N > 0
    assert g.Ne > 0


# --- Temporal ---

@pytest.mark.code
def test_temporal_graph(tmp_path):
    from lrgsglib.graphs.nx import TemporalGraphNX
    g = TemporalGraphNX(
        n_nodes=10, n_snapshots=5,
    )
    assert g is not None


# ===================================================================
# Dynamics instantiation tests
# ===================================================================


@pytest.mark.code
def test_ising_dynamics_instantiation(tmp_path):
    """IsingDynamics can be created with sg= and T= params."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import IsingDynamics

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    ising = IsingDynamics(sg=lat, T=2.0, steps=10, runlang="py", seed=42)
    assert ising.T == 2.0
    assert ising.steps == 10


@pytest.mark.code
def test_contact_process_ei_instantiation(tmp_path):
    """ContactProcessEI can be created with sg= and gamma= params."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import ContactProcessEI

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    cp = ContactProcessEI(sg=lat, gamma=1.5, steps=10, runlang="py", seed=42)
    assert cp.gamma is not None


@pytest.mark.code
def test_voter_model_instantiation(tmp_path):
    """VoterModel can be created with sg= param."""
    from lrgsglib.graphs.nx import Lattice2DNX
    from lrgsglib.statsys import VoterModel

    lat = Lattice2DNX(
        side1=4, geo="sqr", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    voter = VoterModel(sg=lat, steps=10, runlang="py", seed=42)
    assert voter.steps == 10

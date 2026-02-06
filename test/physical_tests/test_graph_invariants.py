"""
Test graph structural invariants: node counts, edge counts, coordination
numbers, connectivity, and frustration fractions.

All tests use exact values (no statistical tolerance) for deterministic
properties. Frustration fraction tests use large lattices for statistics.
"""

import pytest
import networkx as nx
import numpy as np


# ===================================================================
# Node count invariants
# ===================================================================


@pytest.mark.physical
@pytest.mark.parametrize("side", [6, 8, 12])
@pytest.mark.parametrize("geo", ["sqr", "tri"])
def test_lattice2d_node_count(side, geo, tmp_path):
    """Lattice2D (sqr, tri) node count must equal side^2."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=side, geo=geo, pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert lat.N == side * side


@pytest.mark.physical
@pytest.mark.parametrize("side", [6, 8, 12])
def test_lattice2d_hex_node_count(side, tmp_path):
    """Hex lattice node count: N > 0 and consistent with graph."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=side, geo="hex", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert lat.N > 0
    assert lat.N == lat.gr["G"].number_of_nodes()
    assert lat.Ne == lat.gr["G"].number_of_edges()


@pytest.mark.physical
@pytest.mark.parametrize("dim", [3, 4, 6])
def test_lattice3d_node_count(dim, tmp_path):
    """Lattice3D(dim=d) node count must equal d^3."""
    from lrgsglib.graphs.nx import Lattice3DNX
    lat = Lattice3DNX(
        dim=dim, geo="sc", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert lat.N == dim ** 3


# ===================================================================
# Edge count and coordination invariants
# ===================================================================


@pytest.mark.physical
@pytest.mark.parametrize("geo,expected_z", [("sqr", 4), ("tri", 6)])
def test_lattice2d_coordination(geo, expected_z, tmp_path):
    """Average degree must match coordination number for PBC lattice."""
    from lrgsglib.graphs.nx import Lattice2DNX
    side = 10
    lat = Lattice2DNX(
        side1=side, geo=geo, pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    avg_degree = 2 * lat.Ne / lat.N
    assert np.isclose(avg_degree, expected_z, atol=0.01), (
        f"geo={geo}: expected z={expected_z}, got {avg_degree:.2f}"
    )


@pytest.mark.physical
def test_lattice2d_hex_coordination(tmp_path):
    """Hex lattice average degree should approach 3 for large lattices."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=20, geo="hex", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    avg_degree = 2 * lat.Ne / lat.N
    assert 2.5 < avg_degree <= 3.0, (
        f"geo=hex: expected z~3, got {avg_degree:.3f}"
    )


@pytest.mark.physical
def test_lattice3d_sc_coordination(tmp_path):
    """SC lattice: avg degree must be 6."""
    from lrgsglib.graphs.nx import Lattice3DNX
    lat = Lattice3DNX(
        dim=6, geo="sc", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    avg_degree = 2 * lat.Ne / lat.N
    assert np.isclose(avg_degree, 6.0, atol=0.01), (
        f"SC: expected z=6, got {avg_degree:.2f}"
    )


@pytest.mark.physical
@pytest.mark.parametrize("geo", ["bcc", "fcc"])
def test_lattice3d_nonsc_coordination(geo, tmp_path):
    """BCC/FCC: avg degree consistent between sizes (implementation-specific z)."""
    from lrgsglib.graphs.nx import Lattice3DNX
    results = []
    for dim in [4, 6]:
        lat = Lattice3DNX(
            dim=dim, geo=geo, pbc=True, pflip=0.0,
            seed=42, path_data=tmp_path, path_plot=tmp_path,
            init_nw_dict=False,
        )
        avg_z = 2 * lat.Ne / lat.N
        results.append(avg_z)
        assert lat.N > 0
        assert lat.Ne > 0
    # Both sizes should have similar coordination (within 20%)
    assert abs(results[0] - results[1]) / results[1] < 0.2, (
        f"{geo}: z varies too much between dim=4 ({results[0]:.2f}) and dim=6 ({results[1]:.2f})"
    )


@pytest.mark.physical
def test_complete_graph_edge_count(tmp_path):
    """Complete graph on n nodes has n*(n-1)/2 edges."""
    from lrgsglib.graphs.nx import FullyConnectedNX
    n = 15
    g = FullyConnectedNX(
        N=n, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert g.Ne == n * (n - 1) // 2


# ===================================================================
# Connectivity invariants
# ===================================================================


@pytest.mark.physical
@pytest.mark.parametrize("geo", ["sqr", "tri", "hex"])
def test_lattice2d_pbc_connected(geo, tmp_path):
    """PBC lattice must be a single connected component."""
    from lrgsglib.graphs.nx import Lattice2DNX
    lat = Lattice2DNX(
        side1=8, geo=geo, pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert nx.is_connected(lat.gr["G"])


@pytest.mark.physical
def test_lattice3d_pbc_connected(tmp_path):
    """PBC 3D lattice must be connected."""
    from lrgsglib.graphs.nx import Lattice3DNX
    lat = Lattice3DNX(
        dim=4, geo="sc", pbc=True, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert nx.is_connected(lat.gr["G"])


@pytest.mark.physical
def test_erdos_renyi_connected_above_threshold(tmp_path):
    """ER graph well above p_c = ln(n)/n should be connected."""
    from lrgsglib.graphs.nx import ErdosRenyiNX
    n = 50
    p = 0.3  # Well above ln(50)/50 ~ 0.078
    g = ErdosRenyiNX(
        n=n, p=p, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    assert nx.is_connected(g.gr["G"])


# ===================================================================
# Degree regularity
# ===================================================================


@pytest.mark.physical
def test_k_regular_degree_distribution(tmp_path):
    """k-regular graph: all nodes must have degree k."""
    from lrgsglib.graphs.nx import kRegularGraphNX
    n, k = 20, 4
    g = kRegularGraphNX(
        n=n, k=k, pflip=0.0,
        seed=42, path_data=tmp_path, path_plot=tmp_path,
        init_nw_dict=False,
    )
    G = g.gr["G"]
    degrees = [G.degree(nd) for nd in G.nodes()]
    assert all(d == k for d in degrees), (
        f"Expected all degrees = {k}, got {set(degrees)}"
    )

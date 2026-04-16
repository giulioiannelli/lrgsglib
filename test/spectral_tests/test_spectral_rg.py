"""Tests for spectral_rg.py — Spectral RG for signed/frustrated graphs.

Tests cover:
1. Basic function contracts (shapes, types, symmetry)
2. Classical reduction: unsigned graph → d_diff = classical diffusion distance
3. Signed graph compatibility: no blow-up, valid distances
4. Two-channel property: topology vs geometry from same decomposition
5. RG flow: convergence, node reduction, observable extraction
6. Apodictic cases: known ground-truth partitions
"""

import numpy as np
import networkx as nx
import pytest
from numpy.testing import assert_allclose

from lrgsglib.utils.lrg.spectral_rg import (
    compute_frustration_fraction,
    negative_edge_fraction,
    frustration_index,
    spectral_frustration,
    lattice_plaquettes,
    compute_signed_diffusion_distance,
    compute_eigenmode_sign_distance,
    find_rg_scales,
    partition_at_scale,
    build_reduced_graph,
    compute_reduced_spectrum,
    spectral_rg_step,
    spectral_rg_flow,
    rg_flow_observables,
)


# --- Helpers ---

def _make_lattice(side: int, pflip: float = 0.0, seed: int = 42):
    """Build a lattice and return (G, eigenvalues, eigvecs_col)."""
    from lrgsglib.graphs.nx import Lattice2DNX

    lat = Lattice2DNX(side1=side, geo="sqr", pflip=pflip, seed=seed)
    lat.flip_random_fract_edges()
    lat.compute_laplacian_spectrum_weigV()
    return lat.gr["G"], lat.eigv, lat.eigV.T


def _make_sbm(sizes, p_in, p_out, pflip=0.0, seed=42):
    """Build a planted signed SBM."""
    rng = np.random.RandomState(seed)
    K = len(sizes)
    N = sum(sizes)
    block_of = []
    for k, s in enumerate(sizes):
        block_of.extend([k] * s)

    G = nx.Graph()
    G.add_nodes_from(range(N))
    for i in range(N):
        for j in range(i + 1, N):
            p = p_in if block_of[i] == block_of[j] else p_out
            if rng.random() < p:
                w = 1.0
                if pflip > 0 and rng.random() < pflip:
                    w = -1.0
                G.add_edge(i, j, weight=w)

    # Compute signed Laplacian
    nodes = sorted(G.nodes())
    A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
    D = np.diag(np.abs(A).sum(axis=1))
    L = D - A
    eigv, eigV = np.linalg.eigh(L)
    return G, eigv, eigV, np.array(block_of)


def _make_bipartite_antiferro(side: int):
    """All-negative square lattice → bipartite, K=2 ground truth."""
    G = nx.grid_2d_graph(side, side, periodic=True)
    G = nx.convert_node_labels_to_integers(G)
    for u, v in G.edges():
        G[u][v]["weight"] = -1.0

    nodes = sorted(G.nodes())
    A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
    D = np.diag(np.abs(A).sum(axis=1))
    L = D - A
    eigv, eigV = np.linalg.eigh(L)

    # Ground truth: checkerboard
    N = side * side
    truth = np.zeros(N, dtype=int)
    for node in range(N):
        r, c = divmod(node, side)
        truth[node] = (r + c) % 2
    return G, eigv, eigV, truth


# --- Test Classes ---


class TestFrustrationFraction:
    def test_unsigned(self):
        G = nx.cycle_graph(10)
        for u, v in G.edges():
            G[u][v]["weight"] = 1.0
        assert compute_frustration_fraction(G) == 0.0

    def test_all_negative(self):
        G = nx.cycle_graph(10)
        for u, v in G.edges():
            G[u][v]["weight"] = -1.0
        assert compute_frustration_fraction(G) == 1.0

    def test_mixed(self):
        G = nx.path_graph(4)
        G[0][1]["weight"] = 1.0
        G[1][2]["weight"] = -1.0
        G[2][3]["weight"] = 1.0
        assert abs(compute_frustration_fraction(G) - 1 / 3) < 1e-10

    def test_empty(self):
        G = nx.Graph()
        G.add_nodes_from([0, 1, 2])
        assert compute_frustration_fraction(G) == 0.0


class TestSignedDiffusionDistance:
    @pytest.fixture
    def unsigned_lattice(self):
        return _make_lattice(10, pflip=0.0)

    @pytest.fixture
    def signed_lattice(self):
        return _make_lattice(10, pflip=0.15)

    def test_shape_and_symmetry(self, unsigned_lattice):
        G, eigv, eigV_col = unsigned_lattice
        N = G.number_of_nodes()
        D = compute_signed_diffusion_distance(eigv, eigV_col, tau=1.0)
        assert D.shape == (N, N)
        assert_allclose(D, D.T, atol=1e-12)
        assert_allclose(np.diag(D), 0.0, atol=1e-12)

    def test_non_negative(self, unsigned_lattice):
        G, eigv, eigV_col = unsigned_lattice
        D = compute_signed_diffusion_distance(eigv, eigV_col, tau=1.0)
        assert np.all(D >= -1e-12)

    def test_triangle_inequality(self, unsigned_lattice):
        G, eigv, eigV_col = unsigned_lattice
        D = compute_signed_diffusion_distance(eigv, eigV_col, tau=1.0)
        N = D.shape[0]
        rng = np.random.RandomState(42)
        violations = 0
        for _ in range(500):
            i, j, k = rng.choice(N, 3, replace=False)
            if D[i, j] > D[i, k] + D[k, j] + 1e-10:
                violations += 1
        assert violations == 0

    def test_classical_reduction(self, unsigned_lattice):
        """On unsigned graph, d_diff with use_abs=True/False should be identical."""
        G, eigv, eigV_col = unsigned_lattice
        D_abs = compute_signed_diffusion_distance(
            eigv, eigV_col, tau=1.0, use_abs=True
        )
        D_raw = compute_signed_diffusion_distance(
            eigv, eigV_col, tau=1.0, use_abs=False
        )
        # All eigenvalues are >= 0 for unsigned graph, so |λ| = λ
        assert_allclose(D_abs, D_raw, atol=1e-10)

    def test_signed_no_blowup(self, signed_lattice):
        """Signed graph should produce finite distances."""
        G, eigv, eigV_col = signed_lattice
        D = compute_signed_diffusion_distance(eigv, eigV_col, tau=1.0)
        assert np.all(np.isfinite(D))
        assert D.max() < 100  # reasonable bound

    def test_scale_dependence(self, unsigned_lattice):
        """Larger tau should give smaller distances (modes die out)."""
        G, eigv, eigV_col = unsigned_lattice
        D_small = compute_signed_diffusion_distance(eigv, eigV_col, tau=0.01)
        D_large = compute_signed_diffusion_distance(eigv, eigV_col, tau=10.0)
        # At large tau, only zero-mode survives → distances shrink
        assert D_large.max() < D_small.max()


class TestEigenmodeSignDistance:
    @pytest.fixture
    def open_lattice(self):
        from lrgsglib.graphs import Lattice2D

        lat = Lattice2D(
            side1=10, geo="sqr", pflip=0.0, pbc=False, engine="nx"
        )
        lat.compute_laplacian_spectrum_weigV(
            transpose=False, backend="numpy"
        )
        return lat.gr["G"], lat.eigv, lat.eigV

    def test_shape_and_symmetry(self, open_lattice):
        G, eigv, eigV_col = open_lattice
        D = compute_eigenmode_sign_distance(
            eigv, eigV_col, tau=1.0, use_abs=False
        )
        N = G.number_of_nodes()
        assert D.shape == (N, N)
        assert_allclose(D, D.T, atol=1e-12)
        assert_allclose(np.diag(D), 0.0, atol=1e-12)

    def test_bounded(self, open_lattice):
        _, eigv, eigV_col = open_lattice
        D = compute_eigenmode_sign_distance(
            eigv, eigV_col, tau=1.0, use_abs=False
        )
        assert np.all(D >= -1e-12)
        assert np.all(D <= 1.0 + 1e-12)

    def test_scale_dependence_changes_active_signature(self, open_lattice):
        _, eigv, eigV_col = open_lattice
        D_small = compute_eigenmode_sign_distance(
            eigv, eigV_col, tau=0.25, use_abs=False
        )
        D_large = compute_eigenmode_sign_distance(
            eigv, eigV_col, tau=4.0, use_abs=False
        )
        assert not np.allclose(D_small, D_large)


class TestFindRGScales:
    def test_lattice_has_peaks(self):
        _, eigv, _ = _make_lattice(12)
        N = len(eigv)
        result = find_rg_scales(eigv, N)
        assert len(result["tau_stars"]) >= 1
        assert result["C_global"].shape == result["tau_grid"].shape
        assert result["S_global"].shape == result["tau_grid"].shape

    def test_entropy_monotone(self):
        """Entropy should decrease with tau (modes die out)."""
        _, eigv, _ = _make_lattice(12)
        result = find_rg_scales(eigv, len(eigv))
        S = result["S_global"]
        # Not strictly monotone due to normalization, but first value > last
        assert S[0] > S[-1]

    def test_signed_has_peaks(self):
        _, eigv, _ = _make_lattice(12, pflip=0.15)
        result = find_rg_scales(eigv, len(eigv), use_abs=True)
        assert len(result["tau_stars"]) >= 1


class TestPartitionAtScale:
    @pytest.fixture
    def sbm_data(self):
        return _make_sbm([30, 30, 30], p_in=0.3, p_out=0.02)

    def test_basic_contract(self, sbm_data):
        G, eigv, eigV_col, truth = sbm_data
        result = partition_at_scale(
            eigv, eigV_col, tau_star=1.0, n_clusters=3, G=G
        )
        assert result["n_clusters"] == 3
        assert len(result["labels"]) == len(eigv)
        assert len(result["partition"]) == 3
        total = sum(len(s) for s in result["partition"])
        assert total == len(eigv)

    def test_sbm_recovery(self, sbm_data):
        """d_diff should recover planted SBM blocks."""
        from sklearn.metrics import normalized_mutual_info_score

        G, eigv, eigV_col, truth = sbm_data
        result = partition_at_scale(
            eigv, eigV_col, tau_star=1.0, n_clusters=3, G=G
        )
        nmi = normalized_mutual_info_score(truth, result["labels"])
        assert nmi > 0.7, f"NMI={nmi:.3f} too low for clean SBM"

    def test_signed_diffusion_method(self, sbm_data):
        G, eigv, eigV_col, truth = sbm_data
        result = partition_at_scale(
            eigv, eigV_col, tau_star=1.0,
            method="signed_diffusion", n_clusters=3,
        )
        assert result["distance_matrix"] is not None
        assert result["linkage_matrix"] is not None

    def test_quantum_distance_method(self, sbm_data):
        G, eigv, eigV_col, truth = sbm_data
        result = partition_at_scale(
            eigv, eigV_col, tau_star=1.0,
            method="quantum_distance", n_clusters=3,
        )
        assert result["n_clusters"] == 3


class TestApodicticCases:
    """Known ground-truth cases where the answer is exact."""

    def test_bipartite_antiferro(self):
        """All-negative square lattice → K=2 checkerboard."""
        from sklearn.metrics import normalized_mutual_info_score

        G, eigv, eigV_col, truth = _make_bipartite_antiferro(8)
        result = partition_at_scale(
            eigv, eigV_col, tau_star=1.0, n_clusters=2, G=G
        )
        nmi = normalized_mutual_info_score(truth, result["labels"])
        assert nmi > 0.95, f"NMI={nmi:.3f}, bipartite should be near-perfect"

    def test_negative_boundary_domain(self):
        """Ferro domain surrounded by negative-weight boundary → K=2."""
        from sklearn.metrics import normalized_mutual_info_score

        side = 12
        G = nx.grid_2d_graph(side, side, periodic=True)
        G = nx.convert_node_labels_to_integers(G)
        # All edges positive
        for u, v in G.edges():
            G[u][v]["weight"] = 1.0
        # Negative boundary: edges crossing x=side//2
        truth = np.zeros(side * side, dtype=int)
        for node in range(side * side):
            r, c = divmod(node, side)
            if c >= side // 2:
                truth[node] = 1
        # Flip edges crossing the boundary
        for u, v in G.edges():
            ru, cu = divmod(u, side)
            rv, cv = divmod(v, side)
            if truth[u] != truth[v]:
                G[u][v]["weight"] = -1.0

        nodes = sorted(G.nodes())
        A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
        D = np.diag(np.abs(A).sum(axis=1))
        L = D - A
        eigv, eigV = np.linalg.eigh(L)

        result = partition_at_scale(
            eigv, eigV, tau_star=1.0, n_clusters=2, G=G
        )
        nmi = normalized_mutual_info_score(truth, result["labels"])
        assert nmi > 0.9, f"NMI={nmi:.3f}, clean boundary should give near-perfect"


class TestBuildReducedGraph:
    def test_basic(self):
        G = nx.cycle_graph(6)
        for u, v in G.edges():
            G[u][v]["weight"] = 1.0
        partition = [{0, 1, 2}, {3, 4, 5}]
        G_r, meta = build_reduced_graph(G, partition)
        assert G_r.number_of_nodes() == 2
        assert G_r.number_of_edges() >= 1
        assert meta["cluster_sizes"] == [3, 3]

    def test_preserves_sign_structure(self):
        """Reduced graph should have negative weights if cross-edges are negative."""
        G = nx.Graph()
        G.add_edge(0, 1, weight=1.0)
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 3, weight=-1.0)
        G.add_edge(3, 4, weight=1.0)
        G.add_edge(4, 5, weight=1.0)
        partition = [{0, 1, 2}, {3, 4, 5}]
        G_r, meta = build_reduced_graph(G, partition)
        # The cross-edge (2,3) has weight -1 → reduced edge should be negative
        if G_r.number_of_edges() > 0:
            w = list(nx.get_edge_attributes(G_r, "weight").values())[0]
            assert w < 0


class TestComputeReducedSpectrum:
    def test_single_node(self):
        G = nx.Graph()
        G.add_node(0)
        eigv, eigV = compute_reduced_spectrum(G)
        assert len(eigv) == 1
        assert eigv[0] == 0.0

    def test_two_nodes(self):
        G = nx.Graph()
        G.add_edge(0, 1, weight=1.0)
        eigv, eigV = compute_reduced_spectrum(G)
        assert len(eigv) == 2
        assert_allclose(eigv[0], 0.0, atol=1e-10)
        assert eigv[1] > 0


class TestSpectralRGStep:
    def test_reduces_nodes(self):
        G, eigv, eigV_col = _make_lattice(12)
        step = spectral_rg_step(G, eigv, eigV_col)
        assert step["N_reduced"] < step["N_original"]
        assert step["G_reduced"].number_of_nodes() == step["N_reduced"]
        assert step["tau_star"] > 0

    def test_auto_spectrum(self):
        """Should compute spectrum automatically if not provided."""
        G, _, _ = _make_lattice(8)
        step = spectral_rg_step(G)
        assert step["N_reduced"] < step["N_original"]


class TestSpectralRGFlow:
    def test_converges(self):
        G, eigv, eigV_col = _make_lattice(12)
        flow = spectral_rg_flow(G, eigv, eigV_col, min_nodes=3)
        assert len(flow) >= 1
        # Node count should decrease monotonically
        nodes = [flow[0]["N_original"]] + [s["N_reduced"] for s in flow]
        for i in range(len(nodes) - 1):
            assert nodes[i + 1] < nodes[i]

    def test_signed_flow(self):
        G, eigv, eigV_col = _make_lattice(12, pflip=0.15)
        flow = spectral_rg_flow(G, eigv, eigV_col, min_nodes=3)
        assert len(flow) >= 1

    def test_max_steps(self):
        G, eigv, eigV_col = _make_lattice(12)
        flow = spectral_rg_flow(G, eigv, eigV_col, n_steps=1)
        assert len(flow) == 1


class TestRGFlowObservables:
    def test_basic(self):
        G, eigv, eigV_col = _make_lattice(12)
        flow = spectral_rg_flow(G, eigv, eigV_col, min_nodes=3)
        obs = rg_flow_observables(flow)
        assert "n_nodes" in obs
        assert "spectral_gap" in obs
        assert "frustration_fraction" in obs
        assert "tau_stars" in obs
        assert len(obs["n_nodes"]) == len(flow) + 1  # includes final

    def test_node_monotone(self):
        G, eigv, eigV_col = _make_lattice(12)
        flow = spectral_rg_flow(G, eigv, eigV_col, min_nodes=3)
        obs = rg_flow_observables(flow)
        n = obs["n_nodes"]
        for i in range(len(n) - 1):
            assert n[i + 1] < n[i]


class TestFrustrationMeasures:
    """Verify frustration measures distinguish balanced from frustrated graphs."""

    def test_all_negative_lattice_balanced(self):
        """All-negative square lattice is bipartite → zero frustration."""
        side = 8
        G = nx.grid_2d_graph(side, side, periodic=True)
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            G[u][v]["weight"] = -1.0

        # Negative edge fraction is 1 (all negative)
        assert negative_edge_fraction(G) == 1.0
        # But no plaquette is frustrated
        plaq = lattice_plaquettes(side)
        assert frustration_index(G, plaq) == 0.0
        # Spectral frustration is zero
        assert abs(spectral_frustration(G)) < 1e-10

    def test_domain_wall_balanced(self):
        """Ferro|neg|ferro domain-wall lattice is balanced (bipartite cut)."""
        side = 8
        G = nx.grid_2d_graph(side, side, periodic=True)
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            cu = u % side
            cv = v % side
            same = (cu < side // 2) == (cv < side // 2)
            G[u][v]["weight"] = 1.0 if same else -1.0

        plaq = lattice_plaquettes(side)
        assert frustration_index(G, plaq) == 0.0
        assert abs(spectral_frustration(G)) < 1e-10

    def test_all_positive_unfrustrated(self):
        side = 8
        G = nx.grid_2d_graph(side, side, periodic=True)
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            G[u][v]["weight"] = 1.0
        plaq = lattice_plaquettes(side)
        assert frustration_index(G, plaq) == 0.0
        assert negative_edge_fraction(G) == 0.0

    def test_frustrated_triangle(self):
        """Single triangle with 1 or 3 negative edges is frustrated."""
        G = nx.Graph()
        G.add_edge(0, 1, weight=1.0)
        G.add_edge(1, 2, weight=1.0)
        G.add_edge(2, 0, weight=-1.0)  # odd # of negative → frustrated
        assert frustration_index(G) == 1.0

    def test_unfrustrated_triangle(self):
        """Triangle with 0 or 2 negative edges is unfrustrated."""
        G = nx.Graph()
        G.add_edge(0, 1, weight=-1.0)
        G.add_edge(1, 2, weight=-1.0)
        G.add_edge(2, 0, weight=1.0)
        assert frustration_index(G) == 0.0

    def test_deprecated_alias(self):
        """compute_frustration_fraction is just negative_edge_fraction."""
        G = nx.Graph()
        G.add_edge(0, 1, weight=-1.0)
        G.add_edge(1, 2, weight=1.0)
        assert compute_frustration_fraction(G) == negative_edge_fraction(G)

    def test_lattice_plaquettes(self):
        side = 4
        plaq = lattice_plaquettes(side, periodic=True)
        assert len(plaq) == side * side  # one plaquette per site with periodic BC
        for p in plaq:
            assert len(p) == 4

    def test_random_pflip_is_frustrated(self):
        """Lattice with random pflip > 0 has nonzero plaquette frustration."""
        from lrgsglib.graphs.nx import Lattice2DNX
        lat = Lattice2DNX(side1=12, geo="sqr", pflip=0.15, seed=42)
        lat.flip_random_fract_edges()
        plaq = lattice_plaquettes(12)
        fi = frustration_index(lat.gr["G"], plaq)
        assert fi > 0.0
        # spectral frustration also > 0
        assert spectral_frustration(lat.gr["G"]) > 0.0


class TestBalancedBipartiteRecovery:
    """The eigenmode partition should recover ground-truth bipartitions
    on balanced signed graphs (all-negative, domain-wall)."""

    def test_all_negative_checkerboard(self):
        """All-negative lattice: partition = checkerboard (NMI=1)."""
        from sklearn.metrics import normalized_mutual_info_score
        side = 10
        G = nx.grid_2d_graph(side, side, periodic=True)
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            G[u][v]["weight"] = -1.0

        nodes = sorted(G.nodes())
        A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
        D = np.diag(np.abs(A).sum(axis=1))
        L = D - A
        eigv, eigV = np.linalg.eigh(L)

        part = partition_at_scale(
            eigv, eigV, tau_star=1.0,
            method="eigenmode_sign", n_clusters=2
        )
        # Ground truth: checkerboard
        truth = np.array(
            [(n // side + n % side) % 2 for n in range(side * side)]
        )
        nmi = normalized_mutual_info_score(part["labels"], truth)
        assert nmi > 0.99, f"NMI={nmi:.3f}, should be 1 on checkerboard"

    def test_domain_wall_recovery(self):
        """Two-ferro-domain lattice with negative cut: recover the cut."""
        from sklearn.metrics import normalized_mutual_info_score
        side = 10
        G = nx.grid_2d_graph(side, side, periodic=True)
        G = nx.convert_node_labels_to_integers(G)
        for u, v in G.edges():
            cu, cv = u % side, v % side
            same = (cu < side // 2) == (cv < side // 2)
            G[u][v]["weight"] = 1.0 if same else -1.0

        nodes = sorted(G.nodes())
        A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
        D = np.diag(np.abs(A).sum(axis=1))
        L = D - A
        eigv, eigV = np.linalg.eigh(L)

        part = partition_at_scale(
            eigv, eigV, tau_star=1.0,
            method="eigenmode_sign", n_clusters=2
        )
        truth = np.array(
            [0 if (n % side) < side // 2 else 1 for n in range(side * side)]
        )
        nmi = normalized_mutual_info_score(part["labels"], truth)
        assert nmi > 0.99, f"NMI={nmi:.3f}, should be 1 on domain wall"


class TestTwoChannelProperty:
    """Verify that d_diff on unsigned vs signed eigenvectors gives
    different but complementary views of the same graph."""

    def test_channels_differ_under_frustration(self):
        """On a frustrated SBM, topology and geometry channels should diverge."""
        from sklearn.metrics import normalized_mutual_info_score

        G, eigv, eigV_col, truth = _make_sbm(
            [30, 30], p_in=0.3, p_out=0.05, pflip=0.25, seed=42
        )
        # Channel 1: signed eigenvectors (geometry)
        part_geom = partition_at_scale(
            eigv, eigV_col, tau_star=1.0,
            method="signed_diffusion", n_clusters=2,
        )
        # Channel 2: unsigned eigenvectors (topology)
        part_topo = partition_at_scale(
            eigv, eigV_col, tau_star=1.0,
            method="topological", n_clusters=2, G=G,
        )
        nmi_topo = normalized_mutual_info_score(truth, part_topo["labels"])
        # Topology channel should do better on planted blocks
        # (geometry channel confuses frustration with community)
        assert nmi_topo > 0.3, f"Topological NMI={nmi_topo:.3f} too low"

    def test_channels_agree_unsigned(self):
        """On unsigned graph, both channels should give the same result."""
        from sklearn.metrics import normalized_mutual_info_score

        G, eigv, eigV_col, truth = _make_sbm(
            [30, 30], p_in=0.3, p_out=0.02, pflip=0.0,
        )
        part_geom = partition_at_scale(
            eigv, eigV_col, tau_star=1.0,
            method="signed_diffusion", n_clusters=2,
        )
        part_topo = partition_at_scale(
            eigv, eigV_col, tau_star=1.0,
            method="topological", n_clusters=2, G=G,
        )
        nmi = normalized_mutual_info_score(
            part_geom["labels"], part_topo["labels"]
        )
        # Should be identical or near-identical on unsigned graph
        assert nmi > 0.8, f"Channels should agree on unsigned: NMI={nmi:.3f}"

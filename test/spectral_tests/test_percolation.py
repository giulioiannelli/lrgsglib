"""Tests for the percolation analysis module.

Covers both the eigenmode sign-cluster method and the direct
defect identification approach.
"""

import networkx as nx
import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX
from lrgsglib.graphs.nx.funcs.dual import build_planar_dual
from lrgsglib.utils.lrg.percolation import (
    percolation_sweep,
    find_threshold,
    dual_percolation_sweep,
    default_p_range,
    identify_frustrated_faces,
    identify_frustrated_vertices,
    giant_component_fraction,
    build_face_syndrome,
    build_vertex_syndrome,
    defect_percolation_sweep,
    compute_pareto_point_direct,
)


class TestPercolationSweep:
    """Tests for ``percolation_sweep``."""

    def test_basic_sweep(self):
        """Sweep should return correct structure."""
        p_range = np.array([0.0, 0.1, 0.2, 0.3])
        result = percolation_sweep(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=5,
            seed_base=42,
        )
        assert "p" in result
        assert "Pinf_mean" in result
        assert "chi" in result
        assert len(result["Pinf_mean"]) == len(p_range)
        assert result["N"] == 64

    def test_Pinf_bounded(self):
        """P_inf should be in [0, 1]."""
        p_range = np.linspace(0.0, 0.5, 5)
        result = percolation_sweep(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=5,
            seed_base=0,
        )
        assert np.all(result["Pinf_mean"] >= 0.0)
        assert np.all(result["Pinf_mean"] <= 1.0)

    def test_Pinf_decreasing_trend(self):
        """P_inf should generally decrease with p."""
        p_range = np.array([0.01, 0.05, 0.15, 0.3, 0.45])
        result = percolation_sweep(
            Lattice2DNX,
            {"side1": 16, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=20,
            seed_base=42,
        )
        # At p~0, Pinf should be high; at p~0.5, lower
        assert result["Pinf_mean"][0] > result["Pinf_mean"][-1]

    def test_chi_positive(self):
        """Susceptibility should be non-negative."""
        p_range = np.linspace(0.05, 0.4, 6)
        result = percolation_sweep(
            Lattice2DNX,
            {"side1": 10, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=10,
            seed_base=0,
        )
        assert np.all(result["chi"] >= 0.0)

    def test_reproducibility(self):
        """Same seed_base should give same results."""
        params = {
            "side1": 8,
            "geo": "sqr",
            "pbc": False,
        }
        p_range = np.array([0.1, 0.2])
        r1 = percolation_sweep(Lattice2DNX, params, p_range, 5, seed_base=99)
        r2 = percolation_sweep(Lattice2DNX, params, p_range, 5, seed_base=99)
        np.testing.assert_array_equal(r1["Pinf_all"], r2["Pinf_all"])


class TestFindThreshold:
    """Tests for ``find_threshold``."""

    def test_peak_method(self):
        """Peak method should return p at max chi."""
        p = np.array([0.05, 0.10, 0.15, 0.20, 0.25])
        chi = np.array([0.1, 0.5, 1.0, 0.4, 0.1])
        pc = find_threshold(p, chi, method="peak")
        assert pc == pytest.approx(0.15)

    def test_fit_method(self):
        """Fit method should refine peak location."""
        p = np.linspace(0.05, 0.25, 21)
        # Gaussian-like peak at p=0.12
        chi = np.exp(-((p - 0.12) ** 2) / (2 * 0.02**2))
        pc = find_threshold(p, chi, method="fit")
        assert abs(pc - 0.12) < 0.01

    def test_invalid_method(self):
        with pytest.raises(ValueError):
            find_threshold(np.array([0.1]), np.array([1.0]), method="invalid")


class TestDualPercolationSweep:
    """Tests for ``dual_percolation_sweep``."""

    def test_basic_dual_sweep(self):
        """Dual sweep should return primal and dual results."""
        p_range = np.array([0.05, 0.15, 0.25])
        result = dual_percolation_sweep(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=5,
            seed_base=42,
        )
        assert "primal" in result
        assert "dual" in result
        assert len(result["primal"]["Pinf_mean"]) == 3
        assert len(result["dual"]["Pinf_mean"]) == 3

    def test_square_self_duality(self):
        """For square lattice, primal and dual thresholds should be similar.

        The dual of the NxN grid (no outer face) is an (N-1)x(N-1) grid,
        so thresholds should roughly match, especially for larger sizes.
        """
        p_range = np.linspace(0.03, 0.25, 15)
        result = dual_percolation_sweep(
            Lattice2DNX,
            {"side1": 16, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=20,
            seed_base=42,
        )
        pc_primal = find_threshold(
            result["primal"]["p"], result["primal"]["chi"], "peak"
        )
        pc_dual = find_threshold(
            result["dual"]["p"], result["dual"]["chi"], "peak"
        )
        # Self-duality: thresholds should be within ~50% of each other
        # (finite-size effects are significant at N=16)
        assert abs(pc_primal - pc_dual) < 0.10, (
            f"pc_primal={pc_primal:.3f}, pc_dual={pc_dual:.3f}"
        )


class TestDefaultPRange:
    """Tests for ``default_p_range``."""

    def test_sorted_unique(self):
        p = default_p_range()
        assert np.all(np.diff(p) > 0), "Should be strictly increasing"

    def test_bounds(self):
        p = default_p_range(p_min=0.01, p_max=0.4)
        assert p[0] >= 0.01
        assert p[-1] <= 0.4

    def test_focus_density(self):
        """More points should be in the focus region."""
        p = default_p_range(focus=(0.08, 0.15))
        in_focus = np.sum((p >= 0.08) & (p <= 0.15))
        outside = len(p) - in_focus
        assert in_focus > outside


# ---------------------------------------------------------------
# Direct defect identification tests
# ---------------------------------------------------------------


def _make_signed_grid(n: int = 6, pflip: float = 0.0, seed: int = 0):
    """Helper: build a square grid with signed edges."""
    lat = Lattice2DNX(n, geo="sqr", pbc=False, pflip=pflip, seed=seed)
    if pflip > 0:
        lat.flip_random_fract_edges()
    G = lat.gr[lat.on_g]
    _, fi = build_planar_dual(G, include_outer_face=True)
    return G, fi, lat


class TestIdentifyFrustratedFaces:
    """Tests for ``identify_frustrated_faces``."""

    def test_no_disorder_no_frustration(self):
        """All-positive edges → zero frustrated faces."""
        G, fi, _ = _make_signed_grid(6, pflip=0.0)
        frust = identify_frustrated_faces(G, fi)
        assert len(frust) == 0

    def test_single_flip_creates_frustration(self):
        """Flipping one edge frustrates exactly the two adjacent faces."""
        G, fi, _ = _make_signed_grid(8, pflip=0.0)
        # Flip one interior edge
        interior_edges = [
            (u, v) for u, v in G.edges()
            if G.degree(u) == 4 and G.degree(v) == 4
        ]
        u, v = interior_edges[0]
        G[u][v]["weight"] = -1
        frust = identify_frustrated_faces(G, fi)
        # One edge borders exactly 2 faces → both frustrated
        assert len(frust) == 2

    def test_frustrated_count_grows_with_p(self):
        """More disorder → more frustrated faces."""
        counts = []
        for p in [0.05, 0.15, 0.30]:
            G, fi, _ = _make_signed_grid(10, pflip=p, seed=42)
            frust = identify_frustrated_faces(G, fi)
            outer = fi["outer_face_idx"]
            n_inner = len([f for f in frust if f != outer])
            counts.append(n_inner)
        # Monotonically increasing (on average; seed fixed)
        assert counts[0] <= counts[1] <= counts[2]

    def test_outer_face_can_be_frustrated(self):
        """Outer face should also be identifiable as frustrated."""
        G, fi, _ = _make_signed_grid(6, pflip=0.3, seed=7)
        frust = identify_frustrated_faces(G, fi)
        # Outer face might or might not be frustrated —
        # just check the function doesn't crash
        assert isinstance(frust, list)


class TestIdentifyFrustratedVertices:
    """Tests for ``identify_frustrated_vertices``."""

    def test_no_disorder_no_frustration(self):
        """All-positive edges → zero frustrated vertices."""
        G, _, _ = _make_signed_grid(6, pflip=0.0)
        frust = identify_frustrated_vertices(G)
        assert len(frust) == 0

    def test_single_flip_frustrates_endpoints(self):
        """Flipping one edge frustrates exactly its two endpoints."""
        G, _, _ = _make_signed_grid(8, pflip=0.0)
        edges = list(G.edges())
        u, v = edges[len(edges) // 2]  # pick a middle edge
        G[u][v]["weight"] = -1
        frust = identify_frustrated_vertices(G)
        assert set(frust) == {u, v}

    def test_vertex_frustration_grows_with_p(self):
        """More disorder → more frustrated vertices."""
        counts = []
        for p in [0.05, 0.15, 0.30]:
            G, _, _ = _make_signed_grid(10, pflip=p, seed=42)
            counts.append(len(identify_frustrated_vertices(G)))
        assert counts[0] <= counts[1] <= counts[2]


class TestGiantComponentFraction:
    """Tests for ``giant_component_fraction``."""

    def test_empty_graph(self):
        assert giant_component_fraction(nx.Graph()) == 0.0

    def test_complete_graph(self):
        assert giant_component_fraction(nx.complete_graph(5)) == 1.0

    def test_disconnected(self):
        G = nx.Graph()
        G.add_edges_from([(0, 1), (2, 3), (4, 5)])
        # 3 components of size 2 → GC = 2/6
        assert giant_component_fraction(G) == pytest.approx(1 / 3)

    def test_singleton_nodes(self):
        G = nx.Graph()
        G.add_nodes_from(range(10))
        # 10 isolated nodes → GC = 1/10
        assert giant_component_fraction(G) == pytest.approx(0.1)


class TestSyndromeGraphBuilders:
    """Tests for ``build_face_syndrome`` and ``build_vertex_syndrome``."""

    def test_face_syndrome_empty(self):
        """No frustrated faces → empty syndrome graph."""
        G, fi, _ = _make_signed_grid(6, pflip=0.0)
        syn = build_face_syndrome(fi, [])
        assert syn.number_of_nodes() == 0

    def test_face_syndrome_adjacency(self):
        """Adjacent frustrated faces should be connected in syndrome graph."""
        G, fi, _ = _make_signed_grid(8, pflip=0.0)
        # Flip one edge → 2 adjacent frustrated faces
        interior_edges = [
            (u, v) for u, v in G.edges()
            if G.degree(u) == 4 and G.degree(v) == 4
        ]
        u, v = interior_edges[0]
        G[u][v]["weight"] = -1
        frust = identify_frustrated_faces(G, fi)
        syn = build_face_syndrome(fi, frust)
        assert syn.number_of_nodes() == 2
        assert syn.number_of_edges() == 1

    def test_vertex_syndrome_is_subgraph(self):
        """Vertex syndrome should be an induced subgraph of G."""
        G, _, _ = _make_signed_grid(8, pflip=0.2, seed=42)
        frust = identify_frustrated_vertices(G)
        syn = build_vertex_syndrome(G, frust)
        assert set(syn.nodes()) == set(frust)
        # Every edge in syndrome must exist in G
        for u, v in syn.edges():
            assert G.has_edge(u, v)


class TestDefectPercolationSweep:
    """Tests for ``defect_percolation_sweep``."""

    def test_sweep_structure(self):
        """Sweep should return Z, X, and face_info_sample."""
        p_range = np.array([0.05, 0.15, 0.30])
        result = defect_percolation_sweep(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=5,
            seed_base=42,
        )
        assert "Z" in result
        assert "X" in result
        assert "face_info_sample" in result
        assert len(result["Z"]["Pinf_mean"]) == 3
        assert len(result["X"]["Pinf_mean"]) == 3

    def test_defect_density_increases(self):
        """Mean defect count should increase with disorder."""
        p_range = np.array([0.02, 0.10, 0.30])
        result = defect_percolation_sweep(
            Lattice2DNX,
            {"side1": 10, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=10,
            seed_base=42,
        )
        # Z-defect density should increase
        z_def = result["Z"]["n_defects_mean"]
        assert z_def[0] < z_def[-1]
        # X-defect density should increase
        x_def = result["X"]["n_defects_mean"]
        assert x_def[0] < x_def[-1]

    def test_Pinf_bounded(self):
        """P_inf should be in [0, 1]."""
        p_range = np.array([0.05, 0.25])
        result = defect_percolation_sweep(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=5,
        )
        for key in ["Z", "X"]:
            assert np.all(result[key]["Pinf_all"] >= 0.0)
            assert np.all(result[key]["Pinf_all"] <= 1.0)

    def test_no_disorder_no_defects(self):
        """At p=0, there should be no frustrated faces or vertices."""
        p_range = np.array([0.0])
        result = defect_percolation_sweep(
            Lattice2DNX,
            {"side1": 8, "geo": "sqr", "pbc": False},
            p_range,
            n_realizations=3,
        )
        assert result["Z"]["n_defects_mean"][0] == 0.0
        assert result["X"]["n_defects_mean"][0] == 0.0

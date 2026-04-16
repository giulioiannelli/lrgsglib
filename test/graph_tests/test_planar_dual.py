"""Tests for the planar dual graph constructor."""

import networkx as nx
import numpy as np
import pytest

from lrgsglib.graphs.nx.funcs.dual import build_planar_dual, build_signed_dual
from lrgsglib.graphs.nx import Lattice2DNX, SignedGraphNX
from lrgsglib.graphs.nx.SierpinskiNX.SierpinskiNX import SierpinskiNX


# ---------------------------------------------------------------------------
# build_planar_dual
# ---------------------------------------------------------------------------

class TestBuildPlanarDual:
    """Tests for ``build_planar_dual``."""

    def test_euler_formula_square_grid(self):
        """V - E + F = 2 for a square grid."""
        for n in [3, 5, 8]:
            G = nx.grid_2d_graph(n, n)
            G = nx.convert_node_labels_to_integers(G)
            nx.set_edge_attributes(G, 1, "weight")
            _, info = build_planar_dual(G)
            V, E, F = G.number_of_nodes(), G.number_of_edges(), len(info["faces"])
            assert V - E + F == 2, f"Euler failed for {n}x{n} grid"

    def test_euler_formula_hexagonal(self):
        """V - E + F = 2 for a hexagonal lattice."""
        H = nx.hexagonal_lattice_graph(4, 4)
        H = nx.convert_node_labels_to_integers(H)
        nx.set_edge_attributes(H, 1, "weight")
        _, info = build_planar_dual(H)
        V, E, F = H.number_of_nodes(), H.number_of_edges(), len(info["faces"])
        assert V - E + F == 2

    def test_euler_formula_triangular(self):
        """V - E + F = 2 for a triangular lattice."""
        T = nx.triangular_lattice_graph(5, 5)
        T = nx.convert_node_labels_to_integers(T)
        nx.set_edge_attributes(T, 1, "weight")
        _, info = build_planar_dual(T)
        V, E, F = T.number_of_nodes(), T.number_of_edges(), len(info["faces"])
        assert V - E + F == 2

    def test_self_duality_square_grid(self):
        """Dual of NxN grid (no outer face) is (N-1)x(N-1) grid."""
        for n in [4, 6, 8]:
            G = nx.grid_2d_graph(n, n)
            G = nx.convert_node_labels_to_integers(G)
            nx.set_edge_attributes(G, 1, "weight")
            dual, _ = build_planar_dual(G, include_outer_face=False)

            expected_v = (n - 1) ** 2
            expected_e = 2 * (n - 1) * (n - 2)
            assert dual.number_of_nodes() == expected_v
            assert dual.number_of_edges() == expected_e

    def test_outer_face_identification(self):
        """Outer face should be the one with the largest perimeter."""
        G = nx.grid_2d_graph(5, 5)
        G = nx.convert_node_labels_to_integers(G)
        nx.set_edge_attributes(G, 1, "weight")
        _, info = build_planar_dual(G)

        outer_idx = info["outer_face_idx"]
        outer_size = info["face_sizes"][outer_idx]
        assert outer_size == max(info["face_sizes"])

    def test_exclude_outer_face(self):
        """Excluding outer face reduces dual node count by 1."""
        G = nx.grid_2d_graph(4, 4)
        G = nx.convert_node_labels_to_integers(G)
        nx.set_edge_attributes(G, 1, "weight")

        dual_with, _ = build_planar_dual(G, include_outer_face=True)
        dual_without, _ = build_planar_dual(G, include_outer_face=False)
        assert dual_with.number_of_nodes() - dual_without.number_of_nodes() == 1

    def test_sign_transfer(self):
        """Dual edges inherit the sign of the crossed primal edge."""
        G = nx.grid_2d_graph(4, 4)
        G = nx.convert_node_labels_to_integers(G)
        nx.set_edge_attributes(G, 1, "weight")

        # Flip specific edges
        G[0][1]["weight"] = -1
        G[4][5]["weight"] = -1
        G[8][9]["weight"] = -1

        dual, _ = build_planar_dual(G)
        neg_primal = sum(
            1 for _, _, d in G.edges(data=True) if d.get("weight", 1) < 0
        )
        neg_dual = sum(
            1 for _, _, d in dual.edges(data=True) if d.get("weight", 1) < 0
        )
        # Each negative primal edge creates exactly one negative dual edge
        assert neg_dual == neg_primal

    def test_nonplanar_raises(self):
        """Non-planar graph should raise."""
        K5 = nx.complete_graph(5)
        with pytest.raises(nx.NetworkXException, match="not planar"):
            build_planar_dual(K5)

    def test_sierpinski_carpet_multiscale_faces(self):
        """Carpet gen 2 should have faces at multiple scales."""
        carpet = SierpinskiNX(2, variant="carpet", pflip=0.0)
        G = carpet.gr[carpet.on_g]
        _, info = build_planar_dual(G)

        face_sizes = set(info["face_sizes"])
        # Carpet gen 2 has unit squares (4), gen-1 holes (8),
        # center hole (16), and outer face (32)
        assert 4 in face_sizes, "Missing unit-square faces"
        assert len(face_sizes) > 1, "Should have multi-scale faces"

        # Euler check
        V, E, F = G.number_of_nodes(), G.number_of_edges(), len(info["faces"])
        assert V - E + F == 2

    def test_dual_edge_count_matches_primal(self):
        """Including outer face: dual edges <= primal edges.

        Equality holds when no two faces share more than one edge.
        For the interior-only dual, this is always the case.
        """
        G = nx.grid_2d_graph(6, 6)
        G = nx.convert_node_labels_to_integers(G)
        nx.set_edge_attributes(G, 1, "weight")
        dual, _ = build_planar_dual(G, include_outer_face=False)
        # Interior faces of NxN grid: each pair shares exactly 1 edge
        assert dual.number_of_edges() <= G.number_of_edges()


# ---------------------------------------------------------------------------
# build_signed_dual
# ---------------------------------------------------------------------------

class TestBuildSignedDual:
    """Tests for ``build_signed_dual``."""

    def test_basic_construction(self):
        """Signed dual should be a valid SignedGraphNX."""
        lat = Lattice2DNX(8, geo="sqr", pflip=0.15, pbc=False, seed=42)
        lat.flip_random_fract_edges()

        dual_sg = build_signed_dual(lat, include_outer_face=False)
        assert isinstance(dual_sg, SignedGraphNX)
        assert dual_sg.N > 0
        assert dual_sg.Ne > 0

    def test_sign_consistency(self):
        """Dual should inherit negative-edge count from primal."""
        lat = Lattice2DNX(6, geo="sqr", pflip=0.2, pbc=False, seed=7)
        lat.flip_random_fract_edges()

        G = lat.gr[lat.on_g]
        neg_primal = sum(
            1 for _, _, d in G.edges(data=True) if d.get("weight", 1) < 0
        )

        dual_sg = build_signed_dual(lat, include_outer_face=False)
        neg_dual = len(dual_sg.fleset[dual_sg.on_g])
        # Every negative primal edge that borders 2 interior faces
        # creates one negative dual edge. Some boundary neg edges may
        # be lost when excluding the outer face.
        assert neg_dual <= neg_primal
        assert neg_dual > 0 or neg_primal == 0

    def test_spectral_pipeline(self):
        """Full spectral pipeline should work on dual."""
        lat = Lattice2DNX(8, geo="sqr", pflip=0.1, pbc=False, seed=99)
        lat.flip_random_fract_edges()

        dual_sg = build_signed_dual(lat, include_outer_face=False)
        dual_sg.compute_k_eigvV(k=1, backend="scipy")
        dual_sg.compute_pinf(which=0)

        assert 0.0 <= dual_sg.Pinf <= 1.0

    def test_sierpinski_carpet_dual(self):
        """Dual of Sierpinski carpet should work."""
        carpet = SierpinskiNX(2, variant="carpet", pflip=0.1, seed=42)
        carpet.flip_random_fract_edges()

        dual_sg = build_signed_dual(carpet, include_outer_face=False)
        assert dual_sg.N > 0
        dual_sg.compute_k_eigvV(k=1, backend="scipy")
        dual_sg.compute_pinf(which=0)
        assert 0.0 <= dual_sg.Pinf <= 1.0

    def test_provenance_metadata(self):
        """Dual should store face info and reference to primal."""
        lat = Lattice2DNX(6, geo="sqr", pflip=0.0, pbc=False)
        dual_sg = build_signed_dual(lat, include_outer_face=False)

        assert hasattr(dual_sg, "_dual_face_info")
        assert "faces" in dual_sg._dual_face_info
        assert hasattr(dual_sg, "_dual_of")
        assert dual_sg._dual_of is lat

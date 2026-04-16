"""Tests for Kagome and tri-hexagonal lattice geometries."""

import networkx as nx
import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX
from lrgsglib.graphs.nx.funcs.dual import build_planar_dual, build_signed_dual


class TestKagome:
    """Tests for ``Lattice2DNX(geo='kgm')``."""

    def test_construction(self):
        lat = Lattice2DNX(6, geo="kgm", pbc=False)
        assert lat.N > 0
        assert lat.Ne > 0
        assert lat.z == 4

    def test_coordination(self):
        """Bulk coordination should approach 4."""
        lat = Lattice2DNX(12, geo="kgm", pbc=False)
        avg_deg = np.mean([d for _, d in lat.gr[lat.on_g].degree()])
        assert 3.5 < avg_deg < 4.0

    def test_planarity(self):
        lat = Lattice2DNX(8, geo="kgm", pbc=False)
        assert nx.is_planar(lat.gr[lat.on_g])

    def test_spectral_pipeline(self):
        lat = Lattice2DNX(8, geo="kgm", pbc=False, pflip=0.1, seed=42)
        lat.flip_random_fract_edges()
        lat.compute_k_eigvV(k=1, backend="scipy")
        lat.compute_pinf(which=0)
        assert 0.0 <= lat.Pinf <= 1.0

    def test_dual(self):
        lat = Lattice2DNX(8, geo="kgm", pbc=False, pflip=0.1, seed=42)
        lat.flip_random_fract_edges()
        dual = build_signed_dual(lat, include_outer_face=False)
        assert dual.N > 0
        dual.compute_k_eigvV(k=1, backend="scipy")
        dual.compute_pinf(which=0)
        assert 0.0 <= dual.Pinf <= 1.0

    def test_euler(self):
        lat = Lattice2DNX(8, geo="kgm", pbc=False)
        G = lat.gr[lat.on_g]
        _, info = build_planar_dual(G)
        V, E, F = G.number_of_nodes(), G.number_of_edges(), len(info["faces"])
        assert V - E + F == 2


class TestTriHexa:
    """Tests for ``Lattice2DNX(geo='tri_hex')``."""

    def test_construction(self):
        lat = Lattice2DNX(6, geo="tri_hex", pbc=False)
        assert lat.N > 0
        assert lat.Ne > 0
        assert lat.z == 3

    def test_coordination(self):
        """Bulk coordination should approach 3."""
        lat = Lattice2DNX(12, geo="tri_hex", pbc=False)
        avg_deg = np.mean([d for _, d in lat.gr[lat.on_g].degree()])
        assert 2.5 < avg_deg < 3.0

    def test_planarity(self):
        lat = Lattice2DNX(8, geo="tri_hex", pbc=False)
        assert nx.is_planar(lat.gr[lat.on_g])

    def test_spectral_pipeline(self):
        lat = Lattice2DNX(8, geo="tri_hex", pbc=False, pflip=0.1, seed=42)
        lat.flip_random_fract_edges()
        lat.compute_k_eigvV(k=1, backend="scipy")
        lat.compute_pinf(which=0)
        assert 0.0 <= lat.Pinf <= 1.0

    def test_dual(self):
        lat = Lattice2DNX(8, geo="tri_hex", pbc=False, pflip=0.1, seed=42)
        lat.flip_random_fract_edges()
        dual = build_signed_dual(lat, include_outer_face=False)
        assert dual.N > 0
        dual.compute_k_eigvV(k=1, backend="scipy")
        dual.compute_pinf(which=0)
        assert 0.0 <= dual.Pinf <= 1.0

    def test_pbc_raises(self):
        with pytest.raises(NotImplementedError):
            Lattice2DNX(6, geo="tri_hex", pbc=True)

    def test_euler(self):
        lat = Lattice2DNX(8, geo="tri_hex", pbc=False)
        G = lat.gr[lat.on_g]
        _, info = build_planar_dual(G)
        V, E, F = G.number_of_nodes(), G.number_of_edges(), len(info["faces"])
        assert V - E + F == 2

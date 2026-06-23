"""Native graph-tool 2D lattice geometries and their parity with NetworkX.

``Lattice2DGT`` builds every geometry through the same NetworkX-free generators
in ``_generators.py`` (square, triangular, hexagonal, rhomb-octagonal, kagome,
tri-hexagonal) plus small-world rewiring via ``prew``.  These tests assert that
the GT engine reproduces the *exact* topology of the matching ``Lattice2DNX``
lattice -- same node/edge count, degree sequence and adjacency spectrum -- and
resolves sides / periodic defaults identically, so the two engines are
interchangeable through the ``Lattice2D`` factory.
"""

from collections import Counter

import numpy as np
import pytest

pytest.importorskip("graph_tool")

from lrgsglib.graphs import Lattice2D
from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT

# (geo, periodic) cases with exact NX<->GT topological parity.
# tri_hex has no PBC implementation (matches NX), so only the open case.
PARITY_CASES = [
    ("sqr", False),
    ("sqr", True),
    ("tri", False),
    ("tri", True),
    ("hex", False),
    ("hex", True),
    ("oct_sqr", False),
    ("oct_sqr", True),
    ("kgm", False),
    ("kgm", True),
    ("tri_hex", False),
]

# Expected bulk coordination per geometry.
Z_BY_GEO = {"sqr": 4, "tri": 6, "hex": 3, "oct_sqr": 3, "kgm": 4, "tri_hex": 3}

# Even side keeps periodic hex/kagome valid (even rescaled sides / even n).
PARITY_SIDE = 12


def _spectrum(adjacency):
    """Sorted, rounded adjacency eigenvalues (isomorphism invariant)."""
    return np.round(np.sort(np.linalg.eigvalsh(np.asarray(adjacency))), 6)


def _degree_hist_nx(lat):
    return dict(sorted(Counter(d for _, d in lat.gr["G"].degree()).items()))


def _degree_hist_gt(lat):
    degs = (int(v.out_degree()) for v in lat.G.vertices())
    return dict(sorted(Counter(degs).items()))


def _make_pair(geo, periodic, side1=PARITY_SIDE, side2=None):
    """Matched GT/NX lattices from the *same* side spec (both resolve sides)."""
    gt_lat = Lattice2DGT(side1=side1, side2=side2, geo=geo, periodic=periodic)
    nx_kw = {} if side2 is None else {"side2": side2}
    nx_lat = Lattice2DNX(side1=side1, geo=geo, pbc=periodic, **nx_kw)
    return gt_lat, nx_lat


class TestNativeGeometryParity:
    """GT native generators must match NX topology exactly."""

    @pytest.mark.parametrize("geo,periodic", PARITY_CASES)
    def test_node_and_edge_count(self, geo, periodic):
        gt_lat, nx_lat = _make_pair(geo, periodic)
        assert gt_lat.N == nx_lat.gr["G"].number_of_nodes()
        assert gt_lat.num_edges == nx_lat.gr["G"].number_of_edges()

    @pytest.mark.parametrize("geo,periodic", PARITY_CASES)
    def test_degree_sequence(self, geo, periodic):
        gt_lat, nx_lat = _make_pair(geo, periodic)
        assert _degree_hist_gt(gt_lat) == _degree_hist_nx(nx_lat)

    @pytest.mark.parametrize("geo,periodic", PARITY_CASES)
    def test_adjacency_spectrum(self, geo, periodic):
        import networkx as nx

        gt_lat, nx_lat = _make_pair(geo, periodic)
        gt_spec = _spectrum(gt_lat.get_adjacency_matrix())
        nx_spec = _spectrum(nx.to_numpy_array(nx_lat.gr["G"], weight=None))
        np.testing.assert_allclose(gt_spec, nx_spec, atol=1e-6)

    @pytest.mark.parametrize(
        "geo", ["sqr", "tri", "hex", "oct_sqr", "kgm", "tri_hex"]
    )
    def test_rectangular_parity(self, geo):
        """Non-square lattices match too (explicit side2 on both engines)."""
        import networkx as nx

        gt_lat, nx_lat = _make_pair(geo, periodic=False, side1=9, side2=6)
        assert gt_lat.N == nx_lat.gr["G"].number_of_nodes()
        np.testing.assert_allclose(
            _spectrum(gt_lat.get_adjacency_matrix()),
            _spectrum(nx.to_numpy_array(nx_lat.gr["G"], weight=None)),
            atol=1e-6,
        )


class TestSideAndPbcResolution:
    """Side resolution and periodic default must match Lattice2DNX."""

    def test_default_periodic_matches_L2D_PBC(self):
        from lrgsglib.config.const import L2D_PBC

        assert Lattice2DGT(side1=10, geo="sqr").periodic == L2D_PBC

    def test_hex_side_rescale(self):
        """hex rescales side1 -> adjust_to_even(side1/sqrt(3)), side2=side1."""
        gt_lat = Lattice2DGT(side1=12, geo="hex")  # no explicit side2
        nx_lat = Lattice2DNX(side1=12, geo="hex")
        assert (gt_lat.side1, gt_lat.side2) == (nx_lat.side1, nx_lat.side2)
        assert gt_lat.N == nx_lat.gr["G"].number_of_nodes()

    def test_side_swap(self):
        """Larger side becomes side1, matching Lattice2DNX.__init_side__."""
        lat = Lattice2DGT(side1=5, side2=9, geo="sqr")
        assert (lat.side1, lat.side2) == (9, 5)

    def test_default_pbc_engine_agreement(self):
        """No periodic arg -> identical graph on both engines."""
        import networkx as nx

        for geo in ("sqr", "tri", "hex", "oct_sqr", "kgm"):
            gt_lat = Lattice2D(side1=12, geo=geo, engine="gt")
            nx_lat = Lattice2D(side1=12, geo=geo, engine="nx")
            assert gt_lat.N == nx_lat.N
            np.testing.assert_allclose(
                _spectrum(gt_lat.get_adjacency_matrix()),
                _spectrum(nx.to_numpy_array(nx_lat.gr["G"], weight=None)),
                atol=1e-6,
            )


class TestConstruction:
    """Basic sanity of the native GT generators."""

    @pytest.mark.parametrize("geo", ["sqr", "tri", "hex", "oct_sqr", "kgm"])
    def test_builds_nonempty(self, geo):
        lat = Lattice2DGT(side1=6, geo=geo)
        assert lat.N > 0
        assert lat.num_edges > 0

    def test_tri_hex_builds_open(self):
        lat = Lattice2DGT(side1=6, geo="tri_hex", periodic=False)
        assert lat.N > 0 and lat.num_edges > 0

    def test_oct_sqr_node_multiplier(self):
        """Rhomb-octagonal has 4 nodes per square-lattice site."""
        lat = Lattice2DGT(side1=10, side2=7, geo="oct_sqr")
        assert lat.N == 4 * 10 * 7

    @pytest.mark.parametrize(
        "geo,per",
        [
            ("sqr", True),
            ("tri", True),
            ("hex", True),
            ("oct_sqr", True),
            ("kgm", True),
            ("tri_hex", False),
        ],
    )
    def test_bulk_coordination(self, geo, per):
        lat = Lattice2DGT(side1=14, geo=geo, periodic=per)
        avg_deg = 2.0 * lat.num_edges / lat.N
        assert avg_deg <= Z_BY_GEO[geo] + 1e-9
        # periodic lattices reach (sqr/tri) or closely approach the bulk z
        assert avg_deg > Z_BY_GEO[geo] - 1.0

    @pytest.mark.parametrize(
        "geo,per",
        [
            ("sqr", True),
            ("tri", True),
            ("hex", True),
            ("oct_sqr", True),
            ("kgm", True),
            ("tri_hex", False),
        ],
    )
    def test_spectral_pipeline(self, geo, per):
        """Sign flipping + signed Laplacian works on every geometry."""
        lat = Lattice2DGT(side1=10, geo=geo, periodic=per, pflip=0.1, seed=42)
        lat.flip_random_fract_edges()
        assert lat.Ne_n > 0
        ev = np.linalg.eigvalsh(np.asarray(lat.get_signed_laplacian()))
        assert ev[0] > -1e-9  # signed Laplacian is PSD


class TestSmallWorld:
    """``prew`` rewiring (sqr_sw / tri_sw style) on the GT engine."""

    def test_edge_count_preserved(self):
        base = Lattice2DGT(side1=20, geo="sqr")
        sw = Lattice2DGT(side1=20, geo="sqr", prew=0.2, seed=1)
        assert sw.num_edges == base.num_edges

    def test_reproducible(self):
        a = Lattice2DGT(side1=16, geo="tri", prew=0.15, seed=7)
        b = Lattice2DGT(side1=16, geo="tri", prew=0.15, seed=7)
        assert a.num_edges == b.num_edges
        np.testing.assert_array_equal(
            _spectrum(a.get_adjacency_matrix()),
            _spectrum(b.get_adjacency_matrix()),
        )

    def test_rewiring_changes_topology(self):
        base = Lattice2DGT(side1=20, geo="sqr")
        sw = Lattice2DGT(side1=20, geo="sqr", prew=0.3, seed=3)
        assert not np.allclose(
            _spectrum(base.get_adjacency_matrix()),
            _spectrum(sw.get_adjacency_matrix()),
        )


class TestGeoAliasing:
    """Short and full geometry names are both accepted."""

    def test_geometries_set(self):
        for geo in ("sqr", "tri", "hex", "oct_sqr", "kgm", "tri_hex"):
            assert geo in Lattice2DGT.GEOMETRIES

    @pytest.mark.parametrize(
        "full,short",
        [
            ("octagonal_sqr", "oct_sqr"),
            ("kagome", "kgm"),
            ("squared", "sqr"),
            ("triangular", "tri"),
            ("hexagonal", "hex"),
        ],
    )
    def test_full_name_alias(self, full, short):
        lat = Lattice2DGT(side1=8, geo=full)
        assert lat.geo == short

    def test_tri_hexagonal_full_name(self):
        lat = Lattice2DGT(side1=8, geo="tri_hexagonal", periodic=False)
        assert lat.geo == "tri_hex"

    def test_invalid_geo_raises(self):
        with pytest.raises(ValueError):
            Lattice2DGT(side1=8, geo="not_a_lattice")

    def test_tri_hex_periodic_not_implemented(self):
        """GT refuses periodic tri-hexagonal, matching Lattice2DNX."""
        with pytest.raises(NotImplementedError):
            Lattice2DGT(side1=8, geo="tri_hex", periodic=True)


class TestFactoryInterchange:
    """Unified factory yields identical topology across engines (matched PBC)."""

    @pytest.mark.parametrize("geo,periodic", PARITY_CASES)
    def test_factory_parity(self, geo, periodic):
        import networkx as nx

        nx_lat = Lattice2D(
            side1=PARITY_SIDE, geo=geo, engine="nx", pbc=periodic
        )
        gt_lat = Lattice2D(
            side1=PARITY_SIDE, geo=geo, engine="gt", periodic=periodic
        )
        assert nx_lat.N == gt_lat.N
        np.testing.assert_allclose(
            _spectrum(gt_lat.get_adjacency_matrix()),
            _spectrum(nx.to_numpy_array(nx_lat.gr["G"], weight=None)),
            atol=1e-6,
        )

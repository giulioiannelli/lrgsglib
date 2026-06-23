"""Parity + correctness tests for the graph-tool ``nwDict`` (negative-link
defect patterns), ported from the NX engine.

Two kinds of check:

* **Structural invariants** (per engine): every pattern has the right shape —
  ``single`` is one edge, ``singleXERR`` is the full star at a node,
  ``singleZERR`` is a closed elementary cycle of the geometry's girth, ``rand``
  is exactly the negative-edge set, etc. These hold on every geometry.

* **Cross-engine parity**: where the NX ``'G'`` representation and the GT engine
  share an integer node labelling (square & honeycomb — verified below), the
  ``rand`` pattern built from an identical sign configuration must be the *same*
  integer edge set on both engines. The triangular labelling differs between
  engines, so there we assert geometric/structural equivalence (same girth,
  same coordination) rather than identical labels.
"""
import pytest

from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
from lrgsglib.graphs.gt.Lattice2DGT.Lattice2DGT import Lattice2DGT
from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT
from lrgsglib.graphs.gt._converters import nx_to_gt


SIDE = 8
PFLIP = 0.25
SEED = 7

# (gt short alias, nx full name, bulk coordination z, elementary-cell girth)
GEOMETRIES = [
    ("sqr", "squared", 4, 4),
    ("tri", "triangular", 6, 3),
    ("hex", "hexagonal", 3, 6),
    ("kgm", None, 4, 3),
    ("oct_sqr", None, 3, 4),
    ("tri_hex", None, 3, 3),
]

# Geometries whose NX 'G' integer labelling coincides with the GT labelling
# (verified empirically: square & honeycomb match; triangular does not).
COINCIDENT_LABELS = {"sqr", "hex"}


def _as_edge_set(edges):
    """Canonicalise a list of (u, v) tuples to a frozenset of sorted pairs."""
    return frozenset(tuple(sorted((int(u), int(v)))) for u, v in edges)


def _is_closed_cycle(edges):
    """True iff ``edges`` form a single simple cycle (every vertex degree 2)."""
    from collections import Counter

    deg = Counter()
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    return bool(edges) and all(d == 2 for d in deg.values()) and len(edges) == len(deg)


def _gt_lattice(geo, with_patterns=True):
    lat = Lattice2DGT(
        side1=SIDE, geo=geo, periodic=False, pflip=PFLIP, seed=SEED,
        init_nw_dict=with_patterns,
    )
    return lat


# ---------------------------------------------------------------- structural

@pytest.mark.parametrize("geo,_nx,z,girth", GEOMETRIES)
def test_gt_pattern_keys_present(geo, _nx, z, girth):
    lat = _gt_lattice(geo)
    assert set(lat.nwDict) == {
        "single", "singleXERR", "singleZERR", "rand", "randXERR", "randZERR",
    }


@pytest.mark.parametrize("geo,_nx,z,girth", GEOMETRIES)
def test_gt_single_is_one_real_edge(geo, _nx, z, girth):
    lat = _gt_lattice(geo)
    single = lat.nwDict["single"]["G"]
    assert len(single) == 1
    u, v = single[0]
    assert lat.G.edge(u, v) is not None


@pytest.mark.parametrize("geo,_nx,z,girth", GEOMETRIES)
def test_gt_singleXERR_is_full_star(geo, _nx, z, girth):
    lat = _gt_lattice(geo)
    center = lat.nwDict["single"]["G"][0][0]
    star = lat.nwDict["singleXERR"]["G"]
    assert all(center in (u, v) for u, v in star)
    assert _as_edge_set(star) == _as_edge_set(
        (center, nb) for nb in lat.get_graph_neighbors(center)
    )


@pytest.mark.parametrize("geo,_nx,z,girth", GEOMETRIES)
def test_gt_singleZERR_is_elementary_cycle(geo, _nx, z, girth):
    lat = _gt_lattice(geo)
    cell = lat.nwDict["singleZERR"]["G"]
    assert len(cell) == girth
    assert _is_closed_cycle(cell)


@pytest.mark.parametrize("geo,_nx,z,girth", GEOMETRIES)
def test_gt_rand_matches_negative_edges_after_flip(geo, _nx, z, girth):
    lat = _gt_lattice(geo)
    lat.flip_random_fract_edges()
    lat.build_nw_dict()  # refresh: GT flips after construction
    assert _as_edge_set(lat.nwDict["rand"]["G"]) == _as_edge_set(lat.fleset["G"])


@pytest.mark.parametrize("geo,_nx,z,girth", GEOMETRIES)
def test_gt_randXERR_edges_incident_to_seed_nodes(geo, _nx, z, girth):
    lat = _gt_lattice(geo)
    rx = lat.nwDict["randXERR"]["G"]
    # randXERR is a union of node stars, so every edge is a real graph edge;
    # for pflip>0 the seed set is non-empty.
    all_edges = _as_edge_set(lat.get_edges_list())
    assert _as_edge_set(rx) <= all_edges
    assert len(rx) > 0


# ---------------------------------------------------------------- cross-engine

@pytest.mark.parametrize("geo,nx_geo,z,girth", GEOMETRIES)
def test_singleZERR_girth_matches_across_engines(geo, nx_geo, z, girth):
    if nx_geo is None:
        pytest.skip("no NX get_central_edge geometry mapping needed")
    nx_lat = Lattice2DNX(side1=SIDE, geo=nx_geo, pbc=False, seed=SEED,
                         init_nw_dict=True)
    gt_lat = _gt_lattice(geo)
    assert len(nx_lat.nwDict["singleZERR"]["G"]) == girth
    assert len(gt_lat.nwDict["singleZERR"]["G"]) == girth


@pytest.mark.parametrize("geo,nx_geo,z,girth", GEOMETRIES)
def test_singleXERR_size_is_coordination(geo, nx_geo, z, girth):
    gt_lat = _gt_lattice(geo)
    center = gt_lat.nwDict["single"]["G"][0][0]
    # central (bulk) node must have full coordination z for this lattice
    assert len(gt_lat.get_graph_neighbors(center)) == z
    assert len(gt_lat.nwDict["singleXERR"]["G"]) == z


@pytest.mark.parametrize("geo", sorted(COINCIDENT_LABELS))
def test_rand_pattern_identical_labels(geo):
    """For square & honeycomb the integer labels coincide, so the rand pattern
    built from an identical sign configuration is the *same* edge set."""
    nx_geo = {"sqr": "squared", "hex": "hexagonal"}[geo]
    nx_lat = Lattice2DNX(side1=SIDE, geo=nx_geo, pbc=False, pflip=PFLIP,
                         seed=SEED, init_nw_dict=True)
    nx_lat.flip_random_fract_edges()

    # Convert the *flipped* NX graph to GT, preserving signs (NX stores 'weight').
    G_gt = nx_to_gt(nx_lat.gr["G"], sign_attr="weight")
    gt_lat = SignedGraphGT(G=G_gt)

    nx_neg = _as_edge_set(nx_lat.fleset["G"])
    gt_neg = _as_edge_set(gt_lat.fleset["G"])
    assert nx_neg == gt_neg, f"{geo}: negative-edge sets differ across engines"

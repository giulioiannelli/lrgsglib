"""Canonical / oriented ZERR-cell convention (PLAN-Disorder-Creator §5a, §6.1-6.3).

The ZERR ("elementary cell") support flips the smallest cycle through each seed
node. Two guarantees this pins:

1. **Engine-independence (canonical layer 1).** ``cell_edges`` is a deterministic
   function of the node labels, so the NetworkX and graph-tool engines pick the
   *same* cell for every node — on every geometry whose two engines share a node
   labelling (square, hexagonal). The previous implementation tie-broke by
   neighbour-iteration order and so diverged across engines (and even aliased a
   square-lattice node onto a neighbour's face).

2. **Node-distinctness (oriented layer 2).** On the square lattice each node owns
   its ``+x/+y`` north-east plaquette, a bijection in the bulk, so distinct
   ``randZERR`` seeds flip distinct cells instead of silently colliding.

Plus the §5b / §6.3 regression: structured supports must build through the shared
``NwContainer`` on the families whose bespoke containers used to choke — the
coordinate-less ``oct_sqr`` / ``kgm`` geometries (a ``KeyError``), the 3D lattice
(``KeyError`` on the loop ``*ZERR`` patterns it never built), and the DGM graph (a
``NameError`` from an unimported ``random``).

Note: the ``triangular`` generators label nodes differently on the two engines
(same graph up to isomorphism, different integer labels), so NX==GT cannot hold
there for *any* label-based quantity — that is a generator concern, out of scope
here, and deliberately not asserted.
"""
import warnings
from collections import Counter

import pytest

from lrgsglib.graphs import Lattice2D, Lattice3D, DGMgraph, Disorder
from lrgsglib.graphs._shared._nw_geometry import elementary_cell_edges

pytest.importorskip("graph_tool.all")

ENGINES = ["nx", "gt"]
CELL_SIZE = {"sqr": 4, "tri": 3, "hex": 6}


def _norm(edges):
    """Cell as an order-independent frozenset of sorted ``(u, v)`` tuples."""
    return frozenset(tuple(sorted(e)) for e in edges)


def _is_simple_cycle(edges):
    deg = Counter()
    for u, v in edges:
        deg[u] += 1
        deg[v] += 1
    return bool(edges) and all(d == 2 for d in deg.values()) and len(edges) == len(deg)


def _lattice(geo, engine, side=8, pbc=True):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Lattice2D(side1=side, geo=geo, engine=engine, seed=1, pbc=pbc)


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("geo,size", CELL_SIZE.items(), ids=list(CELL_SIZE))
def test_cell_is_simple_cycle_of_expected_size(engine, geo, size):
    """Every interior cell is a single simple cycle of the geometry's face size."""
    lat = _lattice(geo, engine)
    seen = False
    for v in range(lat.N):
        cell = [tuple(sorted(e)) for e in lat.cell_edges(v, "G")]
        if not cell:
            continue  # boundary node on no face (non-pbc) — allowed
        seen = True
        assert _is_simple_cycle(cell), f"{geo}/{engine} node {v}: not a simple cycle"
        assert len(cell) == size, f"{geo}/{engine} node {v}: size {len(cell)} != {size}"
    assert seen, f"{geo}/{engine}: no cells found"


@pytest.mark.parametrize("geo", ["sqr", "hex"])
@pytest.mark.parametrize("pbc", [True, False])
def test_cell_cross_engine_identical(geo, pbc):
    """sqr/hex share a labelling across engines → identical cell for every node."""
    nx_l, gt_l = _lattice(geo, "nx", pbc=pbc), _lattice(geo, "gt", pbc=pbc)
    assert nx_l.N == gt_l.N
    for v in range(nx_l.N):
        assert _norm(nx_l.cell_edges(v, "G")) == _norm(gt_l.cell_edges(v, "G")), (
            f"{geo} pbc={pbc} node {v}: NX and GT picked different cells"
        )


@pytest.mark.parametrize("engine", ENGINES)
def test_square_torus_cells_node_distinct(engine):
    """On the square torus the oriented picker is a bijection: N distinct cells."""
    lat = _lattice("sqr", engine, side=8, pbc=True)
    cells = {_norm(lat.cell_edges(v, "G")) for v in range(lat.N)}
    assert len(cells) == lat.N, (
        f"sqr/{engine}: only {len(cells)}/{lat.N} distinct cells (defects collide)"
    )


@pytest.mark.parametrize("engine", ENGINES)
def test_elementary_cell_edges_deterministic(engine):
    """The geometry-free baseline returns the same cell on repeated calls."""
    lat = _lattice("sqr", engine)
    for v in range(0, lat.N, 5):
        a = _norm(elementary_cell_edges(lat, v, "G"))
        b = _norm(elementary_cell_edges(lat, v, "G"))
        assert a == b


# §5b: structured supports on the coordinate-less geometries must not KeyError.
OCT_KGM_SUPPORTS = ["randXERR", "randZERR", "singleXERR", "singleZERR", "single"]


@pytest.mark.parametrize("geo", ["oct_sqr", "kgm"])
@pytest.mark.parametrize("support", OCT_KGM_SUPPORTS)
def test_structured_support_oct_kgm_no_keyerror(geo, support):
    """oct_sqr/kgm structured disorder builds + applies (was a KeyError)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = Lattice2D(
            side1=6, geo=geo, engine="nx", seed=7,
            disorder=Disorder(support, pflip=0.1),
        )
    assert support in g.nwDict, f"{geo}: {support!r} not built"
    assert g.Ne_n > 0, f"{geo}: {support!r} flipped no edge"


# §6.3 (Shape-A): the Lattice3DNX / DGMgraphNX bespoke containers were retired onto
# the shared NwContainer, so those families gain the full structured-support set — and
# the cascade reaches Lattice3DGT / DGMgraphGT, which subclass the NX impls. Before:
# the 3D lattice raised KeyError on the loop ``*ZERR`` patterns it never built; DGM
# raised NameError on an unimported ``random``. Both engines must now build every one.
STRUCTURED_SUPPORTS = ["single", "singleXERR", "randXERR", "singleZERR", "randZERR"]


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("support", STRUCTURED_SUPPORTS)
def test_structured_support_lattice3d(engine, support):
    """3D-lattice structured disorder builds + flips (``*ZERR`` was a KeyError)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = Lattice3D(
            dim=4, engine=engine, seed=7,
            disorder=Disorder(support, pflip=0.15),
        )
    assert support in g.nwDict, f"L3D/{engine}: {support!r} not built"
    assert g.Ne_n > 0, f"L3D/{engine}: {support!r} flipped no edge"


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("support", STRUCTURED_SUPPORTS)
def test_structured_support_dgm(engine, support):
    """DGM structured disorder builds + flips (every support was a NameError)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = DGMgraph(
            n=4, engine=engine, seed=7,
            disorder=Disorder(support, pflip=0.15),
        )
    assert support in g.nwDict, f"DGM/{engine}: {support!r} not built"
    assert g.Ne_n > 0, f"DGM/{engine}: {support!r} flipped no edge"

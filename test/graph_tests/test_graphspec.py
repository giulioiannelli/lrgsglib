"""Tests for the engine-free graph spec layer (``graphs/_shared/_spec.py``).

Validates the single-contract infrastructure: a generator emits a ``GraphSpec``,
which materializes into either engine and converts between them, all producing
topologically identical graphs.
"""
import numpy as np
import networkx as nx
import pytest

from lrgsglib.graphs._shared._spec import (
    GraphSpec, materialize_nx, materialize_gt, spec_from_nx, spec_from_gt,
    spec_from_graph, convert, build, detect_engine,
)
from lrgsglib.graphs._engine import GraphEngine

gt = pytest.importorskip("graph_tool.all")
gts = pytest.importorskip("graph_tool.spectral")

SIDE = 8


# --- fixtures: a small reusable spec ---------------------------------------- #
def _hex_spec() -> GraphSpec:
    edges = np.array([[0, 1], [1, 2], [0, 3], [1, 4], [2, 5], [3, 4], [4, 5]])
    pos = np.array([[0, 0], [1, 0], [2, 0], [0, 1], [1, 1], [2, 1]], float)
    return GraphSpec(n=6, edges=edges, pos=pos)


def _edgeset(spec: GraphSpec):
    return {tuple(sorted(map(int, e))) for e in spec.edges}


def _spectrum_binarized(adj: np.ndarray) -> np.ndarray:
    return np.sort(np.linalg.eigvalsh((np.abs(adj) > 0).astype(float)))


# --- validation ------------------------------------------------------------- #
def test_spec_validation_rejects_out_of_range_edge():
    with pytest.raises(ValueError):
        GraphSpec(n=3, edges=np.array([[0, 5]]))


def test_spec_validation_rejects_pos_length_mismatch():
    with pytest.raises(ValueError):
        GraphSpec(n=3, edges=np.array([[0, 1]]), pos=np.zeros((2, 2)))


def test_from_genresult_bridge():
    spec = GraphSpec.from_genresult((4, [(0, 1), (1, 2), (2, 3)], None))
    assert spec.n == 4 and spec.num_edges == 3 and spec.pos is None


# --- round-trip both engines ------------------------------------------------ #
def test_roundtrip_nx():
    spec = _hex_spec()
    G = materialize_nx(spec)
    assert G.number_of_nodes() == spec.n
    back = spec_from_nx(G)
    assert _edgeset(back) == _edgeset(spec)
    assert np.allclose(back.pos, spec.pos)


def test_roundtrip_gt():
    spec = _hex_spec()
    G = materialize_gt(spec)
    assert int(G.num_vertices()) == spec.n
    back = spec_from_gt(G)
    assert _edgeset(back) == _edgeset(spec)
    assert np.allclose(back.pos, spec.pos)


def test_detect_engine():
    spec = _hex_spec()
    assert detect_engine(materialize_nx(spec)) is GraphEngine.NETWORKX
    assert detect_engine(materialize_gt(spec)) is GraphEngine.GRAPHTOOL


# --- conversion ------------------------------------------------------------- #
def test_convert_both_directions():
    spec = _hex_spec()
    Gnx, Ggt = materialize_nx(spec), materialize_gt(spec)
    assert _edgeset(spec_from_graph(convert(Gnx, "gt"))) == _edgeset(spec)
    assert _edgeset(spec_from_graph(convert(Ggt, "nx"))) == _edgeset(spec)


# --- the core promise: one generator, both engines, == factory -------------- #
@pytest.mark.parametrize("geo", ["sqr", "tri", "kgm"])
def test_native_generator_through_spec_matches_factory(geo):
    from lrgsglib.graphs.gt.Lattice2DGT import _generators as g2d
    from lrgsglib.graphs import Lattice2D

    fn = {"sqr": g2d.squared, "tri": g2d.triangular, "kgm": g2d.kagome}[geo]
    spec = GraphSpec.from_genresult(fn(SIDE, SIDE, False))

    Anx = nx.to_numpy_array(materialize_nx(spec), nodelist=range(spec.n))
    Agt = gts.adjacency(materialize_gt(spec)).toarray()
    # both engines agree on topology
    assert np.allclose(_spectrum_binarized(Anx), _spectrum_binarized(Agt))

    # and it matches what the factory builds today (topology, sign-agnostic)
    fac = Lattice2D(SIDE, geo=geo, engine="gt", periodic=False)
    Afac = gts.adjacency(fac.gr["G"]).toarray()
    assert np.allclose(_spectrum_binarized(Agt), _spectrum_binarized(Afac))


# --- optional engine-native fast path --------------------------------------- #
def test_build_default_and_fast_path():
    spec = _hex_spec()

    def gen():
        return spec
    G = build(gen, "gt")
    assert int(G.num_vertices()) == spec.n

    sentinel = {}

    def gen_fast():
        return spec
    def native():
        sentinel["used"] = True
        return materialize_gt(spec)
    gen_fast.gt_native = native

    build(gen_fast, "gt")
    assert sentinel.get("used") is True          # fast path taken for gt
    sentinel.clear()
    build(gen_fast, "nx")
    assert "used" not in sentinel                # fast path NOT taken for nx

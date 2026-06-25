"""Universal disorder tools across the whole graph family, both engines.

Proves the two guarantees:
  1. Every graph gets the geometry-free nwDict patterns (rand, randXERR) via
     the engine-neutral default ``NwContainer`` — not just lattices/ER.
  2. ``prew`` (rewiring) and ``pdil`` (dilution) are base-level tools available
     to every graph, on both engines, reachable through the factories.
"""
import warnings

import pytest

from lrgsglib.graphs import (
    BarabasiAlbert, ErdosRenyi, HolmeKim, WattsStrogatz, Lattice2D,
)

pytest.importorskip("graph_tool.all")

ENGINES = ["nx", "gt"]

# (factory, kwargs) for a representative spread: scale-free, random, clustered,
# small-world, and a geometric lattice (which also has its own richer container).
FAMILY = [
    ("BarabasiAlbert", BarabasiAlbert, dict(n=400, m=3)),
    ("ErdosRenyi", ErdosRenyi, dict(n=400, p=0.02)),
    ("HolmeKim", HolmeKim, dict(n=400, m=3, p=0.3)),
    ("WattsStrogatz", WattsStrogatz, dict(n=400, k=4, p=0.1)),
    ("Lattice2D", Lattice2D, dict(side1=16, geo="sqr")),
]

SEED = 1
PFLIP = 0.2


def _num_edges(g):
    G = g.gr["G"]
    return G.number_of_edges() if hasattr(G, "number_of_edges") else int(G.num_edges())


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("name,factory,kw", FAMILY, ids=[f[0] for f in FAMILY])
def test_nwdict_universal(engine, name, factory, kw):
    """init_nw_dict=True yields at least rand + randXERR for every graph."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = factory(**kw, pflip=PFLIP, seed=SEED, engine=engine,
                    init_nw_dict=True)
    keys = set(g.nwDict.keys())
    assert {"rand", "randXERR"} <= keys, f"{name}/{engine} missing geom-free keys"


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("name,factory,kw", FAMILY, ids=[f[0] for f in FAMILY])
def test_prew_preserves_edge_count(engine, name, factory, kw):
    """Rewiring keeps the edge count (it moves endpoints, doesn't add/remove)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = factory(**kw, seed=SEED, engine=engine)
        rewired = factory(**kw, seed=SEED, engine=engine, prew=0.15)
    assert _num_edges(rewired) == _num_edges(base), f"{name}/{engine} edge count moved"


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("name,factory,kw", FAMILY, ids=[f[0] for f in FAMILY])
def test_pdil_reduces_edges(engine, name, factory, kw):
    """Dilution removes edges (and may shrink N via giant-component extraction)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        base = factory(**kw, seed=SEED, engine=engine)
        diluted = factory(**kw, seed=SEED, engine=engine, pdil=0.2)
    assert _num_edges(diluted) < _num_edges(base), f"{name}/{engine} not diluted"

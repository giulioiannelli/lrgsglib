"""Tests for the unified disorder model (PLAN-Disorder-Model).

A disorder realization = support x coupling-law, expressed via the engine-neutral
:class:`Disorder` spec and applied at construction on BOTH engines. These tests
pin the agreed contract:

- default ('rand' + 'flip') auto-flips a pflip fraction at construction;
- pflip=0 leaves a positive graph (backward compatible);
- disorder=None defers (old manual-flip workflow);
- redundant flip_random_fract_edges() is idempotent (no double-flip);
- 'all' support flips every edge; structured supports route through nwDict;
- distributional laws assign continuous couplings (NX; GT is Stage 3);
- a graph handed in pre-signed is not re-flipped.
"""
import numpy as np
import pytest

from lrgsglib.graphs import (
    BarabasiAlbert,
    ErdosRenyi,
    HolmeKim,
    Lattice2D,
    WattsStrogatz,
    Disorder,
)

ENGINES = ["nx", "gt"]

# (factory, kwargs) covering scale-free, random, small-world, geometric.
FAMILIES = [
    ("BarabasiAlbert", BarabasiAlbert, dict(n=400, m=3)),
    ("ErdosRenyi", ErdosRenyi, dict(n=400, p=0.02)),
    ("HolmeKim", HolmeKim, dict(n=400, m=3, p=0.3)),
    ("WattsStrogatz", WattsStrogatz, dict(n=400, k=4, p=0.1)),
    ("Lattice2D", Lattice2D, dict(side1=16, geo="sqr")),
]


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("name,factory,kwargs", FAMILIES, ids=[f[0] for f in FAMILIES])
def test_default_autoflip(name, factory, kwargs, engine):
    """Default disorder ('rand'+'flip') realizes ~pflip*Ne negatives at ctor."""
    pflip = 0.2
    g = factory(**kwargs, pflip=pflip, engine=engine, seed=7)
    assert g.Ne_n > 0, f"{name}/{engine}: pflip>0 must auto-flip at construction"
    # within a loose tolerance of the requested fraction
    assert abs(g.Ne_n / g.Ne - pflip) < 0.05


@pytest.mark.parametrize("engine", ENGINES)
def test_pflip_zero_is_positive(engine):
    g = Lattice2D(side1=16, geo="sqr", pflip=0.0, engine=engine, seed=7)
    assert g.Ne_n == 0


@pytest.mark.parametrize("engine", ENGINES)
def test_disorder_none_defers(engine):
    g = BarabasiAlbert(n=400, m=3, pflip=0.2, engine=engine, seed=7, disorder=None)
    assert g.Ne_n == 0, "disorder=None must defer sign realization"
    g.flip_random_fract_edges()
    assert g.Ne_n > 0, "manual flip must realize the deferred selection"


@pytest.mark.parametrize("engine", ENGINES)
def test_redundant_flip_is_idempotent(engine):
    """After ctor auto-flip, extra no-arg flips must not change the sign set."""
    g = Lattice2D(side1=16, geo="sqr", pflip=0.2, engine=engine, seed=7)
    n0 = g.Ne_n
    g.flip_random_fract_edges()
    g.flip_random_fract_edges()
    assert g.Ne_n == n0


@pytest.mark.parametrize("engine", ENGINES)
def test_support_all_flips_every_edge(engine):
    g = Lattice2D(side1=12, geo="sqr", engine=engine, seed=7, disorder=Disorder("all"))
    assert g.Ne_n == g.Ne


@pytest.mark.parametrize("engine", ENGINES)
def test_structured_support_routes_through_nwdict(engine):
    """A structured support builds nwDict at ctor and applies the pattern."""
    g = BarabasiAlbert(
        n=400, m=3, engine=engine, seed=7, disorder=Disorder("randXERR", pflip=0.05)
    )
    assert "randXERR" in g.nwDict
    assert g.Ne_n > 0


# Central/cell patterns: previously lattice-only (needed a geometric central
# edge + the build_single/build_zerr container flags). PLAN-Disorder-Model B
# makes them universal via a base hub ``get_central_edge`` and on-demand build
# flags, so they must now build + apply on EVERY family.
CENTRAL_SUPPORTS = ["single", "singleXERR", "singleZERR", "randZERR"]


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("support", CENTRAL_SUPPORTS)
@pytest.mark.parametrize("name,factory,kwargs", FAMILIES, ids=[f[0] for f in FAMILIES])
def test_structured_supports_universal(name, factory, kwargs, support, engine):
    """Central/cell nwDict supports build + flip on every family, not just lattices."""
    g = factory(
        **kwargs, engine=engine, seed=7, disorder=Disorder(support, pflip=0.05)
    )
    assert support in g.nwDict, f"{name}/{engine}: {support!r} must be built on demand"
    # All sampled families have cycles through their hub, so even the cell
    # (ZERR) patterns flip at least one edge.
    assert g.Ne_n > 0, f"{name}/{engine}: {support!r} must flip at least one edge"


@pytest.mark.parametrize("engine", ENGINES)
def test_hub_central_edge_is_highest_degree(engine):
    """The geometry-free base get_central_edge returns a hub-incident edge."""
    g = BarabasiAlbert(n=300, m=3, engine=engine, seed=7)
    u, v = g.get_central_edge()
    degs = {n: len(g.get_graph_neighbors(n)) for n in g.get_nodes_list()}
    assert degs[u] == max(degs.values()), "central node must be a highest-degree hub"
    assert v in g.get_graph_neighbors(u), "central edge must be incident to the hub"


def test_string_shorthand_equals_default():
    g1 = Lattice2D(side1=12, geo="sqr", pflip=0.2, engine="nx", seed=7)
    g2 = Lattice2D(side1=12, geo="sqr", pflip=0.2, engine="nx", seed=7, disorder="rand")
    assert g1.Ne_n == g2.Ne_n


def test_bad_coupling_law_rejected():
    with pytest.raises(ValueError):
        Disorder("rand", law="does-not-exist")


def _offdiag_couplings(g, engine):
    """Off-diagonal signed couplings of a graph, engine-agnostic."""
    A = g.get_signed_adjacency() if engine == "gt" else g.adj.toarray()
    iu = np.triu_indices_from(A, k=1)
    vals = A[iu]
    return vals[vals != 0]


@pytest.mark.parametrize("engine", ENGINES)
def test_distributional_gaussian(engine):
    """Distributional couplings: continuous edge weights on BOTH engines.

    NX uses its float ``weight`` attribute directly; GT writes a ``weight``
    (double) property that the signed adjacency reads (Stage 3).
    """
    sigma = 2.0
    g = Lattice2D(
        side1=20,
        geo="sqr",
        engine=engine,
        seed=7,
        disorder=Disorder("all", law="gaussian", params={"mu": 0.0, "sigma": sigma}),
    )
    w = _offdiag_couplings(g, engine)
    assert w.shape[0] == g.Ne
    assert abs(w.mean()) < 0.3 * sigma
    assert abs(w.std() - sigma) < 0.3 * sigma
    assert (w < 0).any() and (w > 0).any()
    # continuous (not just ±1)
    assert len(np.unique(np.round(w, 6))) > 10


@pytest.mark.parametrize("engine", ENGINES)
def test_distributional_signed_laplacian_psd(engine):
    """The signed Laplacian stays symmetric and PSD with continuous couplings."""
    g = Lattice2D(
        side1=12,
        geo="sqr",
        engine=engine,
        seed=7,
        disorder=Disorder("all", law="gaussian", params={"mu": 0.0, "sigma": 1.0}),
    )
    L = g.get_signed_laplacian()
    L = np.asarray(L.todense()) if hasattr(L, "todense") else np.asarray(L)
    assert np.allclose(L, L.T, atol=1e-9)
    ev = np.linalg.eigvalsh(L)
    assert ev.min() > -1e-8


def test_presigned_graph_not_reflipped():
    """Converting a signed NX graph to GT must not layer extra random flips."""
    from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
    from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT
    from lrgsglib.graphs.gt._converters import nx_to_gt

    nx_sg = Lattice2DNX(side1=8, geo="sqr", pflip=0.2, seed=11)
    n_neg = nx_sg.Ne_n
    G_gt = nx_to_gt(nx_sg.gr["G"], sign_attr="weight")
    gt_sg = SignedGraphGT(G=G_gt, pflip=0.2, seed=11)
    assert gt_sg.Ne_n == n_neg, "pre-signed graph must be respected, not re-flipped"

"""Registered supports (§6.4) + disorder composition (§6.8).

Two additive extensions of the ``Disorder`` model:

1. **``register_support``** — a third support class that resolves straight to
   edges and *bypasses* ``nwDict`` (registered names are not in
   ``STRUCTURED_SUPPORTS``), mirroring ``register_coupling``. It composes with any
   coupling law and needs no ``NwContainer`` machinery.

2. **``CompositeDisorder``** — combine several ``Disorder`` specs under one
   ``mode``: ``overlay`` (apply each in turn), ``mixture`` (partition the seed
   budget by ``weights`` — the "half-X/half-Z" case), and the edge-set algebra
   ``union`` / ``intersection`` / ``difference`` / ``symmetric_difference``
   (also reachable via the ``|`` / ``&`` / ``-`` / ``^`` / ``+`` operators).

Parity note: the ``mixture`` planner draws seeds from an explicit seeded
``Generator``, so it is NX==GT identical. The set-algebra modes reuse the
``nwDict`` random-seed patterns, which already differ across engines (a
pre-existing RNG asymmetry — NX seeds global ``random``, GT does not), so only
the *engine-local* set-algebra identities are asserted there, not cross-engine
equality.
"""
import warnings

import pytest

from lrgsglib.graphs import (
    Lattice2D,
    ErdosRenyi,
    Disorder,
    CompositeDisorder,
    register_support,
    registered_supports,
)

pytest.importorskip("graph_tool.all")

ENGINES = ["nx", "gt"]


def _lat(engine, disorder, side=24, pflip=0.0, seed=5):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return Lattice2D(
            side1=side, geo="sqr", engine=engine, seed=seed,
            pflip=pflip, disorder=disorder,
        )


# --------------------------------------------------------------------------
# Coercion + spec shape (pure; no graph)
# --------------------------------------------------------------------------
def test_dict_shorthand_is_mixture():
    from lrgsglib.graphs._shared._disorder import as_disorder

    c = as_disorder({"randXERR": 0.5, "randZERR": 0.5}, pflip=0.1)
    assert isinstance(c, CompositeDisorder)
    assert c.mode == "mixture"
    assert c.weights == [0.5, 0.5]
    assert [x.support for x in c.components] == ["randXERR", "randZERR"]
    assert c.pflip == pytest.approx(0.1)  # budget from the top-level pflip


def test_list_and_weighted_list():
    from lrgsglib.graphs._shared._disorder import as_disorder

    assert as_disorder(["randXERR", "randZERR"]).mode == "overlay"
    w = as_disorder([("randXERR", 0.3), ("randZERR", 0.7)], pflip=0.1)
    assert w.mode == "mixture" and w.weights == pytest.approx([0.3, 0.7])


def test_operator_sugar_and_flattening():
    A, B = Disorder("randXERR", pflip=0.05), Disorder("randZERR", pflip=0.05)
    assert (A | B).mode == "union"
    assert (A & B).mode == "intersection"
    assert (A - B).mode == "difference"
    assert (A ^ B).mode == "symmetric_difference"
    assert (A + B).mode == "overlay"
    assert len((A | B | B).components) == 3          # flat, not nested
    assert (A | B).is_structured and (A | B).is_flip
    assert (A | B).build_flags() == (False, True)     # randZERR needs build_zerr


def test_invalid_composites_rejected():
    A, B = Disorder("randXERR", pflip=0.05), Disorder("randZERR", pflip=0.05)
    with pytest.raises(ValueError):
        CompositeDisorder([A], mode="bogus")
    with pytest.raises(ValueError):
        CompositeDisorder([A, B, A], mode="difference")   # binary only
    with pytest.raises(ValueError):
        CompositeDisorder([], mode="overlay")             # empty


# --------------------------------------------------------------------------
# register_support (§6.4)
# --------------------------------------------------------------------------
def test_hubxerr_registered():
    assert "hubXERR" in registered_supports()
    assert Disorder("hubXERR").is_registered_support
    assert not Disorder("hubXERR").is_structured   # bypasses nwDict


@pytest.mark.parametrize("engine", ENGINES)
def test_registered_support_builds_and_flips(engine):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        g = ErdosRenyi(n=200, p=0.03, engine=engine, seed=7,
                       disorder=Disorder("hubXERR"))
        g3 = ErdosRenyi(n=200, p=0.03, engine=engine, seed=7,
                        disorder=Disorder("hubXERR", support_params={"k": 3}))
    assert g.Ne_n > 0, f"{engine}: hubXERR flipped nothing"
    assert g3.Ne_n >= g.Ne_n, f"{engine}: k=3 should touch >= edges than k=1"
    assert not hasattr(g, "nwDict") or "hubXERR" not in getattr(g, "nwDict", {})


@pytest.mark.parametrize("engine", ENGINES)
def test_custom_registered_support(engine):
    @register_support("_test_node0_star")
    def _builder(sg, pflip, rng, on_g, **params):
        from lrgsglib.graphs._shared._nw_geometry import star_edges

        n0 = sorted(sg.get_nodes_list(on_g), key=str)[0]
        return [(int(u), int(v)) for (u, v) in star_edges(sg, n0, on_g)]

    g = _lat(engine, Disorder("_test_node0_star"))
    assert g.Ne_n > 0


# --------------------------------------------------------------------------
# Mixture (§6.8) — the "half-X / half-Z" headline
# --------------------------------------------------------------------------
@pytest.mark.parametrize("engine", ENGINES)
def test_mixture_builds(engine):
    g = _lat(engine, {"randXERR": 0.5, "randZERR": 0.5}, pflip=0.10)
    assert g.Ne_n > 0
    assert "randZERR" in g.nwDict  # build_zerr forced on for the ZERR arm


def test_mixture_is_cross_engine_identical():
    """The seeded mixture planner is deterministic AND NX==GT."""
    kw = dict(disorder={"randXERR": 0.5, "randZERR": 0.5}, pflip=0.10)
    assert _lat("nx", **kw).Ne_n == _lat("gt", **kw).Ne_n


@pytest.mark.parametrize("engine", ENGINES)
def test_mixture_deterministic(engine):
    kw = dict(disorder={"randXERR": 0.4, "randZERR": 0.6}, pflip=0.12)
    assert _lat(engine, **kw).Ne_n == _lat(engine, **kw).Ne_n


@pytest.mark.parametrize("engine", ENGINES)
def test_mixture_weight_extremes(engine):
    all_x = _lat(engine, {"randXERR": 1.0, "randZERR": 0.0}, pflip=0.10).Ne_n
    all_z = _lat(engine, {"randXERR": 0.0, "randZERR": 1.0}, pflip=0.10).Ne_n
    assert all_x > 0 and all_z > 0


def test_mixture_rejects_non_randfamily():
    with pytest.raises(ValueError):
        _lat("nx", {"single": 0.5, "randZERR": 0.5}, pflip=0.1)


@pytest.mark.parametrize("engine", ENGINES)
def test_mixture_discrete_plus_continuous(engine):
    """A flip arm + a gaussian arm: the flips must survive (GT force_weight)."""
    g = _lat(
        engine,
        [(Disorder("randXERR", law="flip"), 0.5),
         (Disorder("randZERR", law="gaussian", params={"sigma": 1.5}), 0.5)],
        pflip=0.10,
    )
    assert g.Ne_n > 0


# --------------------------------------------------------------------------
# Set algebra (§6.8) — verified per-build (self-consistent; GT nwDict random
# patterns are not reproducible across separate builds, so we never compare
# Ne_n across builds — only against the SAME graph's realized patterns).
# --------------------------------------------------------------------------
def _keyset(pattern):
    tup = pattern["G"] if isinstance(pattern, dict) else pattern
    return {tuple(sorted((int(u), int(v)))) for (u, v) in tup}


@pytest.mark.parametrize("engine", ENGINES)
@pytest.mark.parametrize("op", ["union", "intersection", "difference", "symmetric_difference"])
def test_set_algebra_realizes_operation(engine, op):
    A = Disorder("randXERR", pflip=0.06)
    B = Disorder("randZERR", pflip=0.06)
    spec = {"union": A | B, "intersection": A & B,
            "difference": A - B, "symmetric_difference": A ^ B}[op]
    g = _lat(engine, spec)
    X, Z = _keyset(g.nwDict["randXERR"]), _keyset(g.nwDict["randZERR"])
    expected = {"union": X | Z, "intersection": X & Z,
                "difference": X - Z, "symmetric_difference": X ^ Z}[op]
    assert g.Ne_n == len(expected), f"{engine}/{op}: Ne_n != |{op}(X, Z)|"


@pytest.mark.parametrize("engine", ENGINES)
def test_overlay_realizes_union(engine):
    """Overlay of two flip supports realizes exactly the union of their edges."""
    A = Disorder("randXERR", pflip=0.06)
    B = Disorder("randZERR", pflip=0.06)
    g = _lat(engine, [A, B])
    X, Z = _keyset(g.nwDict["randXERR"]), _keyset(g.nwDict["randZERR"])
    assert g.Ne_n == len(X | Z)


# --------------------------------------------------------------------------
# Backward compatibility
# --------------------------------------------------------------------------
@pytest.mark.parametrize("engine", ENGINES)
def test_single_disorder_unchanged(engine):
    """A plain single-support Disorder still realizes as before."""
    plain = _lat(engine, Disorder("randXERR", pflip=0.06)).Ne_n
    assert plain > 0
    # pflip=0 default => unsigned
    assert _lat(engine, None).Ne_n == 0

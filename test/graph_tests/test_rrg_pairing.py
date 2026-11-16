"""
Exact pairing sampler of k-regular graphs (graphs/_shared/_regular_pairing.py).

- Uniformity on 6 nodes: the sampler hits all 70 labelled cubic graphs, at
  frequencies a chi-square cannot tell from equal. At DRAWS_PER_GRAPH this
  rules out only gross errors (the native NX generator passes it too).
  Exactness rests on the (k!)^n counting argument in _regular_pairing.py;
  the short-cycle test is the one that tells NX from uniform.
- Short cycles: at N=1000 the mean counts of 3- and 4-cycles match the
  uniform-ensemble limits (k-1)^j / (2j). These limits are a literature
  result (Bollobas 1980; Wormald 1999, "Models of random regular graphs"),
  not in .agents/ref; this test checks them numerically. The native NX
  generator was measured to miss the triangle limit (see the
  RRG_SAMPLER_NATIVE comment in config/progargs/defs/RandomRegular.py);
  that measurement is not a test here.
- Invariants and seeding: the graph classes draw from ``pairing_stream``,
  a child of their seed.
- Cost estimate: ``expected_pairing_draws`` is inf where the exponential
  would overflow (k >= 54), never an OverflowError.
- Engine adapters: ``kRegularGraph(method='pairing')`` gives the same edges
  on NX and GT; the native NX draw is reproducible from its seed.
"""

from __future__ import annotations

import math
from collections import Counter

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.stats import chi2

from lrgsglib.graphs import kRegularGraph
from lrgsglib.graphs._shared._regular_pairing import (
    canonical_simple_edges,
    expected_pairing_draws,
    pairing_stream,
    random_regular_edges,
)

N_CUBIC_6 = 70  # labelled 3-regular graphs on 6 nodes (OEIS A002829)
DRAWS_PER_GRAPH = 200
CHI2_P_MIN = 1e-4
CYCLE_N, CYCLE_K, CYCLE_DRAWS = 1000, 3, 600
CYCLE_Z = 5.0  # tolerance in standard errors of the Poisson mean
SEED = 20260930
ADAPTER_N, ADAPTER_K, ADAPTER_SEED = 60, 3, 11
ENGINES = ["nx", "gt"]


def _degrees(edges: np.ndarray, n: int) -> np.ndarray:
    return np.bincount(edges.ravel(), minlength=n)


def _cycle_counts(edges: np.ndarray, n: int, k: int) -> tuple[float, float]:
    """(C3, C4) of a simple k-regular graph from traces of A^3 and A^4."""
    ones = np.ones(len(edges))
    A = sp.csr_matrix((ones, (edges[:, 0], edges[:, 1])), shape=(n, n))
    A = A + A.T
    A2 = A @ A
    c3 = A2.multiply(A).sum() / 6
    trace_a4 = A2.multiply(A2).sum()
    c4 = (trace_a4 - 2 * n * k**2 + n * k) / 8
    return float(c3), float(c4)


def _kregular(engine: str, **kwargs):
    """``kRegularGraph`` on ``engine``; a skip if graph-tool is missing."""
    if engine == "gt":
        pytest.importorskip("graph_tool")
    return kRegularGraph(engine=engine, **kwargs)


def _edges_and_weights(graph) -> tuple[np.ndarray, np.ndarray]:
    """Canonical edges (u < v, sorted) and the weights of a library graph."""
    if hasattr(graph, "get_edge_arrays"):  # GT bulk accessor
        src, dst, weights = graph.get_edge_arrays()
        pairs = np.column_stack((src, dst))
    else:
        data = list(graph.gr["G"].edges(data="weight"))
        pairs = np.array([(u, v) for u, v, _ in data])
        weights = np.array([w for _, _, w in data])
    return canonical_simple_edges(pairs, graph.N), np.asarray(weights)


# ---------------------------------------------------------------------------
# Exactness and invariants
# ---------------------------------------------------------------------------
def test_uniform_over_the_70_labelled_cubic_graphs_on_6_nodes():
    rng = np.random.default_rng(SEED)
    draws = N_CUBIC_6 * DRAWS_PER_GRAPH
    counts = Counter(
        random_regular_edges(6, 3, rng).tobytes() for _ in range(draws)
    )
    assert len(counts) == N_CUBIC_6
    observed = np.array(list(counts.values()))
    stat = np.sum((observed - DRAWS_PER_GRAPH) ** 2) / DRAWS_PER_GRAPH
    assert chi2.sf(stat, N_CUBIC_6 - 1) > CHI2_P_MIN


@pytest.mark.parametrize("n,k", [(10, 3), (50, 4), (9, 2), (12, 1), (5, 0)])
def test_simple_regular_sorted(n, k):
    edges = random_regular_edges(n, k, SEED)
    assert edges.shape == (n * k // 2, 2) and edges.dtype == np.int64
    if k:
        assert np.all(_degrees(edges, n) == k)
    assert np.all(edges[:, 0] < edges[:, 1])
    key = edges[:, 0] * n + edges[:, 1]
    assert np.all(np.diff(key) > 0)  # sorted and no repeated edge


def test_same_seed_same_graph_and_generator_is_consumed():
    assert np.array_equal(
        random_regular_edges(100, 3, 7), random_regular_edges(100, 3, 7)
    )
    rng = np.random.default_rng(7)
    first, second = random_regular_edges(100, 3, rng), random_regular_edges(
        100, 3, rng
    )
    assert not np.array_equal(first, second)


def test_pairing_stream_is_a_child_of_the_seed():
    """The graph classes seed their flips with the bare seed, so the
    topology draw must come from a different stream."""
    assert pairing_stream(None) is None
    child = random_regular_edges(ADAPTER_N, ADAPTER_K, pairing_stream(SEED))
    again = random_regular_edges(ADAPTER_N, ADAPTER_K, pairing_stream(SEED))
    bare = random_regular_edges(ADAPTER_N, ADAPTER_K, SEED)
    assert np.array_equal(child, again) and not np.array_equal(child, bare)


@pytest.mark.parametrize("n,k", [(5, 3), (4, 4), (6, -1), (0, 0)])
def test_impossible_parameters_raise(n, k):
    with pytest.raises(ValueError):
        random_regular_edges(n, k, SEED)


def test_attempt_cap_raises():
    # k = n - 1 forces the complete graph, which one draw almost never hits
    with pytest.raises(RuntimeError, match="draws"):
        random_regular_edges(8, 7, SEED, max_attempts=1)


def test_expected_draws_saturate_at_inf():
    """exp((k^2 - 1) / 4) exceeds the largest float from k = 54 on."""
    assert math.isfinite(expected_pairing_draws(53))
    assert expected_pairing_draws(54) == math.inf


def test_canonical_simple_edges_rejects_loops_and_multi_edges():
    assert canonical_simple_edges(np.array([[1, 1], [0, 2]]), 3) is None
    assert canonical_simple_edges(np.array([[0, 1], [1, 0]]), 3) is None
    edges = canonical_simple_edges(np.array([[2, 0], [1, 0]]), 3)
    assert edges.tolist() == [[0, 1], [0, 2]]
    with pytest.raises(ValueError):
        canonical_simple_edges(np.array([[0, 3]]), 3)


# ---------------------------------------------------------------------------
# Short cycles: uniform-ensemble limits (k-1)^j / (2j), a literature result
# not in .agents/ref, checked numerically here
# ---------------------------------------------------------------------------
@pytest.mark.physical
def test_short_cycle_means_match_the_uniform_ensemble():
    rng = np.random.default_rng(SEED)
    counts = np.array(
        [
            _cycle_counts(
                random_regular_edges(CYCLE_N, CYCLE_K, rng), CYCLE_N, CYCLE_K
            )
            for _ in range(CYCLE_DRAWS)
        ]
    )
    for j, column in ((3, 0), (4, 1)):
        mean = (CYCLE_K - 1) ** j / (2 * j)
        stderr = np.sqrt(mean / CYCLE_DRAWS)
        assert abs(counts[:, column].mean() - mean) < CYCLE_Z * stderr


# ---------------------------------------------------------------------------
# Engine adapters: kRegularGraph(method=...)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("engine", ENGINES)
def test_pairing_gives_the_same_graph_on_both_engines(engine):
    graph = _kregular(
        engine, n=ADAPTER_N, k=ADAPTER_K, seed=ADAPTER_SEED, method="pairing"
    )
    edges, weights = _edges_and_weights(graph)
    expected = random_regular_edges(
        ADAPTER_N, ADAPTER_K, pairing_stream(ADAPTER_SEED)
    )
    assert np.array_equal(edges, expected)  # the same reference on NX and GT
    assert np.all(weights == 1)


def test_seeded_nx_native_is_reproducible():
    def draw(seed: int) -> np.ndarray:
        graph = kRegularGraph(ADAPTER_N, ADAPTER_K, seed=seed, engine="nx")
        edges, _ = _edges_and_weights(graph)
        assert np.all(_degrees(edges, ADAPTER_N) == ADAPTER_K)
        return edges

    first = draw(ADAPTER_SEED)
    assert np.array_equal(first, draw(ADAPTER_SEED))
    assert not np.array_equal(first, draw(ADAPTER_SEED + 1))


@pytest.mark.parametrize("engine", ENGINES)
def test_syshape_path_names_a_non_default_method(engine):
    shape = f"N={ADAPTER_N}_k={ADAPTER_K}"
    for method, path in (
        ("native", shape),
        ("pairing", f"{shape}_method=pairing"),
    ):
        graph = _kregular(
            engine, n=ADAPTER_N, k=ADAPTER_K, seed=ADAPTER_SEED, method=method
        )
        assert graph.syshapePth == path


@pytest.mark.parametrize("engine", ENGINES)
def test_unknown_method_is_rejected(engine):
    with pytest.raises(ValueError, match="method"):
        _kregular(engine, n=10, k=3, seed=1, method="bogus")

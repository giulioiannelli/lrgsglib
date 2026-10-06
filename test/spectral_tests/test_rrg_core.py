"""
TransCluster v2 core (kernels/TransClusterCore) and RRG topology (kernels/RRG).

Checks:
- exact count rule and seed tree (integer keys only, rep range);
- nested disorder, and units that do not depend on their companions;
- record rows against a small-N reference: P(S) rebuilds the finite-cluster
  sums, the gap fields match a dense solve, n_neg = cnt, the nev column;
- p = 0 limit: lambda_0 = 0 and P_inf = 1;
- P_inf equals the library NX ``compute_pinf`` (tie rows skipped);
- SET-X is gauge-equivalent to S on the X-X edges;
- for mu = 0, the G family does not depend on sigma;
- majority and tie rules, neither of which depends on the solver's sign;
  dense vs Lanczos;
- solver attempts: matvec cap, larger-basis retry, gap fallback, solver
  errors, degenerate eigenpairs (NaN, inf, x = 0) and wrong pairs (residual
  test, also on tiny couplings) as a status, the last failure wins; bad
  solver knobs, control values, G coupling scales and couplings (also
  overflowing degree sums) are refused;
- RRG draws are connected, k-regular and reproducible on both engines, a
  disconnected draw is redrawn, bad generator output fails loudly;
- sampler guard (pairing cost) and fail-loud engine check;
- paper anchor (slow, skipped under --quick): on the square lattice at
  N = 4096 the chi peak and P_inf match iannelli2025topological, Fig. 2(c).
"""

from __future__ import annotations

import copy
import dataclasses
import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
import scipy.linalg as sla
from scipy.sparse.linalg import ArpackError

ROOT_DIR = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT_DIR / "src"

# Add src/ to path so `from kernels... import ...` works
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import kernels.RRG as rrg  # noqa: E402
import kernels.TransClusterCore as core  # noqa: E402
from kernels.RRG import draw_topology, resolve_sampler  # noqa: E402
from kernels.TransClusterCore import (  # noqa: E402
    Family,
    GroundState,
    SolverConfig,
    canonical_edges,
    check_ctrl,
    component_labels,
    couplings,
    draw_disorder,
    edge_fingerprint,
    engine_seed,
    exact_count,
    majority_signs,
    make_record,
    measure_unit,
    orient_tie,
    sign_clusters,
    signed_laplacian,
    solve_attempts,
    stream_rng,
    stream_seed,
    sweep_rep,
)
from lrgsglib.config.progargs.defs.RandomRegular import (  # noqa: E402
    RRG_MAX_ATTEMPTS,
    RRG_PAIRING_K_MAX,
    RRG_PAIRING_MAX_EXPECTED_DRAWS,
)
from lrgsglib.config.progargs.defs.TransCluster import (  # noqa: E402
    DEFAULT_TC_NEV,
    TC_EIG_RETRY_NCV_FACTOR,
    TC_ENSEMBLE_BERNOULLI,
    TC_FAMILY_GAUSS,
    TC_FAMILY_SIGNS,
    TC_FAMILY_XNODES,
    TC_G_SCALE_MAX,
    TC_G_SCALE_MIN,
    TC_GAP_NEV,
    TC_GERSHGORIN_FACTOR,
    TC_NO_COUNT,
    TC_NUNREL_EPS_FACTOR,
    TC_SEED_SPAN,
    TC_STATUS_ERROR,
    TC_STATUS_NOCONV,
    TC_STATUS_NONFINITE,
    TC_STATUS_OK,
    TC_STREAM_V0,
)
from lrgsglib.graphs._shared._regular_pairing import (  # noqa: E402
    expected_pairing_draws,
    pairing_stream,
    random_regular_edges,
)

ROOT = 12345
N_SMALL, N_LANCZOS, K = 120, 400, 3
P_GRID = (0.0, 0.05, 0.1, 0.15, 0.2, 0.3)
GAP_EVERY = 2  # rep 0 is a gap rep, rep 1 is not
LIB_N, LIB_SEEDS, LIB_P = 200, range(1, 5), (0.05, 0.12, 0.2, 0.35)
GAUGE_TOL = 1e-12
EIG_TOL = 1e-10
# Solver attempts at N_LANCZOS: a basis this small needs several times the
# matvecs of twice its size, so the larger-basis retry has room to rescue.
SMALL_NCV = 6
STARVED_MV = 5  # far below what any attempt needs
TINY_SCALE = 1e-26  # coupling scale where ARPACK's own test is absolute
GAP_FIELDS = ("lam1", "gap", "dk", "nunrel_dk", "e1")
# Paper anchor: square lattice, S defects (iannelli2025topological.pdf,
# Fig. 2(c) and p. 4). N = 4096 is one of the figure's sizes. The paper does
# not state the boundary conditions: PBC is inferred, as it matches the
# digitised N = 4096 curve and open BC does not.
ANCHOR_L, ANCHOR_REPS = 64, 200
ANCHOR_P = np.round(np.arange(0.07, 0.135, 0.01), 2)
ANCHOR_PC, ANCHOR_PC_TOL = 0.10, 0.02  # "p_c ~ 0.10(1)"; tolerance 2 x error
ANCHOR_PINF = {0.10: 0.731, 0.12: 0.602}  # N = 4096 curve, read off the figure
ANCHOR_PINF_TOL = 0.03  # figure reading +-0.005, sampling error ~0.006
PEAK_FIT_HALF_WIDTH = 2  # parabola through the argmax and 2 points each side
# N = 8 exercises the redraw path: the only disconnected cubic graph on 8
# nodes is two disjoint K4 (each component is cubic, so it has >= 4 nodes,
# and K4 is the only cubic graph on 4). REDRAW_REP was found by scanning
# reps: its attempt-0 pairing draw is two K4.
TINY_N, REDRAW_ROOT, REDRAW_REP = 8, 77, 21
K4 = np.array(list(itertools.combinations(range(4), 2)))
TWO_K4 = np.vstack((K4, K4 + 4))
CUBE = np.array(
    [(u, u ^ b) for u in range(TINY_N) for b in (1, 2, 4) if u < u ^ b]
)
# The cube with edge (6, 7) moved to (5, 6): M = 12 but degrees 4 and 2
REWIRED = np.vstack((CUBE[~np.all(CUBE == (6, 7), axis=1)], [(5, 6)]))


def _rrg(n: int, seed: int = 1) -> np.ndarray:
    return random_regular_edges(n, K, seed)


def _unit(edges, n, family, ctrl, rep=0, solver=SolverConfig()):
    """Fields of one unit, measured the way sweep_rep does."""
    draw = draw_disorder(family, n, len(edges), ROOT, (n, K), rep)
    J, disorder = couplings(family, draw, edges, ctrl)
    v0 = stream_rng(ROOT, (n, K), rep, TC_STREAM_V0).standard_normal(n)
    fields, ps = measure_unit(edges, J, n, v0, solver=solver)
    return J, {**disorder, **fields}, ps


def _lanczos_unit(p: float = 0.2):
    """Edges, S couplings at ``p`` and v0 of rep 0 at N_LANCZOS."""
    edges = _rrg(N_LANCZOS)
    draw = draw_disorder(
        Family("S"), N_LANCZOS, len(edges), ROOT, (N_LANCZOS, K), 0
    )
    J, _ = couplings(Family("S"), draw, edges, p)
    v0 = stream_rng(ROOT, (N_LANCZOS, K), 0, TC_STREAM_V0).standard_normal(
        N_LANCZOS
    )
    return edges, J, v0


def _solve(edges, J, v0, nev, **solver):
    """Fields of one Lanczos unit with the given solver knobs."""
    cfg = SolverConfig(dense_max_n=0, **solver)
    return measure_unit(edges, J, N_LANCZOS, v0, nev, cfg)[0]


def _path(n: int) -> np.ndarray:
    """Edges of the path 0 - 1 - ... - (n - 1)."""
    return np.column_stack((np.arange(n - 1), np.arange(1, n)))


def _engine(engine: str) -> str:
    """``engine``, or a skip if it is graph-tool and that is missing."""
    if engine == "gt":
        pytest.importorskip("graph_tool")
    return engine


# ---------------------------------------------------------------------------
# Count rule and seed tree
# ---------------------------------------------------------------------------
def test_exact_count_rounds_half_up_on_the_decimal_value():
    assert exact_count(0.29, 1500) == 435  # int(0.29 * 1500) == 434
    assert exact_count(0.5, 3) == 2
    assert exact_count(0.1, 3) == 0
    assert exact_count(0.0, 7) == 0 and exact_count(1.0, 7) == 7
    with pytest.raises(ValueError):
        exact_count(1.5, 10)


def test_seed_tree_is_injective_bounded_and_deterministic():
    seeds = [engine_seed(ROOT, (1000, K), rep) for rep in range(1000)]
    assert len(set(seeds)) == len(seeds)
    assert all(1 <= s <= TC_SEED_SPAN for s in seeds)
    assert seeds[5] == engine_seed(ROOT, (1000, K), 5)
    assert engine_seed(ROOT, (1000, K), 5, attempt=1) != seeds[5]
    with pytest.raises(ValueError, match="spawn-key"):
        stream_seed(ROOT, (2**32, K), 0, TC_STREAM_V0)


def test_seed_keys_are_integers_and_reps_lie_in_range():
    key = (N_SMALL, K)
    numpy_key = (np.int64(N_SMALL), np.int64(K))
    assert engine_seed(ROOT, numpy_key, np.int64(3)) == engine_seed(
        ROOT, key, 3
    )
    with pytest.raises(TypeError):  # a float is refused, never truncated
        engine_seed(ROOT, key, 3.0)
    with pytest.raises(TypeError):
        stream_seed(ROOT, (N_SMALL + 0.5, K), 0, TC_STREAM_V0)
    with pytest.raises(TypeError):
        engine_seed(float(ROOT), key, 0)
    for rep in (-1, TC_SEED_SPAN):
        with pytest.raises(ValueError, match="rep must lie"):
            engine_seed(ROOT, key, rep)
        with pytest.raises(ValueError, match="rep must lie"):
            stream_seed(ROOT, key, rep, TC_STREAM_V0)
    assert 1 <= engine_seed(ROOT, key, TC_SEED_SPAN - 1) <= TC_SEED_SPAN


def test_edge_fingerprint_is_canonical():
    edges = _rrg(N_SMALL)
    assert edge_fingerprint(edges) == edge_fingerprint(edges.copy())
    assert edge_fingerprint(edges) != edge_fingerprint(_rrg(N_SMALL, 2))


# ---------------------------------------------------------------------------
# Nested disorder and tiling invariance
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "family",
    [
        Family("S"),
        Family("S", "bernoulli"),
        Family("X"),
        Family("X", "bernoulli"),
    ],
)
def test_disorder_is_nested_in_the_control_value(family):
    edges = _rrg(N_SMALL)
    draw = draw_disorder(family, N_SMALL, len(edges), ROOT, (N_SMALL, K), 3)
    previous = np.zeros(len(edges), dtype=bool)
    for ctrl in P_GRID:
        J, fields = couplings(family, draw, edges, ctrl)
        negative = J < 0
        assert np.all(negative[previous])  # superset of the previous set
        assert fields["n_neg"] == negative.sum()
        previous = negative
    canonical = family.ensemble == "canonical"
    assert (fields["cnt"] != TC_NO_COUNT) == canonical


def test_unit_does_not_depend_on_its_companions():
    edges = _rrg(N_LANCZOS)
    kwargs = dict(root=ROOT, graph_key=(N_LANCZOS, K), rep=4, gap_every=2)
    together = [
        rec
        for rec, _ in sweep_rep(
            edges, N_LANCZOS, Family("S"), P_GRID, range(len(P_GRID)), **kwargs
        )
    ]
    for ictrl in (5, 2):
        ((alone, _),) = sweep_rep(
            edges, N_LANCZOS, Family("S"), P_GRID, [ictrl], **kwargs
        )
        for name in alone.dtype.names:
            if name != "t_eig":
                np.testing.assert_array_equal(
                    alone[name], together[ictrl][name]
                )


# ---------------------------------------------------------------------------
# Record rows against a small-N reference
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n", [N_SMALL, N_LANCZOS])  # dense and Lanczos
def test_record_rows_match_the_small_n_reference(n):
    edges, family = _rrg(n), Family("S")
    draw = draw_disorder(family, n, len(edges), ROOT, (n, K), 0)
    rows = sweep_rep(
        edges,
        n,
        family,
        P_GRID,
        range(len(P_GRID)),
        root=ROOT,
        graph_key=(n, K),
        rep=0,
        gap_every=GAP_EVERY,
    )
    for rec, (sizes, counts) in rows:
        assert rec["status"] == TC_STATUS_OK and rec["nev"] == TC_GAP_NEV
        assert rec["n_neg"] == rec["cnt"]  # exact count rule
        # P(S) excludes the giant and rebuilds the finite-cluster fields
        assert np.sum(sizes * counts) == rec["sum_s_p"]
        assert np.sum(sizes**2 * counts) == rec["sum_s2_p"]
        assert counts.sum() == rec["ncl_p"] - 1
        assert max(sizes, default=0) == rec["s2p"] <= rec["s1p"]
        assert rec["n_pos"] == rec["s1p"] + rec["sum_s_p"] >= n / 2
        assert rec["pinf"] == rec["s1p"] / n
        # gap fields against a dense solve of the same Laplacian
        J, _ = couplings(family, draw, edges, rec["ctrl"])
        lap, _ = signed_laplacian(edges, J, n)
        lam = sla.eigvalsh(lap.toarray())[:TC_GAP_NEV]
        np.testing.assert_allclose(
            [rec["lam0"], rec["lam1"]], lam, rtol=0, atol=EIG_TOL
        )
        assert rec["gap"] == rec["lam1"] - rec["lam0"] > 0
        assert rec["dk"] == pytest.approx(rec["resid"] / rec["gap"])


def test_only_gap_reps_solve_the_second_eigenpair():
    edges = _rrg(N_SMALL)
    kwargs = dict(root=ROOT, graph_key=(N_SMALL, K), gap_every=GAP_EVERY)
    for rep, gap_rep in ((0, True), (1, False)):
        ((rec, _),) = sweep_rep(
            edges, N_SMALL, Family("S"), P_GRID, [3], rep=rep, **kwargs
        )
        assert rec["nev"] == (TC_GAP_NEV if gap_rep else DEFAULT_TC_NEV)
        for name in ("lam1", "gap", "dk", "nunrel_dk", "e1"):
            assert np.isnan(rec[name]) != gap_rep


# ---------------------------------------------------------------------------
# Physics limits and library parity
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("n", [N_SMALL, N_LANCZOS])
def test_zero_disorder_is_the_ferromagnetic_ground_state(n):
    _, fields, ps = _unit(_rrg(n), n, Family("S"), 0.0)
    assert fields["status"] == TC_STATUS_OK
    assert abs(fields["lam0"]) < EIG_TOL
    assert fields["pinf"] == 1.0 and fields["n_neg"] == 0
    assert fields["e0"] == -1.0 and fields["tie"] == 0
    assert ps[0].size == 0


@pytest.mark.physical
def test_pinf_equals_library_nx_compute_pinf():
    from lrgsglib.graphs import kRegularGraph

    checked = 0
    for seed in LIB_SEEDS:
        # seed= also seeds the native NX generator: a reproducible graph
        reference = kRegularGraph(LIB_N, K, pflip=0.0, seed=seed, engine="nx")
        n, edges = canonical_edges(reference)
        for p in LIB_P:
            J, fields, _ = _unit(edges, n, Family("S"), p, rep=seed)
            if fields["tie"]:
                continue
            graph = copy.deepcopy(reference)  # compute_pinf caches clusters
            graph.flip_sel_edges([tuple(e) for e in edges[J < 0]])
            graph.compute_k_eigvV(k=1, which="SA")
            graph.compute_pinf()
            assert graph.Pinf == pytest.approx(fields["pinf"], abs=0)
            checked += 1
    assert checked >= len(LIB_SEEDS) * len(LIB_P) - 1


def _chi_peak(ctrl: np.ndarray, chi: np.ndarray) -> float:
    """Vertex of the least-squares parabola around the argmax of ``chi``."""
    i = int(np.argmax(chi))
    lo = max(0, i - PEAK_FIT_HALF_WIDTH)
    window = slice(lo, i + PEAK_FIT_HALF_WIDTH + 1)
    a, b, _ = np.polyfit(ctrl[window], chi[window], 2)
    return -b / (2 * a)


@pytest.mark.slow
@pytest.mark.physical
def test_square_lattice_reproduces_the_published_transition(quick_mode):
    if quick_mode:
        pytest.skip("paper anchor: 1-3 min of eigensolves")
    from lrgsglib.graphs import Lattice2D

    lattice = Lattice2D(ANCHOR_L, geo="sqr", pbc=True, pflip=0.0, engine="nx")
    n, edges = canonical_edges(lattice)
    pinf = np.full((ANCHOR_REPS, ANCHOR_P.size), np.nan)
    for rep in range(ANCHOR_REPS):
        for rec, _ in sweep_rep(
            edges,
            n,
            Family("S"),
            ANCHOR_P,
            range(ANCHOR_P.size),
            root=ROOT,
            graph_key=(ANCHOR_L,),
            rep=rep,
            gap_every=0,
        ):
            pinf[rep, rec["ictrl"]] = rec["pinf"]
    chi = np.sqrt(n) * pinf.std(axis=0, ddof=1)  # chi = sqrt(N) sigma (p. 4)
    assert abs(_chi_peak(ANCHOR_P, chi) - ANCHOR_PC) <= ANCHOR_PC_TOL
    for p, published in ANCHOR_PINF.items():
        mean = pinf[:, np.isclose(ANCHOR_P, p)].mean()
        assert mean == pytest.approx(published, abs=ANCHOR_PINF_TOL)


def test_set_x_is_gauge_equivalent_to_s_on_xx_edges():
    edges = _rrg(N_SMALL)
    family = Family("X")
    draw = draw_disorder(family, N_SMALL, len(edges), ROOT, (N_SMALL, K), 0)
    J, _ = couplings(family, draw, edges, 0.3)
    is_x = np.zeros(N_SMALL, dtype=bool)
    is_x[draw.order[: exact_count(0.3, N_SMALL)]] = True
    xx_only = np.where(is_x[edges[:, 0]] & is_x[edges[:, 1]], -1.0, 1.0)
    spectra = [
        sla.eigvalsh(signed_laplacian(edges, j, N_SMALL)[0].toarray())
        for j in (J, xx_only)
    ]
    np.testing.assert_allclose(spectra[0], spectra[1], atol=GAUGE_TOL)


def test_gauss_family_with_zero_mean_ignores_sigma():
    edges = _rrg(N_LANCZOS)
    family = Family("G", mu=0.0)
    (_, weak, _), (_, strong, _) = (
        _unit(edges, N_LANCZOS, family, sigma) for sigma in (0.3, 3.0)
    )
    assert weak["pinf"] == strong["pinf"] and weak["s1n"] == strong["s1n"]
    assert strong["lam0"] == pytest.approx(10 * weak["lam0"], rel=1e-8)
    assert weak["cnt"] == TC_NO_COUNT and weak["min_absJ"] > 0


# ---------------------------------------------------------------------------
# Sign conventions and solver paths
# ---------------------------------------------------------------------------
MAJORITY = {  # x -> signs; zeros count for neither side and map to +1
    "negative majority": ((-1.0, -2.0, 0.0, 3.0, -0.5), [1, 1, 1, -1, 1]),
    "not N / 2": ((-1.0, 0.0, 2.0, -3.0), [1, 1, -1, 1]),
    "zeros": ((-1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0), [1, 1, 1, 1, 1, -1, -1]),
}


@pytest.mark.parametrize("x, signs", MAJORITY.values(), ids=MAJORITY)
def test_majority_flip_is_the_compute_k_eigvV_convention(x, signs):
    """Flip when #(x < 0) > #(x > 0), then zeros to +1, as compute_k_eigvV
    does; x and -x give the same signs."""
    for sign in (1.0, -1.0):
        got, tie = majority_signs(sign * np.array(x))
        assert got.tolist() == signs and not tie


def test_the_tie_flag_compares_the_two_signs():
    signs, tie = majority_signs(np.array([-1.0, 0.0, 2.0]))
    assert signs.tolist() == [-1, 1, 1] and tie


def test_tie_is_oriented_toward_the_larger_giant():
    path = np.array([[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
    signs = np.array([1, -1, 1, -1, -1, 1], dtype=np.int8)  # 3 vs 3
    clusters = sign_clusters(path, signs, 6)
    assert clusters.largest(1) == 1 and clusters.largest(-1) == 2
    oriented, flipped = orient_tie(signs, clusters, tie=True)
    assert flipped.largest(1) == 2 and oriented.tolist() == (-signs).tolist()
    kept, same = orient_tie(signs, clusters, tie=False)
    assert kept is signs and same is clusters


TIE_ORIENTATION = {  # tie on a path: signs -> oriented signs
    # the larger giant (3 > 2) goes positive, though its second cluster is
    # the smaller (1 < 2)
    "giant": ([-1, -1, -1, 1, 1, -1, 1, 1], [1, 1, 1, -1, -1, 1, -1, -1]),
    # equal giants (3), so the larger second cluster (2 > 1) goes positive
    "second clusters": (
        [1, 1, 1, -1, -1, -1, 1, -1, -1, 1],
        [-1, -1, -1, 1, 1, 1, -1, 1, 1, -1],
    ),
    # equal size lists ([2] and [2]), so node 0 goes positive
    "node 0": ([-1, -1, 1, 1], [1, 1, -1, -1]),
}


@pytest.mark.parametrize(
    "signs, oriented", TIE_ORIENTATION.values(), ids=TIE_ORIENTATION
)
def test_tie_orientation_does_not_depend_on_the_solver_sign(signs, oriented):
    """The rule reads only the partition, which x and -x share."""
    signs = np.array(signs, dtype=np.int8)
    path = _path(signs.size)
    for s in (signs, -signs):
        got, clusters = orient_tie(s, sign_clusters(path, s, s.size), tie=True)
        assert got.tolist() == oriented
        assert clusters.ranked(1) == sign_clusters(path, got, s.size).ranked(1)


def test_a_tie_unit_does_not_depend_on_the_solver_sign(monkeypatch):
    """x and -x give the same record and P(S) at a tie with equal giants
    (the "second clusters" case above)."""
    signs = np.array(TIE_ORIENTATION["second clusters"][0], dtype=np.float64)
    n = signs.size
    x = signs * np.linspace(1.0, 2.0, n)  # no zeros
    x /= np.linalg.norm(x)
    measured = []
    for sign in (1.0, -1.0):
        state = GroundState(
            np.zeros(1), sign * x[:, None], TC_STATUS_OK, 1, 0.0
        )
        monkeypatch.setattr(core, "ground_state", lambda *a, gs=state: gs)
        measured.append(measure_unit(_path(n), np.ones(n - 1), n, np.ones(n)))
    (fields, ps), (flipped_fields, flipped_ps) = measured
    row, flipped = make_record(fields), make_record(flipped_fields)
    assert row["tie"] == 1 and row["s1p"] == row["s1n"] and row["s2p"] == 2
    for name in row.dtype.names:
        np.testing.assert_array_equal(row[name], flipped[name])
    for a, b in zip(ps, flipped_ps):
        np.testing.assert_array_equal(a, b)


def test_lanczos_agrees_with_dense():
    edges = _rrg(N_LANCZOS)
    dense, lanczos = (
        _unit(edges, N_LANCZOS, Family("S"), 0.12, solver=SolverConfig(**cfg))[
            1
        ]
        for cfg in ({"dense_max_n": N_LANCZOS}, {"dense_max_n": 0})
    )
    assert dense["n_mv"] == -1 and lanczos["n_mv"] > 0
    assert lanczos["lam0"] == pytest.approx(dense["lam0"], abs=EIG_TOL)
    assert lanczos["pinf"] == dense["pinf"]
    assert lanczos["resid"] < EIG_TOL


def test_solve_attempts_end_with_the_nev_1_ladder():
    retry = 40 * TC_EIG_RETRY_NCV_FACTOR
    assert solve_attempts(1, 40) == [(1, 40), (1, retry)]
    assert solve_attempts(2, 40) == [(2, 40), (1, 40), (1, retry)]


@pytest.mark.parametrize("nev", [DEFAULT_TC_NEV, TC_GAP_NEV])
def test_each_attempt_stops_at_its_matvec_cap(nev):
    """No convergence is a status, not an exception, and the cost of a
    unit is bounded exactly: max_mv matvecs per attempt."""
    edges, J, v0 = _lanczos_unit()
    measured = measure_unit(
        edges, J, N_LANCZOS, v0, nev, SolverConfig(max_mv=STARVED_MV)
    )
    fields, ps = measured
    n_attempts = len(solve_attempts(nev, SolverConfig().ncv))
    assert fields["status"] == TC_STATUS_NOCONV and fields["nev"] == nev
    assert fields["n_mv"] == STARVED_MV * n_attempts
    row = make_record(fields)
    assert np.isnan(row["pinf"]) and np.isnan(row["lam0"]) and row["tie"] == -1
    assert ps[0].size == 0


def test_the_larger_basis_retry_rescues_a_starved_attempt():
    edges, J, v0 = _lanczos_unit()
    large_ncv = SMALL_NCV * TC_EIG_RETRY_NCV_FACTOR
    small = _solve(edges, J, v0, 1, ncv=SMALL_NCV)["n_mv"]
    large = _solve(edges, J, v0, 1, ncv=large_ncv)
    assert large["n_mv"] < small  # the premise of this test
    cap = (large["n_mv"] + small) // 2
    rescued = _solve(edges, J, v0, 1, ncv=SMALL_NCV, max_mv=cap)
    assert rescued["status"] == TC_STATUS_OK
    assert rescued["n_mv"] == cap + large["n_mv"]
    assert rescued["lam0"] == large["lam0"] and rescued["pinf"] == large["pinf"]


def test_a_failed_gap_solve_still_measures_lambda_0():
    """The nev = 2 attempt runs out of matvecs. The nev = 1 attempt after it
    is the one any other unit makes, so lambda_0 and P_inf are the same."""
    edges, J, v0 = _lanczos_unit()
    one = _solve(edges, J, v0, DEFAULT_TC_NEV)
    two = _solve(edges, J, v0, TC_GAP_NEV)
    assert one["n_mv"] < two["n_mv"] and two["nev"] == TC_GAP_NEV
    cap = (one["n_mv"] + two["n_mv"]) // 2
    gap = _solve(edges, J, v0, TC_GAP_NEV, max_mv=cap)
    assert gap["status"] == TC_STATUS_OK and gap["nev"] == DEFAULT_TC_NEV
    assert gap["n_mv"] == cap + one["n_mv"]
    for name in ("lam0", "resid", "pinf", "s1p", "sum_s_p", "tie", "e0"):
        assert gap[name] == one[name]
    row = make_record(gap)
    assert all(np.isnan(row[name]) for name in GAP_FIELDS)


@pytest.mark.parametrize("n_failing", [1, 3])
def test_arpack_errors_become_a_status(n_failing, monkeypatch):
    """ARPACK info 3 ("no shifts could be applied") is a plain ArpackError;
    the larger-basis retry is its cure, and if every attempt raises the
    unit gets TC_STATUS_ERROR."""
    edges, J, v0 = _lanczos_unit()
    calls = []

    def flaky(*args, **kwargs):
        calls.append((kwargs["k"], kwargs["ncv"]))
        if len(calls) <= n_failing:
            raise ArpackError(3)
        return real_eigsh(*args, **kwargs)

    real_eigsh = core.eigsh
    monkeypatch.setattr(core, "eigsh", flaky)
    fields = _solve(edges, J, v0, DEFAULT_TC_NEV)
    attempts = solve_attempts(DEFAULT_TC_NEV, SolverConfig().ncv)
    assert calls == attempts[: min(n_failing + 1, len(attempts))]
    rescued = n_failing < len(attempts)
    assert fields["status"] == (TC_STATUS_OK if rescued else TC_STATUS_ERROR)


def test_dense_solver_errors_become_a_status(monkeypatch):
    def broken(*args, **kwargs):
        raise np.linalg.LinAlgError("eigh did not converge")

    monkeypatch.setattr(core.sla, "eigh", broken)
    _, fields, ps = _unit(_rrg(N_SMALL), N_SMALL, Family("S"), 0.2)
    assert fields["status"] == TC_STATUS_ERROR and fields["n_mv"] == -1
    assert np.isnan(make_record(fields)["pinf"]) and ps[0].size == 0


def _scripted_eigsh(monkeypatch, plan: list[str]) -> list[tuple[int, int]]:
    """Make ``core.eigsh`` play ``plan``, one outcome per call (the last one
    repeats): "solve" runs ARPACK, "err" raises its info 3, "nan" returns
    NaN pairs, "inf" an infinite eigenvalue with unit vectors (only the
    value check catches it), "zero" x = 0, as ARPACK did at sigma = 1e307,
    "zero_last" unit vectors but the last, x = 0 (only a per-vector norm
    check catches it), and "wrong" unit vectors that are not eigenvectors
    (only the residual test catches them). Returns the (k, ncv) of each
    call."""
    real_eigsh = core.eigsh
    calls = []

    def scripted(op, **kwargs):
        k = kwargs["k"]
        calls.append((k, kwargs["ncv"]))
        outcome = plan[min(len(calls), len(plan)) - 1]
        if outcome == "solve":
            return real_eigsh(op, **kwargs)
        if outcome == "err":
            raise ArpackError(3)
        if outcome == "inf":
            return np.full(k, np.inf), np.eye(op.shape[0], k)
        if outcome == "zero_last":
            vecs = np.eye(op.shape[0], k)
            vecs[:, -1] = 0.0
            return np.arange(k, dtype=float), vecs
        if outcome == "wrong":
            return np.arange(k, dtype=float), np.eye(op.shape[0], k)
        fill = np.nan if outcome == "nan" else 0.0
        return np.full(k, fill), np.full((op.shape[0], k), fill)

    monkeypatch.setattr(core, "eigsh", scripted)
    return calls


@pytest.mark.parametrize("outcome", ["nan", "inf", "zero"])
def test_degenerate_eigenpairs_become_a_status(outcome, monkeypatch):
    """An attempt that returns a non-finite eigenvalue or vector, or a vector
    of norm != 1, has failed; a unit where all do is NONFINITE, with no order
    parameter."""
    edges, J, v0 = _lanczos_unit()
    calls = _scripted_eigsh(monkeypatch, [outcome])
    fields, ps = measure_unit(edges, J, N_LANCZOS, v0)
    assert calls == solve_attempts(DEFAULT_TC_NEV, SolverConfig().ncv)
    assert fields["status"] == TC_STATUS_NONFINITE
    assert fields["nev"] == DEFAULT_TC_NEV
    row = make_record(fields)
    assert np.isnan(row["pinf"]) and np.isnan(row["lam0"]) and ps[0].size == 0


@pytest.mark.parametrize(
    "nev, outcome",
    [
        (DEFAULT_TC_NEV, "zero"),
        (TC_GAP_NEV, "zero"),
        (TC_GAP_NEV, "zero_last"),
        (DEFAULT_TC_NEV, "wrong"),
        (TC_GAP_NEV, "wrong"),
    ],
)
def test_a_degenerate_attempt_hands_over_to_the_next(nev, outcome, monkeypatch):
    """The unit is then bitwise the clean solve of the next attempt; one
    degenerate vector of a gap solve is enough to fail it, and so are
    wrong pairs."""
    edges, J, v0 = _lanczos_unit()
    attempts = solve_attempts(nev, SolverConfig().ncv)
    k, ncv = attempts[1]
    clean = make_record(_solve(edges, J, v0, k, ncv=ncv))
    calls = _scripted_eigsh(monkeypatch, [outcome, "solve"])
    row = make_record(_solve(edges, J, v0, nev))
    assert calls == attempts[:2] and row["status"] == TC_STATUS_OK
    for name in row.dtype.names:
        if name != "t_eig":
            np.testing.assert_array_equal(row[name], clean[name])


LAST_FAILURE = {  # outcomes of the attempts at a starved cap -> unit status
    "cap, then error": (["solve", "err"], TC_STATUS_ERROR),
    "error, then cap": (["err", "solve"], TC_STATUS_NOCONV),
    "NaN, then cap": (["nan", "solve"], TC_STATUS_NOCONV),
    "error, then wrong pairs": (["err", "wrong"], TC_STATUS_NOCONV),
}


@pytest.mark.parametrize(
    "plan, status", LAST_FAILURE.values(), ids=LAST_FAILURE
)
def test_the_last_failed_attempt_sets_the_status(plan, status, monkeypatch):
    """A starved "solve" stops at its matvec cap, i.e. NOCONV."""
    edges, J, v0 = _lanczos_unit()
    calls = _scripted_eigsh(monkeypatch, plan)
    fields = _solve(edges, J, v0, DEFAULT_TC_NEV, max_mv=STARVED_MV)
    assert calls == solve_attempts(DEFAULT_TC_NEV, SolverConfig().ncv)
    assert fields["status"] == status
    assert fields["n_mv"] == STARVED_MV * plan.count("solve")


def test_a_wrong_pair_on_tiny_couplings_is_not_converged():
    """At this scale ARPACK returns a wrong pair (its P_inf is off) as
    converged, and the residual test turns it into NOCONV. check_ctrl keeps
    the program's couplings far from it; measure_unit guards itself."""
    edges, J, v0 = _lanczos_unit()
    fields = _solve(edges, TINY_SCALE * J, v0, DEFAULT_TC_NEV)
    assert fields["status"] == TC_STATUS_NOCONV and "pinf" not in fields


def test_a_degenerate_dense_solve_becomes_a_status(monkeypatch):
    def zero_vectors(a, subset_by_index):
        nev = subset_by_index[1] + 1
        return np.zeros(nev), np.zeros((len(a), nev))

    monkeypatch.setattr(core.sla, "eigh", zero_vectors)
    _, fields, ps = _unit(_rrg(N_SMALL), N_SMALL, Family("S"), 0.2)
    assert fields["status"] == TC_STATUS_NONFINITE and fields["n_mv"] == -1
    assert np.isnan(make_record(fields)["pinf"]) and ps[0].size == 0


BAD_SOLVER = {
    "NaN tol": ("tol", np.nan),
    "infinite tol": ("tol", np.inf),
    "negative tol": ("tol", -1.0),
    "basis of a gap solve": ("ncv", TC_GAP_NEV),
    "no matvecs": ("max_mv", 0),
    "negative dense_max_n": ("dense_max_n", -1),
}


@pytest.mark.parametrize("knob, value", BAD_SOLVER.values(), ids=BAD_SOLVER)
def test_solver_config_refuses_bad_knobs(knob, value):
    """At construction, not inside the eigsh call of every unit."""
    with pytest.raises(ValueError, match=knob):
        SolverConfig(**{knob: value})


def test_solver_config_accepts_its_bounds_and_drops_negative_zero():
    """cfg.json holds the knobs; JSON prints -0.0, though -0.0 == 0.0."""
    SolverConfig(tol=0.0, ncv=TC_GAP_NEV + 1, max_mv=1, dense_max_n=0)
    signed, zero = (SolverConfig(tol=tol) for tol in (-0.0, 0.0))
    assert json.dumps(dataclasses.asdict(signed)) == json.dumps(
        dataclasses.asdict(zero)
    )


BAD_CTRL = {
    "bernoulli S, NaN": (TC_FAMILY_SIGNS, np.nan),
    "bernoulli S, below 0": (TC_FAMILY_SIGNS, -0.1),
    "bernoulli X, above 1": (TC_FAMILY_XNODES, 1.5),
    "G, negative sigma": (TC_FAMILY_GAUSS, -0.5),
    "G, NaN sigma": (TC_FAMILY_GAUSS, np.nan),
    "G, infinite sigma": (TC_FAMILY_GAUSS, np.inf),
}


@pytest.mark.parametrize("kind, ctrl", BAD_CTRL.values(), ids=BAD_CTRL)
def test_couplings_refuse_bad_control_values(kind, ctrl):
    """Every ensemble fails like the canonical one (whose count rule raises)."""
    ensemble = (
        {} if kind == TC_FAMILY_GAUSS else {"ensemble": TC_ENSEMBLE_BERNOULLI}
    )
    family = Family(kind, **ensemble)
    edges = _rrg(N_SMALL)
    draw = draw_disorder(family, N_SMALL, len(edges), ROOT, (N_SMALL, K), 0)
    with pytest.raises(ValueError, match="control values must lie in"):
        couplings(family, draw, edges, ctrl)


def test_vanishing_or_non_finite_couplings_are_refused():
    edges = _rrg(N_SMALL)
    zero = Family(TC_FAMILY_GAUSS, mu=0.0)
    draw = draw_disorder(zero, N_SMALL, len(edges), ROOT, (N_SMALL, K), 0)
    with pytest.raises(ValueError, match="make every coupling vanish"):
        couplings(zero, draw, edges, 0.0)
    for bad_mu in (np.nan, np.inf):
        with pytest.raises(ValueError, match="mu must be finite"):
            Family(TC_FAMILY_GAUSS, mu=bad_mu)
    v0 = np.ones(N_SMALL)
    for value, message in ((0.0, "couplings vanish"), (np.nan, "finite")):
        J = np.full(len(edges), value)
        with pytest.raises(ValueError, match=message):
            measure_unit(edges, J, N_SMALL, v0)


def test_couplings_whose_degree_sums_overflow_are_refused():
    edges = _rrg(N_SMALL)
    J = np.full(len(edges), np.finfo(np.float64).max)  # finite; K * J is not
    with pytest.raises(ValueError, match="degree sums"):
        measure_unit(edges, J, N_SMALL, np.ones(N_SMALL))


BELOW_G_SCALE = np.nextafter(TC_G_SCALE_MIN, 0.0)  # one ulp outside the window
ABOVE_G_SCALE = np.nextafter(TC_G_SCALE_MAX, np.inf)
BAD_SCALE = {  # G (mu, sigma) whose scale max(|mu|, sigma) is just outside
    "sigma below, mu 0": (0.0, BELOW_G_SCALE),
    "sigma above, mu 0": (0.0, ABOVE_G_SCALE),
    "mu below, sigma 0": (BELOW_G_SCALE, 0.0),
    "mu above": (ABOVE_G_SCALE, 1.0),
    "negative mu above": (-ABOVE_G_SCALE, 1.0),
}


@pytest.mark.parametrize("mu, sigma", BAD_SCALE.values(), ids=BAD_SCALE)
def test_gauss_couplings_refuse_a_scale_outside_the_window(mu, sigma):
    with pytest.raises(ValueError, match="coupling scale"):
        check_ctrl(Family(TC_FAMILY_GAUSS, mu=mu), sigma)


def test_gauss_scale_window_keeps_its_bounds_and_normal_units():
    for mu, sigma in (
        (0.0, TC_G_SCALE_MIN),
        (0.0, TC_G_SCALE_MAX),
        (-TC_G_SCALE_MIN, 0.0),
        (TC_G_SCALE_MAX, 0.0),
    ):
        check_ctrl(Family(TC_FAMILY_GAUSS, mu=mu), sigma)
    edges, family = _rrg(N_LANCZOS), Family(TC_FAMILY_GAUSS, mu=1.0)
    J, fields, _ = _unit(edges, N_LANCZOS, family, 0.5)
    draw = draw_disorder(family, N_LANCZOS, len(edges), ROOT, (N_LANCZOS, K), 0)
    assert np.array_equal(J, 1.0 + 0.5 * draw.values)
    lap, _ = signed_laplacian(edges, J, N_LANCZOS)
    assert fields["status"] == TC_STATUS_OK
    assert fields["lam0"] == pytest.approx(
        sla.eigvalsh(lap.toarray())[0], abs=EIG_TOL
    )


@pytest.mark.parametrize("scale", [TC_G_SCALE_MIN, TC_G_SCALE_MAX])
def test_gauss_units_do_not_depend_on_the_coupling_scale(scale):
    edges = _rrg(N_LANCZOS)
    unit, scaled = (
        _unit(edges, N_LANCZOS, Family(TC_FAMILY_GAUSS, mu=s), 0.5 * s)[1]
        for s in (1.0, scale)
    )
    assert scaled["status"] == unit["status"] == TC_STATUS_OK
    assert scaled["lam0"] == pytest.approx(scale * unit["lam0"], rel=EIG_TOL)
    for field in ("nunrel_eps", "pinf", "s1n", "tie"):
        assert scaled[field] == unit[field], field


def test_resolution_floor_is_c_eps_norm_in_units_of_the_mean_coupling():
    """SI Sec. IV.A floor for couplings +-1, carried to any coupling scale."""
    triangle = np.array([[0, 1], [1, 2], [0, 2]])
    shape = np.array([4.0, -4.0, 1.0])  # D = (5, 8, 5), mean |J| = 3
    norm_in_units = TC_GERSHGORIN_FACTOR * 8.0 / 3.0
    floor = TC_NUNREL_EPS_FACTOR * np.finfo(np.float64).eps * norm_in_units
    # a max |J| unit (4) would flag neither, no unit or a signed mean both
    near, far = 0.9 * floor, 2 * floor
    x = np.array([near, far, np.sqrt(1 - near**2 - far**2)])
    for scale in (1.0, TC_G_SCALE_MIN, TC_G_SCALE_MAX):
        J = scale * shape
        lap, dbar = signed_laplacian(triangle, J, len(x))
        diag = core._eigvec_diagnostics(lap, 0.0, x, dbar, J)
        assert diag["nunrel_eps"] == 1, scale


def test_record_rejects_unknown_fields():
    with pytest.raises(KeyError, match="unknown record field"):
        make_record({"pinff": 0.5})


# ---------------------------------------------------------------------------
# RRG topology: sampler, engine and draws
# ---------------------------------------------------------------------------
def test_validate_rrg_rejects_degenerate_ensembles():
    from kernels.RRG import validate_rrg

    for n, k in ((100, 2), (5, 3), (4, 4)):
        with pytest.raises(ValueError):
            validate_rrg(n, k)
    validate_rrg(100, 3)


def test_auto_sampler_is_pairing_for_small_k():
    from kernels.RRG import identity_engine

    assert resolve_sampler("auto", 3, "gt") == "pairing"
    assert resolve_sampler("auto", 5, "gt") == "native"
    with pytest.warns(UserWarning, match="not uniform"):
        assert resolve_sampler("native", 3, "nx") == "native"
    assert identity_engine("pairing", "gt") is None
    assert identity_engine("native", "gt") == "gt"
    with pytest.raises(ValueError, match="sampler must be one of"):
        resolve_sampler("bogus", 3, "gt")


def test_explicit_pairing_warns_above_auto_k_and_refuses_prohibitive_k():
    # 4 < k <= 6 is slow but feasible; k >= 7 exceeds the expected-draw cap,
    # and from k = 54 the estimate is inf (it would overflow a float)
    assert expected_pairing_draws(6) <= RRG_PAIRING_MAX_EXPECTED_DRAWS
    assert expected_pairing_draws(7) > RRG_PAIRING_MAX_EXPECTED_DRAWS
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert resolve_sampler("pairing", RRG_PAIRING_K_MAX, "gt") == "pairing"
    for k in (5, 6):
        with pytest.warns(UserWarning, match="draws per graph"):
            assert resolve_sampler("pairing", k, "gt") == "pairing"
    for k in (7, 8, 54):
        with pytest.raises(ValueError, match="use --sampler native"):
            resolve_sampler("pairing", k, "gt")


def test_require_engine_fails_loudly_instead_of_falling_back(monkeypatch):
    assert rrg.require_engine("nx") == "nx"
    with pytest.raises(ValueError, match="unavailable"):
        rrg.require_engine("bogus")

    def missing(graph_type, engine, fallback=True):
        assert not fallback  # a fallback would silently run NX
        raise ImportError("backend not importable (simulated)")

    monkeypatch.setattr(rrg, "get_implementation", missing)
    with pytest.raises(ValueError, match="'gt' is unavailable"):
        rrg.require_engine("gt")


@pytest.mark.parametrize(
    "engine,method",
    [("gt", "pairing"), ("gt", "native"), ("nx", "pairing"), ("nx", "native")],
)
def test_draw_topology_is_regular_connected_and_reproducible(engine, method):
    draw = dict(root=ROOT, engine=_engine(engine), method=method)
    first = draw_topology(N_SMALL, K, rep=0, **draw)
    again = draw_topology(N_SMALL, K, rep=0, **draw)
    other = draw_topology(N_SMALL, K, rep=1, **draw)
    degrees = np.bincount(first.edges.ravel(), minlength=N_SMALL)
    assert np.all(degrees == K) and len(first.edges) == N_SMALL * K // 2
    assert component_labels(first.edges, N_SMALL).max() == 0  # connected
    assert first.seed == engine_seed(ROOT, (N_SMALL, K), 0, first.attempt)
    assert first.fingerprint == again.fingerprint
    assert first.fingerprint != other.fingerprint
    if method == "pairing":  # the graph classes draw from pairing_stream
        expected = random_regular_edges(N_SMALL, K, pairing_stream(first.seed))
        assert np.array_equal(first.edges, expected)


def test_pairing_topology_does_not_depend_on_the_engine():
    """The run identity drops the engine for pairing, so both must agree."""
    engines = ("nx", _engine("gt"))
    draws = [
        draw_topology(N_SMALL, K, root=ROOT, rep=7, engine=e, method="pairing")
        for e in engines
    ]
    assert np.array_equal(draws[0].edges, draws[1].edges)
    assert draws[0].fingerprint == draws[1].fingerprint


@pytest.mark.parametrize("engine", ["nx", "gt"])
def test_a_real_disconnected_draw_is_redrawn(engine):
    seed0 = engine_seed(REDRAW_ROOT, (TINY_N, K), REDRAW_REP)
    attempt0 = random_regular_edges(TINY_N, K, pairing_stream(seed0))
    assert component_labels(attempt0, TINY_N).max() == 1  # two K4
    topo = draw_topology(
        TINY_N,
        K,
        root=REDRAW_ROOT,
        rep=REDRAW_REP,
        engine=_engine(engine),
        method="pairing",
    )
    assert topo.attempt == 1
    assert topo.seed == engine_seed(REDRAW_ROOT, (TINY_N, K), REDRAW_REP, 1)
    assert component_labels(topo.edges, TINY_N).max() == 0


class _StubGraph:
    """Just enough of a library graph for ``canonical_edges``."""

    def __init__(self, edges: np.ndarray, n_nodes: int = TINY_N) -> None:
        self.N, self.edges = n_nodes, np.asarray(edges)

    def get_edge_arrays(self):
        return self.edges[:, 0], self.edges[:, 1], np.ones(len(self.edges))


def _stub_generator(monkeypatch, *graphs: _StubGraph) -> list[int]:
    """Make the factory return ``graphs[a]`` at attempt ``a`` (the last one
    once they run out). Returns the list of seeds it is called with."""
    seeds: list[int] = []

    def factory(n_nodes, k, **kwargs):
        seeds.append(kwargs["seed"])
        return graphs[min(len(seeds), len(graphs)) - 1]

    monkeypatch.setattr(rrg, "kRegularGraph", factory)
    return seeds


def _draw_tiny(rep: int = 5):
    return draw_topology(
        TINY_N, K, root=ROOT, rep=rep, engine="nx", method="pairing"
    )


def test_disconnected_draw_is_redrawn_with_the_next_attempt_seed(monkeypatch):
    seeds = _stub_generator(monkeypatch, _StubGraph(TWO_K4), _StubGraph(CUBE))
    topo = _draw_tiny()
    assert topo.attempt == 1 and topo.seed == seeds[1]
    assert seeds == [engine_seed(ROOT, (TINY_N, K), 5, a) for a in (0, 1)]
    assert np.array_equal(topo.edges, CUBE)  # CUBE is already canonical


def test_always_disconnected_draws_raise_after_the_attempt_cap(monkeypatch):
    seeds = _stub_generator(monkeypatch, _StubGraph(TWO_K4))
    with pytest.raises(RuntimeError, match="no connected"):
        _draw_tiny()
    assert len(seeds) == RRG_MAX_ATTEMPTS


@pytest.mark.parametrize(
    "graph,match",
    [
        (_StubGraph(CUBE, n_nodes=TINY_N + 1), "9 nodes"),
        (_StubGraph(CUBE[1:]), "M=11"),
        (_StubGraph(REWIRED), "non-regular"),
    ],
    ids=["node-count", "edge-count", "degrees"],
)
def test_bad_generator_output_fails_loudly(monkeypatch, graph, match):
    seeds = _stub_generator(monkeypatch, graph)
    with pytest.raises(RuntimeError, match=match):
        _draw_tiny()
    assert len(seeds) == 1  # not retried: a broken generator is not bad luck

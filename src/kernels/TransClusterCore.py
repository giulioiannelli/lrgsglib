"""Graph-agnostic core of the TransCluster v2 programs.

A *unit* is one disorder realization (``rep``) at one control value
(``ictrl``). Each unit goes through these steps:

1. canonical edges ``E`` to nested couplings ``J``;
2. the signed Laplacian ``L = D - A`` with ``D_ii = sum_j |A_ij|``
   (iannelli2025topological, Eq. 1);
3. its ground state ``(lambda_0, x)``;
4. signs with the majority set positive;
5. same-sign clusters, where ``P_inf = n_inf / N`` (ibid., p. 4);
6. one record row plus the sparse cluster-size histogram ``P(S)``.

Everything works on plain arrays. ``canonical_edges`` is the only function
that touches a graph object, and it only reads it.

The disorder never goes through ``Disorder`` or ``flip_sel_edges``. Each
rep draws its variates once and thresholds them at every control value. For
one rep, a larger control value therefore always gives a superset of
defects (nested disorder), whatever the tiling or task order.
"""

from __future__ import annotations

import hashlib
import math
import operator
import time
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Iterator, Optional, Sequence, Union

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
from numpy.typing import NDArray
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import (
    ArpackError,
    ArpackNoConvergence,
    LinearOperator,
    eigsh,
)

from lrgsglib.config.const import EXACT_FLOAT_FMT, SG_REPR
from lrgsglib.config.progargs.defs.TransCluster import (
    DEFAULT_TC_EIG_TOL,
    DEFAULT_TC_ENSEMBLE,
    DEFAULT_TC_MU,
    DEFAULT_TC_NEV,
    TC_CTRL_RANGE,
    TC_DENSE_MAX_N,
    TC_EIG_MAX_MV,
    TC_EIG_NCV,
    TC_EIG_NORM_TOL,
    TC_EIG_RESID_FACTOR,
    TC_EIG_RETRY_NCV_FACTOR,
    TC_ENSEMBLE_CANONICAL,
    TC_ENSEMBLES,
    TC_FAMILIES,
    TC_FAMILY_GAUSS,
    TC_FAMILY_SIGNS,
    TC_FINGERPRINT_BYTES,
    TC_G_SCALE_MAX,
    TC_G_SCALE_MIN,
    TC_GAP_NEV,
    TC_GERSHGORIN_FACTOR,
    TC_KEY_WORD_LIMIT,
    TC_NO_COUNT,
    TC_NUNREL_EPS_FACTOR,
    TC_SEED_SCHEMA,
    TC_SEED_SPAN,
    TC_STATUS_ERROR,
    TC_STATUS_NOCONV,
    TC_STATUS_NONFINITE,
    TC_STATUS_OK,
    TC_STREAM_GAUSS,
    TC_STREAM_SIGNS,
    TC_STREAM_SUB,
    TC_STREAM_TOPO,
    TC_STREAM_V0,
    TC_STREAM_XNODES,
)
from lrgsglib.graphs._shared._regular_pairing import canonical_simple_edges
from lrgsglib.utils.basic.arithmetic import (
    bin_sign,
    flip_to_positive_majority_adapted,
)

__all__ = [
    "EMPTY_PS",
    "RECORD_DTYPE",
    "DisorderDraw",
    "Family",
    "GroundState",
    "SignClusters",
    "SolverConfig",
    "blank_record",
    "bond_energy",
    "canonical_edges",
    "check_ctrl",
    "cluster_observables",
    "component_labels",
    "couplings",
    "draw_disorder",
    "edge_fingerprint",
    "engine_seed",
    "exact_count",
    "ground_state",
    "is_gap_rep",
    "majority_signs",
    "make_record",
    "measure_unit",
    "orient_tie",
    "sign_clusters",
    "signed_laplacian",
    "solve_attempts",
    "stream_rng",
    "stream_seed",
    "sweep_rep",
]

#: Sparse P(S) of one unit: (cluster sizes, counts), giant excluded.
PS = tuple[NDArray[np.int64], NDArray[np.int64]]
EMPTY_PS: PS = (np.zeros(0, np.int64), np.zeros(0, np.int64))

# One row per unit. NaN (floats) or -1 (ints) marks what was not computed:
# gap fields outside gap reps, n_xx/gc_nonx outside X, min_absJ outside G,
# every spectral/order field of a row whose status is not OK.
RECORD_DTYPE = np.dtype(
    [
        # identity
        ("rep", np.int64),
        ("ictrl", np.int64),
        ("ctrl", np.float64),
        ("cnt", np.int64),  # prescribed count m or n_x; TC_NO_COUNT if none
        ("attempt", np.int64),  # topology redraws before this graph
        ("topo_fp", np.uint64),  # edge_fingerprint of the topology
        # disorder
        ("n_neg", np.int64),  # negative couplings
        ("p_eff", np.float64),  # n_neg / M
        ("n_xx", np.float64),  # X: edges with both endpoints X
        ("gc_nonx", np.float64),  # X: giant of the non-X subgraph / N
        ("min_absJ", np.float64),  # G: min |J_e|
        # spectrum
        ("status", np.int64),
        # eigenpairs solved (TC_GAP_NEV on gap reps); the requested count
        # on a row whose every attempt failed
        ("nev", np.int64),
        ("lam0", np.float64),
        ("resid", np.float64),  # ||L x - lambda_0 x||
        ("nunrel_eps", np.float64),  # #|x_i| < c eps ||L|| / mean|J|
        ("min_absx", np.float64),
        ("ipr_n", np.float64),  # N sum_i x_i^4
        ("n_zero", np.float64),
        ("n_mv", np.int64),  # Lanczos matvecs; -1 on the dense path
        ("t_eig", np.float64),  # seconds
        ("lam1", np.float64),
        ("gap", np.float64),
        ("dk", np.float64),  # Davis-Kahan bound resid / gap
        ("nunrel_dk", np.float64),  # #|x_i| < dk
        # order parameter (majority positive; clusters of same-sign edges).
        # ncl_p counts the giant too, so the finite clusters are ncl_p - 1.
        ("tie", np.int64),  # #(x < 0) == #(x > 0)
        ("pinf", np.float64),  # s1p / N
        ("n_pos", np.float64),
        ("s1p", np.float64),
        ("s2p", np.float64),
        ("ncl_p", np.float64),  # positive clusters, giant included
        ("sum_s_p", np.float64),  # sum S over positive clusters, giant out
        ("sum_s2_p", np.float64),  # sum S^2 over positive clusters, giant out
        ("s1n", np.float64),
        ("ncl_n", np.float64),  # negative clusters
        ("e0", np.float64),  # -sum_e J_e s_u s_v / M for sign(x_0)
        ("e1", np.float64),  # the same for sign(x_1), gap reps only
    ]
)


# ---------------------------------------------------------------------------
# Seed tree
# ---------------------------------------------------------------------------
def stream_seed(
    root: int,
    graph_key: Sequence[int],
    rep: int,
    stream: int,
    sub: int = TC_STREAM_SUB,
) -> np.random.SeedSequence:
    """SeedSequence of one random stream of ``rep``.

    The stream never depends on the control value, the tile or the time.
    """
    return _keyed_seed(root, *graph_key, _check_rep(rep), stream, sub)


def stream_rng(
    root: int,
    graph_key: Sequence[int],
    rep: int,
    stream: int,
    sub: int = TC_STREAM_SUB,
) -> np.random.Generator:
    """Generator of one random stream (see ``stream_seed``)."""
    return np.random.default_rng(stream_seed(root, graph_key, rep, stream, sub))


def engine_seed(
    root: int, graph_key: Sequence[int], rep: int, attempt: int = 0
) -> int:
    """Integer seed in ``[1, TC_SEED_SPAN]`` for a graph generator.

    Attempt 0 uses ``1 + (h(root, graph_key) + rep) mod span``. It is
    injective on ``0 <= rep < span``, so two reps never share a topology
    through a seed collision. Redraws (``attempt > 0``) hash the attempt in.
    """
    rep = _check_rep(rep)
    if attempt == 0:
        base = _keyed_seed(root, *graph_key, TC_STREAM_TOPO)
        return 1 + (_first_word(base) + rep) % TC_SEED_SPAN
    word = _first_word(
        stream_seed(root, graph_key, rep, TC_STREAM_TOPO, attempt)
    )
    return 1 + word % TC_SEED_SPAN


def _keyed_seed(root: int, *words: int) -> np.random.SeedSequence:
    """SeedSequence of ``root`` with spawn key ``(TC_SEED_SCHEMA, *words)``.

    Every word must be an integer (a float raises instead of being
    truncated) in ``[0, 2**32)``.
    """
    key = (TC_SEED_SCHEMA, *map(operator.index, words))
    if not all(0 <= w < TC_KEY_WORD_LIMIT for w in key):
        raise ValueError(f"spawn-key words must lie in [0, 2**32): {key}")
    return np.random.SeedSequence(operator.index(root), spawn_key=key)


def _check_rep(rep: int) -> int:
    """``rep`` as an int in ``[0, TC_SEED_SPAN)``, where ``engine_seed`` is
    injective."""
    rep = operator.index(rep)
    if not 0 <= rep < TC_SEED_SPAN:
        raise ValueError(f"rep must lie in [0, {TC_SEED_SPAN}), got {rep}")
    return rep


def _first_word(seq: np.random.SeedSequence) -> int:
    """First 32-bit word of a SeedSequence state."""
    return int(seq.generate_state(1, dtype=np.uint32)[0])


# ---------------------------------------------------------------------------
# Topology (the only place that reads a graph object)
# ---------------------------------------------------------------------------
def canonical_edges(graph: Any) -> tuple[int, NDArray[np.int64]]:
    """``(N, E)`` of a library graph, with ``E`` sorted and ``u < v`` per row.

    Reads only the endpoints. It uses the bulk GT accessor when present and
    the NX edge view of ``gr[SG_REPR]`` otherwise. Raises ``ValueError``
    on a non-simple graph.
    """
    n_nodes = int(graph.N)
    if hasattr(graph, "get_edge_arrays"):
        src, dst, _ = graph.get_edge_arrays()
        pairs = np.column_stack((src, dst))
    else:
        pairs = np.array(list(graph.gr[SG_REPR].edges()), dtype=np.int64)
    edges = canonical_simple_edges(pairs, n_nodes)
    if edges is None:
        raise ValueError("graph is not simple (self-loop or repeated edge)")
    return n_nodes, edges


def edge_fingerprint(edges: NDArray[np.int64]) -> int:
    """64-bit blake2b fingerprint of a canonical edge array."""
    data = np.ascontiguousarray(edges, dtype="<i8").tobytes()
    digest = hashlib.blake2b(data, digest_size=TC_FINGERPRINT_BYTES).digest()
    return int.from_bytes(digest, "little")


def component_labels(
    edges: NDArray[np.int64], n_nodes: int
) -> NDArray[np.int32]:
    """Connected-component label of every node of the graph ``edges``."""
    ones = np.ones(len(edges), dtype=np.float64)
    adj = sp.csr_matrix(
        (ones, (edges[:, 0], edges[:, 1])), shape=(n_nodes, n_nodes)
    )
    return connected_components(adj, directed=False)[1]


def _giant_fraction(
    node_mask: NDArray[np.bool_], edges: NDArray[np.int64]
) -> float:
    """Largest component of the subgraph induced by ``node_mask``, over N."""
    if not node_mask.any():
        return 0.0
    keep = node_mask[edges[:, 0]] & node_mask[edges[:, 1]]
    labels = component_labels(edges[keep], node_mask.size)
    return float(np.bincount(labels[node_mask]).max() / node_mask.size)


# ---------------------------------------------------------------------------
# Disorder: nested per-rep variates -> couplings
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Family:
    """Disorder family: ``S`` edge signs, ``X`` node defects, ``G`` gaussian.

    ``ensemble`` applies to S/X only. ``mu`` is the G mean coupling.
    """

    kind: str
    ensemble: str = DEFAULT_TC_ENSEMBLE
    mu: float = DEFAULT_TC_MU

    def __post_init__(self) -> None:
        if self.kind not in TC_FAMILIES:
            raise ValueError(f"family must be one of {TC_FAMILIES}")
        if self.ensemble not in TC_ENSEMBLES:
            raise ValueError(f"ensemble must be one of {TC_ENSEMBLES}")
        if not math.isfinite(self.mu):
            raise ValueError(f"mu must be finite, got {self.mu}")


@dataclass(frozen=True)
class DisorderDraw:
    """The variates of one rep, thresholded at every control value.

    ``values`` is ``u_e`` (S), ``w_i`` (X) or ``z_e`` (G). ``order`` is the
    stable argsort of ``values`` (canonical S/X only).
    """

    values: NDArray[np.float64]
    order: Optional[NDArray[np.int64]]


def exact_count(value: float, total: int) -> int:
    """Round-half-up count ``floor(value * total + 1/2)``.

    It is exact on the value's ``EXACT_FLOAT_FMT`` decimal: 0.29 * 1500
    gives 435, where ``int()`` gives 434.
    """
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"a count fraction must lie in [0, 1], got {value}")
    exact = Fraction(format(value, EXACT_FLOAT_FMT)) * total
    return math.floor(exact + Fraction(1, 2))


def check_ctrl(family: Family, ctrl: float) -> None:
    """Raise ``ValueError`` unless ``ctrl`` is a control value of ``family``.

    It must be finite and lie in ``TC_CTRL_RANGE``. For G, sigma = 0 at
    mu = 0 is refused (every coupling would vanish), and so is a coupling
    scale ``max(|mu|, sigma)`` outside ``[TC_G_SCALE_MIN, TC_G_SCALE_MAX]``.
    """
    low, high = TC_CTRL_RANGE[family.kind]
    if not (math.isfinite(ctrl) and low <= ctrl <= high):
        raise ValueError(
            f"{family.kind} control values must lie in [{low}, {high}], "
            f"got {ctrl}"
        )
    if family.kind != TC_FAMILY_GAUSS:
        return
    if ctrl == 0 and family.mu == 0:
        raise ValueError("mu = 0 and sigma = 0 make every coupling vanish")
    scale = max(abs(family.mu), ctrl)
    if not TC_G_SCALE_MIN <= scale <= TC_G_SCALE_MAX:
        raise ValueError(
            f"the coupling scale max(|mu|, sigma) = {float(scale)!r} lies "
            f"outside [{TC_G_SCALE_MIN:g}, {TC_G_SCALE_MAX:g}]; "
            "rescale mu and sigma"
        )


def draw_disorder(
    family: Family,
    n_nodes: int,
    n_edges: int,
    root: int,
    graph_key: Sequence[int],
    rep: int,
) -> DisorderDraw:
    """Draw the per-rep variates of ``family`` from its own stream."""
    if family.kind == TC_FAMILY_GAUSS:
        rng = stream_rng(root, graph_key, rep, TC_STREAM_GAUSS)
        return DisorderDraw(rng.standard_normal(n_edges), None)
    if family.kind == TC_FAMILY_SIGNS:
        rng = stream_rng(root, graph_key, rep, TC_STREAM_SIGNS)
        values = rng.random(n_edges)
    else:
        rng = stream_rng(root, graph_key, rep, TC_STREAM_XNODES)
        values = rng.random(n_nodes)
    canonical = family.ensemble == TC_ENSEMBLE_CANONICAL
    order = np.argsort(values, kind="stable") if canonical else None
    return DisorderDraw(values, order)


def _defect_mask(
    draw: DisorderDraw, ctrl: float
) -> tuple[NDArray[np.bool_], int]:
    """Defect mask and prescribed count (``TC_NO_COUNT`` if Bernoulli)."""
    if draw.order is None:
        return draw.values < ctrl, TC_NO_COUNT
    count = exact_count(ctrl, draw.values.size)
    mask = np.zeros(draw.values.size, dtype=bool)
    mask[draw.order[:count]] = True
    return mask, count


def couplings(
    family: Family,
    draw: DisorderDraw,
    edges: NDArray[np.int64],
    ctrl: float,
) -> tuple[NDArray[np.float64], dict[str, Any]]:
    """Couplings ``J_e`` at one control value, plus the disorder fields.

    S: the defect edges get ``-1``. G: ``J = mu + sigma * z`` with
    ``sigma = ctrl``. X: every edge touching an X node is SET to ``-1``
    (iannelli2025topological, p. 10). Derived here, not in the paper: SET-X
    is gauge-equivalent to S on the X-X edges alone
    (``test_set_x_is_gauge_equivalent_to_s_on_xx_edges``), so it is
    balanced (``lambda_0 = 0``) whenever the negative X-X edges alone are,
    e.g. at ``n_xx = 0`` (cf. ibid., p. 4). Toggling the signs instead of
    setting them would be a pure gauge transformation.
    """
    check_ctrl(family, ctrl)
    if family.kind == TC_FAMILY_GAUSS:
        J = family.mu + ctrl * draw.values
        gauss = {"cnt": TC_NO_COUNT, "min_absJ": np.abs(J).min()}
        return J, {**gauss, **_negative_fields(J < 0)}
    mask, count = _defect_mask(draw, ctrl)
    extra: dict[str, Any] = {"cnt": count}
    if family.kind == TC_FAMILY_SIGNS:
        negative = mask
    else:
        at_u, at_v = mask[edges[:, 0]], mask[edges[:, 1]]
        negative = at_u | at_v
        extra["n_xx"] = np.count_nonzero(at_u & at_v)
        extra["gc_nonx"] = _giant_fraction(~mask, edges)
    J = np.where(negative, -1.0, 1.0)
    return J, {**extra, **_negative_fields(negative)}


def _negative_fields(negative: NDArray[np.bool_]) -> dict[str, Any]:
    """``n_neg`` and ``p_eff = n_neg / M`` of a negative-edge mask."""
    n_neg = int(np.count_nonzero(negative))
    return {"n_neg": n_neg, "p_eff": n_neg / negative.size}


# ---------------------------------------------------------------------------
# Spectral step
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SolverConfig:
    """Eigensolver knobs; all of them belong to the run identity.

    ``max_mv`` caps the matvecs of each solve attempt.
    """

    tol: float = DEFAULT_TC_EIG_TOL
    ncv: int = TC_EIG_NCV
    max_mv: int = TC_EIG_MAX_MV
    dense_max_n: int = TC_DENSE_MAX_N

    def __post_init__(self) -> None:
        if not (math.isfinite(self.tol) and self.tol >= 0):
            raise ValueError(f"tol must be finite and >= 0, got {self.tol}")
        object.__setattr__(self, "tol", self.tol + 0.0)  # -0.0 becomes 0.0
        if self.ncv <= TC_GAP_NEV:
            raise ValueError(f"ncv must exceed {TC_GAP_NEV}, got {self.ncv}")
        if self.max_mv < 1:
            raise ValueError(f"max_mv must be >= 1, got {self.max_mv}")
        if self.dense_max_n < 0:
            raise ValueError(
                f"dense_max_n must be >= 0, got {self.dense_max_n}"
            )


@dataclass(frozen=True)
class GroundState:
    """Lowest eigenpairs (ascending) with solver status and cost."""

    vals: NDArray[np.float64]
    vecs: NDArray[np.float64]
    status: int
    n_mv: int
    t_eig: float


class _MatvecCapReached(Exception):
    """A solve attempt used up its matvecs."""


class _CountingOperator(LinearOperator):
    """``L @ x`` that counts its calls and stops an attempt at ``limit``."""

    def __init__(self, matrix: sp.csr_matrix) -> None:
        super().__init__(dtype=matrix.dtype, shape=matrix.shape)
        self.matrix = matrix
        self.n_mv = 0
        self.limit = math.inf

    def _matvec(self, x: NDArray[np.float64]) -> NDArray[np.float64]:
        if self.n_mv >= self.limit:
            raise _MatvecCapReached
        self.n_mv += 1
        return self.matrix @ x


def signed_laplacian(
    edges: NDArray[np.int64], J: NDArray[np.float64], n_nodes: int
) -> tuple[sp.csr_matrix, NDArray[np.float64]]:
    """``L = D - A`` as a float64 symmetric CSR, together with ``diag(D)``."""
    u, v = edges[:, 0], edges[:, 1]
    absJ = np.abs(J)
    dbar = np.bincount(u, absJ, n_nodes) + np.bincount(v, absJ, n_nodes)
    nodes = np.arange(n_nodes)
    rows = np.concatenate((u, v, nodes))
    cols = np.concatenate((v, u, nodes))
    data = np.concatenate((-J, -J, dbar)).astype(np.float64)
    lap = sp.csr_matrix((data, (rows, cols)), shape=(n_nodes, n_nodes))
    return lap, dbar


def ground_state(
    lap: sp.csr_matrix,
    v0: NDArray[np.float64],
    nev: int = DEFAULT_TC_NEV,
    solver: SolverConfig = SolverConfig(),
) -> GroundState:
    """Lowest ``nev`` eigenpairs of ``lap``, or fewer (see ``_lanczos``).

    Uses dense ``eigh`` up to ``solver.dense_max_n`` and Lanczos above.
    Solver failures become a status code, never an exception; the pairs of
    a failed unit are NaN.
    """
    t0 = time.perf_counter()
    if lap.shape[0] <= solver.dense_max_n:
        vals, vecs, status = _dense(lap, nev)
        n_mv = -1
    else:
        vals, vecs, n_mv, status = _lanczos(lap, v0, nev, solver)
    return GroundState(vals, vecs, status, n_mv, time.perf_counter() - t0)


def _usable(vals: NDArray[np.float64], vecs: NDArray[np.float64]) -> bool:
    """Finite eigenpairs whose vectors have unit norm."""
    if not (np.isfinite(vals).all() and np.isfinite(vecs).all()):
        return False
    norms = np.linalg.norm(vecs, axis=0)
    return bool(np.all(np.abs(norms - 1) <= TC_EIG_NORM_TOL))


def _converged(
    lap: sp.csr_matrix,
    vals: NDArray[np.float64],
    vecs: NDArray[np.float64],
    tol: float,
) -> bool:
    """Whether every residual ``||L x - lambda x||`` is at most
    ``TC_EIG_RESID_FACTOR * max(eps, tol) * ||L||``.

    ARPACK's own test is relative to ``|lambda|`` only down to
    ``eps**(2/3)`` and absolute below, so on tiny couplings it accepts
    pairs that are wrong. This test is relative to the Gershgorin bound
    ``||L|| <= 2 max_i D_ii`` whatever the scale of the couplings.
    """
    eps = np.finfo(np.float64).eps
    norm_bound = TC_GERSHGORIN_FACTOR * lap.diagonal().max()
    bound = TC_EIG_RESID_FACTOR * max(eps, tol) * norm_bound
    resid = np.linalg.norm(lap @ vecs - vecs * vals, axis=0)
    return bool(np.all(resid <= bound))


def _no_eigenpairs(
    n: int, nev: int
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    return np.full(nev, np.nan), np.full((n, nev), np.nan)


def _dense(
    lap: sp.csr_matrix, nev: int
) -> tuple[NDArray[np.float64], NDArray[np.float64], int]:
    """Dense ``eigh`` of the ``nev`` lowest eigenpairs."""
    try:
        vals, vecs = sla.eigh(lap.toarray(), subset_by_index=[0, nev - 1])
    except np.linalg.LinAlgError:
        return *_no_eigenpairs(lap.shape[0], nev), TC_STATUS_ERROR
    if not _usable(vals, vecs):
        return *_no_eigenpairs(lap.shape[0], nev), TC_STATUS_NONFINITE
    return vals, vecs, TC_STATUS_OK


def solve_attempts(nev: int, ncv: int) -> list[tuple[int, int]]:
    """``(nev, ncv)`` of each Lanczos attempt of a unit, in order.

    Every unit ends with the nev = 1 solve and its larger-basis retry (the
    cure for ARPACK's "no shifts could be applied"), so a failed gap solve
    still gives lambda_0 and P_inf, like any other unit.
    """
    gap = [(nev, ncv)] if nev > 1 else []
    return gap + [(1, ncv), (1, ncv * TC_EIG_RETRY_NCV_FACTOR)]


def _lanczos(
    lap: sp.csr_matrix,
    v0: NDArray[np.float64],
    nev: int,
    solver: SolverConfig,
) -> tuple[NDArray[np.float64], NDArray[np.float64], int, int]:
    """ARPACK ``which='SA'`` over ``solve_attempts``; the first to return
    usable eigenpairs (``_usable``) whose residuals pass ``_converged``
    wins, and ``n_mv`` counts them all.

    Each attempt starts from ``v0`` (a failed nev = 1 attempt returns no
    converged vector to restart from) and stops after ``solver.max_mv``
    matvecs. The status of a unit where all fail is that of the last one.
    """
    op = _CountingOperator(lap)
    n = lap.shape[0]
    status = TC_STATUS_NOCONV
    for k, ncv in solve_attempts(nev, solver.ncv):
        op.limit = op.n_mv + solver.max_mv
        try:
            vals, vecs = eigsh(
                op,
                k=k,
                which="SA",
                tol=solver.tol,
                v0=v0,
                ncv=min(ncv, n),
                maxiter=solver.max_mv,  # restarts: never binds before the cap
            )
        except (ArpackNoConvergence, _MatvecCapReached):
            status = TC_STATUS_NOCONV
            continue
        except ArpackError:
            status = TC_STATUS_ERROR
            continue
        if not _usable(vals, vecs):
            status = TC_STATUS_NONFINITE
            continue
        if not _converged(lap, vals, vecs, solver.tol):
            status = TC_STATUS_NOCONV
            continue
        order = np.argsort(vals)
        return vals[order], vecs[:, order], op.n_mv, TC_STATUS_OK
    return *_no_eigenpairs(n, nev), op.n_mv, status


def _eigvec_diagnostics(
    lap: sp.csr_matrix,
    lam: float,
    x: NDArray[np.float64],
    dbar: NDArray[np.float64],
    J: NDArray[np.float64],
) -> dict[str, Any]:
    """Residual, resolution-limit count and localization of an eigenvector.

    Components below ``c * eps * ||L||`` are unreliable
    (iannelli2025topological, SI Sec. IV.A, couplings +-1). ``||L||`` is
    bounded by Gershgorin, ``||L|| <= 2 max_i D_ii``, and taken in units of
    the mean ``|J|``: ``x`` has unit norm, so its floor must not grow with
    the scale of G couplings. For S and X the unit is 1.
    """
    x = x / np.linalg.norm(x)
    absx = np.abs(x)
    norm_bound = TC_GERSHGORIN_FACTOR * dbar.max() / np.abs(J).mean()
    floor = TC_NUNREL_EPS_FACTOR * np.finfo(np.float64).eps * norm_bound
    return {
        "resid": np.linalg.norm(lap @ x - lam * x),
        "nunrel_eps": np.count_nonzero(absx < floor),
        "min_absx": absx.min(),
        "ipr_n": x.size * np.sum(x**4),
        "n_zero": np.count_nonzero(x == 0),
    }


# ---------------------------------------------------------------------------
# Order parameter
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class SignClusters:
    """Sizes and signs of the same-sign clusters."""

    sizes: NDArray[np.int64]
    signs: NDArray[np.int8]

    def largest(self, sign: int) -> int:
        """Size of the largest cluster of ``sign`` (0 if none)."""
        sizes = self.sizes[self.signs == sign]
        return int(sizes.max()) if sizes.size else 0

    def ranked(self, sign: int) -> list[int]:
        """Sizes of the clusters of ``sign``, largest first."""
        return sorted(self.sizes[self.signs == sign].tolist(), reverse=True)


def majority_signs(x: NDArray[np.float64]) -> tuple[NDArray[np.int8], bool]:
    """Signs of ``x`` with the majority set positive, plus the tie flag.

    Flips when ``#(x < 0) > #(x > 0)`` and maps zeros to ``+1``, as
    ``compute_k_eigvV`` does (``flip_to_positive_majority_adapted``, then
    ``bin_sign``). ``x -> -x`` swaps the two counts, so a non-tie row gets
    the same signs whatever sign the eigensolver returns, zeros included.
    The tie flag is ``#(x < 0) == #(x > 0)``.
    """
    signs = bin_sign(flip_to_positive_majority_adapted(x)).astype(np.int8)
    return signs, np.count_nonzero(x < 0) == np.count_nonzero(x > 0)


def sign_clusters(
    edges: NDArray[np.int64], signs: NDArray[np.int8], n_nodes: int
) -> SignClusters:
    """Components of the subgraph of edges with ``s_u == s_v`` (any J)."""
    same = signs[edges[:, 0]] == signs[edges[:, 1]]
    labels = component_labels(edges[same], n_nodes)
    sizes = np.bincount(labels).astype(np.int64)
    cluster_signs = np.zeros(sizes.size, dtype=np.int8)
    cluster_signs[labels] = signs
    return SignClusters(sizes, cluster_signs)


def orient_tie(
    signs: NDArray[np.int8], clusters: SignClusters, tie: bool
) -> tuple[NDArray[np.int8], SignClusters]:
    """At an exact tie, make the sign with the larger clusters positive.

    The cluster sizes of each sign, largest first, are compared
    lexicographically: the larger giant wins, then the larger second
    cluster, and so on; equal lists put node 0 positive. The rule reads
    only the partition, which ``x`` and ``-x`` share, so the orientation
    does not depend on the sign the eigensolver returns. The exception is
    an ``x`` with exact zeros (``n_zero > 0``, tiny N only): the zeros join
    the positive side, so the partition itself follows that sign.

    A tie row then has ``P_inf = max(P+, P-)``. Any other row keeps the
    giant of the majority sign, which need not be the largest cluster
    overall; this is the library ``compute_pinf`` convention.
    ``max(s1p, s1n) / N`` can be rebuilt from the record.
    """
    if not tie:
        return signs, clusters
    neg, pos = clusters.ranked(-1), clusters.ranked(1)
    if neg > pos or (neg == pos and signs[0] < 0):
        return -signs, SignClusters(clusters.sizes, -clusters.signs)
    return signs, clusters


def cluster_observables(
    clusters: SignClusters, n_nodes: int
) -> tuple[dict[str, Any], PS]:
    """Order-parameter fields and the positive-cluster ``P(S)``.

    The largest positive cluster (the giant) is excluded from ``P(S)``.
    """
    pos = np.sort(clusters.sizes[clusters.signs > 0])[::-1]
    neg = clusters.sizes[clusters.signs < 0]
    rest = pos[1:]
    s1p = int(pos[0]) if pos.size else 0
    fields = {
        "pinf": s1p / n_nodes,
        "n_pos": pos.sum(),
        "s1p": s1p,
        "s2p": int(pos[1]) if pos.size > 1 else 0,
        "ncl_p": pos.size,
        "sum_s_p": rest.sum(),
        "sum_s2_p": np.sum(rest.astype(np.float64) ** 2),
        "s1n": clusters.largest(-1),
        "ncl_n": neg.size,
    }
    sizes, counts = np.unique(rest, return_counts=True)
    return fields, (sizes.astype(np.int64), counts.astype(np.int64))


def bond_energy(
    edges: NDArray[np.int64],
    J: NDArray[np.float64],
    signs: NDArray[np.int8],
) -> float:
    """Ising energy per bond, ``-sum_e J_e s_u s_v / M`` (ibid., Eq. 3)."""
    s = signs.astype(np.float64)
    return float(-np.dot(J, s[edges[:, 0]] * s[edges[:, 1]]) / len(J))


# ---------------------------------------------------------------------------
# One unit, one rep
# ---------------------------------------------------------------------------
def blank_record() -> np.ndarray:
    """A record row with every field marked as not computed."""
    rec = np.zeros((), dtype=RECORD_DTYPE)
    for name in RECORD_DTYPE.names:
        kind = RECORD_DTYPE[name].kind
        rec[name] = np.nan if kind == "f" else (-1 if kind == "i" else 0)
    return rec


def make_record(fields: dict[str, Any]) -> np.ndarray:
    """A record row filled from ``fields``; an unknown name raises."""
    rec = blank_record()
    for name, value in fields.items():
        if name not in RECORD_DTYPE.names:
            raise KeyError(f"unknown record field {name!r}")
        rec[name] = value
    return rec


def is_gap_rep(rep: int, gap_every: int) -> bool:
    """Whether ``rep`` also solves ``lambda_1`` (0 disables gap reps)."""
    return gap_every > 0 and rep % gap_every == 0


def measure_unit(
    edges: NDArray[np.int64],
    J: NDArray[np.float64],
    n_nodes: int,
    v0: NDArray[np.float64],
    nev: int = DEFAULT_TC_NEV,
    solver: SolverConfig = SolverConfig(),
) -> tuple[dict[str, Any], PS]:
    """Spectral and order-parameter fields of one unit (couplings given).

    The ``nev`` field counts the eigenpairs solved: 1 on a gap rep whose
    nev = 2 solve failed, whose gap fields are then NaN.
    """
    _check_couplings(J)
    lap, dbar = signed_laplacian(edges, J, n_nodes)
    if not np.isfinite(dbar).all():
        raise ValueError("the degree sums overflow: rescale the couplings")
    gs = ground_state(lap, v0, nev, solver)
    fields: dict[str, Any] = {
        "status": gs.status,
        "nev": gs.vals.size,
        "n_mv": gs.n_mv,
        "t_eig": gs.t_eig,
    }
    if gs.status != TC_STATUS_OK:
        return fields, EMPTY_PS
    lam0, x0 = gs.vals[0], gs.vecs[:, 0]
    fields.update(lam0=lam0, **_eigvec_diagnostics(lap, lam0, x0, dbar, J))
    signs, tie = majority_signs(x0)
    signs, clusters = orient_tie(
        signs, sign_clusters(edges, signs, n_nodes), tie
    )
    observables, ps = cluster_observables(clusters, n_nodes)
    fields.update(observables, tie=int(tie), e0=bond_energy(edges, J, signs))
    if gs.vals.size > 1:
        fields.update(_gap_fields(gs, fields["resid"], edges, J))
    return fields, ps


def _check_couplings(J: NDArray[np.float64]) -> None:
    """Refuse couplings whose ground state is undefined."""
    if not np.isfinite(J).all():
        raise ValueError("the couplings must be finite")
    if not np.any(J):
        raise ValueError("all couplings vanish: every vector is a ground state")


def _gap_fields(
    gs: GroundState,
    resid: float,
    edges: NDArray[np.int64],
    J: NDArray[np.float64],
) -> dict[str, Any]:
    """lambda_1, the gap, the Davis-Kahan bound and the energy of sign(x_1)."""
    gap = gs.vals[1] - gs.vals[0]
    dk = resid / gap if gap > 0 else np.inf
    signs1, _ = majority_signs(gs.vecs[:, 1])
    return {
        "lam1": gs.vals[1],
        "gap": gap,
        "dk": dk,
        "nunrel_dk": np.count_nonzero(np.abs(gs.vecs[:, 0]) < dk),
        "e1": bond_energy(edges, J, signs1),
    }


def sweep_rep(
    edges: NDArray[np.int64],
    n_nodes: int,
    family: Family,
    ctrl_values: Union[Sequence[float], NDArray[np.float64]],
    ictrls: Sequence[int],
    *,
    root: int,
    graph_key: Sequence[int],
    rep: int,
    gap_every: int,
    solver: SolverConfig = SolverConfig(),
    base_fields: Optional[dict[str, Any]] = None,
) -> Iterator[tuple[np.ndarray, PS]]:
    """Yield ``(record, P(S))`` for rep ``rep`` at ``ctrl_values[ictrls]``.

    The disorder variates and the start vector are drawn once, so the
    output of a unit does not depend on which other units run with it.
    """
    draw = draw_disorder(family, n_nodes, len(edges), root, graph_key, rep)
    v0 = stream_rng(root, graph_key, rep, TC_STREAM_V0).standard_normal(n_nodes)
    nev = TC_GAP_NEV if is_gap_rep(rep, gap_every) else DEFAULT_TC_NEV
    for ictrl in ictrls:
        ctrl = float(ctrl_values[ictrl])
        J, disorder = couplings(family, draw, edges, ctrl)
        fields, ps = measure_unit(edges, J, n_nodes, v0, nev, solver)
        identity = {"rep": rep, "ictrl": ictrl, "ctrl": ctrl}
        row = {**(base_fields or {}), **identity, **disorder, **fields}
        yield make_record(row), ps

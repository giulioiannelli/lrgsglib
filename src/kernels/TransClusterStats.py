"""Statistics of TransCluster v2 runs: complete tiles -> ``agg.npz``.

``collect(rundir, n_nodes)`` reads every complete tile of a run directory,
keeps one record per ``(rep, ctrl)`` (a control value shared by two grids
is computed twice), checks that each rep kept one topology, and writes
``<rundir>/agg.npz``:

- ``summary``: one ``SUMMARY_DTYPE`` row per control value;
- ``ps_ptr``, ``ps_size``, ``ps_count``: the P(S) of the finite positive
  clusters (giant excluded) of each row, summed over reps, in CSR form;
- ``meta``: unit counts and the chi-peak location with its paired-bootstrap
  interval and window-edge flags (``META_DTYPE``);
- the run identity, ``cfg.json``, the provenance of the collect and the
  distinct code states (``tile_git``) of the tiles.

Errors about a tile name its file.

The order parameter is P_inf and chi = sqrt(N) sd(P_inf)
(iannelli2025topological.pdf, p.4); ``chi`` uses ddof 1, ``chi0`` the
legacy ddof 0. Records whose eigensolve failed are counted, not used.
``load_agg`` reads the file back; ``legacy_rows`` gives the six columns of
the legacy ordParam files.
"""

from __future__ import annotations

import hashlib
import json
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, NamedTuple, Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from lrgsglib.config.progargs.defs.TransCluster import (
    DEFAULT_TC_BOOT,
    DEFAULT_TC_BOOT_SEED,
    TC_AGG_NAME,
    TC_AGG_SCHEMA,
    TC_BOOT_MIN_REPS,
    TC_BOOT_QUANTILES,
    TC_KURT_MIN_N,
    TC_LEGACY_FMT,
    TC_LEGACY_TXT_NAME,
    TC_PART_EXT,
    TC_PEAK_MIN_CTRLS,
    TC_STATUS_OK,
    TC_TILE_EXT,
    TC_TILES_DIR,
)
from lrgsglib.statsys._io import atomic_write_bytes
from lrgsglib.statsys._naming import CFG_SIDECAR_NAME

from .TransClusterCore import PS, RECORD_DTYPE
from .TransClusterTiles import (
    TileBuffer,
    concat_int64,
    load_tile,
    provenance,
    write_npz,
)

__all__ = [
    "META_DTYPE",
    "SUMMARY_DTYPE",
    "UNIT_FIELDS",
    "VOLATILE_FIELDS",
    "Deduped",
    "RunRecords",
    "bootstrap_peak",
    "check_topologies",
    "chi_peak",
    "collect",
    "complete_tile_paths",
    "cost_table",
    "drop_duplicate_units",
    "jackknife_chi_se",
    "kurtosis_chi_se",
    "legacy_rows",
    "load_agg",
    "load_records",
    "paired_matrix",
    "pearson_kurtosis",
    "pool_ps",
    "summarize",
    "summary_row",
]

#: Fields fixed by (seed, rep, control value): duplicate units must agree.
UNIT_FIELDS = ("cnt", "attempt", "topo_fp", "n_neg")
#: Fields that may differ between two computations of one unit.
VOLATILE_FIELDS = ("ictrl", "t_eig")
_JACKKNIFE_MIN = 3  # a leave-one-out sd (ddof 1) needs two values left
#: bootstrap_peak without enough reps or control values.
_NO_PEAK = {
    "peak": np.nan,
    "peak_lo": np.nan,
    "peak_hi": np.nan,
    "peak_sd": np.nan,
    "peak_at_edge": -1,
    "peak_edge_frac": np.nan,
}

SUMMARY_DTYPE = np.dtype(
    [
        ("ctrl", "<f8"),  # control value: p, q or sigma
        ("n_rep", "<i8"),  # records used (eigensolve converged)
        ("n_bad", "<i8"),  # records whose eigensolve failed
        ("n_tie", "<i8"),  # records with as many + as - signs
        ("pinf", "<f8"),  # <P_inf>
        ("pinf2", "<f8"),  # <P_inf^2>
        ("pinf_se", "<f8"),  # standard error of <P_inf>
        # Pearson kurtosis of P_inf (3 for a gaussian), plug-in m4 / m2^2:
        # biased low at small n_rep, most on heavy tails (kurtosis_chi_se)
        ("kurt", "<f8"),
        ("chi", "<f8"),  # sqrt(N) sd(P_inf), ddof 1
        ("chi0", "<f8"),  # sqrt(N) sd(P_inf), ddof 0 (legacy)
        ("chi_se_jk", "<f8"),  # jackknife standard error of chi
        # chi sqrt((max(kurt - 1, 0) + 2 / (n_rep - 1)) / (4 n_rep))
        ("chi_se_kurt", "<f8"),
        ("n_neg", "<f8"),  # <negative couplings>
        ("p_eff", "<f8"),  # <negative coupling fraction>
        ("lam0", "<f8"),  # <lambda_0>
        ("e0", "<f8"),  # <energy per edge of the ground-state signs>
        ("smean", "<f8"),  # sum s^2 / sum s over finite positive clusters
        ("n_gap", "<i8"),  # records that also solved lambda_1
        ("gap", "<f8"),  # <lambda_1 - lambda_0> over those records
        ("e1", "<f8"),  # <energy per edge of the lambda_1 signs>
        ("t_eig_med", "<f8"),  # median eigensolve time [s]
        ("t_eig_max", "<f8"),  # longest eigensolve [s]
        ("n_mv_med", "<f8"),  # median matrix-vector products (NaN: dense)
    ]
)

META_DTYPE = np.dtype(
    [
        ("n_nodes", "<i8"),
        ("n_units", "<i8"),  # records kept
        ("n_dupes", "<i8"),  # duplicate units dropped
        ("n_dupes_differ", "<i8"),  # of those, differing beyond timing
        ("n_tiles", "<i8"),  # complete tiles read
        ("n_redundant_topo", "<i8"),  # reps minus distinct topologies
        ("peak", "<f8"),  # chi-peak location (paired reps)
        # peak_lo/hi and peak_edge_frac NaN: no usable resample (--boot 0,
        # below TC_BOOT_MIN_REPS paired reps, no peak, or only flat ones)
        ("peak_lo", "<f8"),  # TC_BOOT_QUANTILES of the bootstrap
        ("peak_hi", "<f8"),
        ("peak_sd", "<f8"),  # bootstrap standard deviation
        ("peak_at_edge", "<i8"),  # 1: the chi maximum is a window end
        ("peak_edge_frac", "<f8"),  # bootstrap maxima at a window end
        ("peak_reps", "<i8"),  # reps measured at every point of the window
        ("window_lo", "<f8"),
        ("window_hi", "<f8"),
        ("n_boot", "<i8"),  # resamples requested, not the usable ones
        ("boot_seed", "<i8"),
    ]
)


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RunRecords:
    """The units of one run directory, once each, sorted by (rep, ctrl)."""

    records: NDArray  # RECORD_DTYPE
    ps: list[PS]
    source: NDArray  # tile file name of each record
    identity: str  # sha256 of cfg.json
    cfg: str  # cfg.json text
    n_tiles: int
    n_dupes: int
    n_dupes_differ: int
    git_states: list[Optional[dict[str, Any]]]  # distinct, of the tiles


class Deduped(NamedTuple):
    """The output of ``drop_duplicate_units``."""

    records: NDArray  # one per (rep, ctrl), sorted
    ps: list[PS]
    source: NDArray  # tile file of each record
    n_dupes: int  # duplicates dropped
    n_differ: int  # of those, differing beyond VOLATILE_FIELDS


def complete_tile_paths(rundir: Path) -> list[Path]:
    """Complete tile files of ``rundir`` (checkpoints excluded), sorted."""
    paths = (Path(rundir) / TC_TILES_DIR).glob(f"*{TC_TILE_EXT}")
    return sorted(p for p in paths if not p.name.endswith(TC_PART_EXT))


def load_records(rundir: Path) -> RunRecords:
    """Every unit of the complete tiles of ``rundir``, once each. Tiles
    must be well formed and match the run's ``cfg.json``; errors name the
    tile file."""
    rundir = Path(rundir)
    cfg = (rundir / CFG_SIDECAR_NAME).read_bytes()
    identity = hashlib.sha256(cfg).hexdigest()
    paths = complete_tile_paths(rundir)
    if not paths:
        raise FileNotFoundError(f"no complete tiles in {rundir}")
    tables, ps, source, states = [], [], [], set()
    for path in paths:
        data = load_tile(path, identity)
        tables.append(data["records"])
        ps.extend(TileBuffer.from_arrays(data).ps)
        source.append(np.full(data["records"].size, path.name, dtype=object))
        git = json.loads(str(data["provenance"])).get("git")
        states.add(json.dumps(git, sort_keys=True))
    kept = drop_duplicate_units(
        np.concatenate(tables), ps, np.concatenate(source)
    )
    return RunRecords(
        records=kept.records,
        ps=kept.ps,
        source=kept.source,
        identity=identity,
        cfg=cfg.decode(),
        n_tiles=len(paths),
        n_dupes=kept.n_dupes,
        n_dupes_differ=kept.n_differ,
        git_states=[json.loads(state) for state in sorted(states)],
    )


def drop_duplicate_units(
    records: NDArray, ps: Sequence[PS], source: Optional[NDArray] = None
) -> Deduped:
    """One record per ``(rep, ctrl)``, sorted by (rep, ctrl).

    Among duplicates a converged record beats a failed one, and a
    gap-solved record beats one without the gap, whatever the tile order.
    Duplicates must agree on ``UNIT_FIELDS``; ``source``, the tile of each
    record, names the files of a clash. ``n_differ`` counts the dropped
    copies that differ from the kept record beyond ``VOLATILE_FIELDS``.
    """
    if source is None:
        source = np.full(records.size, "?", dtype=object)
    failed = records["status"] != TC_STATUS_OK
    no_gap = np.isnan(records["gap"])
    order = np.lexsort((no_gap, failed, records["ctrl"], records["rep"]))
    records, source = records[order], source[order]
    ps = [ps[i] for i in order]
    rep, ctrl = records["rep"], records["ctrl"]
    dup = (rep[1:] == rep[:-1]) & (ctrl[1:] == ctrl[:-1])
    for name in UNIT_FIELDS:
        clash = dup & (records[name][1:] != records[name][:-1])
        if clash.any():
            i = int(np.flatnonzero(clash)[0])
            raise ValueError(
                f"rep {rep[i]} at {ctrl[i]}: tiles disagree on {name}"
                + _where(source, i)
            )
    keep = np.concatenate(([True], ~dup))
    drop = ~keep
    kept = np.flatnonzero(keep)[np.cumsum(keep) - 1]  # the kept one of its unit
    same = _same_records(records[drop], records[kept[drop]])
    return Deduped(
        records[keep],
        [p for p, k in zip(ps, keep) if k],
        source[keep],
        int(dup.sum()),
        int(np.count_nonzero(~same)),
    )


def _same_records(a: NDArray, b: NDArray) -> NDArray[np.bool_]:
    """Row-wise ``a == b`` outside ``VOLATILE_FIELDS``; NaN equals NaN."""
    same = np.ones(a.size, dtype=bool)
    for name in RECORD_DTYPE.names:
        if name in VOLATILE_FIELDS:
            continue
        x, y = a[name], b[name]
        equal = x == y
        if x.dtype.kind == "f":
            equal |= np.isnan(x) & np.isnan(y)
        same &= equal
    return same


def _where(source: Optional[NDArray], i: int) -> str:
    """`` (tile a, tile b)`` of records ``i`` and ``i + 1``, if known."""
    if source is None:
        return ""
    return f" ({source[i]}, {source[i + 1]})"


def check_topologies(records: NDArray, source: Optional[NDArray] = None) -> int:
    """Raise unless each rep has one topology in all its records. Return
    the redundant reps, reps minus distinct topologies: 0 unless two reps
    drew the same graph, which only a small N makes likely."""
    rep, fp, attempt = records["rep"], records["topo_fp"], records["attempt"]
    same = rep[1:] == rep[:-1]  # records are sorted by rep
    clash = same & ((fp[1:] != fp[:-1]) | (attempt[1:] != attempt[:-1]))
    if clash.any():
        i = int(np.flatnonzero(clash)[0])
        raise ValueError(
            f"rep {rep[i]} has records of two topologies" + _where(source, i)
        )
    firsts = fp[np.concatenate(([True], ~same))]
    return int(firsts.size - np.unique(firsts).size)


# ---------------------------------------------------------------------------
# Per control value
# ---------------------------------------------------------------------------
def summarize(records: NDArray, n_nodes: int) -> NDArray:
    """One ``SUMMARY_DTYPE`` row per control value."""
    ctrls = np.unique(records["ctrl"])
    out = np.zeros(ctrls.size, dtype=SUMMARY_DTYPE)
    for i, ctrl in enumerate(ctrls):
        row = summary_row(records[records["ctrl"] == ctrl], n_nodes)
        for name, value in row.items():
            out[name][i] = value
    return out


def summary_row(recs: NDArray, n_nodes: int) -> dict[str, float]:
    """Statistics of the records of one control value."""
    ok = recs[recs["status"] == TC_STATUS_OK]
    pinf, n = ok["pinf"], ok.size
    gap = ok[np.isfinite(ok["gap"])]
    chi = _sd(pinf, ddof=1) * np.sqrt(n_nodes)
    kurt = pearson_kurtosis(pinf) if n >= TC_KURT_MIN_N else np.nan
    return {
        "ctrl": recs["ctrl"][0],
        "n_rep": n,
        "n_bad": recs.size - n,
        "n_tie": int(np.sum(ok["tie"] == 1)),
        "pinf": _mean(pinf),
        "pinf2": _mean(pinf**2),
        "pinf_se": _sd(pinf, ddof=1) / np.sqrt(n) if n else np.nan,
        "kurt": kurt,
        "chi": chi,
        "chi0": _sd(pinf, ddof=0) * np.sqrt(n_nodes),
        "chi_se_jk": jackknife_chi_se(pinf, n_nodes),
        "chi_se_kurt": kurtosis_chi_se(chi, kurt, n),
        "n_neg": _mean(ok["n_neg"]),
        "p_eff": _mean(ok["p_eff"]),
        "lam0": _mean(ok["lam0"]),
        "e0": _mean(ok["e0"]),
        "smean": _ratio(ok["sum_s2_p"].sum(), ok["sum_s_p"].sum()),
        "n_gap": gap.size,
        "gap": _mean(gap["gap"]),
        "e1": _mean(gap["e1"]),
        "t_eig_med": _finite_stat(np.median, recs["t_eig"]),
        "t_eig_max": _finite_stat(np.max, recs["t_eig"]),
        "n_mv_med": _finite_stat(np.median, recs["n_mv"][recs["n_mv"] >= 0]),
    }


def _mean(x: NDArray) -> float:
    return float(np.mean(x)) if x.size else np.nan


def _sd(x: NDArray, ddof: int) -> float:
    """Sample sd; exactly 0 for a constant ``x`` (np.std gives round-off,
    1e-16 at P_inf = 0.998)."""
    if x.size <= ddof:
        return np.nan
    return float(np.std(x, ddof=ddof)) if np.ptp(x) > 0 else 0.0


def _ratio(num: float, den: float) -> float:
    return float(num / den) if den else np.nan


def _finite_stat(stat, x: NDArray) -> float:
    x = x[np.isfinite(x)]
    return float(stat(x)) if x.size else np.nan


def pearson_kurtosis(x: NDArray) -> float:
    """m4 / m2^2 of the central sample moments; NaN if ``x`` is constant."""
    if x.size < 2 or np.ptp(x) == 0:
        return np.nan
    y = x - x.mean()
    return float(np.mean(y**4) / np.mean(y**2) ** 2)


def kurtosis_chi_se(chi: float, kurt: float, n: int) -> float:
    """Standard error of chi = sqrt(N) sd(P_inf) from the kurtosis.

    The sd of a sample sd is sd sqrt((kurt - 1 + 2 / (n - 1)) / (4 n)) to
    leading order (delta method on the exact variance of the sample
    variance). NaN below ``TC_KURT_MIN_N`` records or without a kurtosis.

    ``kurt`` is the plug-in ``pearson_kurtosis``, biased low at small ``n``
    and more so on heavy tails, and so is this SE. On the flank of a chi
    peak (RRG k = 3, N = 300, kurtosis 11) it measured about 0.55, 0.66 and
    0.8 of the true sd of chi at n = 10, 20 and 50, ``jackknife_chi_se``
    0.76, 0.83 and 0.9; at the peak (kurtosis 2.2), 0.93 or more.
    """
    if n < TC_KURT_MIN_N:
        return np.nan
    if chi == 0:
        return 0.0
    if not np.isfinite(kurt):
        return np.nan
    spread = max(kurt - 1.0, 0.0)  # b2 >= 1; clamp the round-off
    return float(chi * np.sqrt((spread + 2 / (n - 1)) / (4 * n)))


def jackknife_chi_se(x: NDArray, n_nodes: int) -> float:
    """Jackknife standard error of chi = sqrt(N) sd(x, ddof=1)."""
    n = x.size
    if n < _JACKKNIFE_MIN:
        return np.nan
    y = x - x.mean()
    s_i, q_i = y.sum() - y, (y**2).sum() - y**2
    var_i = np.maximum((q_i - s_i**2 / (n - 1)) / (n - 2), 0.0)
    chi_i = np.sqrt(n_nodes * var_i)
    return float(np.sqrt((n - 1) / n * np.sum((chi_i - chi_i.mean()) ** 2)))


def pool_ps(records: NDArray, ps: Sequence[PS], ctrls: NDArray) -> dict:
    """P(S) of each control value summed over its converged records."""
    ok = records["status"] == TC_STATUS_OK
    ptr, sizes, counts = [0], [], []
    for ctrl in ctrls:
        idx = np.flatnonzero(ok & (records["ctrl"] == ctrl))
        size = concat_int64([ps[i][0] for i in idx])
        uniq, inverse = np.unique(size, return_inverse=True)
        pooled = np.zeros(uniq.size, dtype=np.int64)
        np.add.at(pooled, inverse, concat_int64([ps[i][1] for i in idx]))
        sizes.append(uniq)
        counts.append(pooled)
        ptr.append(ptr[-1] + uniq.size)
    return {
        "ps_ptr": np.array(ptr, dtype=np.int64),
        "ps_size": concat_int64(sizes),
        "ps_count": concat_int64(counts),
    }


# ---------------------------------------------------------------------------
# Chi peak
# ---------------------------------------------------------------------------
def chi_peak(ctrl: NDArray, chi: NDArray) -> float:
    """Vertex of the parabola through the chi maximum and its neighbours
    (the grid point itself at an edge or without curvature)."""
    i = int(np.argmax(chi))
    if i == 0 or i == chi.size - 1:
        return float(ctrl[i])
    (x0, x1, x2), (y0, y1, y2) = ctrl[i - 1 : i + 2], chi[i - 1 : i + 2]
    slope = (y1 - y0) / (x1 - x0)
    curvature = ((y2 - y1) / (x2 - x1) - slope) / (x2 - x0)
    if curvature >= 0:
        return float(x1)
    return float((x0 + x1) / 2 - slope / (2 * curvature))


def paired_matrix(
    records: NDArray, window: tuple[float, float]
) -> tuple[NDArray, NDArray]:
    """``P_inf[rep, j]`` at the control values ``ctrl[j]`` in ``window``,
    for the reps with a converged record at every one of them."""
    lo, hi = window
    ctrl = records["ctrl"]
    keep = (records["status"] == TC_STATUS_OK) & (ctrl >= lo) & (ctrl <= hi)
    ok = records[keep]
    ctrls, col = np.unique(ok["ctrl"], return_inverse=True)
    reps, row = np.unique(ok["rep"], return_inverse=True)
    matrix = np.full((reps.size, ctrls.size), np.nan)
    matrix[row, col] = ok["pinf"]
    return ctrls, matrix[~np.isnan(matrix).any(axis=1)]


def bootstrap_peak(
    records: NDArray,
    n_nodes: int,
    *,
    window: tuple[float, float],
    n_boot: int,
    seed: int,
) -> dict[str, float]:
    """Chi-peak location and its paired-bootstrap spread: whole reps are
    resampled, so the correlation between control values is kept.

    ``peak_lo`` and ``peak_hi`` are the ``TC_BOOT_QUANTILES`` (68 %)
    percentiles of the resampled peaks: the sampling spread of this
    estimator (grid argmax, then the parabola vertex) on the grid that was
    run, so a coarse grid bounds it.

    ``peak_at_edge`` flags a maximum at an end of the window, where the
    location is only a bound; ``peak_edge_frac`` is the bootstrap fraction
    of such maxima. Resamples with a constant chi (one rep drawn every
    time) have no peak and are left out. Below ``TC_BOOT_MIN_REPS`` reps
    there is no bootstrap: the interval fields are NaN, the peak is kept.
    """
    _check_window(window)
    ctrls, matrix = paired_matrix(records, window)
    n = matrix.shape[0]
    if ctrls.size < TC_PEAK_MIN_CTRLS or n < 2:  # sd with ddof 1 needs two reps
        return _NO_PEAK | {"peak_reps": n}
    peak, at_edge = _peak_of(ctrls, matrix, n_nodes)
    if np.isnan(peak):  # chi is constant over the window
        return _NO_PEAK | {"peak_reps": n}
    rng = np.random.default_rng(seed)
    n_boot = n_boot if n >= TC_BOOT_MIN_REPS else 0
    boot = [_peak_of(ctrls, matrix, n_nodes, rng) for _ in range(n_boot)]
    peaks = np.array([b[0] for b in boot], dtype=np.float64)
    edges = np.array([b[1] for b in boot], dtype=bool)
    ok = ~np.isnan(peaks)
    lo = hi = edge_frac = np.nan  # no usable resample (or n_boot = 0)
    if ok.any():
        lo, hi = np.quantile(peaks[ok], TC_BOOT_QUANTILES)
        edge_frac = float(edges[ok].mean())
    return {
        "peak": peak,
        "peak_lo": float(lo),
        "peak_hi": float(hi),
        "peak_sd": _sd(peaks[ok], ddof=1),
        "peak_at_edge": int(at_edge),
        "peak_edge_frac": edge_frac,
        "peak_reps": n,
    }


def _check_window(window: tuple[float, float]) -> None:
    """Raise unless LO < HI: a NaN bound or LO >= HI leaves at most one
    control value in the window, never a peak."""
    lo, hi = window
    if not (lo < hi):
        raise ValueError(f"the chi-peak window needs LO < HI, got ({lo}, {hi})")


def _peak_of(ctrls, matrix, n_nodes, rng=None) -> tuple[float, bool]:
    """Chi peak of ``matrix`` (or of one resample of its rows, the reps)
    and whether its maximum is a window end; NaN if chi is constant."""
    if rng is not None:
        matrix = matrix[rng.integers(0, matrix.shape[0], matrix.shape[0])]
    varies = np.ptp(matrix, axis=0) > 0  # a constant column has chi 0 exactly
    chi = np.sqrt(n_nodes) * np.where(varies, matrix.std(axis=0, ddof=1), 0.0)
    if np.ptp(chi) == 0:
        return np.nan, False
    i = int(np.argmax(chi))
    return chi_peak(ctrls, chi), i in (0, chi.size - 1)


# ---------------------------------------------------------------------------
# Collect, files and reports
# ---------------------------------------------------------------------------
def collect(
    rundir: Path,
    n_nodes: int,
    *,
    window: Optional[tuple[float, float]] = None,
    n_boot: int = DEFAULT_TC_BOOT,
    seed: int = DEFAULT_TC_BOOT_SEED,
    argv: Sequence[str] = (),
    legacy: bool = False,
) -> dict[str, NDArray]:
    """Aggregate ``rundir`` into ``agg.npz`` (and the legacy txt); return
    the arrays written. The ``window`` needs LO < HI, checked before any
    tile is read. A legacy txt left by an earlier collect is refreshed even
    without ``legacy``, so it never lags behind ``agg.npz``."""
    rundir = Path(rundir)
    window = window or (-np.inf, np.inf)
    _check_window(window)
    run = load_records(rundir)
    _check_n_nodes(run.cfg, n_nodes)
    redundant = check_topologies(run.records, run.source)
    _warn_run(rundir, run, redundant)
    summary = summarize(run.records, n_nodes)
    peak = bootstrap_peak(
        run.records, n_nodes, window=window, n_boot=n_boot, seed=seed
    )
    meta = {
        "n_nodes": n_nodes,
        "n_units": run.records.size,
        "n_dupes": run.n_dupes,
        "n_dupes_differ": run.n_dupes_differ,
        "n_tiles": run.n_tiles,
        "n_redundant_topo": redundant,
        **peak,
        "window_lo": window[0],
        "window_hi": window[1],
        "n_boot": n_boot,
        "boot_seed": seed,
    }
    agg = {
        "schema": np.int64(TC_AGG_SCHEMA),
        "identity": np.str_(run.identity),
        "cfg": np.str_(run.cfg),
        "provenance": np.str_(provenance(argv)),
        "tile_git": np.str_(json.dumps(run.git_states, sort_keys=True)),
        "summary": summary,
        **pool_ps(run.records, run.ps, summary["ctrl"]),
        "meta": np.array(tuple(meta[k] for k in META_DTYPE.names), META_DTYPE),
    }
    write_npz(rundir / TC_AGG_NAME, agg)
    _write_legacy_txt(rundir / TC_LEGACY_TXT_NAME, summary, legacy)
    return agg


def _write_legacy_txt(path: Path, summary: NDArray, legacy: bool) -> None:
    """Write the legacy txt if ``legacy``, or refresh it, with a note on
    stderr, if an earlier collect wrote it: it would miss the tiles added
    since. It is never deleted."""
    if not (legacy or path.exists()):
        return
    rows = legacy_rows(summary)
    atomic_write_bytes(path, lambda fh: np.savetxt(fh, rows, fmt=TC_LEGACY_FMT))
    if not legacy:
        print(f"note: refreshed {path} to match {TC_AGG_NAME}", file=sys.stderr)


def _check_n_nodes(cfg: str, n_nodes: int) -> None:
    """Raise unless ``n_nodes`` is the N of ``cfg.json`` (when it has one):
    chi scales with sqrt(N)."""
    recorded = json.loads(cfg).get("graph", {}).get("N")
    if recorded is not None and recorded != n_nodes:
        raise ValueError(
            f"n_nodes = {n_nodes}, but {CFG_SIDECAR_NAME} has N = {recorded}"
        )


def _warn_run(rundir: Path, run: RunRecords, redundant: int) -> None:
    """Warnings about the data behind a collect (they do not stop it). Each
    names the run, so a collect of several runs warns about each."""
    notes = []
    if redundant:
        notes.append(f"{redundant} reps repeat another rep's topology")
    if run.n_dupes_differ:
        notes.append(
            f"{run.n_dupes_differ} duplicate units differ beyond timing "
            "(e.g. another --gap-every); kept the converged, gap-solved ones"
        )
    clean = [g for g in run.git_states if g and not g["dirty"]]
    if len(run.git_states) > 1 or len(clean) != len(run.git_states):
        notes.append(
            "tiles do not come from one clean commit: "
            f"{json.dumps(run.git_states, sort_keys=True)}"
        )
    for note in notes:
        warnings.warn(f"{rundir}: {note}")


def load_agg(rundir: Path) -> dict[str, NDArray]:
    """The arrays of ``<rundir>/agg.npz`` (see the module docstring)."""
    with np.load(Path(rundir) / TC_AGG_NAME, allow_pickle=False) as npz:
        agg = {name: npz[name] for name in npz.files}
    if int(agg["schema"]) != TC_AGG_SCHEMA:
        raise ValueError(f"agg schema {int(agg['schema'])} != {TC_AGG_SCHEMA}")
    return agg


def legacy_rows(summary: NDArray) -> NDArray[np.float64]:
    """Legacy ordParam columns: reps, control value, <n_neg>, <P_inf>,
    <P_inf^2> and sd(P_inf) with ddof 0, as the legacy programs wrote it.

    Unlike the legacy per-p files this is one table sorted by control value.
    Column 1 counts converged records (the legacy programs counted
    attempts); column 2 is sigma for ``-ew normal``.
    """
    sd0 = np.sqrt(np.maximum(summary["pinf2"] - summary["pinf"] ** 2, 0.0))
    columns = ("n_rep", "ctrl", "n_neg", "pinf", "pinf2")
    return np.column_stack([summary[c] for c in columns] + [sd0])


def cost_table(summary: NDArray) -> str:
    """Per control value: reps, failures, P_inf, chi and the solve cost."""
    head = "    ctrl   reps  bad  <P_inf>      chi  t_eig med/max [s]  matvecs"
    lines = [head]
    for row in summary:
        lines.append(
            f"{row['ctrl']:8.4g} {row['n_rep']:6d} {row['n_bad']:4d} "
            f"{row['pinf']:8.4f} {row['chi']:8.4f} "
            f"{row['t_eig_med']:8.3g} {row['t_eig_max']:8.3g} "
            f"{row['n_mv_med']:8.4g}"
        )
    return "\n".join(lines)

"""TransCluster v2 statistics (kernels/TransClusterStats) on synthetic rows."""

import hashlib
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

SRC = Path(__file__).resolve().parents[2] / "src"
sys.path.insert(0, str(SRC))

import kernels.TransClusterStats as tcs  # noqa: E402
import kernels.TransClusterTiles as tiles  # noqa: E402
from kernels.TransClusterCore import RECORD_DTYPE, make_record  # noqa: E402
from kernels.TransClusterTiles import (  # noqa: E402
    Tile,
    TileBuffer,
    TileHeader,
    TileUnreadable,
    tile_paths,
    write_tile,
)
from lrgsglib.config.progargs.defs.TransCluster import (  # noqa: E402
    TC_AGG_NAME,
    TC_BOOT_MIN_REPS,
    TC_KURT_MIN_N,
    TC_LEGACY_TXT_NAME,
    TC_STATUS_NOCONV,
    TC_STATUS_OK,
)
from lrgsglib.statsys._naming import (  # noqa: E402
    cfg_sidecar_bytes,
    publish_cfg_sidecar,
)

N_NODES = 400
CTRLS = (0.0, 0.1, 0.25, 0.4)
N_REPS = 12
PEAK = 0.137  # centre of the synthetic chi bump, vertex of the parabola
IDENTITY = {"graph": {"N": N_NODES}, "model": {"f": "S"}, "run": {"s": 1}}
GAP_EVERY = 0  # the synthetic records solve no lambda_1
GIT = {"commit": "c0ffee", "dirty": False}  # code state of the tile writers


@pytest.fixture(autouse=True)
def _one_clean_commit(monkeypatch):
    """Tiles written from one clean commit, whatever the checkout's state."""
    monkeypatch.setattr(tiles, "git_state", lambda: GIT)


# ---------------------------------------------------------------------------
# Synthetic records
# ---------------------------------------------------------------------------
def _pinf(rep: int, ctrl: float) -> float:
    """Deterministic P_inf whose spread over reps peaks near ``PEAK``."""
    width = 0.01 + 0.05 * np.exp(-(((ctrl - PEAK) / 0.1) ** 2))
    return 1.0 - ctrl + width * np.sin(3.0 * rep + 7.0 * ctrl)


def _record(rep: int, ictrl: int, **fields) -> np.ndarray:
    ctrl = CTRLS[ictrl]
    base = {
        "rep": rep,
        "ictrl": ictrl,
        "ctrl": ctrl,
        "cnt": int(round(100 * ctrl)),
        "attempt": 0,
        "topo_fp": 1000 + rep,
        "n_neg": int(round(100 * ctrl)),
        "p_eff": ctrl,
        "status": TC_STATUS_OK,
        "lam0": ctrl,
        "t_eig": 0.01 * (rep + 1),
        "n_mv": 10 + rep,
        "tie": int(rep == 3),
        "pinf": _pinf(rep, ctrl),
        "sum_s_p": 2.0,
        "sum_s2_p": 6.0,
        "e0": -1.0 + ctrl,
    }
    return make_record({**base, **fields})


def _ps(rep: int, ictrl: int):
    """P(S) of one record: sizes 1 and 2 with counts that depend on it."""
    return (
        np.array([1, 2], dtype=np.int64),
        np.array([rep + 1, ictrl + 1], dtype=np.int64),
    )


def _table(rows) -> np.ndarray:
    return np.array([np.asarray(r)[()] for r in rows], dtype=RECORD_DTYPE)


def _assert_table_equal(a: np.ndarray, b: np.ndarray) -> None:
    """Field-wise equality (``==`` on structured rows treats NaN as unequal)."""
    assert a.dtype == b.dtype
    for name in a.dtype.names:
        np.testing.assert_array_equal(a[name], b[name], name)


def _grid_records():
    rows = [_record(r, i) for r in range(N_REPS) for i in range(len(CTRLS))]
    ps = [_ps(r, i) for r in range(N_REPS) for i in range(len(CTRLS))]
    return _table(rows), ps


def _write_tile(rundir: Path, ghash: str, tile: Tile, rows, grid=CTRLS):
    """One complete tile of the run with ``rows``, (record, P(S)) pairs."""
    identity = hashlib.sha256(cfg_sidecar_bytes(**IDENTITY)).hexdigest()
    header = TileHeader(identity, np.array(grid), GAP_EVERY, ("test",))
    buf = TileBuffer()
    for record, ps in rows:
        buf.add(record, ps)
    path = tile_paths(rundir, ghash, tile).final
    return write_tile(path, buf, header, tile, complete=True)


def _write_run(rundir: Path, split: int = 2, record=_record) -> None:
    """cfg.json plus the synthetic records as ``split`` complete tiles."""
    publish_cfg_sidecar(rundir, **IDENTITY)
    bounds = np.linspace(0, N_REPS, split + 1).astype(int).tolist()
    for r0, r1 in zip(bounds[:-1], bounds[1:]):
        rows = [
            (record(rep, ictrl), _ps(rep, ictrl))
            for rep in range(r0, r1)
            for ictrl in range(len(CTRLS))
        ]
        _write_tile(rundir, "deadbeef", Tile(0, len(CTRLS), r0, r1), rows)


# ---------------------------------------------------------------------------
# Per control value
# ---------------------------------------------------------------------------
def test_summary_chi_both_ddof_and_failed_rows():
    records, _ = _grid_records()
    bad = _record(N_REPS, 1, status=TC_STATUS_NOCONV, pinf=np.nan)
    summary = tcs.summarize(_table([*records, bad]), N_NODES)
    assert summary["ctrl"].tolist() == list(CTRLS)
    for row in summary:
        pinf = np.array([_pinf(r, row["ctrl"]) for r in range(N_REPS)])
        sqrt_n = np.sqrt(N_NODES)
        assert row["n_rep"] == N_REPS
        assert row["pinf"] == pytest.approx(pinf.mean(), rel=1e-12)
        assert row["pinf2"] == pytest.approx(np.mean(pinf**2), rel=1e-12)
        assert row["chi"] == pytest.approx(sqrt_n * pinf.std(ddof=1))
        assert row["chi0"] == pytest.approx(sqrt_n * pinf.std(ddof=0))
        assert row["n_tie"] == 1
        assert row["smean"] == pytest.approx(3.0)
    assert summary["n_bad"].tolist() == [0, 1, 0, 0]


@pytest.mark.parametrize("n_rep", [2, TC_KURT_MIN_N - 1, TC_KURT_MIN_N])
def test_kurtosis_and_its_chi_error_need_tc_kurt_min_n_reps(n_rep):
    rows = [_record(r, i) for r in range(n_rep) for i in range(len(CTRLS))]
    summary = tcs.summarize(_table(rows), N_NODES)
    assert np.all(summary["n_rep"] == n_rep)
    assert np.all(np.isfinite(summary["chi"]))  # ddof 1 needs two reps
    defined = n_rep >= TC_KURT_MIN_N
    for name in ("kurt", "chi_se_kurt"):
        assert np.all(np.isfinite(summary[name]) == defined), name


def test_a_constant_p_inf_has_chi_0_and_no_kurtosis():
    """np.std of twelve copies of 0.998 is about 2e-16, not 0."""
    rows = [_record(r, 1, pinf=0.998) for r in range(N_REPS)]
    row = tcs.summary_row(_table(rows), N_NODES)
    assert row["chi"] == row["chi0"] == row["pinf_se"] == 0.0
    assert row["chi_se_kurt"] == row["chi_se_jk"] == 0.0
    assert np.isnan(row["kurt"])


def test_the_matvec_median_leaves_out_the_dense_path():
    """Dense solves record n_mv = -1: they make no matvecs."""
    rows = [_record(r, 1, n_mv=-1) for r in range(N_REPS)]
    assert np.isnan(tcs.summary_row(_table(rows), N_NODES)["n_mv_med"])


def test_jackknife_matches_leave_one_out():
    x = np.random.default_rng(3).gamma(2.0, size=15)
    chi_i = [np.sqrt(N_NODES) * np.delete(x, i).std(ddof=1) for i in range(15)]
    chi_i = np.array(chi_i)
    brute = np.sqrt(14 / 15 * np.sum((chi_i - chi_i.mean()) ** 2))
    assert tcs.jackknife_chi_se(x, N_NODES) == pytest.approx(brute, rel=1e-10)
    assert np.isnan(tcs.jackknife_chi_se(x[:2], N_NODES))


def test_kurtosis_and_the_standard_error_of_chi():
    x = np.random.default_rng(4).exponential(size=50)
    expected = stats.kurtosis(x, fisher=False, bias=True)
    assert tcs.pearson_kurtosis(x) == pytest.approx(expected, rel=1e-12)
    assert np.isnan(tcs.pearson_kurtosis(np.ones(5)))
    chi, kurt, n = 2.0, 5.0, 100
    se = chi * np.sqrt((kurt - 1 + 2 / (n - 1)) / (4 * n))
    assert tcs.kurtosis_chi_se(chi, kurt, n) == pytest.approx(se)
    floor = chi * np.sqrt(2 / (n - 1) / (4 * n))  # kurt - 1 clamped at 0
    assert tcs.kurtosis_chi_se(chi, 1 - 1e-15, n) == pytest.approx(floor)
    assert np.isnan(tcs.kurtosis_chi_se(chi, kurt, TC_KURT_MIN_N - 1))
    assert tcs.kurtosis_chi_se(0.0, np.nan, n) == 0.0
    assert np.isnan(tcs.kurtosis_chi_se(chi, np.nan, n))


@pytest.mark.parametrize(
    "law, sd, kurt",  # population sd and Pearson kurtosis
    [
        ("normal", 1.0, 3.0),
        ("uniform", np.sqrt(1 / 12), 1.8),
        ("exponential", 1.0, 9.0),
    ],
)
def test_standard_error_of_chi_is_its_spread_over_runs(law, sd, kurt):
    """From the population kurtosis, ``kurtosis_chi_se`` is the sd of chi
    over many runs of ``n`` reps (for the exponential law kurt - 1 is four
    times the gaussian value, so the kurtosis term is tested)."""
    n, runs = 200, 4000
    x = getattr(np.random.default_rng(7), law)(size=(runs, n))
    chi = np.sqrt(N_NODES) * x.std(axis=1, ddof=1)
    se = tcs.kurtosis_chi_se(np.sqrt(N_NODES) * sd, kurt, n)
    assert chi.std(ddof=1) == pytest.approx(se, rel=0.05)


def test_pool_ps_sums_converged_records_only():
    records, ps = _grid_records()
    bad = _record(N_REPS, 0, status=TC_STATUS_NOCONV)
    pooled = tcs.pool_ps(
        _table([*records, bad]), [*ps, _ps(99, 99)], np.array(CTRLS)
    )
    ptr = pooled["ps_ptr"]
    assert ptr.tolist() == [0, 2, 4, 6, 8]
    for i in range(len(CTRLS)):
        size = pooled["ps_size"][ptr[i] : ptr[i + 1]]
        count = pooled["ps_count"][ptr[i] : ptr[i + 1]]
        assert size.tolist() == [1, 2]
        n = N_REPS
        assert count.tolist() == [n * (n + 1) // 2, n * (i + 1)]


# ---------------------------------------------------------------------------
# Loading: duplicates and topologies
# ---------------------------------------------------------------------------
def test_duplicates_keep_the_first_record():
    records, ps = _grid_records()
    twin = _record(2, 1, pinf=-5.0)  # same unit, other measured values
    table = _table([twin, *records])
    kept = tcs.drop_duplicate_units(table, [_ps(7, 7), *ps])
    assert kept.n_dupes == 1 and kept.n_differ == 1
    recs = kept.records
    assert recs.size == records.size and kept.source.size == recs.size
    order = np.lexsort((recs["ctrl"], recs["rep"]))
    np.testing.assert_array_equal(order, np.arange(recs.size))
    (i,) = np.flatnonzero((recs["rep"] == 2) & (recs["ictrl"] == 1))
    assert recs["pinf"][i] == -5.0  # the first one in input order
    assert kept.ps[i][1].tolist() == [8, 8]


def test_duplicates_prefer_converged_then_gap_solved_records():
    """Whatever the input (tile) order. ``n_differ`` counts the twins that
    differ beyond ``VOLATILE_FIELDS``: the timing and the grid index."""
    records, ps = _grid_records()
    failed = _record(2, 1, status=TC_STATUS_NOCONV, pinf=np.nan)
    gapped = _record(3, 2, lam1=0.5, gap=0.25, e1=-0.5)
    regridded = _record(4, 0, t_eig=9.0)  # the same unit on another grid
    regridded["ictrl"] = 7
    table = _table([failed, *records, gapped, regridded])
    source = np.array(["new", *["old"] * records.size, "new", "new"], object)
    ps = [_ps(0, 0), *ps, _ps(0, 0), _ps(0, 0)]
    kept = tcs.drop_duplicate_units(table, ps, source)
    assert kept.n_dupes == 3 and kept.n_differ == 2
    unit = {
        (int(r["rep"]), float(r["ctrl"])): (r, src)
        for r, src in zip(kept.records, kept.source)
    }
    record, src = unit[2, CTRLS[1]]
    assert record["status"] == TC_STATUS_OK and src == "old"
    record, src = unit[3, CTRLS[2]]
    assert record["gap"] == 0.25 and src == "new"
    record, src = unit[4, CTRLS[0]]
    assert record["ictrl"] == 0 and src == "old"  # the first of equals


@pytest.mark.parametrize(
    "copies, n_differ",
    [
        ([{"lam1": 0.5, "gap": 0.25, "e1": -0.5}, {}], 2),  # gap-solved kept
        ([{"pinf": -5.0}, {}], 1),  # the last copy equals the kept first one
    ],
)
def test_each_dropped_copy_is_compared_with_the_kept_one(copies, n_differ):
    """Unit (2, CTRLS[1]) three times: the grid's record, then ``copies``.
    ``n_differ`` counts the dropped copies unlike the kept record, not
    those unlike the copy sorted just before them."""
    records, ps = _grid_records()
    extra = [_record(2, 1, **fields) for fields in copies]
    ps = [*ps, *(_ps(2, 1) for _ in extra)]
    kept = tcs.drop_duplicate_units(_table([*records, *extra]), ps)
    assert kept.n_dupes == 2 and kept.n_differ == n_differ


@pytest.mark.parametrize("field", ["cnt", "attempt", "topo_fp", "n_neg"])
def test_duplicates_must_agree_on_the_unit_fields(field):
    records, ps = _grid_records()
    twin = _record(2, 1, **{field: 77})
    with pytest.raises(ValueError, match=field):
        tcs.drop_duplicate_units(_table([*records, twin]), [*ps, _ps(2, 1)])


def test_topologies_one_per_rep():
    records, _ = _grid_records()
    assert tcs.check_topologies(records) == 0
    shared = records.copy()
    shared["topo_fp"][shared["rep"] == 5] = 1000 + 4  # rep 5 reuses rep 4's
    assert tcs.check_topologies(shared) == 1
    clash = records.copy()
    clash["topo_fp"][(clash["rep"] == 5) & (clash["ictrl"] == 2)] = 1
    with pytest.raises(ValueError, match="rep 5"):
        tcs.check_topologies(clash)


# ---------------------------------------------------------------------------
# Chi peak
# ---------------------------------------------------------------------------
def test_chi_peak_parabola_vertex_on_a_nonuniform_grid():
    ctrl = np.array([0.0, 0.05, 0.12, 0.15, 0.3])
    chi = 4.0 - 50.0 * (ctrl - PEAK) ** 2
    assert tcs.chi_peak(ctrl, chi) == pytest.approx(PEAK, abs=1e-12)
    assert tcs.chi_peak(ctrl, ctrl) == ctrl[-1]  # maximum at the edge
    tied = np.array([0.0, 1.0, 2.0, 2.0, 0.0])  # the vertex lies between
    assert tcs.chi_peak(ctrl, tied) == pytest.approx((ctrl[2] + ctrl[3]) / 2)


def test_paired_matrix_keeps_reps_complete_over_the_window():
    records, _ = _grid_records()
    table = records[~((records["rep"] == 4) & (records["ictrl"] == 2))]
    ctrls, matrix = tcs.paired_matrix(table, (0.05, 0.5))
    assert ctrls.tolist() == list(CTRLS[1:])
    assert matrix.shape == (N_REPS - 1, 3)
    _, whole = tcs.paired_matrix(table, (0.0, 0.1))
    assert whole.shape == (N_REPS, 2)  # rep 4 is complete in this window


def test_bootstrap_is_reproducible_and_brackets_the_peak():
    records, _ = _grid_records()
    kwargs = {"window": (-np.inf, np.inf), "n_boot": 200}
    a = tcs.bootstrap_peak(records, N_NODES, seed=1, **kwargs)
    b = tcs.bootstrap_peak(records, N_NODES, seed=1, **kwargs)
    c = tcs.bootstrap_peak(records, N_NODES, seed=2, **kwargs)
    assert a == b and a != c
    assert a["peak_reps"] == N_REPS
    assert a["peak_lo"] <= PEAK <= a["peak_hi"] and a["peak_sd"] > 0
    assert a["peak_at_edge"] == 0 and a["peak_edge_frac"] == 0
    skipped = tcs.bootstrap_peak(
        records, N_NODES, window=(0, 1), n_boot=0, seed=1
    )
    assert skipped["peak"] == a["peak"] and skipped["peak_at_edge"] == 0
    assert np.isnan(skipped["peak_lo"]) and np.isnan(skipped["peak_sd"])
    assert np.isnan(skipped["peak_edge_frac"])
    narrow = tcs.bootstrap_peak(
        records, N_NODES, window=(0, 0.1), n_boot=9, seed=1
    )
    assert np.isnan(narrow["peak"]) and narrow["peak_reps"] == N_REPS
    assert narrow["peak_at_edge"] == -1


def test_bootstrap_flags_a_maximum_at_a_window_end():
    """The chi maximum is at CTRLS[1]: the low end of (0.05, 0.5), inside
    (0, 0.25). At an end the location is only a bound."""
    records, _ = _grid_records()
    kwargs = {"n_boot": 200, "seed": 1}
    end = tcs.bootstrap_peak(records, N_NODES, window=(0.05, 0.5), **kwargs)
    assert end["peak_at_edge"] == 1 and end["peak"] == CTRLS[1]
    assert end["peak_edge_frac"] > 0.5  # most resamples peak there too
    inner = tcs.bootstrap_peak(records, N_NODES, window=(0, 0.25), **kwargs)
    assert inner["peak_at_edge"] == 0 and CTRLS[0] < inner["peak"] < CTRLS[2]
    assert inner["peak_edge_frac"] < 0.5


@pytest.mark.parametrize("n_rep", [2, TC_BOOT_MIN_REPS - 1, TC_BOOT_MIN_REPS])
def test_the_bootstrap_needs_tc_boot_min_reps(n_rep):
    """Two reps have 3 distinct resamples, and every usable one is the pair
    itself: below the minimum the interval would collapse onto the peak,
    so it is NaN. The peak of the data is kept."""
    records, _ = _grid_records()
    kwargs = {"window": (-np.inf, np.inf), "n_boot": 200, "seed": 1}
    few = tcs.bootstrap_peak(records[records["rep"] < n_rep], N_NODES, **kwargs)
    assert np.isfinite(few["peak"]) and few["peak_reps"] == n_rep
    assert few["peak_at_edge"] in (0, 1)
    interval = [few[key] for key in ("peak_lo", "peak_hi", "peak_sd")]
    defined = n_rep >= TC_BOOT_MIN_REPS
    assert np.isfinite(interval + [few["peak_edge_frac"]]).all() == defined
    assert np.isnan(interval + [few["peak_edge_frac"]]).all() != defined


def test_no_chi_peak_is_flagged_minus_one():
    """One rep has no chi; a chi that is constant over the window has no
    maximum, also when it is 0 (np.std alone would leave round-off with a
    maximum). All give the no-peak result (flag -1, NaN fields)."""
    records, _ = _grid_records()
    kwargs = {"window": (-np.inf, np.inf), "n_boot": 50, "seed": 1}
    one = tcs.bootstrap_peak(records[records["rep"] < 1], N_NODES, **kwargs)
    flat = records.copy()
    flat["pinf"] = 0.5 + 0.01 * flat["rep"]  # the same spread at every ctrl
    const = tcs.bootstrap_peak(flat, N_NODES, **kwargs)
    steps = records.copy()  # one P_inf for every rep at each ctrl
    steps["pinf"] = np.array([0.998, 0.997, 0.993, 0.991])[steps["ictrl"]]
    zero = tcs.bootstrap_peak(steps, N_NODES, **kwargs)
    nan_fields = [key for key in tcs._NO_PEAK if key != "peak_at_edge"]
    for result, n_rep in ((one, 1), (const, N_REPS), (zero, N_REPS)):
        assert result["peak_at_edge"] == -1 and result["peak_reps"] == n_rep
        assert np.isnan([result[key] for key in nan_fields]).all()


@pytest.mark.parametrize(
    "window", [(np.nan, 0.3), (0.1, np.nan), (0.25, 0.25), (0.4, 0.1)]
)
def test_a_chi_peak_window_needs_lo_below_hi(tmp_path, window):
    """A NaN bound or LO >= HI leaves at most one control value, never a
    peak. ``bootstrap_peak`` refuses it, and so does ``collect`` (notebooks
    call it directly) before it reads a tile: this run has none."""
    records, _ = _grid_records()
    with pytest.raises(ValueError, match="LO < HI"):
        tcs.bootstrap_peak(records, N_NODES, window=window, n_boot=0, seed=0)
    publish_cfg_sidecar(tmp_path, **IDENTITY)
    with pytest.raises(ValueError, match="LO < HI"):
        tcs.collect(tmp_path, N_NODES, window=window, n_boot=0, seed=0)


# ---------------------------------------------------------------------------
# Collect and files
# ---------------------------------------------------------------------------
def _shared_graph(rep: int, ictrl: int, **fields) -> np.ndarray:
    """A synthetic record in which rep 5 drew the same graph as rep 4."""
    topo_fp = 1000 + (4 if rep == 5 else rep)
    return _record(rep, ictrl, topo_fp=topo_fp, **fields)


def test_collect_writes_agg_and_the_legacy_file(tmp_path):
    _write_run(tmp_path)
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # one clean commit, nothing to flag
        agg = tcs.collect(tmp_path, N_NODES, n_boot=50, seed=0, legacy=True)
    again = tcs.load_agg(tmp_path)
    assert sorted(agg) == sorted(again)
    _assert_table_equal(agg["summary"], again["summary"])
    meta = again["meta"][()]
    assert meta["n_units"] == N_REPS * len(CTRLS) and meta["n_tiles"] == 2
    assert meta["n_dupes"] == 0 and meta["n_nodes"] == N_NODES
    assert meta["n_dupes_differ"] == 0 and meta["n_redundant_topo"] == 0
    assert meta["peak_at_edge"] == 0 and meta["peak_reps"] == N_REPS
    assert json.loads(str(again["tile_git"])) == [GIT]
    records, _ = _grid_records()
    expected = tcs.summarize(records, N_NODES)
    _assert_table_equal(again["summary"], expected)
    legacy = np.loadtxt(tmp_path / TC_LEGACY_TXT_NAME)
    rows = tcs.legacy_rows(expected)
    np.testing.assert_allclose(legacy, rows, rtol=1e-6, atol=1e-12)
    assert legacy.shape == (len(CTRLS), 6)


def test_collect_refreshes_the_legacy_file_of_an_earlier_collect(
    tmp_path, capsys
):
    """Without ``legacy`` no txt is made, but one that an earlier collect
    wrote is rewritten, with a note: a second grid would leave it behind."""
    _write_run(tmp_path)
    txt = tmp_path / TC_LEGACY_TXT_NAME
    tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)
    assert not txt.exists()
    tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0, legacy=True)
    assert capsys.readouterr().err == ""
    grid = (0.2, 0.3)  # a second grid: 0.3 is a new control value
    rows = [(_record(r, 1, ctrl=grid[1]), _ps(r, 1)) for r in range(N_REPS)]
    _write_tile(tmp_path, "0ther", Tile(1, 2, 0, N_REPS), rows, grid)
    summary = tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)["summary"]
    (note,) = capsys.readouterr().err.splitlines()
    assert "refreshed" in note and str(txt) in note
    expected = tcs.legacy_rows(summary)
    assert expected.shape == (len(CTRLS) + 1, 6)
    np.testing.assert_allclose(np.loadtxt(txt), expected, rtol=1e-6, atol=1e-12)


def test_collect_counts_a_unit_shared_by_two_tiles_once(tmp_path):
    _write_run(tmp_path)
    rows = [(_record(0, 1), _ps(0, 1))]  # rep 0 at CTRLS[1], again
    _write_tile(tmp_path, "0ther", Tile(1, 2, 0, 1), rows)
    meta = tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)["meta"][()]
    assert meta["n_dupes"] == 1 and meta["n_tiles"] == 3
    assert meta["n_dupes_differ"] == 0  # the same record twice


def test_collect_reports_differing_twins_and_repeated_graphs(tmp_path):
    _write_run(tmp_path, record=_shared_graph)
    twin = _shared_graph(0, 1, pinf=0.5)  # unit (0, CTRLS[1]), redone
    _write_tile(tmp_path, "0ther", Tile(1, 2, 0, 1), [(twin, _ps(0, 1))])
    with pytest.warns(UserWarning) as caught:
        meta = tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)["meta"][()]
    assert meta["n_dupes"] == meta["n_dupes_differ"] == 1
    assert meta["n_redundant_topo"] == 1
    text = " ".join(str(w.message) for w in caught)
    assert "differ beyond timing" in text and "repeat another rep" in text


@pytest.mark.parametrize(
    "fields, grid, problem",
    [
        ({"cnt": 77}, CTRLS, "disagree on cnt"),
        ({"ctrl": 0.3, "topo_fp": 1}, (0.2, 0.3), "two topologies"),
    ],
)
def test_collect_names_both_tiles_of_a_clash(tmp_path, fields, grid, problem):
    """Rep 2 again at ictrl 1: the same unit with another count, or a
    point of another grid on another graph."""
    _write_run(tmp_path, split=1)
    rows = [(_record(2, 1, **fields), _ps(2, 1))]
    _write_tile(tmp_path, "0ther", Tile(1, 2, 2, 3), rows, grid)
    names = [path.name for path in tcs.complete_tile_paths(tmp_path)]
    with pytest.raises(ValueError, match=problem) as err:
        tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)
    assert len(names) == 2 and all(name in str(err.value) for name in names)
    assert not (tmp_path / TC_AGG_NAME).exists()


@pytest.mark.parametrize("damage", ["empty", "truncated"])
def test_collect_names_an_unreadable_tile(tmp_path, damage):
    _write_run(tmp_path)
    path = tcs.complete_tile_paths(tmp_path)[0]
    data = path.read_bytes()
    path.write_bytes(b"" if damage == "empty" else data[: len(data) // 2])
    with pytest.raises(TileUnreadable, match=re.escape(path.name)):
        tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)
    assert not (tmp_path / TC_AGG_NAME).exists()


def test_collect_warns_about_each_run(tmp_path, monkeypatch):
    """--collect calls collect once per run, and a warning text is shown
    once per process: each text names its run, so none is dropped."""
    dirty = {"commit": "b0b", "dirty": True}
    monkeypatch.setattr(tiles, "git_state", lambda: dirty)
    runs = [tmp_path / "a", tmp_path / "b"]
    for rundir in runs:
        _write_run(rundir, split=1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default")  # a repeated text is not shown
        for rundir in runs:
            tcs.collect(rundir, N_NODES, n_boot=0, seed=0)
    texts = [str(w.message) for w in caught if "clean commit" in str(w.message)]
    assert len(texts) == len(runs)
    assert all(text.startswith(f"{run}: ") for run, text in zip(runs, texts))


def test_collect_refuses_another_n(tmp_path):
    """chi scales with sqrt(N): a wrong N would rescale it silently."""
    _write_run(tmp_path, split=1)
    with pytest.raises(ValueError, match=f"has N = {N_NODES}"):
        tcs.collect(tmp_path, N_NODES + 1, n_boot=0, seed=0)
    assert not (tmp_path / TC_AGG_NAME).exists()


def test_collect_lists_the_code_states_of_the_tiles(tmp_path, monkeypatch):
    _write_run(tmp_path)
    dirty = {"commit": "b0b", "dirty": True}
    monkeypatch.setattr(tiles, "git_state", lambda: dirty)
    rows = [(_record(0, 1), _ps(0, 1))]
    _write_tile(tmp_path, "0ther", Tile(1, 2, 0, 1), rows)
    with pytest.warns(UserWarning, match="one clean commit"):
        agg = tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)
    tile_git = json.loads(str(agg["tile_git"]))
    assert len(tile_git) == 2 and GIT in tile_git and dirty in tile_git


def test_load_agg_rejects_another_schema(tmp_path, monkeypatch):
    _write_run(tmp_path, split=1)
    tcs.collect(tmp_path, N_NODES, n_boot=0, seed=0)
    monkeypatch.setattr(tcs, "TC_AGG_SCHEMA", tcs.TC_AGG_SCHEMA + 1)
    with pytest.raises(ValueError, match="schema"):
        tcs.load_agg(tmp_path)


def test_collect_needs_complete_tiles(tmp_path):
    publish_cfg_sidecar(tmp_path, **IDENTITY)
    with pytest.raises(FileNotFoundError, match="no complete tiles"):
        tcs.collect(tmp_path, N_NODES)
    assert not (tmp_path / TC_AGG_NAME).exists()


def test_cost_table_has_a_row_per_control_value():
    records, _ = _grid_records()
    table = tcs.cost_table(tcs.summarize(records, N_NODES))
    assert len(table.splitlines()) == 1 + len(CTRLS)

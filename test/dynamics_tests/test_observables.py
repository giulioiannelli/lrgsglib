"""Phase-2 spine: the generic observable codecs (:mod:`statsys._io`) and the
reusable ``Observable`` / ``ObservableSet`` objects (:mod:`statsys._observables`)
that VoterModel (and, in later phases, every dynamics model) persists through.

These pin the round-trip + lazy-view + streaming/batch contract the disk-backed
observables rely on, independent of any particular model.
"""

import numpy as np
import pytest

from lrgsglib.statsys import _io
from lrgsglib.statsys._observables import (
    HistogramSeries,
    ObservableSet,
    PersistMode,
    RowTrajectory,
    ScalarSeries,
)


# ---------------------------------------------------------------------------
# _io codecs
# ---------------------------------------------------------------------------
@pytest.mark.code
def test_io_array_roundtrip(tmp_path):
    vals = [0.0, -1.5, 3.25, 7.0]
    p = _io.save_array(tmp_path / "a.bin", vals, np.float64)
    assert np.array_equal(_io.load_array(p, np.float64), np.asarray(vals))


@pytest.mark.code
def test_io_rows_roundtrip_and_mmap(tmp_path):
    rows = np.array([[1, -1, 1], [-1, -1, 1]], dtype=np.int8)
    p = _io.save_rows(tmp_path / "r.bin", rows, np.int8)
    loaded = _io.load_rows(p, ncols=3, dtype=np.int8)
    assert np.array_equal(loaded, rows)
    mm = _io.load_rows(p, ncols=3, dtype=np.int8, mmap=True)
    assert isinstance(mm, np.memmap) and np.array_equal(mm, rows)


@pytest.mark.code
def test_io_row_stream_roundtrip(tmp_path):
    path = tmp_path / "s.bin"
    fh, tmp = _io.open_row_stream(path)
    rows = [np.array([1, -1], np.int8), np.array([-1, 1], np.int8)]
    for r in rows:
        _io.write_row(fh, r, np.int8)
    assert not path.exists()  # only materializes on close (atomic rename)
    _io.close_row_stream(fh, tmp, path)
    assert path.exists() and not tmp.exists()
    assert np.array_equal(_io.load_rows(path, 2, np.int8), np.asarray(rows))


@pytest.mark.code
def test_io_histogram_roundtrip(tmp_path):
    series = [{1: 4, 2: 1}, {}, {3: 2}]
    p = _io.save_histogram_series(tmp_path / "h.npz", series)
    assert _io.load_histogram_series(p) == series


@pytest.mark.code
def test_io_atomic_write_leaves_no_temp(tmp_path):
    _io.save_array(tmp_path / "a.bin", [1.0], np.float64)
    assert [q.name for q in tmp_path.iterdir()] == ["a.bin"]  # no .tmp sibling


# ---------------------------------------------------------------------------
# ScalarSeries
# ---------------------------------------------------------------------------
@pytest.mark.code
def test_scalar_series_buffer_then_lazy_load(tmp_path):
    obs = ScalarSeries("magn", dtype=np.float64)
    for v in (0.5, -0.5, 1.0):
        obs.append(v)
    assert obs.data == [0.5, -0.5, 1.0]      # RAM buffer
    obs.set_path(tmp_path / "m.bin")
    obs.persist()
    assert obs.persist_mode is PersistMode.BATCH
    # A fresh instance pointed at the file lazily loads it.
    fresh = ScalarSeries("magn", dtype=np.float64)
    fresh.set_path(tmp_path / "m.bin")
    assert fresh.data == [0.5, -0.5, 1.0]


@pytest.mark.code
def test_scalar_series_persist_noop_when_empty(tmp_path):
    obs = ScalarSeries("magn")
    obs.set_path(tmp_path / "m.bin")
    obs.persist()                            # nothing buffered -> no file
    assert not (tmp_path / "m.bin").exists()
    obs.reset()
    assert obs.data == [] and obs.path is None


# ---------------------------------------------------------------------------
# RowTrajectory
# ---------------------------------------------------------------------------
@pytest.mark.code
def test_row_trajectory_stream(tmp_path):
    obs = RowTrajectory("sout", ncols=2, dtype=np.int8)
    obs.set_path(tmp_path / "sout.bin")
    obs.open_stream()
    assert obs.streaming and obs.persist_mode is PersistMode.STREAM
    rows = [np.array([1, -1], np.int8), np.array([-1, -1], np.int8)]
    for r in rows:
        obs.stream_row(r)
    obs.close_stream()
    assert not obs.streaming and obs.nrec == 2
    obs.finalize_batch()                     # cache empty -> just drop it
    assert np.array_equal(np.asarray(obs.data), np.asarray(rows))
    assert isinstance(obs.data, np.memmap)


@pytest.mark.code
def test_row_trajectory_batch_with_subsample(tmp_path):
    obs = RowTrajectory("sout", ncols=2, dtype=np.int8)
    obs.set_path(tmp_path / "sout.bin")
    full = [np.array([i % 2, -(i % 2)], np.int8) for i in range(6)]
    obs.set_rows([r.copy() for r in full])
    obs.finalize_batch(subsample=2)          # keep rows 0,2,4
    assert obs.nrec == 3
    assert np.array_equal(np.asarray(obs.data), np.asarray(full[::2]))


@pytest.mark.code
def test_row_trajectory_cross_instance_reload(tmp_path):
    # The batch-loader / GUI path: a fresh instance pointed at an existing file
    # (cache still empty) must memmap it, like ScalarSeries/HistogramSeries do.
    rows = np.array([[1, -1], [-1, -1], [1, 1]], dtype=np.int8)
    _io.save_rows(tmp_path / "sout.bin", rows, np.int8)
    fresh = RowTrajectory("sout", ncols=2, dtype=np.int8)
    fresh.set_path(tmp_path / "sout.bin")
    assert np.array_equal(np.asarray(fresh.data), rows)


@pytest.mark.code
def test_row_trajectory_zero_byte_file_is_safe(tmp_path):
    # A degenerate zero-row file must not crash data (np.memmap rejects empty).
    (tmp_path / "sout.bin").write_bytes(b"")
    obs = RowTrajectory("sout", ncols=2, dtype=np.int8)
    obs.set_path(tmp_path / "sout.bin")
    obs.mark_on_disk()                       # cache -> None, force the disk branch
    assert obs.data == []


@pytest.mark.code
def test_row_trajectory_in_ram_when_no_path():
    obs = RowTrajectory("sout", ncols=2, dtype=np.int8)
    obs.append_row(np.array([1, 1], np.int8))
    obs.append_row(np.array([-1, 1], np.int8))
    assert len(obs.data) == 2                # stays a RAM list (savedisk=False path)
    obs.reset()
    assert obs.data == [] and obs.nrec == 0


# ---------------------------------------------------------------------------
# HistogramSeries
# ---------------------------------------------------------------------------
@pytest.mark.code
def test_histogram_series_buffer_then_lazy_load(tmp_path):
    obs = HistogramSeries("cldist")
    obs.append({4: 1})
    obs.append({1: 2, 2: 1})
    assert obs.data == [{4: 1}, {1: 2, 2: 1}]
    obs.set_path(tmp_path / "c.npz")
    obs.persist()
    fresh = HistogramSeries("cldist")
    fresh.set_path(tmp_path / "c.npz")
    assert fresh.data == [{4: 1}, {1: 2, 2: 1}]


# ---------------------------------------------------------------------------
# ObservableSet
# ---------------------------------------------------------------------------
@pytest.mark.code
def test_observable_set_access_and_reset(tmp_path):
    oset = ObservableSet([
        ScalarSeries("magn"),
        RowTrajectory("sout", ncols=3),
        HistogramSeries("cldist"),
    ])
    assert oset.names() == ["magn", "sout", "cldist"]
    assert "magn" in oset and "nope" not in oset
    assert len(oset) == 3
    assert {o.name for o in oset} == {"magn", "sout", "cldist"}
    oset["magn"].append(1.0)
    oset["magn"].set_path(tmp_path / "m.bin")
    oset.reset()
    assert oset["magn"].data == [] and oset["magn"].path is None

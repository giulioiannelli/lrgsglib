"""Reusable on-disk observables for dynamics models.

A dynamics model produces *observables* -- time series sampled as it steps (a
magnetization series, a spin trajectory, a cluster-size-distribution series, an
energy series, ...). This module factors the lifecycle those series share --
accumulate in RAM, optionally persist to disk, expose a **lazy** disk-backed
view -- out of any single model, on top of the generic codecs in :mod:`._io`.

Three archetypes cover the cases seen so far:

- :class:`ScalarSeries`    -- a 1-D numeric series (e.g. magnetization, energy).
- :class:`RowTrajectory`   -- a 2-D ``(n_rec, ncols)`` trajectory (e.g. spins),
                              either *streamed* row-by-row during the run or
                              *batch-dumped* from a whole-array backend result.
- :class:`HistogramSeries` -- a list of ``{int: int}`` histograms (e.g. a
                              cluster-size distribution per recorded sweep).

A model owns an :class:`ObservableSet` (``self.observables``) and drives it: it
decides *which* observables are enabled, *where* their files live, and -- for a
:class:`RowTrajectory` -- whether to stream or batch. The objects own the
*mechanics* (buffers, codecs, atomic writes, lazy load/memmap); the model owns
the *policy* (gating flags, subsampling, backend-specific dumps). ``savedisk``
is the model-level master switch: when off, a model simply never assigns paths,
so every observable stays purely in RAM.

The lazy-view contract every archetype honors:

- while the RAM buffer/cache is populated, :attr:`data` returns it verbatim;
- once the series lives on disk *and* the RAM copy has been dropped, :attr:`data`
  loads (or memmaps) the file on demand;
- with neither, :attr:`data` is an empty list.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping, Sequence
from enum import Enum
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import DTypeLike

from . import _io

__all__ = [
    "PersistMode",
    "Observable",
    "ScalarSeries",
    "RowTrajectory",
    "HistogramSeries",
    "ObservableSet",
]


class PersistMode(Enum):
    """How an observable's data is currently held.

    ``RAM`` -- in memory only (no path assigned; ``savedisk=False`` or pre-run).
    ``STREAM`` -- appended row-by-row to an open file as the run steps.
    ``BATCH`` -- accumulated in RAM, written once when the run finalizes.
    """

    RAM = "ram"
    STREAM = "stream"
    BATCH = "batch"


class Observable:
    """Base class: one named time series with an optional on-disk backing.

    Subclasses implement the codec-specific bits (:meth:`persist`, :attr:`data`,
    and how values accumulate); this base owns the shared name/path bookkeeping.
    """

    def __init__(self, name: str) -> None:
        self.name = name
        self.persist_mode = PersistMode.RAM
        self._path: Path | None = None

    @property
    def path(self) -> Path | None:
        """On-disk location of this observable (``None`` => not on disk)."""
        return self._path

    def set_path(self, path) -> None:
        """Point this observable at a file (used as the persistence target and
        the lazy-load source)."""
        self._path = Path(path) if path is not None else None

    def reset(self) -> None:
        """Clear all in-RAM buffers and on-disk-backing state."""
        self.persist_mode = PersistMode.RAM
        self._path = None

    def persist(self) -> None:
        """Write the accumulated RAM buffer to :attr:`path` (no-op if there is no
        path or nothing buffered)."""
        raise NotImplementedError

    @property
    def data(self) -> Any:
        """Lazy view: the RAM buffer if populated, else the disk file, else
        an empty list."""
        raise NotImplementedError


class ScalarSeries(Observable):
    """A 1-D numeric series (raw-binary ``dtype``), e.g. a magnetization series.

    Always batch-persisted: values accumulate in a RAM list and are dumped once
    at finalize. :attr:`data` returns the buffer (a Python ``list``) when
    populated, else lazily loads (and caches) the file.
    """

    def __init__(self, name: str, dtype: DTypeLike = np.float64) -> None:
        super().__init__(name)
        self.dtype = np.dtype(dtype)
        self._buf: list = []

    def reset(self) -> None:
        super().reset()
        self._buf = []

    def append(self, value) -> None:
        """Append one sample to the RAM buffer."""
        self._buf.append(value)

    def set_values(self, values: Iterable | None) -> None:
        """Replace the RAM buffer wholesale (e.g. a backend's returned series)."""
        self._buf = list(values) if values is not None else []

    def persist(self) -> None:
        if self._path is not None and self._buf:
            _io.save_array(self._path, self._buf, self.dtype)
            self.persist_mode = PersistMode.BATCH

    @property
    def data(self) -> list:
        if self._buf:
            return self._buf
        if self._path is not None and Path(self._path).exists():
            self._buf = _io.load_array(self._path, self.dtype).tolist()
        return self._buf


class RowTrajectory(Observable):
    """A 2-D ``(n_rec, ncols)`` trajectory (raw-binary ``dtype``), e.g. spins.

    Supports two persistence routes the backends actually use:

    - **stream** (the Python loop): :meth:`open_stream` then :meth:`stream_row`
      per recorded sweep keeps the full trajectory off-RAM; :meth:`close_stream`
      finalizes the file.
    - **batch** (whole-array backends, e.g. pybind): :meth:`set_rows` caches the
      returned trajectory, then :meth:`finalize_batch` dumps it (optionally
      subsampled) and drops the RAM copy.

    :attr:`data` returns the RAM cache while it is held (a ``list`` of rows), a
    read-only memmap once the trajectory lives on disk, or an empty list.
    """

    def __init__(self, name: str, ncols: int, dtype: DTypeLike = np.int8) -> None:
        super().__init__(name)
        self.ncols = int(ncols)
        self.dtype = np.dtype(dtype)
        # In-RAM trajectory rows. A *non-empty* cache is the live data; an empty
        # cache (``[]`` fresh, or ``None`` after a backend wrote the file) defers
        # ``data`` to the on-disk memmap. ``finalize_batch`` / ``mark_on_disk``
        # set it to ``None`` once the trajectory lives on disk.
        self._cache: list | None = []
        self._fh = None
        self._tmp: Path | None = None
        self._nrec = 0

    def reset(self) -> None:
        super().reset()
        self._cache = []
        self._fh = None
        self._tmp = None
        self._nrec = 0

    def reset_stream(self) -> None:
        """Clear streaming state (handle, temp path, recorded-row count) without
        touching the RAM cache or path -- called at the start of a run."""
        self._fh = None
        self._tmp = None
        self._nrec = 0

    @property
    def nrec(self) -> int:
        """Number of rows actually recorded to disk."""
        return self._nrec

    @property
    def streaming(self) -> bool:
        """True while a row stream is open."""
        return self._fh is not None

    def append_row(self, row) -> None:
        """Append one row to the in-RAM cache (the ``savedisk=False`` path)."""
        if self._cache is None:
            self._cache = []
        self._cache.append(row)

    def set_rows(self, rows: Sequence | None) -> None:
        """Cache a whole trajectory in RAM (e.g. a backend's returned array, or a
        direct ``s_t`` assignment)."""
        self._cache = rows

    def open_stream(self) -> None:
        """Open a row stream on :attr:`path` for row-by-row writing."""
        self._fh, self._tmp = _io.open_row_stream(self._path)
        self.persist_mode = PersistMode.STREAM

    def stream_row(self, row) -> None:
        """Write one row to the open stream and count it."""
        _io.write_row(self._fh, row, self.dtype)
        self._nrec += 1

    def close_stream(self) -> None:
        """Finalize a streamed trajectory (no-op if not streaming)."""
        if self._fh is not None:
            _io.close_row_stream(self._fh, self._tmp, self._path)
            self._fh = None
            self._tmp = None

    def finalize_batch(self, subsample: int = 1) -> None:
        """Dump a RAM-cached trajectory to disk (optionally subsampled) and drop
        the RAM copy so :attr:`data` memmaps the file.

        No-op when streaming already wrote the rows (cache empty); in that case it
        only drops the cache once the file exists.
        """
        if self._path is None:
            return
        if self._cache:
            arr = np.asarray(self._cache, dtype=self.dtype)
            if subsample > 1:
                arr = arr[::subsample]
            _io.save_rows(self._path, arr, self.dtype)
            self._nrec = int(arr.shape[0])
            self.persist_mode = PersistMode.BATCH
        if Path(self._path).exists():
            self._cache = None

    def mark_on_disk(self) -> None:
        """Drop the RAM cache so :attr:`data` memmaps :attr:`path` (used when a
        backend wrote the file itself)."""
        self._cache = None

    def persist(self) -> None:
        """Default persistence == batch dump (alias for :meth:`finalize_batch`)."""
        self.finalize_batch()

    @property
    def data(self):
        # A populated RAM cache wins (the savedisk=False / pre-persist path); an
        # empty cache falls through to the on-disk memmap, so a fresh instance
        # pointed at an existing file (the batch-loader / GUI path) reloads it --
        # consistent with ScalarSeries/HistogramSeries. The size guard avoids
        # ``np.memmap`` raising on a (degenerate) zero-row file.
        if self._cache:
            return self._cache
        p = self._path
        if p is not None and Path(p).exists() and Path(p).stat().st_size > 0:
            return _io.load_rows(p, self.ncols, self.dtype, mmap=True)
        return []


class HistogramSeries(Observable):
    """A list of ``{int: int}`` histograms (compressed ``.npz``), e.g. a
    cluster-size distribution per recorded sweep.

    Batch-persisted: histograms accumulate in a RAM list and are written once at
    finalize. :attr:`data` returns the buffer when populated, else lazily loads
    (and caches) the file.
    """

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self._buf: list[dict[int, int]] = []

    def reset(self) -> None:
        super().reset()
        self._buf = []

    def append(self, hist: Mapping[int, int]) -> None:
        """Append one histogram to the RAM buffer."""
        self._buf.append(dict(hist))

    def set_values(self, series: Iterable[Mapping[int, int]] | None) -> None:
        """Replace the RAM buffer wholesale (e.g. a backend's returned series)."""
        self._buf = list(series) if series is not None else []

    def persist(self) -> None:
        if self._path is not None and self._buf:
            _io.save_histogram_series(self._path, self._buf)
            self.persist_mode = PersistMode.BATCH

    @property
    def data(self) -> list[dict[int, int]]:
        if self._buf:
            return self._buf
        if self._path is not None and Path(self._path).exists():
            self._buf = _io.load_histogram_series(self._path)
        return self._buf


class ObservableSet:
    """An ordered, name-keyed collection of a model's :class:`Observable`\\ s.

    A thin container the model populates once (``self.observables = ObservableSet(
    [...])``) and accesses by name (``self.observables['magn']``). Gives a uniform
    surface for consumers (a GUI, a batch loader) to enumerate and read a model's
    series without knowing the model's internals.
    """

    def __init__(self, observables: Iterable[Observable] = ()) -> None:
        self._by_name: dict[str, Observable] = {}
        for obs in observables:
            self.add(obs)

    def add(self, obs: Observable) -> None:
        """Register an observable (keyed by its ``name``)."""
        self._by_name[obs.name] = obs

    def reset(self) -> None:
        """Reset every observable to its empty, in-RAM state."""
        for obs in self._by_name.values():
            obs.reset()

    def names(self) -> list[str]:
        """Registered observable names, in insertion order."""
        return list(self._by_name)

    def __getitem__(self, name: str) -> Observable:
        return self._by_name[name]

    def __contains__(self, name: object) -> bool:
        return name in self._by_name

    def __iter__(self) -> Iterator[Observable]:
        return iter(self._by_name.values())

    def __len__(self) -> int:
        return len(self._by_name)

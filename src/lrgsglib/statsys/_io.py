"""Generic, model-agnostic codecs for on-disk dynamics observables.

Generalized from ``VoterModel/_voter_io.py`` so *every* statsys model can stream
big simulation output to disk **during** a run and load it back only when needed,
keeping RAM bounded. Three archetypes, each mirroring the C-backend file
convention (``<base>_p=<pflip>[_suffix].<ext>`` under the run's ``dynpath``):

- **array**     -- a 1-D numeric series of any ``dtype`` (e.g. a magnetization or
                   energy series), raw binary.
- **rows**      -- a 2-D ``(n_rec, ncols)`` matrix of any ``dtype`` (e.g. a spin
                   trajectory), row-major raw binary. Streamable row-by-row *or*
                   dumped one-shot, and loadable as a memmap. ``ncols`` is not
                   stored in the stream, so the reader must supply it.
- **histogram** -- a list of ``{int: int}`` histograms (e.g. a cluster-size
                   distribution time series), stored compressed as three flat
                   ``int64`` arrays (``rec``/``size``/``count``) plus ``n_rec`` in
                   a single ``.npz``.

Every writer goes through a temp sibling file with a random suffix and an atomic
``os.replace``, so a crash mid-write never leaves a truncated file and concurrent
runs sharing a name never clobber a half-written one. These are pure functions
(path in, array out); a model's methods are thin wrappers, so the same format is
loadable without a live model instance.
"""

from __future__ import annotations

import io
import os
from collections.abc import Iterable, Mapping
from pathlib import Path

import numpy as np
from numpy.typing import DTypeLike

from ..utils.basic.strings import generate_random_id

__all__ = [
    "tmp_sibling",
    "atomic_write_bytes",
    "save_array",
    "load_array",
    "open_row_stream",
    "write_row",
    "close_row_stream",
    "save_rows",
    "load_rows",
    "save_histogram_series",
    "load_histogram_series",
]


def tmp_sibling(path: Path) -> Path:
    """Hidden, randomly-suffixed temp path in the same directory as ``path``
    (same filesystem -> ``os.replace`` is atomic)."""
    path = Path(path)
    return path.with_name(f".{path.name}.{generate_random_id()}.tmp")


def atomic_write_bytes(path, write_fn) -> Path:
    """Run ``write_fn(file_handle)`` on a temp file, then atomically rename it
    onto ``path``. The temp file is removed on any failure."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = tmp_sibling(path)
    try:
        with open(tmp, "wb") as fh:
            write_fn(fh)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass
    return path


# ---------------------------------------------------------------------------
# array: a 1-D numeric series (any dtype)
# ---------------------------------------------------------------------------
def save_array(path, values, dtype: DTypeLike) -> Path:
    """Atomically write a 1-D numeric series as raw ``dtype`` binary."""
    arr = np.asarray(values, dtype=dtype)
    return atomic_write_bytes(path, arr.tofile)


def load_array(path, dtype: DTypeLike) -> np.ndarray:
    """Load a series written by :func:`save_array`."""
    return np.fromfile(path, dtype=dtype)


# ---------------------------------------------------------------------------
# rows: a 2-D (n_rec, ncols) matrix (any dtype), row-major
# ---------------------------------------------------------------------------
def open_row_stream(path):
    """Open a streaming handle to a temp sibling of ``path``.

    Returns ``(file_handle, tmp_path)``; write rows with :func:`write_row` and
    finalize with :func:`close_row_stream` (which atomically renames the temp
    file onto ``path``). Streaming keeps the full trajectory off-RAM.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = tmp_sibling(path)
    return open(tmp, "wb"), tmp


def write_row(fh, row, dtype: DTypeLike) -> None:
    """Append one ``dtype`` row to an open row stream."""
    np.asarray(row, dtype=dtype).tofile(fh)


def close_row_stream(fh, tmp, path) -> Path:
    """Close the stream and atomically rename the temp file onto ``path``."""
    path = Path(path)
    fh.close()
    os.replace(tmp, path)
    return path


def save_rows(path, rows, dtype: DTypeLike) -> Path:
    """One-shot atomic dump of an ``(n_rec, ncols)`` matrix.

    Used by backends that return the whole trajectory at once; the streaming
    path writes row-by-row instead (see :func:`open_row_stream`).
    """
    arr = np.ascontiguousarray(np.asarray(rows, dtype=dtype))
    return atomic_write_bytes(path, arr.tofile)


def load_rows(path, ncols: int, dtype: DTypeLike, mmap: bool = False) -> np.ndarray:
    """Load a matrix written row-major as ``(n_rec, ncols)`` ``dtype``.

    With ``mmap=True`` the file is memory-mapped (read-only) and reshaped without
    reading it into RAM -- the right choice for a large trajectory.
    """
    if mmap:
        return np.memmap(path, dtype=dtype, mode="r").reshape(-1, ncols)
    return np.fromfile(path, dtype=dtype).reshape(-1, ncols)


# ---------------------------------------------------------------------------
# histogram: a list of {int: int} histograms (compressed .npz)
# ---------------------------------------------------------------------------
def save_histogram_series(path, series: Iterable[Mapping[int, int]]) -> Path:
    """Atomically write a list of ``{key: count}`` histograms.

    Stored as three flat ``int64`` arrays (record index, key, count) plus the
    total record count, so it round-trips exactly and stays compact for the
    typically-small per-record distributions.
    """
    series = list(series)
    recs: list[int] = []
    keys: list[int] = []
    counts: list[int] = []
    for t, dist in enumerate(series):
        for key, count in dist.items():
            recs.append(t)
            keys.append(int(key))
            counts.append(int(count))
    buf = io.BytesIO()
    np.savez_compressed(
        buf,
        rec=np.asarray(recs, dtype=np.int64),
        size=np.asarray(keys, dtype=np.int64),
        count=np.asarray(counts, dtype=np.int64),
        n_rec=np.int64(len(series)),
    )
    data = buf.getvalue()
    return atomic_write_bytes(path, lambda fh: fh.write(data))


def load_histogram_series(path) -> list[dict[int, int]]:
    """Reconstruct the list of ``{key: count}`` histograms from a ``.npz`` file."""
    with np.load(path) as z:
        rec = z["rec"].tolist()
        size = z["size"].tolist()
        count = z["count"].tolist()
        n_rec = int(z["n_rec"])
    out: list[dict[int, int]] = [dict() for _ in range(n_rec)]
    for t, sz, c in zip(rec, size, count):
        out[t][int(sz)] = int(c)
    return out

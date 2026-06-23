"""Pure save/load helpers for the VoterModel on-disk observables.

Mirrors the C-backend file convention -- ``<base>_p=<pflip>[_suffix].<ext>`` under
the run's ``voter/N=<n>/`` directory -- so an analysis program can stream big
simulation output to disk *during* a run and load it back only when needed,
keeping RAM bounded. Three artifacts:

- ``m``      : per-spin magnetization series, ``float64`` raw binary.
- ``sout``   : spin-configuration trajectory, ``int8`` raw binary, row-major
               ``(n_rec, N)`` (``n_rec`` = recorded sweeps; the reshape needs
               ``N``, which is not stored in the raw stream).
- ``cldist`` : cluster-size-distribution time series -- a list of
               ``{component_size: count}`` histograms -- stored compressed as
               three flat ``int64`` arrays (``rec``/``size``/``count``) plus
               ``n_rec``. DISTINCT from Ising's ``outcl{i}`` (see defaults.py).

Every writer goes through a temp file with a random suffix and an atomic
``os.replace``, so a crash mid-write never leaves a truncated file and
concurrent runs sharing a name never clobber a half-written one. These are
free functions (path in, array out); the VoterModel methods are thin wrappers
so the same format is loadable without a live model instance.
"""

from __future__ import annotations

import io
import os
from pathlib import Path

import numpy as np

from ...utils.basic.strings import generate_random_id


def _tmp_sibling(path: Path) -> Path:
    """Hidden, randomly-suffixed temp path in the same directory as ``path``
    (same filesystem -> ``os.replace`` is atomic)."""
    return path.with_name(f".{path.name}.{generate_random_id()}.tmp")


def _atomic_write_bytes(path, write_fn) -> Path:
    """Run ``write_fn(file_handle)`` on a temp file, then atomically rename it
    onto ``path``. The temp file is removed on any failure."""
    path = Path(path)
    tmp = _tmp_sibling(path)
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
# magnetization series (float64)
# ---------------------------------------------------------------------------
def save_magn(path, magn) -> Path:
    """Atomically write the per-spin magnetization series as ``float64``."""
    arr = np.asarray(magn, dtype=np.float64)
    return _atomic_write_bytes(path, arr.tofile)


def load_magn(path) -> np.ndarray:
    """Load a magnetization series written by :func:`save_magn`."""
    return np.fromfile(path, dtype=np.float64)


# ---------------------------------------------------------------------------
# spin-configuration trajectory (int8, row-major (n_rec, N))
# ---------------------------------------------------------------------------
def open_sout_stream(path):
    """Open a streaming handle to a temp sibling of ``path``.

    Returns ``(file_handle, tmp_path)``; write rows with :func:`write_sout_row`
    and finalize with :func:`close_sout_stream` (which atomically renames the
    temp file onto ``path``). Streaming keeps the full trajectory off-RAM.
    """
    path = Path(path)
    tmp = _tmp_sibling(path)
    return open(tmp, "wb"), tmp


def write_sout_row(fh, s) -> None:
    """Append one int8 configuration row to an open sout stream."""
    np.asarray(s, dtype=np.int8).tofile(fh)


def close_sout_stream(fh, tmp, path) -> Path:
    """Close the stream and atomically rename the temp file onto ``path``."""
    path = Path(path)
    fh.close()
    os.replace(tmp, path)
    return path


def save_sout(path, s_t) -> Path:
    """One-shot atomic dump of an ``(n_rec, N)`` int8 trajectory.

    Used by the backends that return the whole trajectory at once (pybind);
    the Python loop streams row-by-row instead (see ``open_sout_stream``).
    """
    arr = np.ascontiguousarray(np.asarray(s_t, dtype=np.int8))
    return _atomic_write_bytes(path, arr.tofile)


def load_sout(path, N: int, mmap: bool = False) -> np.ndarray:
    """Load a trajectory written row-major as ``(n_rec, N)`` int8.

    With ``mmap=True`` the file is memory-mapped (read-only) and reshaped
    without reading it into RAM -- the right choice for a large ``sout``.
    """
    if mmap:
        return np.memmap(path, dtype=np.int8, mode="r").reshape(-1, N)
    return np.fromfile(path, dtype=np.int8).reshape(-1, N)


# ---------------------------------------------------------------------------
# cluster-size-distribution time series (compressed .npz)
# ---------------------------------------------------------------------------
def save_cldist(path, cluster_dist) -> Path:
    """Atomically write a list of ``{size: count}`` histograms.

    Stored as three flat ``int64`` arrays (record index, size, count) plus the
    total record count, so it round-trips exactly and stays compact for the
    typically-small per-sweep distributions.
    """
    recs: list[int] = []
    sizes: list[int] = []
    counts: list[int] = []
    for t, dist in enumerate(cluster_dist):
        for size, count in dist.items():
            recs.append(t)
            sizes.append(int(size))
            counts.append(int(count))
    buf = io.BytesIO()
    np.savez_compressed(
        buf,
        rec=np.asarray(recs, dtype=np.int64),
        size=np.asarray(sizes, dtype=np.int64),
        count=np.asarray(counts, dtype=np.int64),
        n_rec=np.int64(len(cluster_dist)),
    )
    data = buf.getvalue()
    return _atomic_write_bytes(path, lambda fh: fh.write(data))


def load_cldist(path) -> list[dict[int, int]]:
    """Reconstruct the list of ``{size: count}`` histograms from a .npz file."""
    with np.load(path) as z:
        rec = z["rec"].tolist()
        size = z["size"].tolist()
        count = z["count"].tolist()
        n_rec = int(z["n_rec"])
    out: list[dict[int, int]] = [dict() for _ in range(n_rec)]
    for t, sz, c in zip(rec, size, count):
        out[t][int(sz)] = int(c)
    return out

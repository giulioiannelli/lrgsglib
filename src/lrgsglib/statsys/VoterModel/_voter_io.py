"""Back-compat shim for the VoterModel on-disk observable IO.

The codecs now live in the generic, model-agnostic :mod:`lrgsglib.statsys._io`;
these voter-named wrappers bind them to the historical voter dtypes/format so the
old import surface (and any analysis script that loaded ``m_*``/``sout_*``/
``cldist_*`` files) keeps working byte-for-byte. New code should prefer the
``Observable``/``ObservableSet`` objects in :mod:`lrgsglib.statsys._observables`
(VoterModel itself now drives those) or :mod:`lrgsglib.statsys._io` directly.

On-disk format (unchanged):

- ``m``      : per-spin magnetization series, ``float64`` raw binary.
- ``sout``   : spin-configuration trajectory, ``int8`` raw binary, row-major
               ``(n_rec, N)`` (``N`` not stored; supply it to :func:`load_sout`).
- ``cldist`` : cluster-size-distribution time series -- a list of
               ``{component_size: count}`` histograms -- as a compressed ``.npz``.
               DISTINCT from Ising's ``outcl{i}`` (see defaults.py).
"""

from __future__ import annotations

import numpy as np

from .._io import (  # noqa: F401  (re-exported for back-compat)
    atomic_write_bytes as _atomic_write_bytes,
    close_row_stream,
    load_array,
    load_histogram_series,
    load_rows,
    open_row_stream,
    save_array,
    save_histogram_series,
    save_rows,
    tmp_sibling as _tmp_sibling,
    write_row,
)


# -- magnetization series (float64) -----------------------------------------
def save_magn(path, magn):
    """Atomically write the per-spin magnetization series as ``float64``."""
    return save_array(path, magn, np.float64)


def load_magn(path):
    """Load a magnetization series written by :func:`save_magn`."""
    return load_array(path, np.float64)


# -- spin-configuration trajectory (int8, row-major (n_rec, N)) -------------
def open_sout_stream(path):
    """Open a streaming handle for the int8 trajectory (see :func:`write_sout_row`)."""
    return open_row_stream(path)


def write_sout_row(fh, s):
    """Append one int8 configuration row to an open sout stream."""
    write_row(fh, s, np.int8)


def close_sout_stream(fh, tmp, path):
    """Close the stream and atomically rename the temp file onto ``path``."""
    return close_row_stream(fh, tmp, path)


def save_sout(path, s_t):
    """One-shot atomic dump of an ``(n_rec, N)`` int8 trajectory."""
    return save_rows(path, s_t, np.int8)


def load_sout(path, N: int, mmap: bool = False):
    """Load a trajectory written row-major as ``(n_rec, N)`` int8."""
    return load_rows(path, N, np.int8, mmap=mmap)


# -- cluster-size-distribution time series (compressed .npz) ----------------
def save_cldist(path, cluster_dist):
    """Atomically write a list of ``{size: count}`` histograms."""
    return save_histogram_series(path, cluster_dist)


def load_cldist(path):
    """Reconstruct the list of ``{size: count}`` histograms from a .npz file."""
    return load_histogram_series(path)

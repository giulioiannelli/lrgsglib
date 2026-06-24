"""Phase-4 ContactProcess disk-backed observables (density + snapshots).

Mirrors ``test_observables.py`` (VoterModel): the contact process exposes a
``density`` scalar series and a ``snapshots`` spin trajectory through the shared
:class:`ObservableSet`, kept in RAM (``savedisk=False``) or streamed/dumped under
``dynpath`` and lazily disk-backed (``savedisk=True``).

ContactProcessSIR is used for parity assertions because its Python sampler draws
only from CPython ``np.random`` (reseedable), unlike the EI ``@njit`` sweep whose
numba-global RNG is not reset by ``np.random.seed``.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys.ContactProcess import ContactProcessSIR
from lrgsglib.statsys.ContactProcess.defaults import (
    CP_OBS_DENSITY,
    CP_OBS_SNAPSHOTS,
    CP_SNAPSHOT_MAX_BYTES,
)

pytestmark = pytest.mark.code

_SIDE = 8
_N = _SIDE * _SIDE
_STEPS = 25


def _lattice():
    return Lattice2D(side1=_SIDE, side2=_SIDE, geo="squared", pbc=True,
                     init_nw_dict=False)


def _sir(savedisk, dynpath=None, **kw):
    np.random.seed(123)
    cp = ContactProcessSIR(_lattice(), mu=0.6, runlang="py", ic="homogeneous",
                           state_type="binary", seed=123, steps=_STEPS,
                           savedisk=savedisk, **kw)
    if dynpath is not None:
        cp.dynpath = dynpath
        cp._ensure_dynpath_exists()
    cp.init_contact_dynamics()
    cp.run(tqdm_on=False)
    return cp


def test_observable_set_registered():
    cp = ContactProcessSIR(_lattice(), mu=0.5, runlang="py", state_type="binary")
    assert cp.observables is not None
    assert CP_OBS_DENSITY in cp.observables
    assert CP_OBS_SNAPSHOTS in cp.observables


def test_default_run_records_nothing():
    # savedensity defaults off; savedyn defaults off -> both series empty.
    cp = _sir(savedisk=True)
    assert len(cp.density) == 0
    assert len(cp.s_t) == 0


def test_ram_density_and_snapshots():
    cp = _sir(savedisk=False, savedyn=True, savedensity=True)
    assert isinstance(cp.s_t, list)
    assert np.asarray(cp.s_t).shape == (_STEPS, _N)
    assert len(cp.density) == _STEPS
    # density == per-sweep active-site fraction (mean of each snapshot row)
    assert np.allclose(np.asarray(cp.density), np.asarray(cp.s_t).mean(axis=1))


def test_disk_parity_with_ram(tmp_path):
    ram = _sir(savedisk=False, savedyn=True, savedensity=True)
    dsk = _sir(savedisk=True, dynpath=tmp_path, savedyn=True, savedensity=True)
    # Files land under the stable out_suffix-based names.
    written = {p.name for p in tmp_path.glob("*")}
    assert any(n.startswith("rho") for n in written)
    assert any(n.startswith("sout") for n in written)
    # On disk the trajectory is a read-only memmap; values match the RAM run.
    assert np.asarray(dsk.s_t).shape == (_STEPS, _N)
    assert np.array_equal(np.asarray(ram.s_t), np.asarray(dsk.s_t))
    assert np.allclose(np.asarray(ram.density), np.asarray(dsk.density))


def test_savedisk_false_writes_nothing(tmp_path):
    cp = _sir(savedisk=False, dynpath=tmp_path, savedyn=True, savedensity=True)
    assert list(tmp_path.glob("*")) == []
    # data still available in RAM
    assert len(cp.s_t) == _STEPS


def test_snapshot_subsampling(tmp_path):
    every = 5
    cp = _sir(savedisk=True, dynpath=tmp_path, savedyn=True, snapshot_every=every)
    # Streaming keeps every k-th of the _STEPS recorded sweeps.
    expected = -(-_STEPS // every)  # ceil
    assert np.asarray(cp.s_t).shape == (expected, _N)


def test_snapshot_size_guard(tmp_path):
    # A trajectory exceeding the byte cap must raise before writing.
    cp = ContactProcessSIR(_lattice(), mu=0.5, runlang="py", state_type="binary",
                           seed=1, savedyn=True, savedisk=True)
    cp.dynpath = tmp_path
    cp._ensure_dynpath_exists()
    cp.init_contact_dynamics()
    huge = CP_SNAPSHOT_MAX_BYTES // _N + 2
    with pytest.raises(ValueError, match="CP_SNAPSHOT_MAX_BYTES"):
        cp.run(steps=huge, tqdm_on=False)

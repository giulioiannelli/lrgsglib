"""Phase-3 spine: the model-agnostic frame producer (:mod:`statsys._frames`).

Pins the snapshot-cadence + ``on_frame`` hook + ``make_sweep_fn`` contract the
live-view GUI and the contact-process animation share, independent of any
particular model. Uses ``ContactProcessSIR`` (the simplest binary model, which
inherits the default ``BinDynSys.make_sweep_fn``) as the concrete driver.
"""

import numpy as np
import pytest

from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
from lrgsglib.statsys._frames import collect_frames, iter_frames
from lrgsglib.statsys.ContactProcess import ContactProcessSIR


def _make_cp(steps=4, **kw):
    # Note: seed must be truthy -- DynSys.__init_randomness__ treats seed=0 as
    # "unset" (``seed or <time-based>``), so 0 would give a nondeterministic run.
    lattice = Lattice2D(side1=6, side2=6, geo="squared", pbc=True, init_nw_dict=False)
    cp = ContactProcessSIR(
        lattice, mu=0.4, runlang="py", ic="delta",
        state_type="binary", savedyn=False, seed=1, steps=steps, **kw,
    )
    cp.init_contact_dynamics()
    return cp


@pytest.mark.code
def test_iter_frames_count_includes_t0():
    cp = _make_cp()
    frames = list(iter_frames(cp, steps=6, snapshot_every=1))
    assert len(frames) == 1 + 6  # t0 + one per sweep
    assert all(isinstance(f, np.ndarray) for f in frames)


@pytest.mark.code
def test_iter_frames_snapshot_every():
    cp = _make_cp()
    frames = list(iter_frames(cp, steps=6, snapshot_every=2))
    assert len(frames) == 1 + 6 // 2  # t0 + sweeps 2,4,6


@pytest.mark.code
def test_iter_frames_no_t0():
    cp = _make_cp()
    frames = list(iter_frames(cp, steps=5, snapshot_every=1, include_t0=False))
    assert len(frames) == 5


@pytest.mark.code
def test_iter_frames_yields_independent_copies():
    cp = _make_cp()
    frames = list(iter_frames(cp, steps=3, snapshot_every=1))
    # Mutating one frame must not bleed into another or into the live state.
    frames[0][:] = 7
    assert not np.array_equal(frames[0], frames[1])


@pytest.mark.code
def test_on_frame_hook_indices_and_values():
    cp = _make_cp()
    seen = []
    frames = list(
        iter_frames(
            cp, steps=6, snapshot_every=2,
            on_frame=lambda t, s: seen.append((t, s.copy())),
        )
    )
    assert [t for t, _ in seen] == [0, 2, 4, 6]  # t0 then 1-based sweep counts
    assert len(seen) == len(frames)
    for (_, hooked), produced in zip(seen, frames):
        assert np.array_equal(hooked, produced)


@pytest.mark.code
def test_snapshot_every_must_be_positive():
    cp = _make_cp()
    with pytest.raises(ValueError):
        list(iter_frames(cp, steps=3, snapshot_every=0))


@pytest.mark.code
def test_make_sweep_fn_advances_state_in_place():
    cp = _make_cp()
    sweep = cp.make_sweep_fn()
    before = cp.s.copy()
    for _ in range(20):
        sweep()
    # A supercritical-ish contact process from a single seed should move.
    assert not np.array_equal(before, cp.s)


@pytest.mark.code
def test_collect_frames_matches_iter_frames():
    # Build-then-immediately-run each model: the constructor reseeds the global
    # RNG (seed=1), so each run starts from the same state only if no other run
    # has consumed draws in between.
    a = collect_frames(_make_cp(), steps=4, snapshot_every=1)
    b = list(iter_frames(_make_cp(), steps=4, snapshot_every=1))
    assert len(a) == len(b)
    for x, y in zip(a, b):
        assert np.array_equal(x, y)

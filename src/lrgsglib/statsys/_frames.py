"""Model-agnostic frame producer for dynamics models.

A *frame* is a snapshot of a model's configuration (``model.s``) at a Monte
Carlo sweep boundary. This module turns any :class:`~.DynSys.DynSys.DynSys`
that exposes a one-sweep primitive (:meth:`make_sweep_fn`) into a stream of
frames, optionally feeding each one to a :data:`FrameHook` as it is produced.

This generalizes the contact-process-specific
``ContactProcess.frames.iter_contact_process_frames`` (which now delegates its
in-process Python path here) so the live-view GUI and any model can share one
stepping loop. Two consumers are served by the same generator:

- **offline / batch**: ``collect_frames(model, ...)`` materializes the whole
  list (e.g. to hand to a lattice's ``animate`` accessor for a GIF/MP4);
- **live**: pass ``on_frame=lambda t, s: ...`` to stream each snapshot to a
  canvas as the run steps, without buffering the trajectory.

The per-model "advance one sweep" logic lives on the model (``make_sweep_fn``),
not here: binary models inherit the standard shuffle-and-update sweep from
``BinDynSys``; models with a vectorized kernel (e.g. ``ContactProcessEI``)
override it. This loop only owns the snapshot cadence and the hook.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ._solver import FrameHook

if TYPE_CHECKING:
    from .DynSys.DynSys import DynSys

__all__ = ["FrameHook", "iter_frames", "collect_frames"]

Frame = NDArray[np.int8]


def iter_frames(
    model: "DynSys",
    *,
    steps: int | None = None,
    simref: float | None = None,
    snapshot_every: int = 1,
    include_t0: bool = True,
    on_frame: FrameHook | None = None,
    tqdm_on: bool = False,
) -> Iterator[Frame]:
    """Yield configurations of ``model`` at Monte Carlo sweep boundaries.

    The model is advanced one sweep at a time via ``model.make_sweep_fn()`` and
    its state (``model.s``) is snapshotted every ``snapshot_every`` sweeps. Each
    yielded frame is an independent copy, so buffering the stream is safe.

    Snapshots are intentionally kept *outside* the model's own observables /
    ``s_t`` accumulation -- this is a view facility, not a persisted observable.

    Parameters
    ----------
    model
        A dynamics model exposing ``make_sweep_fn()`` (binary models inherit it
        from ``BinDynSys``). The caller is responsible for having initialized the
        model's state (and any backend caches) before iterating.
    steps, simref
        Run horizon in sweeps (Markov time). ``steps`` takes precedence; when
        either is given the model's time controls are (re)resolved. If neither is
        given the model's current ``steps`` is used.
    snapshot_every
        Emit one snapshot every ``snapshot_every`` sweeps (``1`` => every sweep).
    include_t0
        Emit the initial configuration (sweep ``0``) before any sweep runs.
    on_frame
        Optional ``(sweep_index, state_copy) -> None`` hook invoked for every
        emitted frame (the live-view producer). ``sweep_index`` is ``0`` for the
        ``include_t0`` frame and the 1-based sweep count otherwise.
    tqdm_on
        Wrap the sweep loop in a tqdm progress bar.
    """
    if snapshot_every < 1:
        raise ValueError("snapshot_every must be >= 1.")

    if steps is not None or simref is not None:
        model._set_time_controls(steps=steps, simref=simref)

    sweep = model.make_sweep_fn()

    if include_t0:
        frame = np.asarray(model.s).copy()
        if on_frame is not None:
            on_frame(0, frame)
        yield frame

    if tqdm_on:
        import tqdm as _tqdm

        iterator: Iterable[int] = _tqdm.tqdm(range(model.steps))
    else:
        iterator = range(model.steps)

    for t in iterator:
        sweep()
        if (t + 1) % snapshot_every == 0:
            frame = np.asarray(model.s).copy()
            if on_frame is not None:
                on_frame(t + 1, frame)
            yield frame


def collect_frames(
    model: "DynSys",
    *,
    steps: int | None = None,
    simref: float | None = None,
    snapshot_every: int = 1,
    include_t0: bool = True,
    on_frame: FrameHook | None = None,
    tqdm_on: bool = False,
) -> list[Frame]:
    """Materialize :func:`iter_frames` into a list (the batch/offline path)."""
    return list(
        iter_frames(
            model,
            steps=steps,
            simref=simref,
            snapshot_every=snapshot_every,
            include_t0=include_t0,
            on_frame=on_frame,
            tqdm_on=tqdm_on,
        )
    )

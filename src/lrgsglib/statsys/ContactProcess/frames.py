from __future__ import annotations

"""ContactProcess frame collection (state snapshots of a running model).

Steps a ``ContactProcess`` and snapshots its ``+/-1`` configuration at sweep
boundaries. This is a *dynamics* facility -- it produces the state vectors and is
graph-agnostic (works on any graph the process runs on), so it lives with the
model, not with any graph type. Rendering those frames is a separate, structural
concern: hand them to the lattice's own animation methods (bound from
:mod:`lrgsglib.graphs._shared.animation`).

Typical usage from a notebook:

    from lrgsglib.graphs.nx import Lattice2D
    from lrgsglib.statsys.ContactProcess import ContactProcessEI
    from lrgsglib.statsys.ContactProcess.frames import collect_contact_process_frames

    lattice = Lattice2D(side1=64, geo="squared", pbc=True, init_nw_dict=False)
    cp = ContactProcessEI(lattice, gamma=1.2, runlang="py", ic="delta", state_type="binary", seed=0)
    cp.init_contact_dynamics()
    frames = collect_contact_process_frames(cp, steps=200, snapshot_every=1, backend="py").frames
    anim = lattice.animate.states(frames, model=cp, save="cp.gif")  # structural render
"""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from typing import Callable, Literal, TypeVar

import numpy as np
from numpy.typing import NDArray

from . import ContactProcessBase, ContactProcessEI


Frame = NDArray[np.int8]
TContact = TypeVar("TContact", bound=ContactProcessBase)


@dataclass(frozen=True)
class ContactFrames:
    frames: list[Frame]
    contact_process: ContactProcessBase


def _python_sweep_contact_process(
    nodes: NDArray[np.int64],
    dsNstep: Callable[[NDArray[np.int64]], object],
) -> None:
    np.random.shuffle(nodes)
    dsNstep(nodes)


def _python_sweep_contact_process_ei(cp: ContactProcessEI) -> None:
    if cp._lambda_arr is None:
        raise RuntimeError("ContactProcessEI lambda cache is not initialized.")
    if cp._neigh_indices is None or cp._neigh_weights is None or cp._neigh_offsets is None:
        raise RuntimeError("ContactProcessEI neighbour cache is not initialized.")
    if cp._reverse_weights is None:
        raise RuntimeError("ContactProcessEI reverse-weights cache is not initialized.")
    cp._sweep_ei_lambda_cache(
        cp.s,
        cp._lambda_arr,
        cp._neigh_indices,
        cp._neigh_weights,
        cp._neigh_offsets,
        cp._reverse_weights,
        cp.gamma_eff,
        cp._activation_fn,
        cp.N,
    )


def iter_contact_process_frames(
    cp: TContact,
    *,
    steps: int | None = None,
    simref: float | None = None,
    snapshot_every: int = 1,
    include_t0: bool = True,
    tqdm_on: bool = False,
    backend: Literal["auto", "py", "c_step"] = "auto",
    clean_export: bool = True,
) -> Iterator[Frame]:
    """Yield contact-process configurations at MCMC sweep boundaries.

    This intentionally keeps snapshots *outside* of `ContactProcessBase.s_t`.

    Parameters
    ----------
    cp
        A configured contact process. Call `init_contact_dynamics()` before
        iterating if you want control over the initial condition; otherwise
        this function will call it.
    steps, simref
        Run horizon in sweeps (MCMC time). `steps` takes precedence.
    snapshot_every
        Yield one snapshot every `snapshot_every` sweeps.
    include_t0
        Yield the initial configuration before any sweeps.
    tqdm_on
        Wrap the sweep loop in a tqdm progress bar.
    backend
        - `auto`: uses Python stepping when possible; for C backends falls back
          to `c_step`.
        - `py`: force Python stepping (requires `runlang="py"`).
        - `c_step`: run the C backend one sweep at a time to collect snapshots
          (slow; writes temporary exported files).
    clean_export
        Whether to clean exported files after C stepping.
    """

    if snapshot_every < 1:
        raise ValueError("snapshot_every must be >= 1.")

    # Ensure state is initialized (and for EI/python also initializes caches).
    if getattr(cp, "sini", None) is None:
        cp.init_contact_dynamics()
    cp.initialize_run_parameters(steps=steps, simref=simref)

    if include_t0:
        yield cp.s.copy()

    if tqdm_on:
        import tqdm as _tqdm

        iterator: Iterable[int] = _tqdm.tqdm(range(cp.steps))
    else:
        iterator = range(cp.steps)

    runlang_upper = cp.runlang.upper()
    resolved_backend = backend
    if backend == "auto":
        resolved_backend = "c_step" if runlang_upper.startswith("C") else "py"

    if resolved_backend == "py":
        if runlang_upper != "PY":
            raise ValueError("backend='py' requires cp.runlang='py'.")
        if isinstance(cp, ContactProcessEI):
            sweep = lambda: _python_sweep_contact_process_ei(cp)
        else:
            nodes = np.arange(cp.N, dtype=np.int64)
            dsNstep = cp.dsNstep()
            sweep = lambda: _python_sweep_contact_process(nodes, dsNstep)
        for t in iterator:
            sweep()
            if (t + 1) % snapshot_every == 0:
                yield cp.s.copy()
        return

    if resolved_backend == "c_step":
        if not runlang_upper.startswith("C"):
            raise ValueError("backend='c_step' requires a C runlang (e.g. 'C1c').")
        # Run one sweep at a time, re-exporting the current state.
        # This is slow but avoids needing C-side snapshot support.
        cp.initialize_run_parameters(steps=1, simref=None)
        for t in iterator:
            cp.build_cprogram_command()
            cp.export_s_init()
            cp.run_cprogram(verbose=False)
            if clean_export:
                cp._remove_sfout()
            if (t + 1) % snapshot_every == 0:
                yield cp.s.copy()
        if clean_export:
            cp.remove_run_c_files()
            cp.sg.remove_exported_files()
        return

    raise ValueError(f"Unknown backend: {backend}")


def collect_contact_process_frames(
    cp: TContact,
    *,
    steps: int | None = None,
    simref: float | None = None,
    snapshot_every: int = 1,
    include_t0: bool = True,
    tqdm_on: bool = False,
    backend: Literal["auto", "py", "c_step"] = "auto",
    clean_export: bool = True,
) -> ContactFrames:
    frames = list(
        iter_contact_process_frames(
            cp,
            steps=steps,
            simref=simref,
            snapshot_every=snapshot_every,
            include_t0=include_t0,
            tqdm_on=tqdm_on,
            backend=backend,
            clean_export=clean_export,
        )
    )
    return ContactFrames(frames=frames, contact_process=cp)

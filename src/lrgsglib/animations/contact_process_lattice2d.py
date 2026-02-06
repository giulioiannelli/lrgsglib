from __future__ import annotations

"""Animate ContactProcess dynamics on Lattice2D.

Typical usage from a notebook:

    from lrgsglib.graphs.nx import Lattice2D
    from lrgsglib.statsys.ContactProcess import ContactProcessEI
    from lrgsglib.animations import collect_contact_process_frames, animate_contact_process_on_lattice2d

    lattice = Lattice2D(side1=64, geo="squared", pbc=True, init_nw_dict=False)
    cp = ContactProcessEI(lattice, gamma=1.2, runlang="py", ic="delta", state_type="binary", seed=0)
    cp.init_contact_dynamics()
    frames = collect_contact_process_frames(cp, steps=200, snapshot_every=1, backend="py").frames
    anim = animate_contact_process_on_lattice2d(lattice, frames, out_path="cp.gif", fps=20)
"""

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal, TypeVar

import numpy as np
from numpy.typing import NDArray

from ..statsys.ContactProcess import ContactProcessBase, ContactProcessEI


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


def animate_contact_process_on_lattice2d(
    lattice: object,
    frames: Sequence[Frame],
    *,
    interval_ms: int = 50,
    cmap: str = "viridis",
    add_colorbar: bool = True,
    autoscale: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    blit: bool = False,
    out_path: str | Path | None = None,
    fps: int = 20,
    dpi: int = 150,
    writer: str | None = None,
) -> object:
    """Create (and optionally save) an animation of frames on a Lattice2D."""

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    result = lattice.make_animation(
        fig,
        ax,
        frames,
        interval_ms=interval_ms,
        cmap=cmap,
        add_colorbar=add_colorbar,
        autoscale=autoscale,
        vmin=vmin,
        vmax=vmax,
        blit=blit,
    )

    if out_path is not None:
        from ..graphs.nx.Lattice2DNX._animations import save_animation

        save_animation(result.animation, out_path, fps=fps, dpi=dpi, writer=writer)

    return result.animation

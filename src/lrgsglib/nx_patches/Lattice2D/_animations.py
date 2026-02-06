from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from numpy.typing import NDArray


FrameLike = NDArray[np.integer] | NDArray[np.floating]


@dataclass(frozen=True)
class LatticeAnimationResult:
    animation: object
    image: object
    colorbar: object | None


def _to_2d_frame(
    frame: FrameLike,
    *,
    syshape: tuple[int, int],
    n_expected: int,
) -> NDArray[np.floating]:
    arr = np.asarray(frame)
    if arr.ndim == 1:
        if arr.size != n_expected:
            raise ValueError(
                f"Frame has size {arr.size}, expected {n_expected} (or 2D "
                f"shape {syshape})."
            )
        arr = arr.reshape(syshape)
    elif arr.ndim == 2:
        if tuple(arr.shape) != syshape:
            raise ValueError(f"Frame has shape {arr.shape}, expected {syshape}.")
    else:
        raise ValueError("Frames must be 1D state vectors or 2D arrays.")

    return arr.astype(np.float64, copy=False)


def make_lattice2d_animation(
    lattice: object,
    fig: object,
    ax: object,
    frames: Sequence[FrameLike],
    *,
    interval_ms: int = 50,
    cmap: str = "viridis",
    add_colorbar: bool = True,
    autoscale: bool = False,
    vmin: float | None = None,
    vmax: float | None = None,
    blit: bool = False,
) -> LatticeAnimationResult:
    """Animate 2D lattice configurations with matplotlib.

    Parameters
    ----------
    lattice
        A `Lattice2D` instance (only needs `syshape` and `N` attributes).
    fig, ax
        Matplotlib Figure/Axes to draw into.
    frames
        Sequence of configurations. Each frame can be a 1D vector of length
        `lattice.N` (row-major reshaped into `lattice.syshape`) or a 2D array
        matching `lattice.syshape`.
    interval_ms
        Delay between frames in milliseconds.
    cmap
        Matplotlib colormap name.
    add_colorbar
        Add a colorbar to the axis.
    autoscale
        If True, update color limits per frame. If False, fixed `vmin/vmax`
        is used (falling back to the first frame if not provided).
    vmin, vmax
        Fixed color limits when `autoscale=False`.
    blit
        Forwarded to `matplotlib.animation.FuncAnimation`.
    """

    import matplotlib.animation as mpl_animation

    syshape = tuple(getattr(lattice, "syshape"))
    if len(syshape) != 2:
        raise ValueError("lattice.syshape must be a 2-tuple for Lattice2D.")
    n_expected = int(getattr(lattice, "N"))
    if not frames:
        raise ValueError("frames must be non-empty.")

    first = _to_2d_frame(frames[0], syshape=syshape, n_expected=n_expected)
    if not autoscale:
        if vmin is None:
            vmin = float(np.min(first))
        if vmax is None:
            vmax = float(np.max(first))

    im = ax.imshow(
        first,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal",
    )
    cbar = fig.colorbar(im, ax=ax) if add_colorbar else None

    def _update(i: int):
        arr = _to_2d_frame(frames[i], syshape=syshape, n_expected=n_expected)
        im.set_data(arr)
        if autoscale:
            im.set_clim(float(np.min(arr)), float(np.max(arr)))
        return (im,)

    anim = mpl_animation.FuncAnimation(
        fig,
        _update,
        frames=len(frames),
        interval=interval_ms,
        blit=blit,
    )
    return LatticeAnimationResult(animation=anim, image=im, colorbar=cbar)


def save_animation(
    animation: object,
    out_path: str | Path,
    *,
    fps: int = 20,
    dpi: int = 150,
    writer: str | None = None,
) -> Path:
    """Save a matplotlib animation to disk.

    Notes
    -----
    For `.mp4` outputs this typically requires `ffmpeg` to be installed.
    For `.gif` outputs this typically requires `pillow`.
    """

    path = Path(out_path)
    suffix = path.suffix.lower()
    if writer is None:
        if suffix == ".gif":
            writer = "pillow"
        elif suffix == ".mp4":
            writer = "ffmpeg"
        else:
            raise ValueError("Unknown output extension; pass writer explicitly.")

    save = getattr(animation, "save", None)
    if save is None:
        raise TypeError("animation object does not look like a matplotlib animation.")

    save(str(path), fps=fps, dpi=dpi, writer=writer)
    return path


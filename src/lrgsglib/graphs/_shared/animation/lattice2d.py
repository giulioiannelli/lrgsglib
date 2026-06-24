"""Engine-agnostic 2D-lattice animations.

Both ``Lattice2DNX`` and ``Lattice2DGT`` represent a configuration on a 2D
lattice as a length-``N`` vector that reshapes (row-major) into ``syshape``, so
their state animations are *identical* ``imshow`` movies regardless of backend.
This module is the single home for that shared logic; the two lattice classes
merely bind these functions as methods (mirroring ``graphs/_shared/_draw.py``).

It relies only on accessors present on **both** engines -- ``syshape``, ``N``,
``get_neighbors_with_weights`` -- so the same call works whatever the active
engine::

    lat = Lattice2D(side1=64, geo='sqr', pbc=True, engine='nx')  # or 'gt'
    lat.animate.states(states)                       # inline JS-HTML
    lat.animate.states(states, inline=False)         # raw matplotlib animation
    lat.animate.states(vm, save='run.gif')           # frames from vm.s_t; saved
                                                     # under <plot>/<structure>/voter/

Graphs whose 2D rendering genuinely differs between engines (no shared
``imshow`` representation) should *not* use this module: they bind their own
engine-specific renderers, which can still reuse the back-end-neutral
primitives in :mod:`lrgsglib.plotlib.animation._core`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from ....config.const import GIF
from ....plotlib.animation._core import LatticeAnimationResult, render_animation
from ._accessor import _Accessor, frames_and_model
from ._common import resolve_plot_path


FrameLike = NDArray[np.integer] | NDArray[np.floating]


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


# --------------------------------------------------------------------------- #
# Fast save path: build the whole frame stack as (n, H, W, 3) uint8 once and    #
# stream it to ffmpeg/Pillow, skipping matplotlib's per-frame savefig pipeline  #
# (the lattice is a fixed pixel grid, so re-rendering a Figure per frame is pure #
# overhead -- ~70 s -> ~1 s for 3000 frames). See plotlib.animation.video.      #
# --------------------------------------------------------------------------- #
def _fast_output_ok(save_path, inline: bool) -> bool:
    """Whether the fast RGB->video path applies.

    Needs either a save target or inline display (with ``save=None`` and
    ``inline=False`` the caller wants the raw matplotlib animation object). A
    ``.gif`` target uses Pillow (no ffmpeg); anything else needs an ffmpeg binary.
    """
    if save_path is None and not inline:
        return False
    if save_path is not None and Path(save_path).suffix.lower() == GIF:
        return True
    from ....plotlib.animation.video import ffmpeg_available

    return ffmpeg_available()


def _states_to_rgb(frames, *, syshape, cmap, vmin, vmax) -> NDArray[np.uint8]:
    """Vectorised ``(n, H, W, 3)`` uint8 from state frames via a 256-entry LUT."""
    import matplotlib

    arr = np.asarray(frames, dtype=np.float32)
    arr = arr.reshape(arr.shape[0], syshape[0], syshape[1])
    lo = float(np.min(arr)) if vmin is None else float(vmin)
    hi = float(np.max(arr)) if vmax is None else float(vmax)
    span = (hi - lo) or 1.0
    idx = (np.clip((arr - lo) / span, 0.0, 1.0) * 255).astype(np.uint8)
    lut = (
        matplotlib.colormaps[cmap](np.linspace(0.0, 1.0, 256))[:, :3] * 255
    ).astype(np.uint8)
    return lut[idx]


def _cluster_to_rgb(
    frames, masks, *, syshape, pos_color, neg_color, bg_color
) -> NDArray[np.uint8]:
    """Vectorised ``(n, H, W, 3)`` uint8: largest cluster coloured by spin, rest grey."""
    s = np.asarray(frames, dtype=np.int8).reshape(len(frames), -1)
    n, npix = s.shape
    pos8 = np.asarray(np.asarray(pos_color, np.float64) * 255, np.uint8)
    neg8 = np.asarray(np.asarray(neg_color, np.float64) * 255, np.uint8)
    bg8 = np.asarray(np.asarray(bg_color, np.float64) * 255, np.uint8)
    net = (s * masks).sum(axis=1)                      # net spin of the largest cluster
    col = np.where(net[:, None] > 0, pos8, neg8)        # (n, 3) per-frame cluster colour
    rgb = np.broadcast_to(bg8, (n, npix, 3)).copy()
    rgb = np.where(masks[:, :, None], col[:, None, :], rgb).astype(np.uint8)
    return rgb.reshape(n, syshape[0], syshape[1], 3)


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
        A 2D lattice instance (only needs ``syshape`` and ``N`` attributes;
        works for both ``Lattice2DNX`` and ``Lattice2DGT``).
    fig, ax
        Matplotlib Figure/Axes to draw into.
    frames
        Sequence of configurations. Each frame can be a 1D vector of length
        ``lattice.N`` (row-major reshaped into ``lattice.syshape``) or a 2D
        array matching ``lattice.syshape``.
    interval_ms
        Delay between frames in milliseconds.
    cmap
        Matplotlib colormap name.
    add_colorbar
        Add a colorbar to the axis.
    autoscale
        If True, update color limits per frame. If False, fixed ``vmin/vmax``
        is used (falling back to the first frame if not provided).
    vmin, vmax
        Fixed color limits when ``autoscale=False``.
    blit
        Forwarded to ``matplotlib.animation.FuncAnimation``.
    """

    import matplotlib.animation as mpl_animation

    syshape = tuple(getattr(lattice, "syshape"))
    if len(syshape) != 2:
        raise ValueError("lattice.syshape must be a 2-tuple for a 2D lattice.")
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


def make_lattice2d_cluster_animation(
    lattice: object,
    fig: object,
    ax: object,
    frames: Sequence[FrameLike],
    idx: Sequence[NDArray],
    b: Sequence[NDArray],
    *,
    interval_ms: int = 50,
    pos_color: tuple[float, float, float] = (0.84, 0.19, 0.15),
    neg_color: tuple[float, float, float] = (0.13, 0.40, 0.74),
    bg_color: tuple[float, float, float] = (0.82, 0.82, 0.82),
    blit: bool = False,
) -> LatticeAnimationResult:
    """Animate the largest active-edge cluster in colour, the rest greyed out.

    A *cluster* is a connected component of the **active**-edge subgraph
    (``b_ij * s_i * s_j > 0``); see
    :func:`lrgsglib.utils.statsys.cluster_components`. Each frame colours the
    largest such cluster by its (uniform) spin -- ``pos_color`` for ``+1``,
    ``neg_color`` for ``-1`` -- and greys (``bg_color``) everything else.

    Parameters
    ----------
    lattice
        A 2D lattice instance (only needs ``syshape``).
    fig, ax
        Matplotlib Figure/Axes to draw into.
    frames
        Sequence of 1D state vectors (length ``N``, row-major reshaped into
        ``lattice.syshape``) or 2D arrays matching ``lattice.syshape``.
    idx, b
        Ragged per-node neighbour indices and edge signs (see
        :func:`lrgsglib.utils.statsys.edge_sign_arrays`).
    interval_ms
        Delay between frames in milliseconds.
    pos_color, neg_color, bg_color
        RGB triples in ``[0, 1]`` for the +1 cluster, -1 cluster, and background.
    blit
        Forwarded to ``matplotlib.animation.FuncAnimation``.
    """

    import matplotlib.animation as mpl_animation

    from ....utils.statsys import cluster_components

    syshape = tuple(getattr(lattice, "syshape"))
    if len(syshape) != 2:
        raise ValueError("lattice.syshape must be a 2-tuple for a 2D lattice.")
    if not len(frames):
        raise ValueError("frames must be non-empty.")

    pos = np.asarray(pos_color, dtype=np.float64)
    neg = np.asarray(neg_color, dtype=np.float64)
    bg = np.asarray(bg_color, dtype=np.float64)

    def _rgb(state: FrameLike) -> NDArray[np.floating]:
        s = np.asarray(state, dtype=np.int8).ravel()
        label = cluster_components(s, idx, b)
        largest = label == np.bincount(label).argmax()
        rgb = np.tile(bg, (s.size, 1))
        rgb[largest] = pos if int(s[largest][0]) > 0 else neg
        return rgb.reshape(syshape[0], syshape[1], 3)

    im = ax.imshow(_rgb(frames[0]), interpolation="nearest", aspect="equal")

    def _update(i: int):
        im.set_data(_rgb(frames[i]))
        return (im,)

    anim = mpl_animation.FuncAnimation(
        fig, _update, frames=len(frames), interval=interval_ms, blit=blit
    )
    return LatticeAnimationResult(animation=anim, image=im, colorbar=None)


def make_animation(
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
    """Build a configuration animation into a caller-supplied ``fig``/``ax``.

    Thin convenience over :func:`make_lattice2d_animation` for callers that want
    the raw :class:`LatticeAnimationResult` (animation + image + colorbar) and
    will drive saving/embedding themselves. Bound as a method on both lattice
    engines.
    """
    return make_lattice2d_animation(
        lattice,
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


def animate_states(
    lattice: object,
    states,
    *,
    n_frames: int | None = None,
    fps: int = 12,
    cmap: str = "coolwarm",
    vmin: float | None = -1.0,
    vmax: float | None = 1.0,
    add_colorbar: bool = False,
    figsize: tuple[float, float] = (4.0, 4.0),
    inline: bool = True,
    fast: bool = True,
    save=None,
    dpi: int = 150,
    writer: str | None = None,
    model=None,
    subfolder=None,
):
    """Animation of state configurations on a 2D lattice (engine-agnostic).

    ``states`` is a sequence of 1D state vectors (length ``N``) or 2D arrays
    (``syshape``); with ``n_frames`` the sequence is evenly subsampled to that
    many frames.

    Output is controlled by ``inline``: ``True`` (default) returns an
    ``IPython.display.HTML`` JS animation for inline notebook display;
    ``False`` returns the raw ``matplotlib`` animation object (pure matplotlib).
    Pass ``save="movie.gif"`` / ``"movie.mp4"`` to also write a file -- a bare
    filename lands under the plot tree
    (``<path_plot>/<structure>/<subfolder>``, the subfolder derived from
    ``model`` when given), while an explicit/absolute path is used as given.

    With ``fast=True`` (default) and no colorbar, frames are built as one RGB
    stack and streamed straight to ffmpeg (``.mp4``) / Pillow (``.gif``),
    bypassing matplotlib for a ~50x speed-up; ``fast=False`` forces the
    matplotlib ``FuncAnimation`` path (needed for a colorbar or a raw animation
    object).
    """
    from ....utils.basic import subsample

    frames = subsample(list(states), n_frames) if n_frames else list(states)
    save_path = resolve_plot_path(lattice, save, model=model, subfolder=subfolder)

    if fast and not add_colorbar and _fast_output_ok(save_path, inline):
        from ....plotlib.animation.video import render_rgb_stack

        syshape = tuple(getattr(lattice, "syshape"))
        rgb = _states_to_rgb(frames, syshape=syshape, cmap=cmap, vmin=vmin, vmax=vmax)
        return render_rgb_stack(rgb, save_path, fps=fps, inline=inline)

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    if not add_colorbar:
        # Let the image fill the figure -- no surrounding white border.
        ax.set_position([0.0, 0.0, 1.0, 1.0])
    result = make_lattice2d_animation(
        lattice,
        fig,
        ax,
        frames,
        interval_ms=max(1, round(1000 / fps)),
        cmap=cmap,
        add_colorbar=add_colorbar,
        vmin=vmin,
        vmax=vmax,
    )
    return render_animation(
        result,
        fig,
        fps=fps,
        inline=inline,
        save=save_path,
        dpi=dpi,
        writer=writer,
    )


def animate_largest_cluster(
    lattice: object,
    states,
    *,
    cluster_mode: str = "rawspin",
    n_frames: int | None = None,
    fps: int = 12,
    figsize: tuple[float, float] = (4.0, 4.0),
    inline: bool = True,
    fast: bool = True,
    save=None,
    dpi: int = 150,
    writer: str | None = None,
    pos_color: tuple[float, float, float] = (0.84, 0.19, 0.15),
    neg_color: tuple[float, float, float] = (0.13, 0.40, 0.74),
    bg_color: tuple[float, float, float] = (0.82, 0.82, 0.82),
    masks=None,
    model=None,
    subfolder=None,
):
    """Animation of the largest same-sign/satisfied cluster (engine-agnostic).

    Colours the largest active-edge cluster (see
    :func:`lrgsglib.utils.statsys.cluster_components`) by its spin and greys the
    rest. The active-edge convention is read from this lattice's own signed
    adjacency: ``cluster_mode='rawspin'`` (equal spins) or ``'satisfied'``
    (``sign(w)``-aligned). ``n_frames`` subsamples the trajectory.

    Output and ``save`` behave exactly as in :func:`animate_states`: ``inline``
    toggles JS-HTML vs a raw matplotlib animation, and a bare ``save`` filename
    is written under the plot tree (``<path_plot>/<structure>/<subfolder>``).
    """
    from ....utils.basic import subsample
    from ....utils.statsys import signed_neighbor_arrays

    # Ragged signed neighbour arrays from this lattice's own adjacency
    # (same convention as VoterModel._gillespie_neighbors).
    idx, b = signed_neighbor_arrays(lattice, cluster_mode)

    frames = subsample(list(states), n_frames) if n_frames else list(states)
    cluster_masks = masks
    if cluster_masks is not None and n_frames:
        cluster_masks = subsample(list(cluster_masks), n_frames)
    save_path = resolve_plot_path(lattice, save, model=model, subfolder=subfolder)

    if fast and _fast_output_ok(save_path, inline):
        from ....plotlib.animation.video import render_rgb_stack
        from ....utils.statsys import compute_largest_cluster_masks

        syshape = tuple(getattr(lattice, "syshape"))
        fr = np.asarray(frames, dtype=np.int8)
        # Largest-cluster mask per frame: reuse the caller's buffered ``masks`` if
        # given, else compute them in one process-parallel pass (replaces the
        # matplotlib path's per-frame scipy connected_components).
        cl = (
            np.asarray(cluster_masks, dtype=bool)
            if cluster_masks is not None
            else compute_largest_cluster_masks(fr, idx, b)
        )
        rgb = _cluster_to_rgb(
            fr, cl, syshape=syshape,
            pos_color=pos_color, neg_color=neg_color, bg_color=bg_color,
        )
        return render_rgb_stack(rgb, save_path, fps=fps, inline=inline)

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")
    ax.set_position([0.0, 0.0, 1.0, 1.0])
    result = make_lattice2d_cluster_animation(
        lattice,
        fig,
        ax,
        frames,
        idx,
        b,
        interval_ms=max(1, round(1000 / fps)),
        pos_color=pos_color,
        neg_color=neg_color,
        bg_color=bg_color,
    )
    return render_animation(
        result,
        fig,
        fps=fps,
        inline=inline,
        save=save_path,
        dpi=dpi,
        writer=writer,
    )


# === Accessors (bound as ``lat.plot`` / ``lat.animate`` on both engines) ===


class _Lattice2DAnimate(_Accessor):
    """``lat.animate`` for 2D lattices: imshow state/cluster movies.

    The positional ``source`` may be a dynamics object (e.g. a ``VoterModel`` --
    its ``s_t`` trajectory is used as frames and the figure is filed under that
    model's plot subtree) or an explicit sequence of state vectors.
    """

    def states(self, source, *, model=None, subfolder=None, **kw):
        frames, src_model = frames_and_model(source)
        return animate_states(
            self._sg, frames, model=model or src_model, subfolder=subfolder, **kw
        )

    def largest_cluster(self, source, *, model=None, subfolder=None, **kw):
        frames, src_model = frames_and_model(source)
        return animate_largest_cluster(
            self._sg, frames, model=model or src_model, subfolder=subfolder, **kw
        )

    def make(self, fig, ax, frames, **kw):
        """Low-level builder: returns a ``LatticeAnimationResult`` for a
        caller-supplied ``fig``/``ax`` (no saving)."""
        return make_animation(self._sg, fig, ax, frames, **kw)


class _Lattice2DPlot(_Accessor):
    """``lat.plot`` for 2D lattices: static renderings."""

    def draw(self, **kw):
        return self._sg.draw(**kw)

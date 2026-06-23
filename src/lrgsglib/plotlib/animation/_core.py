"""Engine- and graph-agnostic animation primitives.

These helpers are independent of any graph engine (NetworkX / graph-tool) and of
any particular graph type. They handle the two things every animation needs:

* :func:`save_animation` -- write a finished ``matplotlib`` animation to a
  ``.gif`` / ``.mp4`` file.
* :func:`render_animation` -- turn a built animation into the requested output:
  an inline JS-HTML object for notebooks, **or** the raw ``matplotlib``
  animation object (pure matplotlib, for custom handling / display / saving).

This is the bottom layer of the animation stack and lives in ``plotlib`` because
it knows nothing about graphs. Higher layers build their ``FuncAnimation`` and
defer here, so the "inline vs save vs pure-matplotlib" policy lives in exactly
one place:

* structural, graph-type renderers --
  :mod:`lrgsglib.graphs._shared.animation` (e.g. the engine-agnostic 2D-lattice
  ``imshow`` movies), bound as methods on the graph classes;
* engine-specific renderers for graphs that draw differently on NX vs GT, under
  ``lrgsglib.graphs.<engine>.<Class>.animation``;
* dynamics frame collection -- e.g.
  :mod:`lrgsglib.statsys.ContactProcess.frames`.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class LatticeAnimationResult:
    """Bundle returned by the low-level animation builders.

    Attributes
    ----------
    animation
        The ``matplotlib.animation.FuncAnimation`` instance.
    image
        The ``AxesImage`` (or equivalent artist) being updated each frame.
    colorbar
        The colorbar, when one was added; otherwise ``None``.
    """

    animation: object
    image: object
    colorbar: object | None


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
    For ``.mp4`` outputs this typically requires ``ffmpeg`` to be installed.
    For ``.gif`` outputs this typically requires ``pillow``.
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
        raise TypeError(
            "animation object does not look like a matplotlib animation."
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    save(str(path), fps=fps, dpi=dpi, writer=writer)
    return path


def render_animation(
    result: LatticeAnimationResult,
    fig: object,
    *,
    fps: int,
    inline: bool = True,
    save: str | Path | None = None,
    dpi: int = 150,
    writer: str | None = None,
) -> object:
    """Finalise an animation: optionally save, then return it in the asked form.

    Parameters
    ----------
    result
        The built animation bundle.
    fig
        The figure the animation draws into.
    fps
        Frames per second (used both for saving and for the inline player).
    inline
        If ``True`` (default), return an ``IPython.display.HTML`` JS animation
        for notebook display and close ``fig`` (so the static first frame is
        not rendered a second time below the player). If ``False``, return the
        raw ``matplotlib`` animation object (pure matplotlib) and leave ``fig``
        open so the caller can display, embed, or save it.
    save
        When given, also write the animation to this (already-resolved) path
        (``.gif`` -> pillow, ``.mp4`` -> ffmpeg).
    dpi, writer
        Forwarded to :func:`save_animation`.

    Returns
    -------
    object
        ``IPython.display.HTML`` when ``inline`` is ``True``; otherwise the
        ``matplotlib.animation.FuncAnimation`` instance.
    """

    if save is not None:
        save_animation(result.animation, save, fps=fps, dpi=dpi, writer=writer)

    if not inline:
        return result.animation

    import matplotlib.pyplot as plt
    from IPython.display import HTML

    html = result.animation.to_jshtml(fps=fps)
    plt.close(fig)
    return HTML(html)

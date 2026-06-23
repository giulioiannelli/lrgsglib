"""Generic rasterised-frame animation codec (engine- and graph-agnostic).

For renderers that cannot use ``FuncAnimation`` blitting (e.g. 3D voxels), each
frame is rasterised to an image and the stack is assembled into a looping GIF.
This is the rasterised counterpart of the ``FuncAnimation`` path in
:mod:`lrgsglib.plotlib.animation._core`; both keep the "inline vs save" output
policy out of the graph-specific renderers.
"""
from __future__ import annotations

from pathlib import Path


def figure_to_image(fig):
    """Rasterise the current state of a matplotlib figure to a PIL RGB image."""
    from PIL import Image

    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = fig.canvas.buffer_rgba().tobytes()
    return Image.frombytes("RGBA", (w, h), buf).convert("RGB")


def frames_to_gif_html(
    images,
    *,
    fps: int = 12,
    save: str | Path | None = None,
    loop: int = 0,
    disposal: int = 2,
    optimize: bool = True,
):
    """Assemble PIL images into a looping GIF and return it as inline HTML.

    Parameters
    ----------
    images
        Non-empty sequence of PIL ``Image`` frames (e.g. from
        :func:`figure_to_image`).
    fps
        Frames per second.
    save
        When given, also write the GIF bytes to this (already-resolved) path.
    loop, disposal, optimize
        Forwarded to ``PIL.Image.save`` (the GIF encoder).

    Returns
    -------
    IPython.display.HTML
        An ``<img>`` element with the GIF inlined as base64.
    """
    import base64
    import io

    from IPython.display import HTML

    images = list(images)
    if not images:
        raise ValueError("images must be non-empty.")

    duration = max(1, round(1000 / fps))
    buf = io.BytesIO()
    images[0].save(
        buf,
        format="gif",
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=loop,
        disposal=disposal,
        optimize=optimize,
    )
    data = buf.getvalue()
    if save is not None:
        path = Path(save)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    b64 = base64.b64encode(data).decode()
    return HTML(f'<img src="data:image/gif;base64,{b64}"/>')

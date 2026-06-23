"""Generic rasterised-frame animation codec (engine- and graph-agnostic).

For renderers that cannot use ``FuncAnimation`` blitting (e.g. 3D voxels), each
frame is rasterised to an image once; this module then tight-crops the stack,
optionally writes a looping GIF, and wraps the frames in matplotlib's interactive
JS player (play/pause, frame slider, loop) for inline display. The player is
rebuilt from the *already rasterised* frames via cheap 2D ``imshow`` redraws, so
the expensive 3D rendering is never repeated -- and the controls behave exactly
like ``FuncAnimation.to_jshtml`` everywhere (classic Notebook, JupyterLab, VS
Code), unlike a raw ``<script>`` player.
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


def autocrop_frames(images, *, bg=(255, 255, 255), pad=4):
    """Trim the uniform background border off every frame using ONE bounding box
    -- the union of all frames' content -- so the subject keeps a fixed size and
    position across the movie (no per-frame jitter). 3D axes leave wide margins;
    this removes them after the fact, independent of matplotlib's 3D padding."""
    from PIL import Image, ImageChops

    images = list(images)
    ref = Image.new("RGB", images[0].size, bg)
    box = None
    for im in images:
        b = ImageChops.difference(im.convert("RGB"), ref).getbbox()
        if b is None:
            continue
        box = b if box is None else (
            min(box[0], b[0]), min(box[1], b[1]),
            max(box[2], b[2]), max(box[3], b[3]),
        )
    if box is None:
        return images
    w, h = images[0].size
    l, t, r, btm = box
    crop = (max(0, l - pad), max(0, t - pad), min(w, r + pad), min(h, btm + pad))
    return [im.crop(crop) for im in images]


def frames_to_player_html(images, *, fps=12, save=None, loop=0, disposal=2,
                          optimize=True, crop=True):
    """Tight-crop the frames, optionally write a looping GIF, and return them as
    matplotlib's interactive JS player (play/pause, frame slider, loop) for
    inline display.

    Parameters
    ----------
    images
        Non-empty sequence of PIL ``Image`` frames (e.g. from
        :func:`figure_to_image`).
    fps
        Frames per second (sets the default player speed and the GIF timing).
    save
        When given, also write the GIF bytes to this (already-resolved) path.
    loop, disposal, optimize
        Forwarded to the ``PIL`` GIF encoder.
    crop
        Auto-trim the uniform background border (default ``True``).
    """
    import io

    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as manim
    from IPython.display import HTML

    images = autocrop_frames(images) if crop else list(images)
    if not images:
        raise ValueError("images must be non-empty.")

    if save is not None:
        buf = io.BytesIO()
        images[0].save(
            buf, format="gif", save_all=True, append_images=images[1:],
            duration=max(1, round(1000 / fps)), loop=loop, disposal=disposal,
            optimize=optimize,
        )
        path = Path(save)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(buf.getvalue())

    # Rebuild matplotlib's standard JS player from the rasterised frames: a
    # borderless 2D axis sized to the (cropped) frame, one ``imshow`` per frame.
    w, h = images[0].size
    dpi = 100
    fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_axis_off()
    artists = [[ax.imshow(np.asarray(im), animated=True)] for im in images]
    anim = manim.ArtistAnimation(
        fig, artists, interval=max(1, round(1000 / fps)), blit=True
    )
    html = anim.to_jshtml(fps=fps)
    plt.close(fig)
    return HTML(html)

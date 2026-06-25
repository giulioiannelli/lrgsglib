"""Fast RGB-stack -> video (direct ffmpeg pipe), bypassing matplotlib.

The 2D-lattice animations are just a fixed pixel grid recoloured every frame, so
matplotlib's per-frame ``savefig`` pipeline is pure overhead (~70 s for 3000
small frames). Here an ``(n, H, W, 3)`` uint8 stack is streamed once as raw
``rgb24`` to a single ffmpeg process -- ~1 s for the same input.

Discrete spin colours (red / blue / grey) have near-identical luma, so lossy
H.264 (chroma subsampling + DCT quantisation) smears them into many shades. The
fix is **lossless** ``yuv420p`` (``crf 0``) plus an **even integer upscale**: the
upscale aligns chroma one-per-source-pixel *and* gives a viewable frame, and
``yuv420p`` stays a pixel format every player (including the VS Code webview)
can decode -- so the result is an exact 2-3 colour "chessboard".

Saving ``.gif`` instead uses Pillow (palette, also exact) and needs no ffmpeg.
"""
from __future__ import annotations

import base64
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from ...config.const import (
    ANIM_VIDEO_CRF,
    ANIM_VIDEO_PRESET,
    ANIM_VIDEO_TARGET_PX,
    GIF,
    MP4,
)


def ffmpeg_available() -> bool:
    """True if an ``ffmpeg`` binary is on ``PATH`` (needed for ``.mp4`` output)."""
    return shutil.which("ffmpeg") is not None


def even_upscale(long_edge: int, target_px: int = ANIM_VIDEO_TARGET_PX) -> int:
    """Even integer nearest-neighbour upscale taking ``long_edge`` near ``target_px``.

    Even so the enlarged frame keeps even dimensions (required by ``yuv420p``,
    whatever the lattice parity) and so chroma samples land one-per-source-pixel,
    keeping flat colours crisp.
    """
    up = max(2, round(target_px / max(1, int(long_edge))))
    return up + (up & 1)


def encode_rgb_video(
    frames,
    path,
    *,
    fps: int,
    upscale: int = 1,
    crf: int = ANIM_VIDEO_CRF,
    preset: str = ANIM_VIDEO_PRESET,
) -> Path:
    """Encode an ``(n, H, W, 3)`` uint8 stack to an H.264 ``.mp4`` via one ffmpeg pipe.

    ``crf=0`` (default) is lossless -> exact discrete colours. ``upscale`` is an
    integer nearest-neighbour enlargement applied by ffmpeg (see
    :func:`even_upscale`).
    """
    stack = np.ascontiguousarray(frames, dtype=np.uint8)
    if stack.ndim != 4 or stack.shape[-1] != 3:
        raise ValueError("frames must be an (n, H, W, 3) uint8 array.")
    _, h, w, _ = stack.shape
    cmd = [
        "ffmpeg", "-y", "-loglevel", "error",
        "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{w}x{h}",
        "-r", str(int(fps)), "-i", "-",
    ]
    if upscale > 1:
        cmd += ["-vf", f"scale=iw*{upscale}:ih*{upscale}:flags=neighbor"]
    cmd += [
        "-an", "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-preset", preset, "-crf", str(int(crf)), str(path),
    ]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    proc.stdin.write(stack.tobytes())
    proc.stdin.close()
    if proc.wait() != 0:
        raise RuntimeError("ffmpeg video encode failed (is ffmpeg installed?).")
    return Path(path)


def save_rgb_gif(frames, path, *, fps: int, upscale: int = 1) -> Path:
    """Crisp palette ``.gif`` from an RGB stack (Pillow; no ffmpeg required)."""
    from PIL import Image

    stack = np.ascontiguousarray(frames, dtype=np.uint8)
    if upscale > 1:
        stack = np.repeat(np.repeat(stack, upscale, axis=1), upscale, axis=2)
    imgs = [Image.fromarray(f) for f in stack]
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    imgs[0].save(
        str(path),
        save_all=True,
        append_images=imgs[1:],
        duration=max(1, round(1000 / fps)),
        loop=0,
        disposal=1,
    )
    return Path(path)


def video_to_html(path, *, display_px: int = ANIM_VIDEO_TARGET_PX):
    """Inline a saved ``.mp4`` / ``.gif`` as an autoplaying, looping, pixel-crisp element."""
    path = Path(path)
    data = base64.b64encode(path.read_bytes()).decode()
    style = f"image-rendering:pixelated;width:{int(display_px)}px;height:auto"
    from IPython.display import HTML

    if path.suffix.lower() == GIF:
        return HTML(f'<img src="data:image/gif;base64,{data}" style="{style}">')
    return HTML(
        f'<video controls autoplay loop muted style="{style}">'
        f'<source src="data:video/mp4;base64,{data}" type="video/mp4"></video>'
    )


def render_rgb_stack(
    frames,
    save,
    *,
    fps: int,
    inline: bool,
    target_px: int = ANIM_VIDEO_TARGET_PX,
):
    """Save an ``(n, H, W, 3)`` uint8 stack and return it inline or as a ``Path``.

    Routes by extension: ``.gif`` -> Pillow palette (no ffmpeg), anything else
    -> lossless H.264 ``.mp4``. With ``save=None`` and ``inline=True`` a temp
    ``.mp4`` is used purely to build the inline player.
    """
    stack = np.ascontiguousarray(frames, dtype=np.uint8)
    up = even_upscale(max(stack.shape[1], stack.shape[2]), target_px)
    tmp = None
    if save is None:
        tmp = Path(tempfile.mkstemp(suffix=MP4)[1])
        save = tmp
    save = Path(save)
    if save.suffix.lower() == GIF:
        out = save_rgb_gif(stack, save, fps=fps, upscale=up)
    else:
        out = encode_rgb_video(stack, save, fps=fps, upscale=up)
    if not inline:
        return out
    html = video_to_html(out, display_px=target_px)
    if tmp is not None:
        tmp.unlink(missing_ok=True)
    return html

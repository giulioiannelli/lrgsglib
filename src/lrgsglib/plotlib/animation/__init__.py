"""Engine- and graph-agnostic animation primitives.

The generic, bottom layer of the animation stack: codec/output handling that is
independent of any graph engine or graph type. Graph-specific renderers live
under :mod:`lrgsglib.graphs` (structural) and dynamics frame collection under
:mod:`lrgsglib.statsys`; both reuse the primitives exposed here.
"""
from ._core import LatticeAnimationResult, render_animation, save_animation
from .frames import make_animation_fromFrames
from .raster import autocrop_frames, figure_to_image, frames_to_player_html
from .video import (
    encode_rgb_video,
    even_upscale,
    ffmpeg_available,
    render_rgb_stack,
    save_rgb_gif,
    video_to_html,
)

__all__ = [
    "LatticeAnimationResult",
    "save_animation",
    "render_animation",
    "make_animation_fromFrames",
    "figure_to_image",
    "autocrop_frames",
    "frames_to_player_html",
    "ffmpeg_available",
    "encode_rgb_video",
    "save_rgb_gif",
    "even_upscale",
    "render_rgb_stack",
    "video_to_html",
]

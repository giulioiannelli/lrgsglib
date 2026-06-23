"""Engine- and graph-agnostic animation primitives.

The generic, bottom layer of the animation stack: codec/output handling that is
independent of any graph engine or graph type. Graph-specific renderers live
under :mod:`lrgsglib.graphs` (structural) and dynamics frame collection under
:mod:`lrgsglib.statsys`; both reuse the primitives exposed here.
"""
from ._core import LatticeAnimationResult, render_animation, save_animation
from .frames import make_animation_fromFrames
from .raster import figure_to_image, frames_to_gif_html

__all__ = [
    "LatticeAnimationResult",
    "save_animation",
    "render_animation",
    "make_animation_fromFrames",
    "figure_to_image",
    "frames_to_gif_html",
]

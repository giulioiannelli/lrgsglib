"""Kitchen-sink surface: library-wide chains plus raw third-party names.

Order matters (star-import "last one wins"); it mirrors the historical
single-file ``lrgsglib.notebooks`` module.
"""
from IPython.display import clear_output, display, HTML
# Matplotlib surface used by labs/report notebooks. `plotlib.__all__`
# whitelists only lrgsglib helpers, so the raw matplotlib artists / norms /
# locators imported by the plotlib hub (`const_plotlib`) are re-exported
# here explicitly — notebooks import everything from `lrgsglib.notebooks`.
from ..plotlib.const_plotlib import (
    mpl, animation, gridspec, rc_context, cycler, colormaps,
    Axes, Axes3D, Text, Line2D, PolyCollection,
    Circle, Rectangle, Ellipse, PathPatch, ConnectionPatch, RegularPolygon,
    Polygon,
    GridSpec, blended_transform_factory, AxesDivider,
    ScalarMappable, Colorbar, ColorbarBase,
    Colormap, ListedColormap, LinearSegmentedColormap,
    BoundaryNorm, LightSource, Normalize, LogNorm, SymLogNorm,
    ScalarFormatter, MultipleLocator, SymmetricalLogLocator, LogLocator,
    LogFormatterMathtext, FixedLocator, FuncFormatter,
)
from math import floor   # `ceil` comes from utils.basic.arithmetic
import py3Dmol
import plotly.graph_objects as go


from ..shared import *
from ..core import *
from ..plotlib import *
from ..utils import *
from ..utils.ipy import *

import json   # surfaced to labs/notebooks via the canonical import-* pattern

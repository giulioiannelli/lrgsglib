"""Kitchen-sink surface: library-wide chains plus raw third-party names.

Order matters (star-import "last one wins"); it mirrors the historical
single-file ``lrgsglib.notebooks`` module.
"""

import json  # surfaced to labs/notebooks via the canonical import-* pattern
from math import floor  # `ceil` comes from utils.basic.arithmetic

import plotly.graph_objects as go
import py3Dmol
from IPython.display import HTML, clear_output, display

from ..core import *
from ..plotlib import *

# Matplotlib surface used by labs/report notebooks. `plotlib.__all__`
# whitelists only lrgsglib helpers, so the raw matplotlib artists / norms /
# locators imported by the plotlib hub (`const_plotlib`) are re-exported
# here explicitly — notebooks import everything from `lrgsglib.notebooks`.
from ..plotlib.const_plotlib import (
    Axes,
    Axes3D,
    AxesDivider,
    BoundaryNorm,
    Circle,
    Colorbar,
    ColorbarBase,
    Colormap,
    ConnectionPatch,
    Ellipse,
    FixedLocator,
    FuncFormatter,
    GridSpec,
    LightSource,
    Line2D,
    LinearSegmentedColormap,
    ListedColormap,
    LogFormatterMathtext,
    LogLocator,
    LogNorm,
    MultipleLocator,
    Normalize,
    PathPatch,
    PolyCollection,
    Polygon,
    Rectangle,
    RegularPolygon,
    ScalarFormatter,
    ScalarMappable,
    SymLogNorm,
    SymmetricalLogLocator,
    Text,
    animation,
    blended_transform_factory,
    colormaps,
    cycler,
    gridspec,
    mpl,
    rc_context,
)
from ..shared import *
from ..utils import *
from ..utils.ipy import *

#
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
#
from matplotlib import gridspec, rc_context, cycler, colormaps
from matplotlib.axes import Axes
from matplotlib.cm import hsv, twilight, ScalarMappable
from matplotlib.collections import PolyCollection
from matplotlib.colorbar import ColorbarBase, Colorbar
from matplotlib.colors import Colormap, ListedColormap, BoundaryNorm, \
    LightSource, LinearSegmentedColormap, Normalize, SymLogNorm, to_rgb,\
    rgb2hex, to_hex
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle, Ellipse, PathPatch, ConnectionPatch,RegularPolygon, Polygon
from matplotlib.text import Text
from matplotlib.ticker import ScalarFormatter, MultipleLocator, \
    SymmetricalLogLocator, LogLocator,LogFormatterMathtext
from matplotlib.transforms import blended_transform_factory
#
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import AxesDivider
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset, zoomed_inset_axes
from mpl_toolkits.mplot3d import Axes3D
#
from ..utils.basic import linspace
from ..config.const import ColorType
#
twilight_lim_low = 0.2
twilight_lim_high = 0.8
twilight_lim_blu = 0.65
cred, cblu = twilight(twilight_lim_low), twilight(twilight_lim_blu)
restr_twilight_vals = twilight(
    linspace(twilight_lim_low, twilight_lim_high)
)
restr_twilight = LinearSegmentedColormap.from_list(
    "restr_twilight", restr_twilight_vals
)
credcblu = ListedColormap([cred, cblu])
mpl.colormaps.register(name="restr_twilight", cmap=restr_twilight)
mpl.colormaps.register(name='restr_twilight_r', cmap=restr_twilight.reversed())

red_blue = LinearSegmentedColormap.from_list("red_blue", ["red", "blue"])
mpl.colormaps.register(name="red_blue", cmap=red_blue)
mpl.colormaps.register(name="red_blue_r", cmap=red_blue.reversed())

PLT_SL2DSQ_SIDE1 = 7
PLT_SL2DSQ_SIDE2 = 7
PLT_SL2DSQ_UNIDS = 1.
PLT_SL2DSQ_LNEXTL = .75
PLT_SL2DSQ_KWNODE = dict(marker='o', ms=20, mec='k', mfc='w')
PLT_SL2DSQ_KWEXTL = dict(marker='', ls=':', zorder=0)
PLT_SL2DSQ_KWLINE = dict(lw=3)
PLT_SL2DSQ_PEC = cblu
PLT_SL2DSQ_CPEC = cred
PLT_SL2DSQ_MODE = 'rand'
PLT_SL2DSQ_KWTXT = dict(fontsize=24, c='k', ha='center', va='center')
PLT_SL2DSQ_VDCPEC = [(0, 2), (2, 2), (3, 2), (1, 1), (1, 0)]
PLT_SL2DSQ_HDCPEC = [(3, 2), (3, 1), (1, 0), (2, 0)]
PLT_SL2DSQ_PFLIP = 0.2
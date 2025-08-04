import lmfit
import random
#
import numpy as np
#
from typing import Dict
#
from ..utils.basic import symmetric_logarithm_unchecked, symmetric_logarithm,\
    flip_to_positive_majority, unravel_1d_to_2d_nodemap
from .const_plotlib import *
from .color import set_alpha_torgb

__all__ = [
    'scheme_Lattice2DSquared',
]
#
def defects_on_lattice_plot(sizes, lattices, ax, direction: str = 'parallel', 
                            geometry: str = 'squared', cell: str = 'single', 
                            fit_mode: str = 'lmfit'):
    #
    from scipy.optimize import curve_fit
    #
    newLowerBound = None
    xShiftConst = 0
    kwlogfit = dict(marker='', c='red')# label=r'$a\log(x) + b$'
    if direction == 'perpendicular':
        ylabel = r'${\phi(\bar{x}_\perp,\, y)}/{\phi_{\min}}$'
        ax.set_xlabel(r'$y$')
    else:
        ylabel = r'${\phi(x,\, \bar{y}_\parallel)}/{\phi_{\min}}$'
        ax.set_xlabel(r'$x$')
    if geometry == 'squared':
        if cell == 'single':
            if direction == 'parallel':
                xShiftConst = 0
                newLowerBound = .5
                slice_cut = lambda side: np.s_[:, lattices[side].side2//2]
            else:
                xShiftConst = -.5
                newLowerBound = 0
                slice_cut = lambda side: np.s_[lattices[side].side1//2, :]
        if cell == 'singleZERR':
            if direction == 'parallel':
                xShiftConst = 0
                newLowerBound = -.5
                slice_cut = lambda side: np.s_[:, lattices[side].side2//2-1]
            else:
                xShiftConst = 1
                newLowerBound = -.5
                slice_cut = lambda side: np.s_[lattices[side].side1//2, :]

        if cell == 'singleXERR':
            if direction == 'parallel':
                xShiftConst = -.5
                slice_cut = lambda side: np.s_[:, lattices[side].side2//2-1]
            else:
                xShiftConst = +.5
                slice_cut = lambda side: np.s_[lattices[side].side1//2, :]
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV = np.max(eigV) - eigV
                eigV /= np.max(eigV)
                return eigV
        else:
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV = np.max(eigV) - eigV
                eigV /= np.max(eigV)
                return eigV
    elif geometry == 'triangular':
        xShiftConst = 0
        if cell != 'singleXERR':
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV /= np.min(eigV)
                return eigV
        else:
            def func(eigV):
                eigV = flip_to_positive_majority(eigV)
                eigV /= np.max(eigV)
                return eigV
        if direction == 'parallel':
            xShiftConst = -1
            if cell == 'singleZERR':
                xShiftConst = 0
            elif cell == 'singleXERR':
                xShiftConst = -.5
            slice_cut = lambda side: np.s_[lattices[side].side1//2, :]
            if cell != 'singleXERR':
                newLowerBound = .5 
        else:
            xShiftConst = -.5
            slice_cut = lambda side: np.s_[:, lattices[side].side2//2]

    ax.set_ylabel(ylabel, labelpad=10)
    ax.set_xscale('symlog')
    ax.set_yscale('log')
    #
    cmapv = restr_twilight(np.linspace(0, 1, len(sizes)))
    lista = []
    for side, c in zip(sizes[::-1], cmapv):
        kwdict = {"ls": '-',
                'c': c, 
                'marker': 'H', 
                'ms': 10, 
                'mfc': set_alpha_torgb(c, 0.75), 
                'label': fr"$N={side**2}$"} 
        eigV = unravel_1d_to_2d_nodemap(lattices[side].eigV[0], lattices[side].invnode_map, lattices[side].syshape)
        eigen_state = func(eigV)
        phi_plot = eigen_state[slice_cut(side)]
        if side == sizes[-1]:
            phi_plot0 = phi_plot
        # print(np.min(eigen_state), min(phi_plot), slice_cut(side))
        if geometry == 'squared':
            sideop = side
        elif geometry == 'triangular':
            if direction == 'parallel':
                sideop = lattices[side].side2 
            elif direction == 'perpendicular':
                sideop = lattices[side].side1
        x = np.rint(np.linspace(-sideop//2, sideop//2, num=sideop)+xShiftConst)
        ax.plot(x, phi_plot, **kwdict)
        ax.plot(np.argmax(np.abs(x))+2, phi_plot[np.argmax(np.abs(x))]+2, 'or')
        lista.append(x[np.argmax(np.abs(x))+2]+sideop//2)
        print(np.argmax(np.abs(x)), x[np.argmax(np.abs(x)):np.argmax(np.abs(x))+2], )
    #
    if fit_mode and cell != 'singleXERR':
        # idx = (x < -side//3) | (x > side//3)
        x_values = np.linspace(-sizes[-1]//2, sizes[-1]//2, num=sizes[-1])+xShiftConst
        y_values = phi_plot0
        
        x_plot = np.concatenate(
            [
                np.linspace(-sizes[-1]//2, -1, num=100),
                np.linspace(-1, 1, num=2000),
                np.linspace(1, sizes[-1]//2, num=100)
            ]
        )
        if fit_mode == 'lmfit':
            true_params = dict(a=1, b=2, c=0.5, d=0.1)
            # Create a Model with the function and fit to the data
            regressor = lmfit.Model(symmetric_logarithm, nan_policy='omit')
            params = regressor.make_params(a=1, b=1.5, c=1, d=0)
            params['b'].set(min=1e-10)  # Prevent b from being zero or negative
            params['d'].set(min=-np.inf, max=np.inf) 
            result = regressor.fit(y_values, params, x=x_values)
            # y_fit = regressor.eval(params=result.params, x=x_plot)
            popt = np.array([result.params[key].value for key in result.params])
        elif fit_mode == 'scipy':
            weights = np.abs(x)
            popt, _ = curve_fit(symmetric_logarithm, x, phi_plot0, sigma=1/weights, 
                                p0=[1, 1.5, 1, 0])
        y_fit = symmetric_logarithm_unchecked(x_plot, *popt)
        ax.plot(x_plot, y_fit, **kwlogfit)
    # set_new_lower_ybound(ax, newLowerBound)
    #
    kwvlines = {'ls':'--', 'lw':1, 'c':'k'}
    for i in [-1, +1, -lattices[sizes[0]].r_c, +lattices[sizes[0]].r_c]:
        ax.axvline(i, **kwvlines)
    return lista
#
def scheme_Lattice2DSquared(*args, **kwargs) -> None:
    """
    Alias for the function scheme_Lattice2DSquared_vXX.
    """
    return scheme_Lattice2DSquared_v02(*args, **kwargs)
#
def scheme_Lattice2DSquared_v01(
    ax: Axes,
    size: int = 7,
    kwargs_nodes: Dict = None,
    kwargs_extl: Dict = None,
    etxl_len: float = 0.75,
    kwargs_lines: Dict = None,
    pec: str = "blue",
    cpec: str = "red",
) -> None:
    """
    Plot a 2D square lattice with customizable nodes, lines, and external 
    lines, each with their own styling parameters.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object on which the lattice will be plotted.
    size : int, optional
        The size of the lattice (number of nodes in one dimension). Default 
        is 7.
    kwargs_nodes : dict, optional
        Keyword arguments for styling the nodes in the plot. Includes 
        'marker' (shape of the marker), 'ms' (marker size), 'mec' (marker 
        edge color), and 'mfc' (marker face color). Default is a white-filled 
        circle.
    kwargs_extl : dict, optional
        Keyword arguments for styling the external (boundary) lines of the 
        lattice. Includes 'ls' (line style). Default is dotted lines.
    etxl_len : float, optional
        Length of the external lines extending from the lattice boundaries. 
        Default is 0.75.
    kwargs_lines : dict, optional
        Keyword arguments for styling the internal lines of the lattice. 
        Includes 'lw' (line width). Default is a line width of 3.
    pec : str, optional
        Primary edge color for the lattice lines. Default is 'blue'.
    cpec : str, optional
        Complementary edge color for the lattice lines. Default is 'red'.

    Returns
    -------
    None
        This function does not return any value. It modifies the provided 
        axes object in place.

    Notes
    -----
    - The function randomly assigns the primary or complementary edge color 
      to each line in the lattice. This includes both the internal lines and 
      the external boundary lines.
    - The nodes are plotted over the lines for a clear visualization.

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> scheme_Lattice2DSquared(ax, size=5, pec='green', cpec='magenta')
    >>> plt.show()

    This will plot a 5x5 lattice with green and magenta lines, default node 
    style, and default external line style.
    """
    import numpy as np
    import random
    from matplotlib import rc_context

    kwargs_nodes = kwargs_nodes or {"marker": "o", "ms": 20, "mec": "k", "mfc": "w"}
    kwargs_extl = kwargs_extl or {"ls": ":"}
    kwargs_lines = kwargs_lines or {"lw": 3}

    x, y = np.meshgrid(range(size), range(size))

    # Plot each point in the lattice
    for i in range(size):
        for j in range(size):
            kwargs_lines["color"] = random.choice([pec, cpec])
            if kwargs_lines["color"] == pec:
                if i < size - 1:
                    ax.plot(
                        [x[i, j], x[i + 1, j]],
                        [y[i, j], y[i + 1, j]],
                        zorder=1,
                        **kwargs_lines,
                    )  # Vertical
                if j < size - 1:
                    ax.plot(
                        [x[i, j], x[i, j + 1]],
                        [y[i, j], y[i, j + 1]],
                        zorder=1,
                        **kwargs_lines,
                    )  # Horizontal
            else:
                with rc_context({"path.sketch": (5, 15, 1)}):
                    if i < size - 1:
                        ax.plot(
                            [x[i, j], x[i + 1, j]],
                            [y[i, j], y[i + 1, j]],
                            zorder=1,
                            **kwargs_lines,
                        )  # Vertical
                    if j < size - 1:
                        ax.plot(
                            [x[i, j], x[i, j + 1]],
                            [y[i, j], y[i, j + 1]],
                            zorder=1,
                            **kwargs_lines,
                        )  # Horizontal
            ax.plot(x[i, j], y[i, j], zorder=2, **kwargs_nodes)  # Nodes

    # Adding dashed lines on the boundaries
    for i in range(size):
        kwargs_extl["color"] = random.choice([pec, cpec])
        # Left and right boundaries
        ax.plot(
            [x[i, 0], x[i, 0] - etxl_len],
            [y[i, 0], y[i, 0]],
            zorder=0,
            **kwargs_extl
        )
        ax.plot(
            [x[i, -1], x[i, -1] + etxl_len],
            [y[i, -1], y[i, -1]],
            zorder=0,
            **kwargs_extl
        )

        # Top and bottom boundaries
        ax.plot(
            [x[0, i], x[0, i]],
            [y[0, i], y[0, i] - etxl_len],
            zorder=0,
            **kwargs_extl
        )
        ax.plot(
            [x[-1, i], x[-1, i]],
            [y[-1, i], y[-1, i] + etxl_len],
            zorder=0,
            **kwargs_extl
        )

    # Remove axes
    ax.axis("off")

def scheme_Lattice2DSquared_v02(
        ax: Axes, 
        side1: int = PLT_SL2DSQ_SIDE1,
        side2: int = PLT_SL2DSQ_SIDE2,
        mode: str = PLT_SL2DSQ_MODE,
        kwargNodes: dict = PLT_SL2DSQ_KWNODE,
        kwargsExtl: dict = PLT_SL2DSQ_KWEXTL,
        lenExtl: float = PLT_SL2DSQ_LNEXTL,
        kwargsLines: dict = PLT_SL2DSQ_KWLINE,
        pec: ColorType = PLT_SL2DSQ_PEC,
        cpec: ColorType = PLT_SL2DSQ_CPEC,
        pflip: float = PLT_SL2DSQ_PFLIP,
        kwargsTxt: dict = PLT_SL2DSQ_KWTXT
        ) -> None:
    """
        Function to plot a square lattice where the color and style of each link
        can be controlled individually. This includes the ability to create defects
        or randomly alter the appearance of links and nodes within the lattice.

        Parameters
        ----------
        ax : Axes
            The matplotlib axes to plot on.
        side1 : int, optional
            Number of nodes on one side of the square lattice.
        side2 : int, optional
            Number of nodes on the other side of the square lattice.
        mode : str, optional
            The mode of coloring the lattice links. Can be 'rand' for random 
            colors or 'defects' to specify certain links with different colors.
        kwargNodes : dict, optional
            Keyword arguments for plotting the nodes.
        kwargsExtl : dict, optional
            Keyword arguments for plotting the extended links.
        lenExtl : float, optional
            Length of the extended links.
        kwargsLines : dict, optional
            Keyword arguments for plotting the lattice links.
        pec : ColorType, optional
            Primary color for the lattice links.
        cpec : ColorType, optional
            Color for the specified or random links based on the mode.
        pflip : float, optional
            Probability of flipping the color of a link in 'rand' mode.
        kwargsTxt : dict, optional
            Keyword arguments for plotting text labels on the lattice.

        Returns
        -------
        None
            This function does not return a value. It modifies the given Axes object
            in-place, adding a visual representation of a 2D squared lattice.

        Notes
        -----
        The 'defects' mode allows for specifying particular links to have a different
        appearance (e.g., color) to represent defects or special cases in the lattice.
        The 'rand' mode uses the `pflip` parameter to randomly change the appearance of
        links, simulating randomness or disorder within the lattice.

        This function is designed to be flexible, with many aspects of the visualization
        customizable through keyword arguments. This allows for a wide range of visual
        styles and representations to suit different requirements or preferences.
    """
    if side1 != PLT_SL2DSQ_SIDE1 and side2 == PLT_SL2DSQ_SIDE1:
        side2 = side1
    if kwargsExtl == PLT_SL2DSQ_KWEXTL:
        kwargsExtl.update(dict(c=pec))
    x, y = np.meshgrid(np.arange(side1), np.arange(side2), indexing='ij')
    #
    if mode == 'defects':
        def determine_line_color(i, j, direction):
            if direction == "vertical":
                if (i, j) in PLT_SL2DSQ_VDCPEC:
                    return cpec
                return pec
            # Conditions for horizontal lines
            elif direction == "horizontal":
                if (i, j) in PLT_SL2DSQ_HDCPEC:
                    return cpec
                return pec
            else: 
                return pec
    elif mode == 'rand':
        def determine_line_color(*args):
            return random.choices([pec, cpec], [1-pflip, pflip], k=1)[0]
    #
    for i in range(side1):
        for j in range(side2):
            # Vertical lines
            if i < side1 - 1:
                xPt = [x[i, j], x[i+1, j]]
                yPt = [y[i, j], y[i+1, j]]
                lc = determine_line_color(i, j, "vertical")
                ax.plot(xPt, yPt, zorder=1, **kwargsLines, color=lc)
                # plot_line_with_style(ax, xPt, yPt, lc, kwargsLines, cpec)
            # Horizontal lines
            if j < side2 - 1:
                xPt = [x[i, j], x[i, j+1]]
                yPt = [y[i, j], y[i, j+1]]
                lc = determine_line_color(i, j, "horizontal")
                # plot_line_with_style(ax, xPt, yPt, lc, kwargsLines, cpec)
            # Nodes
            ax.plot(x[i, j], y[i, j], zorder=2, **kwargNodes)
    
    for j in range(side2):
        # Extend left from the left nodes
        ax.plot([x[0, j], x[0, j]- lenExtl], [y[0, j], y[0, j]], **kwargsExtl)
        # Extend right from the right nodes
        ax.plot([x[-1, j], x[-1, j]+ lenExtl], 
                [y[-1, j], y[-1, j]], **kwargsExtl)

    for i in range(side1):
        # Extend downwards from the bottom nodes
        ax.plot([x[i, 0], x[i, 0]], [y[i, 0], y[i, 0] - lenExtl], **kwargsExtl)
        # Extend upwards from the top nodes
        ax.plot([x[i, -1], x[i, -1]], 
                [y[i, -1], y[i, -1] + lenExtl], **kwargsExtl)


    if mode == 'defects':
        ax.text(x[0, 2] + .5, y[0, 2] + 0.3, rf'$S$', color=cpec, **kwargsTxt)
        ax.text(x[3, 2] + .5, y[3, 2] + 0.3, r"$X$", color=cpec, **kwargsTxt)
        ax.text(x[2, 0] + .3, y[2, 0] + 0.5, r"$Z$", color=cpec, **kwargsTxt)


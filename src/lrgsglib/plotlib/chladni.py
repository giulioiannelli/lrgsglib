from .const_plotlib import *
from .ax_patches import add_external_border_axis
from .tilings import plot_hex_tiling_from_pos, plot_honeycomb_grid_fast, plot_octagonal_square_tiling_from_pos

from numpy.typing import NDArray

__all__ = [
    'make_plot_chladni_states',
]

def make_plot_chladni_states(
        ax: list[Axes], 
        state_i: NDArray, 
        state_f: NDArray, 
        geo: str = 'square',
        pos: NDArray | None = None,
        linewidth: float = 1.5, 
        color: ColorType = 'black', 
        cmap: Colormap = 'binary',
        shape: tuple[int, int] | None = None,
        **kwargs
) -> None:
    """
    Plot initial and final state arrays with a common colormap and border.

    Parameters
    ----------
    ax : list[Axes]
        List containing two matplotlib axes objects for plotting the initial
        and final states.
    state_i : NDArray
        2D array representing the initial state to be plotted.
    state_f : NDArray
        2D array representing the final state to be plotted.
    linewidth : float, optional
        Width of the border line. Default is 1.5.
    color : ColorType, optional
        Color of the border. Default is 'black'.
    cmap : Colormap, optional
        Colormap to use for both plots. Default is 'binary'.
    **kwargs
        Additional keyword arguments passed to `imshow`.

    Returns
    -------
    None
        This function modifies the provided axes in place.

    Notes
    -----
    The function plots the initial state on ax[0] and the final state on ax[1],
    adds a border to both, and turns off the axis.
    """
    if geo == 'square' or geo == 'sqr' or geo == 'hexagonal' or geo == 'hex':
        state_i = state_i.reshape(shape)
        state_f = state_f.reshape(shape)
    if geo == 'square' or geo == 'sqr':
        ax[0].imshow(state_i, cmap=cmap, **kwargs)
        ax[1].imshow(state_f, cmap=cmap, **kwargs)
    elif geo == 'hexagonal' or geo == 'hex':
        plot_honeycomb_grid_fast(
            state_i,
            ax[0],
            cmap=cmap,
            **kwargs
        )
        plot_honeycomb_grid_fast(
            state_f,
            ax[1],
            cmap=cmap,
            **kwargs
        )
    elif geo == 'triangular' or geo == 'tri':
        plot_hex_tiling_from_pos(
            pos, 
            state_i, 
            ax[0], 
            cmap=cmap, 
            **kwargs
        )
        plot_hex_tiling_from_pos(
            pos, 
            state_f, 
            ax[1], 
            cmap=cmap, 
            **kwargs
        )
    elif geo == 'oct_sqr' or geo == 'octagonal_square':
        plot_octagonal_square_tiling_from_pos(
            pos, 
            state_i, 
            ax[0], 
            cmap=cmap, 
            **kwargs
        )
        plot_octagonal_square_tiling_from_pos(
            pos, 
            state_f, 
            ax[1], 
            cmap=cmap, 
            **kwargs
        )
    else:
        raise ValueError(f"Unsupported geometry: {geo}. Supported geometries \
                         are 'square', 'triangular', and 'hexagonal'.")
    for axx in ax[:2]:
        add_external_border_axis(
            axx, 
            color=color, 
            linewidth=linewidth, 
            axis_off=True, 
            fill_opt=False, 
        )
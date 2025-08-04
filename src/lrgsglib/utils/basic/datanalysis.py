import numpy as np
#
from typing import Tuple
#
from numpy.typing import NDArray
from scipy.interpolate import griddata

__all__ = [
    "interpolate_grid_data"
]

def interpolate_grid_data(
    x: NDArray, 
    y: NDArray, 
    z: NDArray, 
    num_points: int = 1000,
    method: str = 'cubic',
    fill_value: float = np.nan
) -> Tuple[NDArray, NDArray, NDArray]:
    """
    Performs interpolation of z-data over a uniformly spaced grid defined by x and y coordinates.

    Parameters
    ----------
    x : NDArray
        Meshgrid of x-coordinate values.
    y : NDArray
        Meshgrid of y-coordinate values.
    z : NDArray
        Matrix of computed z-values corresponding to each (x, y) pair.
    num_points : int, optional
        Number of points along each axis for the new grid. Defaults to 1000.
    method : str, optional
        Interpolation method to be used. Options are 'linear', 'nearest', and 'cubic'. Defaults to 'cubic'.
    fill_value : float, optional
        Value used to fill in for requested points outside of the convex hull of the input points. Defaults to np.nan.

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray]
        grid_x : NDArray
            Meshgrid of interpolated x-coordinate values.
        grid_y : NDArray
            Meshgrid of interpolated y-coordinate values.
        z_new : NDArray
            Interpolated z-values on the new grid.
    """
    points = np.column_stack((x.ravel(), y.ravel()))
    grid_x, grid_y = np.meshgrid(
        np.linspace(x.min(), x.max(), num_points),
        np.linspace(y.min(), y.max(), num_points)
    )
    z_new = griddata(
        points, z.ravel(), (grid_x, grid_y), method=method, fill_value=fill_value
    )
    return grid_x, grid_y, z_new
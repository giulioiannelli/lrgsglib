from typing import Any, Tuple
#
from .const_plotlib import *
#
__all__ = [
    'imshow_colorbar_caxdivider',
]
#
def imshow_colorbar_caxdivider(
    mappable: Any,
    ax: Axes,
    position: str = "right",
    size: str = "5%",
    pad: float = 0.05,
    orientation: str = "vertical",
    axis_dim: int = 2,
    **kwargs
) -> Tuple[AxesDivider, Axes, Any]:
    """
    Display a colorbar in a specified position and orientation relative to a 
    given axis. Supports both 2D and 3D axes.

    Parameters
    ----------
    mappable : ScalarMappable
        The image or plot to which the colorbar applies.
    ax : Axes
        The axis where the image is displayed.
    position : {'left', 'right', 'top', 'bottom'}, optional
        Position of the colorbar relative to `ax`. Default is "right".
    size : str or float, optional
        Size of the colorbar relative to `ax` (e.g., "5%") or in points. 
        Default is "5%".
    pad : float, optional
        Padding between `ax` and the colorbar in relative units. Default is 
        0.05.
    orientation : {'vertical', 'horizontal'}, optional
        Orientation of the colorbar. Default is "vertical".
    axis_dim : {2, 3}, optional
        Dimension of the axis. Use 2 for 2D axes and 3 for 3D axes. Default is 2.
    **kwargs : dict, optional
        Additional keyword arguments passed to `plt.colorbar`.

    Returns
    -------
    divider : AxesDivider
        Divider used to position the colorbar.
    cax : Axes
        Axis containing the colorbar.
    clb : Colorbar
        The created Colorbar object.

    Notes
    -----
    - This function uses `mpl_toolkits.axes_grid1.make_axes_locatable` to 
      create a divider for the axis and append a new axis for the colorbar.
    - The `size` parameter can be specified as a percentage (e.g., "5%") or 
      as an absolute size in points.
    - The `pad` parameter controls the spacing between the axis and the 
      colorbar.
    - For 3D axes, a 2D axis is appended to host the colorbar to avoid 
      projection errors.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from mpl_toolkits.axes_grid1 import make_axes_locatable
    >>> fig, ax = plt.subplots()
    >>> data = np.random.rand(10, 10)
    >>> im = ax.imshow(data, cmap="viridis")
    >>> imshow_colorbar_caxdivider(im, ax, position="right", size="5%", pad=0.1)
    >>> plt.show()

    For 3D axes:
    >>> from mpl_toolkits.mplot3d import Axes3D
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> scatter = ax.scatter([1, 2, 3], [4, 5, 6], [7, 8, 9], c=[10, 20, 30])
    >>> imshow_colorbar_caxdivider(scatter, ax, axis_dim=3)
    >>> plt.show()
    """
    import matplotlib.axes as _maxes

    divider = make_axes_locatable(ax)
    if axis_dim == 3:
        cax = divider.append_axes(
            position, size=size, pad=pad, axes_class=_maxes.Axes
        )
        clb = ax.figure.colorbar(
            mappable, cax=cax, orientation=orientation, **kwargs
        )
    else:
        cax = divider.append_axes(position, size=size, pad=pad)
        clb = plt.colorbar(
            mappable=mappable, cax=cax, orientation=orientation, **kwargs
        )

    return divider, cax, clb
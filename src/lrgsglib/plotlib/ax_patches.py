from typing import Any, Optional, Dict
#
from .const_plotlib import *
#
__all__ = [
    'add_external_border_axis',
    'add_rectangle_patch_to_axis',
]
#
def add_external_border_axis(
    axx: Axes, 
    color: str = 'black', 
    linewidth: float = 2.0, 
    axis_off: bool = False, 
    fill_opt: bool = False, 
    tranform: Optional[Any] = None, 
    **more_kwargs: Dict[str, Any]
) -> None:
    """
    Add an external border to the given axis.

    Parameters
    ----------
    axx : Axes
        The axis to which the border will be added.
    color : str, optional
        The color of the border. Default is 'black'.
    linewidth : float, optional
        The width of the border line. Default is 2.0.
    axis_off : bool, optional
        If True, the axis will be turned off. Default is False.
    fill_opt : bool, optional
        If True, the border will be filled. Default is False.
    tranform : Optional[Any], optional
        The transformation to apply to the border. If None, defaults to 
        `axx.transAxes`. Default is None.
    **more_kwargs : Dict[str, Any]
        Additional keyword arguments to pass to the Rectangle.

    Returns
    -------
    None
    """
    add_rectangle_patch_to_axis(axx, (0, 0, 1, 1), color, linewidth, axis_off, 
                                fill_opt, tranform, **more_kwargs)

def add_rectangle_patch_to_axis(
    axx: Axes,
    rect_coords: tuple[float, float, float, float],
    color: str = 'black',
    linewidth: float = 2.0,
    axis_off: bool = False,
    fill_opt: bool = False,
    tranform: Optional[Any] = None,
    **more_kwargs: Dict[str, Any]
) -> None:
    """
    Add a rectangle patch to the given axis.

    Parameters
    ----------
    axx : Axes
        The axis to which the rectangle will be added.
    rect_coords : tuple[float, float, float, float]
        The coordinates of the rectangle as (x, y, width, height).
    color : str, optional
        The color of the rectangle. Default is 'black'.
    linewidth : float, optional
        The width of the rectangle border line. Default is 2.0.
    axis_off : bool, optional
        If True, the axis will be turned off. Default is False.
    fill_opt : bool, optional
        If True, the rectangle will be filled. Default is False.
    tranform : Optional[Any], optional
        The transformation to apply to the rectangle. If None, defaults to 
        `axx.transAxes`. Default is None.
    **more_kwargs : Dict[str, Any]
        Additional keyword arguments to pass to the Rectangle.

    Returns
    -------
    None
    """
    rect_kwargs = {
        'transform': axx.transAxes if tranform is None else tranform,
        'fill': fill_opt,
        'edgecolor': color,
        'linewidth': linewidth,
    }
    rectangle = Rectangle(rect_coords[:2], rect_coords[2], rect_coords[3], **rect_kwargs, **more_kwargs)
    axx.add_patch(rectangle)
    if axis_off:
        axx.axis('off')
#
def set_ax_ratio_1_withlim(
    ax: Axes,
    margin: float = 0.0,
    maintain_center: bool = True
) -> None:
    """
    Adjust the axis limits to ensure a 1:1 aspect ratio, optionally adding 
    a margin and maintaining the center of the plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to adjust.
    margin : float, optional
        Additional margin to add around the plot, as a fraction of the 
        largest range. Default is 0.0 (no margin).
    maintain_center : bool, optional
        If True, the plot will remain centered around its original center. 
        If False, the limits will be adjusted to fit the data tightly. 
        Default is True.

    Returns
    -------
    None
        This function modifies the axis limits in place.

    Notes
    -----
    - The function calculates the largest range (x or y) and adjusts the 
      limits to ensure a square plot.
    - If `margin` is provided, it is added to the largest range before 
      setting the limits.
    - If `maintain_center` is False, the plot will not be centered, and the 
      limits will be adjusted to fit the data tightly.

    Examples
    --------
    >>> fig, ax = plt.subplots()
    >>> ax.plot([0, 1], [0, 2])
    >>> set_ax_ratio_1_withlim(ax, margin=0.1)
    >>> plt.show()
    """
    # Calculate the ranges and centers
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    x_center = (x_max + x_min) / 2
    y_center = (y_max + y_min) / 2
    x_range = x_max - x_min
    y_range = y_max - y_min

    # Determine the largest range and add margin
    max_range = max(x_range, y_range) / 2
    max_range *= (1 + margin)

    if maintain_center:
        # Set the new limits to ensure the plot is square and centered
        ax.set_xlim(x_center - max_range, x_center + max_range)
        ax.set_ylim(y_center - max_range, y_center + max_range)
    else:
        # Adjust limits to fit the data tightly
        ax.set_xlim(x_min, x_min + 2 * max_range)
        ax.set_ylim(y_min, y_min + 2 * max_range)
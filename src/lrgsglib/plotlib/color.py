from .const_plotlib import *
#
from typing import Optional, Sequence
#
__all__ = [
    'convert_to_RGB',
    'get_opposite_color',
    'set_alpha_tocolor',
    'set_alpha_torgb',
    'set_color_cycle',
]
#

def convert_to_RGB(color: ColorType) -> tuple[int, int, int]:
    """
    Convert a given color to an RGB tuple of three integers.

    Parameters
    ----------
    color : ColorType
        The input color, which can be a string (e.g., 'red', '#FF0000'),
        a tuple of three integers (R, G, B) in the range 0-255, or a tuple
        of three floats (R, G, B) in the range 0.0-1.0.

    Returns
    -------
    tuple[int, int, int]
        The RGB representation of the color as a tuple of three integers
        in the range 0-255.

    Examples
    --------
    >>> convert_to_rgb('red')
    (255, 0, 0)
    >>> convert_to_rgb((0.5, 0.5, 0.5))
    (128, 128, 128)
    >>> convert_to_rgb((255, 0, 0))
    (255, 0, 0)
    """
    if isinstance(color, str):
        rgb = to_rgb(color)
        return tuple(int(c * 255) for c in rgb)
    elif isinstance(color, tuple) and len(color) == 3:
        if all(isinstance(c, int) for c in color):
            return color
        elif all(isinstance(c, float) for c in color):
            return tuple(int(c * 255) for c in color)
    raise ValueError(
        "Invalid color format. Must be a string, a tuple of three integers, "
        "or a tuple of three floats."
    )
#
def convert_to_rgb(color: ColorType) -> tuple[float, float, float]:
    """
    Convert a given color to an RGB tuple of three floats in the range [0, 1].

    Parameters
    ----------
    color : ColorType
        The input color, which can be a string (e.g., 'red', '#FF0000'),
        a tuple of three integers (R, G, B) in the range 0-255, or a tuple
        of three floats (R, G, B) in the range 0.0-1.0.

    Returns
    -------
    tuple[float, float, float]
        The RGB representation of the color as a tuple of three floats
        in the range [0.0, 1.0].

    Examples
    --------
    >>> convert_to_rgb('red')
    (1.0, 0.0, 0.0)
    >>> convert_to_rgb((128, 128, 128))
    (0.5019607843137255, 0.5019607843137255, 0.5019607843137255)
    >>> convert_to_rgb((0.5, 0.5, 0.5))
    (0.5, 0.5, 0.5)
    """
    rgb = convert_to_RGB(color)
    return tuple(c / 255 for c in rgb)
#

def get_opposite_color(
    color: ColorType, 
    col_type: str = 'hex'
) -> ColorType:
    """
    Calculate the opposite color for a given color in various formats.

    Parameters
    ----------
    color : ColorType
        The input color, which can be:
        - A hexadecimal string (e.g., '#6496c8')
        - A named color (e.g., 'red')
        - An RGB tuple with values in the range [0, 1] (e.g., (0.1, 0.5, 0.8))
        - An RGB tuple with values in the range [0, 255] (e.g., (100, 150, 200))
        - An RGBA tuple with values in the range [0, 1] or [0, 255] (e.g.,
          (0.1, 0.5, 0.8, 1.0))
    col_type : str, optional
        The desired output format of the opposite color. Options are:
        - 'hex': Returns the color as a hexadecimal string (default)
        - 'rgb': Returns the color as an RGB tuple with values in [0, 1]

    Returns
    -------
    ColorType
        The opposite color represented in the specified format. If the input
        color includes an alpha channel, it will be preserved.

    Examples
    --------
    >>> get_opposite_color((0.1, 0.5, 0.8))
    '#e66a33'
    >>> get_opposite_color('#6496c8', col_type='rgb')
    (0.6078431372549019, 0.4117647058823529, 0.21568627450980393)
    >>> get_opposite_color((0.1, 0.5, 0.8, 0.7))
    (0.9, 0.5, 0.2, 0.7)
    """
    from numpy import array

    # Convert the input color to an RGB(A) tuple (values between 0 and 1)
    try:
        if isinstance(color, tuple):
            if len(color) == 3:
                if max(color) > 1:  # Assuming the values are in [0, 255]
                    rgb = tuple(c / 255 for c in color)
                else:
                    rgb = color
                alpha = None
            elif len(color) == 4:  # Preserve alpha channel for RGBA
                if max(color) > 1:  # Assuming the values are in [0, 255]
                    rgb = tuple(c / 255 for c in color[:3])
                    alpha = color[3]
                else:
                    rgb = color[:3]
                    alpha = color[3]
            else:
                raise ValueError(
                    "Unsupported tuple length for color. Must be 3 (RGB) or "
                    "4 (RGBA)."
                )
        else:
            rgb = to_rgb(color)
            alpha = None
    except ValueError:
        raise ValueError(
            "Unsupported color format. Please provide a valid color."
        )

    # Convert RGB to the opposite color
    rgb_array = array(rgb)
    opposite_rgb = 1 - rgb_array  # Find the opposite for each channel

    # Return the opposite color in the desired output format
    if col_type == 'hex':
        return (
            to_hex(opposite_rgb)
            if alpha is None
            else to_hex(opposite_rgb) + f"{int(alpha * 255):02x}"
        )
    elif col_type == 'rgb':
        return (
            tuple(opposite_rgb)
            if alpha is None
            else (*opposite_rgb, alpha)
        )
    else:
        raise ValueError("Unsupported output format. Please use 'hex' or 'rgb'.")
#

def set_alpha_tocolor(
    color: ColorType, 
    alpha: float = 0.5,
) -> tuple[int, int, int, float]:
    """
    Convert a color to RGB, then set the alpha (transparency) channel.

    Parameters
    ----------
    color : ColorType
        The input color, which can be a string (e.g., 'red', '#FF0000'),
        a tuple of three integers (R, G, B) in the range 0-255, or a tuple
        of three floats (R, G, B) in the range 0.0-1.0.
    alpha : float, optional
        The alpha (transparency) value to set for the color. Must be
        between 0.0 (fully transparent) and 1.0 (fully opaque). Default is 0.5.

    Returns
    -------
    tuple[int, int, int, float]
        A new RGBA color tuple in the format (R, G, B, alpha).

    Examples
    --------
    >>> set_alpha_tocolor('red', 0.2)
    (255, 0, 0, 0.2)
    >>> set_alpha_tocolor((0.5, 0.5, 0.5), 0.8)
    (128, 128, 128, 0.8)
    """
    rgb = convert_to_rgb(color)
    return set_alpha_torgb(rgb, alpha)
#

def set_alpha_torgb(
    rgbcol: tuple[int, int, int], 
    alpha: float = 0.5
) -> tuple[int, int, int, float]:
    """
    Set the alpha (transparency) channel of an RGB color tuple.

    Parameters
    ----------
    rgbcol : tuple[int, int, int]
        A tuple representing an RGB color in the format (R, G, B), where
        R, G, B are integers between 0 and 255.
    alpha : float, optional
        The alpha (transparency) value to set for the color. Must be
        between 0.0 (fully transparent) and 1.0 (fully opaque). Default is 0.5.

    Returns
    -------
    tuple[int, int, int, float]
        A new RGBA color tuple in the format (R, G, B, alpha).

    Examples
    --------
    >>> set_alpha_torgb((255, 0, 0), 0.2)
    (255, 0, 0, 0.2)
    """
    return (rgbcol[0], rgbcol[1], rgbcol[2], alpha)
#
def set_color_cycle(
    arr: Sequence,
    ax: Optional[Axes] = None,
    my_cmap: Optional[Colormap] = None
) -> None:
    """
    Sets the color cycle of the given Axes based on the provided colormap and 
    array length.

    Parameters
    ----------
    arr : Sequence
        An array-like object that determines the number of distinct colors 
        needed. The length of `arr` dictates how many colors will be 
        generated from the colormap.
    ax : matplotlib.axes.Axes, optional
        The Matplotlib Axes object to which the color cycle will be applied. 
        This is where the color properties will be set, affecting subsequent 
        plot elements added to this Axes. If not provided, the current axes 
        (`plt.gca()`) will be used.
    my_cmap : matplotlib.colors.Colormap, optional
        A Matplotlib colormap instance used to map normalized values to 
        colors. This colormap defines the range and variation of colors in 
        the cycle. If not provided, defaults to the custom colormap 
        `'restr_twilight'`.

    Returns
    -------
    None
        This function does not return any value. It modifies the `ax` object 
        in place by setting its color property cycle.

    Raises
    ------
    ValueError
        If `arr` is empty, as at least one color is required to set the color 
        cycle.
    LookupError
        If the default colormap `'restr_twilight'` is not found and `my_cmap` 
        is not provided.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> from matplotlib import cm
    >>> arr = [1, 2, 3, 4, 5]
    >>> x = np.linspace(0, 10, 100)
    >>> fig, ax = plt.subplots()
    >>> set_color_cycle(arr, ax=ax)
    >>> for i, val in enumerate(arr):
    >>>     y = np.sin(x + i)
    >>>     ax.plot(x, y, label=f'Line {i+1}')
    >>> ax.legend()
    >>> plt.show()

    Notes
    -----
    - Ensure that the length of `arr` corresponds to the number of distinct 
      elements you plan to plot to avoid color repetition.
    - The colormap (`my_cmap`) can be any Matplotlib colormap. You can create 
      custom colormaps or use predefined ones like `'viridis'`, `'plasma'`, 
      `'inferno'`, `'magma'`, etc.
    - If using the default `'restr_twilight'` colormap, ensure it is 
      registered with Matplotlib using `plt.register_cmap()` before calling 
      this function.
    """
    if not arr:
        raise ValueError(
            "The input array 'arr' must contain at least one element to set "
            "the color cycle."
        )

    # Use the current axes if none is provided
    ax = ax or plt.gca()

    # Use the default colormap if none is provided
    if my_cmap is None:
        try:
            my_cmap = plt.get_cmap('restr_twilight')
        except ValueError:
            raise LookupError(
                "The default colormap 'restr_twilight' is not found. Please "
                "provide a valid colormap."
            )

    # Generate evenly spaced values between 0 and 1 based on the length of arr
    normalized_values = linspace(0.0, 1.0, len(arr))

    # Map the normalized values to colors using the provided colormap
    colors = my_cmap(normalized_values)

    # Apply the color cycle to the Axes
    ax.set_prop_cycle(cycler(color=colors))
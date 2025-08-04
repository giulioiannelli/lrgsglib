from .const_plotlib import *
from .color import convert_to_RGB, to_hex
#
__all__ = [
    'create_custom_colormap',
]
#
def create_custom_colormap(
    c1: ColorType = "#0000ff", 
    c2: ColorType = "#fc0303", 
    nc: int = 0
) -> LinearSegmentedColormap:
    """
    Create a custom colormap transitioning between two specified colors.

    Parameters
    ----------
    c1 : ColorType, optional
        The starting color, which can be a string (e.g., 'red', '#FF0000'),
        a tuple of three integers (R, G, B) in the range 0-255, or a tuple
        of three floats (R, G, B) in the range 0.0-1.0. Default is "#0000ff".
    c2 : ColorType, optional
        The ending color, which can be a string (e.g., 'red', '#FF0000'),
        a tuple of three integers (R, G, B) in the range 0-255, or a tuple
        of three floats (R, G, B) in the range 0.0-1.0. Default is "#fc0303".
    nc : int, optional
        The number of discrete colors in the colormap. If 0, the colormap
        will be continuous. Default is 0.

    Returns
    -------
    LinearSegmentedColormap
        A custom colormap transitioning from `c1` to `c2`.

    Examples
    --------
    >>> custom_cmap = create_custom_colormap(c1="red", c2="blue", nc=10)
    >>> custom_cmap = create_custom_colormap(c1=(255, 0, 0), c2=(0, 0, 255))
    """
    start_color = to_hex(convert_to_RGB(c1))
    end_color = to_hex(convert_to_RGB(c2))
    colors = [start_color, end_color]
    nocol = dict(N=nc) if nc else dict()

    cmap = LinearSegmentedColormap.from_list(
        "custom_colormap", colors, **nocol
    )
    return cmap
#

def generate_maxpercdiff_colormap(
    number_of_distinct_colors: int = 80, 
    number_of_shades: int = 7
) -> ListedColormap:
    """
    Generate a perceptually distinct colormap using a saw-tooth pattern in 
    the HSV color space.

    Parameters
    ----------
    number_of_distinct_colors : int, optional
        The number of distinct colors in the colormap. Default is 80.
    number_of_shades : int, optional
        The number of shades for each distinct color. Default is 7.

    Returns
    -------
    ListedColormap
        A ListedColormap object representing the generated colormap.

    Notes
    -----
    - The function generates a perceptually distinct colormap using a 
      saw-tooth pattern in the HSV color space. The number of distinct 
      colors can be customized by providing the `number_of_distinct_colors` 
      parameter.
    - The HSV colormap is cyclic, which is leveraged to create cyclic color 
      variations.
    - The lower half of the colormap is modified to make colors towards the 
      beginning darker, while the upper half is adjusted to make colors 
      towards the end brighter.

    Examples
    --------
    >>> colormap = generate_maxpercdiff_colormap(100)
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> data = np.random.rand(10, 10)
    >>> plt.imshow(data, cmap=colormap)
    >>> plt.colorbar()
    >>> plt.show()
    """
    import numpy as np

    # Calculate the total number of colors, ensuring it is a multiple of the 
    # number of shades
    total_colors = int(
        np.ceil(number_of_distinct_colors / number_of_shades) * number_of_shades
    )

    # Create an array with uniformly distributed floats in the range [0, 1)
    linear_nums = np.arange(total_colors) / total_colors

    # Reshape the array into a matrix with `number_of_shades` rows
    shade_matrix = linear_nums.reshape(
        number_of_shades, total_colors // number_of_shades
    )

    # Transpose the matrix to create a saw-tooth pattern
    sawtooth_matrix = shade_matrix.T

    # Flatten the matrix into a single array
    sawtooth_array = sawtooth_matrix.reshape(-1)

    # Generate the initial HSV colormap using the saw-tooth pattern
    initial_colormap = hsv(sawtooth_array)

    # Modify the lower half of the colormap to make colors darker
    lower_half = total_colors // 2
    for i in range(3):
        initial_colormap[:lower_half, i] *= np.linspace(0.2, 1, lower_half)

    # Modify the upper half of the colormap to make colors brighter
    for i in range(3):
        for j in range(
            total_colors // number_of_shades - lower_half // number_of_shades
        ):
            start = lower_half + j * number_of_shades
            end = start + number_of_shades
            modifier = (
                (1 - initial_colormap[start:end, i]) * 
                (j / (total_colors // number_of_shades - 
                      lower_half // number_of_shades))
            )
            initial_colormap[start:end, i] += modifier

    return ListedColormap(initial_colormap)
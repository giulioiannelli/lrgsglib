from numpy.typing import NDArray
from typing import Dict, Tuple, Iterable

__all__ = [
    "zoom_into_array",
    "shift_with_wrap",
    "unravel_1d_to_2d_nodemap"
]

def zoom_into_array(array: NDArray, x: int, y: int, width: int, height: int) -> NDArray:
    """
    Extract a subarray from a 2D array based on a given coordinate and dimensions.

    Parameters
    ----------
    array : NDArray
        The input 2D array.
    x : int
        The x-coordinate of the center of the zoom.
    y : int
        The y-coordinate of the center of the zoom.
    width : int
        The width of the zoomed area.
    height : int
        The height of the zoomed area.

    Returns
    -------
    NDArray
        The zoomed 2D subarray.

    Notes
    -----
    - If the zoomed area exceeds the boundaries of the array, it will be clipped.
    - The function assumes `array` is a 2D NumPy array.

    Examples
    --------
    >>> array = np.arange(100).reshape(10, 10)
    >>> zoom_into_array(array, x=5, y=5, width=3, height=3)
    array([[44, 45, 46],
           [54, 55, 56],
           [64, 65, 66]])
    """
    x_start = max(0, x - width // 2)
    x_end = min(array.shape[0], x + width // 2 + 1)
    y_start = max(0, y - height // 2)
    y_end = min(array.shape[1], y + height // 2 + 1)
    return array[x_start:x_end, y_start:y_end]

def shift_with_wrap(image: NDArray, shift_right: int, shift_down: int) -> NDArray:
    """
    Shift a 2D image with wrap-around at the edges.

    Parameters
    ----------
    image : NDArray
        A 2D NumPy array representing the image to be shifted.
    shift_right : int
        Number of pixels to shift the image to the right.
    shift_down : int
        Number of pixels to shift the image downward.

    Returns
    -------
    NDArray
        The shifted image with wrap-around applied at the edges.

    Notes
    -----
    - The function ensures that the shifts are within the bounds of the image dimensions by using modulo operations.
    - The shifting is performed using `np.roll` for efficient computation.

    Examples
    --------
    >>> import numpy as np
    >>> image = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    >>> shift_with_wrap(image, shift_right=1, shift_down=1)
    array([[9, 7, 8],
           [3, 1, 2],
           [6, 4, 5]])
    """
    from numpy import roll
    # Ensure shifts are within the bounds of the image dimensions
    shift_right %= image.shape[1]
    shift_down %= image.shape[0]

    # Perform the shift
    shifted_image = roll(image, shift=shift_down, axis=0)  # Shift down
    shifted_image = roll(shifted_image, shift=shift_right, axis=1)  # Then shift right

    return shifted_image

def unravel_1d_to_2d_nodemap(arr1d: NDArray, imap: Dict[int, Tuple[int, int]], dims: Tuple[int, int] = None) -> NDArray:
    """
    Transforms a 1D array into a 2D array based on a given index mapping and optional dimensions.

    Parameters:
    ---------------
    arr1d : NDArray
        The 1D input numpy array to be transformed.
    imap : Dict[int, Tuple[int, int]]
        A dictionary where keys are indices in the 1D array and values are (row, col) positions in the 2D output.
    dims : Tuple[int, int], optional
        The dimensions of the 2D output array. If not provided, it's assumed to be a square array.

    Returns:
    ---------------
    NDArray
        The 2D array obtained from rearranging `arr1d` according to `imap`.
    """
    from numpy import empty, sqrt
    if not dims:
        side = int(sqrt(len(arr1d)))
        dims = (side, side)
    arr_2d = empty(dims, dtype=arr1d.dtype)
    for idx_1d, (row, col) in imap.items():
        arr_2d[row, col] = arr1d[idx_1d]  # Adjusted indexing to [row, col] for clarity
    return arr_2d

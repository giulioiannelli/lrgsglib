from .common import *
#
def adjust_to_even(x: float) -> int:
    """
    Round a real number to the nearest even integer.

    If the value is exactly halfway between two even integers, the function
    rounds up to the larger even integer.

    Parameters
    ----------
    x : float
        The number to round.

    Returns
    -------
    int
        The nearest even integer to `x`.

    Examples
    --------
    >>> import numpy as np
    >>> adjust_to_even(128 * np.sqrt(3))
    222
    >>> adjust_to_even(5.5)
    6
    >>> adjust_to_even(2.1)
    2
    """
    # Compute the largest even integer ≤ x
    lower_even = int(x) - (int(x) % 2)
    # The next even integer above lower_even
    upper_even = lower_even + 2
    # Choose the closer one; if exactly halfway, pick upper_even
    return lower_even if (x - lower_even) < (upper_even - x) else upper_even

# review the commenting and typing
def sign_with_threshold(arr: np.ndarray, threshold: float = 1e-17) -> np.ndarray:
    """
    Apply the sign function to an array with a threshold for zero.

    Parameters
    ----------
    arr : np.ndarray
        The input array to which the sign function will be applied.
    threshold : float, optional
        Values with absolute value below this threshold are set to 0. Default is 1e-17.

    Returns
    -------
    np.ndarray
        An array where each entry is the sign of the input, with values below the threshold set to 0.

    Notes
    -----
    - The function modifies the input array in-place by setting values below the threshold to 0.
    - The sign of each element is then computed using `np.sign`.

    Examples
    --------
    >>> import numpy as np
    >>> arr = np.array([-0.5, 0, 0.5, 1e-18])
    >>> sign_with_threshold(arr)
    array([-1.,  0.,  1.,  0.])
    """
    # Create a mask for values below the threshold
    mask = np.abs(arr) < threshold
    # Apply the mask to set these values to 0
    arr[mask] = 0
    # Use np.sign on the modified array
    return np.sign(arr)


def bin_sign(arr: Iterable) -> NDArray:
    """
    Regularizes and binarizes an array by setting all zeros to +1 and taking 
    the sign of each element.
    
    Parameters
    ----------
    arr : Iterable
        Input iterable.

    Returns
    -------
    NDArray
        Binarized array with -1, +1 values.

    Examples
    --------
    >>> import numpy as np
    >>> arr = [-3, 0, 2, -1, 0]
    >>> bin_sign(arr)
    array([-1,  1,  1, -1,  1])

    Notes
    -----
    - This function replaces zero elements with +1 and converts all non-zero 
      elements to their respective signs.
    - It can be used to transform continuous data into a binary representation.
    """
    arr = np.asarray(arr)
    return np.sign(np.where(arr == 0, 1, arr))

def flip_to_positive_majority(arr):
    """ 
    Flips the elements of an array to ensure a majority of positive components.

    Given a numerical array, this function checks if the majority of its 
    components are negative. If so, it multiplies every element by -1, 
    effectively flipping the signs of all components to ensure a majority of 
    positive values. This operation is intended for arrays where elements can be
    distinctly categorized as positive or negative (including zero as a 
    non-negative value).

    Parameters:
    -----------
    arr : array_like
        The input array containing numerical data. This array can be a list, 
        tuple, or any array-like object convertible to a NumPy array. The 
        function is optimized for NumPy arrays for performance reasons.

    Returns:
    --------
    numpy.ndarray
        An array of the same shape and type as the input, with elements flipped
        if the original array had a majority of negative components. If the 
        input array already had a majority of positive components, it is 
        returned unchanged.

    Notes:
    ------
    The function utilizes NumPy for efficient computation, especially for large
    arrays. The determination of "majority" is based purely on the count of 
    positive vs. negative elements, without weighting by magnitude.

    Examples:
    ---------
    >>> import numpy as np
    >>> arr = np.array([1, -2, -3, 4, -5])
    >>> flip_to_positive_majority(arr)
    array([ -1,  2,  3,  -4,  5])

    >>> arr = np.array([-1, 2, 3])
    >>> flip_to_positive_majority(arr)
    array([-1,  2,  3])
    """
    if np.sum(arr < 0) > len(arr) / 2:
        # Flip all elements by multiplying by -1
        arr = arr * -1
    return arr

def flip_to_positive_majority_adapted(arr: NDArray) -> NDArray:
    """
    Flips the elements of an array to ensure a majority of positive components.
    
    This function assesses the balance of positive and negative values in the 
    given numerical array. If the count of negative values exceeds the count of
    positive values, all elements of the array are multiplied by -1. This 
    operation ensures that the majority of elements in the transformed array
    are positive. The function is particularly useful in data processing
    scenarios where the sign of data points affects subsequent analysis or
    visualization.

    Parameters:
    -----------
    arr : ndarray
        The input array containing numerical data. This array can be a list, 
        tuple, or any array-like object convertible to a NumPy array. The 
        function is optimized for NumPy arrays for performance reasons.

    Returns:
    --------
    ndarray
        An array of the same shape and type as the input, with elements flipped
        if the original array had a majority of negative components. If the 
        input array already had a majority of positive components, it is 
        returned unchanged.

    Examples:
    ---------
    >>> import numpy as np
    >>> arr = np.array([1, -2, -3, 4, -5])
    >>> flip_to_positive_majority_adapted(arr)
    array([-1,  2,  3, -4, 5])

    >>> arr = np.array([-1, 2, 3])
    >>> flip_to_positive_majority_adapted(arr)
    array([-1, 2, 3])
    """
    num_negatives = np.sum(arr < 0)
    num_positives = np.sum(arr > 0)

    if num_negatives > num_positives:
        arr = arr * -1
    return arr

def is_in_range(number, range_start, range_end):
    """
    Checks if a given number is within a specified range, inclusive.

    Parameters:
    -----------
    number : int or float
        The number to check if it lies within the given range.
    range_start : int or float
        The starting value of the range (inclusive).
    range_end : int or float
        The ending value of the range (inclusive).

    Returns:
    --------
    bool
        Returns True if the number is within the range specified by
        range_start and range_end, inclusive; otherwise, False.

    Notes:
    ------
    This function uses Python's ability to chain comparison operators,
    making the check concise and efficient. It is versatile enough to handle
    both integer and floating-point numbers.

    Example:
    --------
    >>> is_in_range(5, 1, 10)
    True
    >>> is_in_range(15, 1, 10)
    False
    """
    return range_start <= number <= range_end
#
def ceil(x: float) -> int:
    """
    Compute the ceiling of a real number.

    This function returns the smallest integer greater than or equal to the
    given number. If the number is already an integer, it is returned as is.

    Parameters
    ----------
    x : float
        The number for which the ceiling value is to be computed.

    Returns
    -------
    int
        The ceiling value of the input number.

    Examples
    --------
    >>> ceil(3.2)
    4
    >>> ceil(-3.2)
    -3
    >>> ceil(5.0)
    5
    """
    int_x = int(x)
    return int_x if x == int_x else int_x + 1 if x > 0 else int_x
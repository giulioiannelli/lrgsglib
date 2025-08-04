import numpy as np
from fractions import Fraction
from numpy.typing import NDArray
from typing import Union, Sequence, Any, Type, List
from ...config.const import *
#
__all__ = [
    'elements_within_eta_numpy',
    'dtype_numerical_precision',
    'linspace',
    'round_sigfig_n',
    'symmetric_logarithm_unchecked',
    'symmetric_logarithm',
    'to_fraction',
    'width_interval',
    'is_int',
    'is_positive_int',
]
#
def elements_within_eta_numpy(
    array: Union[np.ndarray, Sequence[float]], 
    eta: float
) -> np.ndarray:
    """
    Return elements that are within a threshold eta from the minimum value.

    Parameters
    ----------
    array : array-like
        Input array of numerical values. It will be converted to a NumPy array.
    eta : float
        Threshold value. Elements whose difference from the minimum value is 
        less than or equal to eta are returned.

    Returns
    -------
    np.ndarray
        Array of elements that are within eta of the minimum value.
    """
    array = np.array(array)
    min_val = np.min(array)
    mask = (array - min_val) <= eta
    filtered_elements = array[mask]
    return filtered_elements
#
def dtype_numerical_precision(dtype: Type[np.floating] = float) -> float:
    """
    Returns the smallest positive number that can be represented in 
    floating-point arithmetic for the specified data type.
    
    This value, known as machine epsilon, provides a measure of the numerical 
    precision or the difference between 1 and the smallest floating point 
    number greater than 1 for the given data type.
    
    Args:
        dtype (Type[np.floating]): The floating-point data type to get the 
                                 precision for (e.g., np.float32, np.float64).
                                 Defaults to float (typically np.float64).
    
    Returns:
        float: The machine epsilon for the given data type.
    """
    return np.finfo(dtype).eps
#
def linspace(
    start: float, 
    stop: float, 
    num: int = 50, 
    endpoint: bool = True
) -> List[float]:
    """
    Generate `num` evenly spaced samples from `start` to `stop`.

    Parameters
    ----------
    start : float
        The starting value of the sequence.
    stop : float
        The end value of the sequence.
    num : int, optional
        Number of samples to generate (default is 50). Must be >= 1.
    endpoint : bool, optional
        If True, `stop` is the last sample; otherwise it is not included 
        (default True).

    Returns
    -------
    List[float]
        List of evenly spaced samples.

    Raises
    ------
    ValueError
        If `num < 1`.

    Example
    -------
    >>> linspace(0.0, 1.0, num=5)
    [0.0, 0.25, 0.5, 0.75, 1.0]
    >>> linspace(0.0, 1.0, num=5, endpoint=False)
    [0.0, 0.2, 0.4, 0.6, 0.8]
    """
    if num < 1:
        raise ValueError(f"num must be >= 1, got {num}")
    if num == 1:
        return [stop] if endpoint else [start]
    steps = num - 1 if endpoint else num
    delta = (stop - start) / steps
    return [start + i * delta for i in range(num)]
#
def round_sigfig_n(num, n: int = 1):
    """
    Round a number or array of numbers to a specified number of significant 
    figures.

    Parameters:
    -----------
    num : float or array-like
        The number or array of numbers to be rounded to 'n' significant 
        figures.

    n : int, optional
        The desired number of significant figures (default is 1).
        Must be within the range [1, 15].

    Returns:
    --------
    float or numpy.ndarray
        The rounded number or array of numbers with 'n' significant figures.

    Raises:
    -------
    ValueError
        If 'n' is not within the valid range [1, 15].

    Notes:
    ------
    - The function calculates the exponent required to obtain 'n' significant 
      figures based on the absolute value of 'num'.
    - It then applies the rounding operation to 'num' with the calculated 
      exponent to achieve the desired number of significant figures.
    - If 'num' is an array-like object, it processes each element separately.

    Example:
    --------
    num = 123.456789
    n = 3
    rounded_num = round_sigfig_n(num, n)
    # The result is 123.0, rounded to 3 significant figures.

    num_array = [123.456789, 0.00123456789]
    rounded_array = round_sigfig_n(num_array, n)
    # The result is [123.0, 0.00123], rounded to 3 significant figures.
    """
    if n not in range(1, DEFAULT_MAX_DIGITS_ROUND_SIGFIG):
        raise ValueError("Significant figures number not in [1, 15].")
    expn = -np.floor(np.log10(np.abs(num))).astype('int')
    if hasattr(num, "__len__"):
        rr = np.array([np.round(nn, ee+n-1) for nn,ee in zip(num, expn)])
    else:
        rr = np.round(num, expn+n-1)
    return rr
#
def symmetric_logarithm_unchecked(
    x: Union[NDArray[np.floating], float],
    a: float,
    b: float,
    c: float,
    d: float
) -> Union[NDArray[np.floating], float]:
    """
    Compute a symmetric logarithmic function without safety adjustments.

    This applies
        a * log( b * (|x| - d) ) + c
    directly, which may produce invalid values if |x| ≤ d.

    Parameters
    ----------
    x : array-like or float
        Input values.
    a : float
        Scaling factor applied to the logarithm.
    b : float
        Multiplicative factor inside the log.
    c : float
        Additive constant.
    d : float
        Offset subtracted from |x| before taking the log.

    Returns
    -------
    array-like or float
        Result of the expression; may contain NaN or -inf where |x| ≤ d.
    """
    return a * np.log(b * (np.abs(x) - d)) + c
#
def symmetric_logarithm(
    x: Union[NDArray[np.floating], float],
    a: float,
    b: float,
    c: float,
    d: float,
    tol: float = 1e-10
) -> Union[NDArray[np.floating], float]:
    """
    Compute a symmetric logarithmic function with safety tolerance.

    Values of x whose magnitude is within `d + tol` are clamped to ±(d + tol)
    before applying the unchecked version, preventing invalid log arguments.

    Parameters
    ----------
    x : array-like or float
        Input values.
    a : float
        Scaling factor applied to the logarithm.
    b : float
        Multiplicative factor inside the log.
    c : float
        Additive constant.
    d : float
        Offset subtracted from |x| before taking the log.
    tol : float, optional
        Small positive tolerance to ensure |x| - d ≥ tol (default: 1e-10).

    Returns
    -------
    array-like or float
        Safe log-transformed values.
    """
    abs_x = np.abs(x)
    safe_val = np.where(abs_x > d + tol, abs_x, d + tol)
    return symmetric_logarithm_unchecked(safe_val, a, b, c, d)
#
def to_fraction(data: Any) -> Any:
    """
    Recursively convert numeric data to `Fraction`, preserving structure.

    - NumPy arrays → object‐dtype arrays of `Fraction`.
    - Sequences (e.g. lists, tuples) → lists of `Fraction`.
    - Floats → `Fraction.from_float(...).limit_denominator()`.
    - Other scalars (int, Decimal, str, numpy scalar, etc.) → `Fraction(data)`.

    Parameters
    ----------
    data
        A numeric scalar, sequence of numerics, or NumPy array.

    Returns
    -------
    Fraction or list of Fraction or ndarray of object
        Same structure as `data`, but all numeric elements converted to 
        `Fraction`.
    """
    # Handle NumPy arrays in one go
    if isinstance(data, np.ndarray):
        return np.vectorize(to_fraction, otypes=[object])(data)

    # Handle any non‐string sequence
    if isinstance(data, Sequence) and not isinstance(data, (str, bytes)):
        return [to_fraction(item) for item in data]

    # Floats get rational approximation
    if isinstance(data, float):
        return Fraction.from_float(data).limit_denominator()

    # Fallback covers ints, Decimals, numpy scalars, strings, etc.
    return Fraction(data)
#
def width_interval(a, b):
    """
    Calculate the width of an interval between two values.

    Parameters:
    -----------
    a : float or numeric
        The first value defining the interval.

    b : float or numeric
        The second value defining the interval.

    Returns:
    --------
    float or numeric
        The width of the interval, which is the absolute difference between 
        'a' and 'b'.

    Example:
    --------
    a = 5
    b = 8
    interval_width = width_interval(a, b)
    # The result is 3, representing the width of the interval [5, 8].
    """
    return np.abs(a - b)

def is_int(x: Any) -> bool:
    """
    Check if the input is an instance of GeneralInteger.

    Parameters
    ----------
    x : Any
        The input value to check.

    Returns
    -------
    bool
        True if the input is an instance of GeneralInteger, otherwise False.

    Example
    -------
    >>> is_int(5)
    True

    >>> is_int(3.14)
    False

    >>> is_int("string")
    False

    >>> is_int(np.int32(10))
    True
    """
    return isinstance(x, GeneralInteger)

def is_positive_int(x: Any) -> bool:
    """
    Check if the input is a positive integer.

    Parameters
    ----------
    x : Any
        The input value to check.

    Returns
    -------
    bool
        True if the input is an instance of GeneralInteger and greater than 0,
        otherwise False.

    Example
    -------
    >>> is_positive_int(5)
    True

    >>> is_positive_int(-3)
    False

    >>> is_positive_int("string")
    False
    """
    return is_int(x) and x > 0
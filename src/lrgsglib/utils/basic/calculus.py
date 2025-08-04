import numpy as np
from numpy.typing import NDArray
#
__all__ = [
    "line",
    "dv",
    "model_1r"
]
# 
def line(x, a, b):
    """
    Calculate the values of a straight line equation for given 'x' values.

    Parameters:
    -----------
    x : float or array-like
        The input values of 'x' for which the corresponding 'y' values are calculated.

    a : float
        The slope of the straight line.

    b : float
        The y-intercept of the straight line.

    Returns:
    --------
    float or numpy.ndarray
        The calculated 'y' values for the given 'x' values based on the straight line equation.

    Notes:
    ------
    - The straight line equation is 'y = a * x + b', where 'a' is the slope and 'b' is the y-intercept.
    - The function computes 'y' values for single 'x' values or arrays of 'x' values.

    Example:
    --------
    x_value = 2.0
    slope = 1.5
    y_intercept = 2.0
    y_result = line(x_value, slope, y_intercept)
    # The result is 5.0, which is the 'y' value for 'x' = 2.0 in the given line equation.

    """
    return a * x + b

#
def dv(f_x: NDArray, x: NDArray = None) -> NDArray:
    """
    Compute the derivative of an array with respect to another array.

    Parameters:
    -----------
    f_x : ndarray
        A NumPy array of N dimensions where the outermost dimension is the one
        where the derivative is computed using the `np.diff` method.

    x : ndarray, optional
        An array with the same dimensions as the outermost axis of 'f_x'.
        If not provided, it is assumed to be the range [0, len(f_x)].

    Returns:
    --------
    df_dx : ndarray
        The computed derivative of 'f_x' with respect to 'x'.

    Notes:
    ------
    - This function calculates the derivative of 'f_x' with respect to 'x' using
      the finite difference method. It computes the derivative along the outermost axis.
    - If 'x' is not provided, it is assumed to be the range [0, len(f_x)].

    Example:
    --------
    f_x = np.array([[1, 2, 3], [4, 5, 6]])
    x = np.array([0, 1, 2])
    df_dx = dv(f_x, x)
    # The result 'df_dx' contains the computed derivatives along the outermost axis.

    """
    if x is None:
        x = np.linspace(0, f_x.shape[-1], num=f_x.shape[-1])
    df_dx = np.diff(f_x, axis=-1) / np.diff(x)
    return df_dx

#
def model_1r(x, C_1, C_2):
    """
    Compute the one-over-r model with offset.

    Parameters:
    -----------
    x : float or array-like
        The input values of 'x' for which the model is evaluated.

    C_1 : float
        The scaling parameter for the 1/r term.

    C_2 : float
        The offset parameter.

    Returns:
    --------
    float or numpy.ndarray
        The calculated values based on the model C_2 * (1 - C_1 / (|x| * C_2)).

    Notes:
    ------
    - This model is commonly used for fitting percolation defect profiles.
    - The model has the form: C_2 * (1 - C_1 / (|x| * C_2)) = C_2 - C_1 / |x|
    - The absolute value of x is used to handle both positive and negative inputs.

    Example:
    --------
    x_values = np.array([1, 2, 3, 4, 5])
    C_1 = 2.0
    C_2 = 1.0
    result = model_1r(x_values, C_1, C_2)
    # Result will be the model evaluated at each x value
    """
    return C_2 * (1 - C_1 / (np.abs(x) * C_2))
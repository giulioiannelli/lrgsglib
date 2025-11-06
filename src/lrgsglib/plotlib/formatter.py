"""
Axis formatter utilities for matplotlib plotting.

This module provides custom formatters for matplotlib axes, particularly
for logarithmic and scientific notation formatting.
"""

from typing import Union, Optional
import numpy as np
from numpy.typing import NDArray


__all__ = ['log_latex_formatter']


def log_latex_formatter(
    x: Union[float, int],
    pos: Optional[int] = None,
    precision: int = 0,
    threshold: float = 0.01
) -> str:
    """
    Format a number for logarithmic axis display in matplotlib.
    
    This formatter converts numbers to LaTeX-formatted scientific notation,
    with special handling for zero and powers of 10. It's designed to be used
    with matplotlib's FuncFormatter for axis tick labels.
    
    Parameters
    ----------
    x : float or int
        The value to format.
    pos : int, optional
        The tick position (required by matplotlib FuncFormatter interface,
        but not used in this implementation). Default is None.
    precision : int, optional
        Number of decimal places to display for the coefficient when it's
        not close to 1. Default is 0.
    threshold : float, optional
        Threshold for considering a coefficient close to 1.0. If the absolute
        difference between the coefficient and 1.0 is less than this value,
        only the power of 10 is displayed. Default is 0.01.
    
    Returns
    -------
    str
        LaTeX-formatted string representation of the number.
        
    Examples
    --------
    >>> log_formatter(0, None)
    '$0$'
    
    >>> log_formatter(1000, None)
    '$10^{3}$'
    
    >>> log_formatter(2500, None)
    '$3\\times 10^{3}$'
    
    >>> log_formatter(0.001, None)
    '$10^{-3}$'
    
    Notes
    -----
    The function returns:
    - '$0$' for zero
    - '$10^{exp}$' for powers of 10 (when coefficient ≈ 1)
    - '$coeff\\times 10^{exp}$' for other values
    
    This function is thread-safe and has no side effects.
    
    See Also
    --------
    matplotlib.ticker.FuncFormatter : Matplotlib's function-based formatter
    matplotlib.ticker.LogFormatter : Matplotlib's built-in log formatter
    """
    # Handle zero case
    if x == 0:
        return r'$0$'
    
    # Calculate exponent and coefficient
    exponent = int(np.floor(np.log10(abs(x))))
    coefficient = x / (10 ** exponent)
    
    # Format based on coefficient value
    if abs(coefficient - 1.0) < threshold:
        # Display only power of 10
        return rf'$10^{{{exponent}}}$'
    else:
        # Display coefficient × 10^exponent
        return rf'${coefficient:.{precision}f}\times 10^{{{exponent}}}$'

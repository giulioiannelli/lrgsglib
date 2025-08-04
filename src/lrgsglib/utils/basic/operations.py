from .common import *
#
def no_operation(*args: Any, **kwargs: Any) -> None:
    """
    A no-op (no operation) function that accepts any arguments and does nothing.

    Parameters
    ----------
    *args : Any
        Positional arguments (ignored).
    **kwargs : Any
        Keyword arguments (ignored).

    Returns
    -------
    None
    """
    return None

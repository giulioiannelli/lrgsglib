from .common import *

def verbose_print(
    message: Any,
    verbose: bool = False,
    **print_kwargs: Any
) -> None:
    """
    Print a message only if verbose mode is enabled.

    Parameters
    ----------
    message : Any
        The object to be printed.
    verbose : bool, optional
        If True, the message is printed (default is False).
    **print_kwargs : Any
        Additional keyword arguments forwarded to `print()`.
    """
    if verbose:
        print(message, **print_kwargs)
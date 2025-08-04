from .common import *
#
def compose(
    f: Callable[..., Any],
    g: Callable[..., Any],
    g_args: Tuple[Any, ...] = (),
    g_kwargs: Optional[Dict[str, Any]] = None
) -> Callable[..., Any]:
    """
    Return a function that applies f to its arguments and then applies g to f’s result.

    Parameters
    ----------
    f : Callable[..., Any]
        First function to apply; accepts any arguments.
    g : Callable[..., Any]
        Second function to apply; its first argument is the output of `f`,
        followed by any `g_args` and `g_kwargs`.
    g_args : Tuple[Any, ...], optional
        Positional arguments to append when calling `g` (default is ()).
    g_kwargs : Dict[str, Any], optional
        Keyword arguments to pass to `g` (default is None, treated as {}).

    Returns
    -------
    Callable[..., Any]
        A function `h` such that `h(*args, **kwargs)` returns
        `g(f(*args, **kwargs), *g_args, **g_kwargs)`.
    """
    g_kwargs = g_kwargs or {}

    def composed(*args: Any, **kwargs: Any) -> Any:
        result = f(*args, **kwargs)
        return g(result, *g_args, **g_kwargs)

    return composed
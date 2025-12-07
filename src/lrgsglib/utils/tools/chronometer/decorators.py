from functools import wraps
from .Chronometer import Chronometer

def time_function_accumulate(auto_log: bool=False):
    """
    Decorator to time a function and accumulate its duration.

    Always accumulates timing data. The auto_log parameter controls whether
    to print timing info immediately after each call.

    Parameters
    ----------
    auto_log : bool
        If True, log each call's duration immediately. If False, data is
        still accumulated and can be printed later with -pc flag.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Always create chronometer to accumulate data
            chrono = Chronometer(func.__name__, auto_log=auto_log)
            result = func(*args, **kwargs)
            chrono.end()
            return result
        return wrapper
    return decorator


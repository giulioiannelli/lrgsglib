from functools import wraps
from .Chronometer import Chronometer

def time_function_accumulate(auto_log: bool=False):
    """
    Decorator to time a function and accumulate its duration.

    Parameters
    ----------
    auto_log : bool
        If True, log each call's duration.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            log_enabled = auto_log or Chronometer.enabled
            if not log_enabled:
                return func(*args, **kwargs)
            chrono = Chronometer(func.__name__, auto_log=log_enabled)
            result = func(*args, **kwargs)
            chrono.end()
            return result
        return wrapper
    return decorator


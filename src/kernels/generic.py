import logging
import os
import numpy as np
from pathlib import Path
from typing import Optional, Any
from lrgsglib.loglib import setup_custom_logger

__all__ = [
    "setup_logger",
    "initialize_custom_logger",
    "resolve_backend",
    "resolve_float_type",
]

# Valid backends and float types
_VALID_BACKENDS = {"numpy", "cupy"}
_VALID_FLOAT_TYPES = {"float32": np.float32, "float64": np.float64}

def resolve_backend(backend: str) -> str:
    """Validate and normalise the backend identifier."""
    backend_norm = backend.lower()
    if backend_norm not in _VALID_BACKENDS:
        raise ValueError(f"Unsupported backend '{backend}'. Expected one of {_VALID_BACKENDS}.")
    if backend_norm == "cupy":
        try:
            import cupy  # noqa: F401
        except ImportError as exc:  # pragma: no cover - informative failure path
            raise RuntimeError(
                "Backend 'cupy' requested but CuPy is not available. Install cupy or choose 'numpy'."
            ) from exc
    return backend_norm

def resolve_float_type(float_type: str) -> type:
    """Map a float type string to the corresponding NumPy dtype."""
    try:
        return _VALID_FLOAT_TYPES[float_type]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported float type '{float_type}'. Expected one of {tuple(_VALID_FLOAT_TYPES)}."
        ) from exc
#
def setup_logger(loglogger: Optional[Any]) -> Any:
    """
    Setup the logger. If no logger is provided, returns a null logger.
    
    Parameters
    ----------
    loglogger : Optional[Any]
        A logger instance or None.

    Returns
    -------
    Any
        A valid logger instance.
    """
    if loglogger is None:
        loglogger = logging.getLogger("null")
        loglogger.addHandler(logging.NullHandler())
    return loglogger
#
def initialize_custom_logger(args_dict: dict, args_dict_keys: dict, out_suffix: str = None, level: int = logging.INFO) -> logging.Logger:
    """
    Initialize a custom logger.

    Parameters:
    -----------
    args_dict : dict
        Dictionary of parsed arguments.
    args_dict_keys : dict
        Dictionary of arguments.
    out_suffix : str
        Output suffix for the logger name.
    level : int, optional
        Logging level to be passed to the logger. Defaults to logging.INFO.

    Returns:
    --------
    logging.Logger
        Configured logger instance.
    """
    __thisfile__ = Path(__file__).stem
    if out_suffix is None:
        out_suffix = args_dict.get('out_suffix', '')
    argstr = '_'.join([str(k)+str(args_dict[k]) for k in args_dict_keys.keys()] + [out_suffix])
    loggername = __thisfile__ + argstr + f'-{os.getpid()}'
    return setup_custom_logger(loggername, level=level)
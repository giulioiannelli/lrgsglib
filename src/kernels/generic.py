import logging
import os
#
from pathlib import Path
from typing import Optional, Any
#
from lrgsglib.loglib import setup_custom_logger
#
__all__ = ["setup_logger", "initialize_custom_logger"]
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
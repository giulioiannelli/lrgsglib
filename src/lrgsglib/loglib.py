"""Legacy logging module - DEPRECATED.

Use `lrgsglib.utils.tools.loglib` instead for the new opt-in logging system.

This module is kept for backward compatibility but will be removed in a future
version. The new loglib provides:

- Disabled by default (no logs unless enabled)
- Session-specific naming
- Size-bounded rotating files (10MB limit)
- WARNING level by default

Migration example::

    # Old (deprecated):
    from lrgsglib.loglib import setup_custom_logger
    logger = setup_custom_logger(__name__)

    # New (recommended):
    from lrgsglib.utils.tools.loglib import get_logger, enable_logging
    logger = get_logger(__name__)
    enable_logging(session_id="my_session")  # Only when you want logs
"""

import logging
import os
import warnings
from pathlib import Path
from .config.const import DEFAULT_LOG_DIR, LOG


def setup_custom_logger(
    name: str,
    level: int = logging.INFO,
    log_file: Path = None
) -> logging.Logger:
    """
    Configure a logger with a FileHandler whose filename is based on `name`.

    Parameters
    ----------
    name : str
        Name of the logger (typically use `__name__`).
    level : int
        Logging level (default: logging.INFO).
    log_file : Path, optional
        Explicit path to log file. If None, it will be `DEFAULT_LOG_DIR / <name>.log`.

    Returns
    -------
    logger : logging.Logger
        The configured logger instance.

    .. deprecated::
        Use `lrgsglib.utils.tools.loglib.get_logger()` and `enable_logging()` instead.
    """
    warnings.warn(
        "setup_custom_logger() is deprecated. "
        "Use lrgsglib.utils.tools.loglib.get_logger() and enable_logging() instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    log_file = Path(log_file or os.environ.get("MYLIB_LOG_FILE") or Path(DEFAULT_LOG_DIR) / (name+LOG))
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(level)

    if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_file)
               for h in logger.handlers):
        fh = logging.FileHandler(log_file, encoding="utf-8")
        fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        fh.setFormatter(logging.Formatter(fmt))
        logger.addHandler(fh)

    return logger
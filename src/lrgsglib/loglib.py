import logging
import os
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
    """
    log_file = Path(log_file or os.environ.get("MYLIB_LOG_FILE") or DEFAULT_LOG_DIR / (name+LOG))
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
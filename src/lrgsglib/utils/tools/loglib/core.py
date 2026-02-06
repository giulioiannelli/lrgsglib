"""Core logging functionality for lrgsglib.

Provides opt-in logging with sensible defaults:
- Disabled by default (NullHandler)
- WARNING level when enabled
- Size-bounded rotating files
- Session-specific naming
"""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .handlers import SessionRotatingHandler, SizeAwareNullHandler

# Module-level state
_logging_enabled: bool = False
_active_handlers: list[logging.Handler] = []
_root_logger_name: str = "lrgsglib"


@dataclass
class LogConfig:
    """Configuration for lrgsglib logging.

    Attributes
    ----------
    max_bytes : int
        Maximum log file size before rotation (default: 10MB).
    backup_count : int
        Number of backup files to keep (default: 3).
    default_level : int
        Default logging level when enabled (default: WARNING).
    log_dir : Path
        Directory for log files.
    """

    max_bytes: int = 10 * 1024 * 1024  # 10 MB
    backup_count: int = 3
    default_level: int = logging.WARNING
    log_dir: Path = Path(".log")

    @classmethod
    def from_env(cls) -> "LogConfig":
        """Create config from environment variables.

        Environment variables:
        - LRGSGLIB_LOG_LEVEL: DEBUG, INFO, WARNING, ERROR, CRITICAL
        - LRGSGLIB_LOG_DIR: Path to log directory
        - LRGSGLIB_LOG_MAX_BYTES: Maximum file size
        - LRGSGLIB_LOG_BACKUP_COUNT: Number of backups
        """
        config = cls()

        # Log level from environment
        level_str = os.environ.get("LRGSGLIB_LOG_LEVEL", "").upper()
        level_map = {
            "DEBUG": logging.DEBUG,
            "INFO": logging.INFO,
            "WARNING": logging.WARNING,
            "ERROR": logging.ERROR,
            "CRITICAL": logging.CRITICAL,
        }
        if level_str in level_map:
            config.default_level = level_map[level_str]

        # Log directory from environment (prefer LRGSG_LOG if set)
        log_dir = os.environ.get("LRGSGLIB_LOG_DIR") or os.environ.get("LRGSG_LOG")
        if log_dir:
            config.log_dir = Path(log_dir)

        # Size limits from environment
        if max_bytes := os.environ.get("LRGSGLIB_LOG_MAX_BYTES"):
            try:
                config.max_bytes = int(max_bytes)
            except ValueError:
                pass

        if backup_count := os.environ.get("LRGSGLIB_LOG_BACKUP_COUNT"):
            try:
                config.backup_count = int(backup_count)
            except ValueError:
                pass

        return config


# Default config (can be modified before enable_logging)
_config = LogConfig.from_env()


def get_logger(name: str) -> logging.Logger:
    """Get a logger for the given module name.

    Returns a logger configured for lrgsglib. By default, all loggers
    use NullHandler (silent). Call enable_logging() to activate file logging.

    Parameters
    ----------
    name : str
        Logger name, typically __name__ from the calling module.

    Returns
    -------
    logging.Logger
        Configured logger instance.

    Examples
    --------
    >>> logger = get_logger(__name__)
    >>> logger.debug("This is silently discarded by default")
    >>>
    >>> # To enable logging:
    >>> enable_logging(session_id="my_session")
    >>> logger.warning("This goes to file")
    """
    # Ensure name is under lrgsglib namespace
    if not name.startswith(_root_logger_name):
        name = f"{_root_logger_name}.{name}"

    logger = logging.getLogger(name)

    # Add NullHandler if no handlers (prevents "No handler" warnings)
    if not logger.handlers and not logger.parent.handlers:
        logger.addHandler(SizeAwareNullHandler())

    return logger


def enable_logging(
    session_id: Optional[str] = None,
    level: Optional[int] = None,
    console: bool = False,
    log_dir: Optional[Path] = None,
) -> None:
    """Enable file logging for lrgsglib.

    Activates logging to a rotating file with the given session ID.
    If no session_id is provided, uses "default" with timestamp.

    Parameters
    ----------
    session_id : str, optional
        Session identifier for log filename. If None, uses "default".
    level : int, optional
        Logging level (default: WARNING or from LRGSGLIB_LOG_LEVEL env).
    console : bool
        Whether to also log to console/stderr (default: False).
    log_dir : Path, optional
        Directory for log files (default: from config or LRGSG_LOG env).

    Examples
    --------
    >>> enable_logging(session_id="mcg_batch_001")
    # Creates: .log/lrgsglib_mcg_batch_001.log

    >>> enable_logging(level=logging.DEBUG, console=True)
    # Also prints to console at DEBUG level
    """
    global _logging_enabled, _active_handlers

    # Clean up any existing handlers first
    disable_logging()

    # Use defaults from config
    actual_level = level if level is not None else _config.default_level
    actual_log_dir = log_dir if log_dir is not None else _config.log_dir

    # Generate session ID if not provided
    if session_id is None:
        from datetime import datetime
        session_id = f"default_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    # Get root logger
    root_logger = logging.getLogger(_root_logger_name)
    root_logger.setLevel(actual_level)

    # Remove any existing NullHandlers
    for handler in root_logger.handlers[:]:
        if isinstance(handler, (logging.NullHandler, SizeAwareNullHandler)):
            root_logger.removeHandler(handler)

    # Add rotating file handler
    file_handler = SessionRotatingHandler(
        log_dir=actual_log_dir,
        session_id=session_id,
        max_bytes=_config.max_bytes,
        backup_count=_config.backup_count,
    )
    file_handler.setLevel(actual_level)
    root_logger.addHandler(file_handler)
    _active_handlers.append(file_handler)

    # Optionally add console handler
    if console:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(actual_level)
        fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        console_handler.setFormatter(logging.Formatter(fmt))
        root_logger.addHandler(console_handler)
        _active_handlers.append(console_handler)

    _logging_enabled = True


def disable_logging() -> None:
    """Disable file logging and clean up handlers.

    Removes all active handlers added by enable_logging() and restores
    the silent NullHandler behavior.
    """
    global _logging_enabled, _active_handlers

    root_logger = logging.getLogger(_root_logger_name)

    # Remove active handlers
    for handler in _active_handlers:
        try:
            handler.close()
        except Exception:
            pass
        try:
            root_logger.removeHandler(handler)
        except Exception:
            pass

    _active_handlers.clear()

    # Ensure NullHandler is present
    if not any(isinstance(h, logging.NullHandler) for h in root_logger.handlers):
        root_logger.addHandler(SizeAwareNullHandler())

    _logging_enabled = False


def set_level(level: int) -> None:
    """Set the logging level for all lrgsglib loggers.

    Parameters
    ----------
    level : int
        Logging level (e.g., logging.DEBUG, logging.INFO).
    """
    root_logger = logging.getLogger(_root_logger_name)
    root_logger.setLevel(level)

    for handler in _active_handlers:
        handler.setLevel(level)


def is_logging_enabled() -> bool:
    """Check if logging is currently enabled.

    Returns
    -------
    bool
        True if enable_logging() has been called and not disabled.
    """
    return _logging_enabled


def get_config() -> LogConfig:
    """Get the current logging configuration.

    Returns
    -------
    LogConfig
        Current configuration dataclass.
    """
    return _config


def set_config(config: LogConfig) -> None:
    """Set the logging configuration.

    Should be called before enable_logging() to take effect.

    Parameters
    ----------
    config : LogConfig
        New configuration to use.
    """
    global _config
    _config = config

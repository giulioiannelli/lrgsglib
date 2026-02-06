"""Lightweight modular logging for lrgsglib.

This module provides opt-in logging with sensible defaults:

- **Disabled by default** - No logs created unless explicitly enabled
- **Session naming** - Logs identify which execution they belong to
- **Size bounded** - RotatingFileHandler with 10MB limit, 3 backups
- **WARNING level** - Only problems logged unless DEBUG requested

Basic Usage
-----------
>>> from lrgsglib.utils.tools.loglib import get_logger
>>> logger = get_logger(__name__)
>>> logger.warning("This is silently discarded by default")

Enable Logging
--------------
>>> from lrgsglib.utils.tools.loglib import enable_logging
>>> enable_logging(session_id="mcg_batch_001")
# Creates: .log/lrgsglib_mcg_batch_001.log

With Graph Context
------------------
>>> from lrgsglib.utils.tools.loglib import session_from_graph, enable_logging
>>> from lrgsglib import MultiplicativeCascadeGraph
>>> mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6)
>>> ctx = session_from_graph(mc, prefix="mcg")
>>> enable_logging(session_id=ctx.generate_id())
# Creates: .log/lrgsglib_mcg_N=4096_p1=0.8_p2=0.6_20260202_143052.log

Debug Mode
----------
>>> import logging
>>> enable_logging(level=logging.DEBUG, console=True)

Environment Variables
---------------------
- LRGSGLIB_LOG_LEVEL: DEBUG, INFO, WARNING, ERROR, CRITICAL
- LRGSGLIB_LOG_DIR: Path to log directory (falls back to LRGSG_LOG)
- LRGSGLIB_LOG_MAX_BYTES: Maximum file size before rotation
- LRGSGLIB_LOG_BACKUP_COUNT: Number of backup files to keep
"""

from .core import (
    get_logger,
    enable_logging,
    disable_logging,
    set_level,
    is_logging_enabled,
    get_config,
    set_config,
    LogConfig,
)
from .context import (
    SessionContext,
    session_from_graph,
    session_from_dynamics,
)
from .handlers import (
    SessionRotatingHandler,
    SizeAwareNullHandler,
)

__all__ = [
    # Core functions
    "get_logger",
    "enable_logging",
    "disable_logging",
    "set_level",
    "is_logging_enabled",
    "get_config",
    "set_config",
    "LogConfig",
    # Context utilities
    "SessionContext",
    "session_from_graph",
    "session_from_dynamics",
    # Handlers (for advanced use)
    "SessionRotatingHandler",
    "SizeAwareNullHandler",
]

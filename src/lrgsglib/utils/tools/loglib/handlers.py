"""Custom logging handlers for lrgsglib.

Provides rotating file handlers with session naming and size-aware null handlers
for efficient silent mode.
"""

from __future__ import annotations

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Optional


class SessionRotatingHandler(RotatingFileHandler):
    """Rotating file handler with session-based naming.

    Creates log files with the pattern: lrgsglib_{session_id}.log
    Rotates when file exceeds max_bytes, keeping backup_count old files.

    Parameters
    ----------
    log_dir : Path
        Directory for log files.
    session_id : str
        Session identifier for the log filename.
    max_bytes : int
        Maximum size in bytes before rotation (default: 10MB).
    backup_count : int
        Number of backup files to keep (default: 3).
    encoding : str
        File encoding (default: utf-8).

    Examples
    --------
    >>> handler = SessionRotatingHandler(
    ...     log_dir=Path(".log"),
    ...     session_id="mcg_batch_001",
    ...     max_bytes=10*1024*1024,
    ...     backup_count=3
    ... )
    # Creates: .log/lrgsglib_mcg_batch_001.log
    """

    def __init__(
        self,
        log_dir: Path,
        session_id: str,
        max_bytes: int = 10 * 1024 * 1024,
        backup_count: int = 3,
        encoding: str = "utf-8",
    ) -> None:
        log_dir = Path(log_dir)
        log_dir.mkdir(parents=True, exist_ok=True)

        filename = log_dir / f"lrgsglib_{session_id}.log"

        super().__init__(
            filename=str(filename),
            maxBytes=max_bytes,
            backupCount=backup_count,
            encoding=encoding,
        )

        self.session_id = session_id

        # Set default format
        fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        self.setFormatter(logging.Formatter(fmt))


class SizeAwareNullHandler(logging.NullHandler):
    """Null handler that tracks suppressed message counts.

    Useful for debugging when logging is disabled - can report how many
    messages would have been logged.

    Attributes
    ----------
    suppressed_count : int
        Number of messages suppressed.
    suppressed_by_level : dict[int, int]
        Count of suppressed messages by logging level.

    Examples
    --------
    >>> handler = SizeAwareNullHandler()
    >>> logger.addHandler(handler)
    >>> # ... run code that logs ...
    >>> print(f"Suppressed {handler.suppressed_count} messages")
    """

    def __init__(self) -> None:
        super().__init__()
        self.suppressed_count: int = 0
        self.suppressed_by_level: dict[int, int] = {}

    def handle(self, record: logging.LogRecord) -> bool:
        """Count suppressed records."""
        self.suppressed_count += 1
        level = record.levelno
        self.suppressed_by_level[level] = self.suppressed_by_level.get(level, 0) + 1
        return super().handle(record)

    def reset_counts(self) -> None:
        """Reset suppression counters."""
        self.suppressed_count = 0
        self.suppressed_by_level.clear()

    def get_summary(self) -> str:
        """Return summary of suppressed messages."""
        if self.suppressed_count == 0:
            return "No messages suppressed"

        level_names = {
            logging.DEBUG: "DEBUG",
            logging.INFO: "INFO",
            logging.WARNING: "WARNING",
            logging.ERROR: "ERROR",
            logging.CRITICAL: "CRITICAL",
        }

        parts = [f"Suppressed {self.suppressed_count} messages:"]
        for level, count in sorted(self.suppressed_by_level.items()):
            name = level_names.get(level, f"LEVEL_{level}")
            parts.append(f"  {name}: {count}")

        return "\n".join(parts)

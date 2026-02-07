"""Tests for lrgsglib.utils.tools.loglib module.

Tests the opt-in logging system including:
- Default silent behavior (NullHandler)
- enable_logging / disable_logging
- Session naming and context generation
- Rotating file handler
- Log level control
"""

import logging
import os
import tempfile
from pathlib import Path

import pytest


class TestGetLogger:
    """Tests for get_logger function."""

    def test_returns_logger_instance(self):
        """get_logger should return a logging.Logger instance."""
        from lrgsglib.utils.tools.loglib import get_logger

        logger = get_logger("test_module")
        assert isinstance(logger, logging.Logger)

    def test_logger_name_prefixed(self):
        """Logger name should be prefixed with 'lrgsglib'."""
        from lrgsglib.utils.tools.loglib import get_logger

        logger = get_logger("my_module")
        assert logger.name.startswith("lrgsglib")

    def test_already_prefixed_name_not_doubled(self):
        """Logger name already prefixed should not be doubled."""
        from lrgsglib.utils.tools.loglib import get_logger

        logger = get_logger("lrgsglib.some_module")
        assert logger.name == "lrgsglib.some_module"
        assert not logger.name.startswith("lrgsglib.lrgsglib")

    def test_default_handler_is_null(self):
        """By default, logger should have NullHandler (silent)."""
        from lrgsglib.utils.tools.loglib import get_logger, disable_logging

        # Ensure logging is disabled
        disable_logging()

        logger = get_logger("test_null_handler")
        # Should have at least one handler, and root should have NullHandler
        root = logging.getLogger("lrgsglib")
        has_null = any(
            isinstance(h, logging.NullHandler) for h in root.handlers
        )
        assert has_null or len(root.handlers) == 0


class TestEnableDisableLogging:
    """Tests for enable_logging and disable_logging functions."""

    def test_enable_creates_log_file(self):
        """enable_logging should create a log file in the specified directory."""
        from lrgsglib.utils.tools.loglib import enable_logging, disable_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            enable_logging(session_id="test_session", log_dir=Path(tmpdir))

            # Check that log file was created
            log_files = list(Path(tmpdir).glob("*.log"))
            assert len(log_files) == 1
            assert "test_session" in log_files[0].name

            disable_logging()

    def test_disable_stops_file_logging(self):
        """disable_logging should stop writing to files."""
        from lrgsglib.utils.tools.loglib import (
            enable_logging,
            disable_logging,
            get_logger,
            is_logging_enabled,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            enable_logging(session_id="test_disable", log_dir=Path(tmpdir))
            assert is_logging_enabled()

            disable_logging()
            assert not is_logging_enabled()

    def test_console_handler_added_when_requested(self):
        """enable_logging with console=True should add StreamHandler."""
        from lrgsglib.utils.tools.loglib import enable_logging, disable_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            enable_logging(
                session_id="test_console",
                log_dir=Path(tmpdir),
                console=True,
            )

            root = logging.getLogger("lrgsglib")
            has_stream = any(
                isinstance(h, logging.StreamHandler)
                and not isinstance(h, logging.FileHandler)
                for h in root.handlers
            )
            assert has_stream

            disable_logging()

    def test_level_can_be_set(self):
        """enable_logging should respect the level parameter."""
        from lrgsglib.utils.tools.loglib import enable_logging, disable_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            enable_logging(
                session_id="test_level",
                log_dir=Path(tmpdir),
                level=logging.DEBUG,
            )

            root = logging.getLogger("lrgsglib")
            assert root.level == logging.DEBUG

            disable_logging()


class TestSessionContext:
    """Tests for SessionContext and session_from_graph."""

    def test_session_context_generates_id(self):
        """SessionContext should generate a valid session ID."""
        from lrgsglib.utils.tools.loglib import SessionContext

        ctx = SessionContext(
            prefix="mcg",
            params={"N": 4096, "p1": 0.8},
            include_timestamp=False,
        )
        session_id = ctx.generate_id()

        assert session_id.startswith("mcg")
        assert "N=4096" in session_id
        assert "p1=0.8" in session_id

    def test_session_context_with_timestamp(self):
        """SessionContext with timestamp should include date."""
        from lrgsglib.utils.tools.loglib import SessionContext

        ctx = SessionContext(
            prefix="test",
            params={},
            include_timestamp=True,
        )
        session_id = ctx.generate_id()

        # Should contain a timestamp pattern (YYYYMMDD_HHMMSS)
        import re
        assert re.search(r"\d{8}_\d{6}", session_id)

    def test_session_context_add_param_chainable(self):
        """add_param should return self for chaining."""
        from lrgsglib.utils.tools.loglib import SessionContext

        ctx = SessionContext(prefix="test", include_timestamp=False)
        result = ctx.add_param("key1", "val1").add_param("key2", "val2")

        assert result is ctx
        assert ctx.params["key1"] == "val1"
        assert ctx.params["key2"] == "val2"


class TestSessionFromGraph:
    """Tests for session_from_graph function."""

    def test_extracts_basic_params(self):
        """session_from_graph should extract N and other parameters."""
        from lrgsglib.utils.tools.loglib import session_from_graph

        # Create a mock graph object
        class MockGraph:
            N = 100
            pflip = 0.3
            p1 = 0.8

        ctx = session_from_graph(MockGraph(), prefix="mock", include_timestamp=False)
        session_id = ctx.generate_id()

        assert "N=100" in session_id
        assert "pflip=0.3" in session_id
        assert "p1=0.8" in session_id

    def test_skips_default_pflip(self):
        """session_from_graph should skip pflip=0.0 (default)."""
        from lrgsglib.utils.tools.loglib import session_from_graph

        class MockGraph:
            N = 100
            pflip = 0.0

        ctx = session_from_graph(MockGraph(), prefix="mock", include_timestamp=False)
        session_id = ctx.generate_id()

        assert "pflip" not in session_id


class TestRotatingHandler:
    """Tests for SessionRotatingHandler."""

    def test_creates_file_with_session_name(self):
        """Handler should create file with lrgsglib_{session_id}.log pattern."""
        from lrgsglib.utils.tools.loglib.handlers import SessionRotatingHandler

        with tempfile.TemporaryDirectory() as tmpdir:
            handler = SessionRotatingHandler(
                log_dir=Path(tmpdir),
                session_id="rotation_test",
            )
            handler.close()

            log_files = list(Path(tmpdir).glob("*.log"))
            assert len(log_files) == 1
            assert log_files[0].name == "lrgsglib_rotation_test.log"

    def test_has_formatter(self):
        """Handler should have a formatter set."""
        from lrgsglib.utils.tools.loglib.handlers import SessionRotatingHandler

        with tempfile.TemporaryDirectory() as tmpdir:
            handler = SessionRotatingHandler(
                log_dir=Path(tmpdir),
                session_id="formatter_test",
            )

            assert handler.formatter is not None
            handler.close()


class TestSizeAwareNullHandler:
    """Tests for SizeAwareNullHandler."""

    def test_counts_suppressed_messages(self):
        """NullHandler should track suppressed message count."""
        from lrgsglib.utils.tools.loglib.handlers import SizeAwareNullHandler

        handler = SizeAwareNullHandler()
        logger = logging.getLogger("test_null_count")
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)

        logger.debug("Message 1")
        logger.info("Message 2")
        logger.warning("Message 3")

        assert handler.suppressed_count == 3
        assert handler.suppressed_by_level[logging.DEBUG] == 1
        assert handler.suppressed_by_level[logging.INFO] == 1
        assert handler.suppressed_by_level[logging.WARNING] == 1

        # Cleanup
        logger.removeHandler(handler)

    def test_reset_counts(self):
        """reset_counts should clear all counters."""
        from lrgsglib.utils.tools.loglib.handlers import SizeAwareNullHandler

        handler = SizeAwareNullHandler()
        handler.suppressed_count = 10
        handler.suppressed_by_level[logging.DEBUG] = 5

        handler.reset_counts()

        assert handler.suppressed_count == 0
        assert len(handler.suppressed_by_level) == 0


class TestLogConfig:
    """Tests for LogConfig dataclass."""

    def test_default_values(self):
        """LogConfig should have sensible defaults."""
        from lrgsglib.utils.tools.loglib import LogConfig

        config = LogConfig()

        assert config.max_bytes == 10 * 1024 * 1024  # 10 MB
        assert config.backup_count == 3
        assert config.default_level == logging.WARNING

    def test_from_env_reads_level(self):
        """LogConfig.from_env should read LRGSGLIB_LOG_LEVEL."""
        from lrgsglib.utils.tools.loglib import LogConfig

        # Save original value
        original = os.environ.get("LRGSGLIB_LOG_LEVEL")

        try:
            os.environ["LRGSGLIB_LOG_LEVEL"] = "DEBUG"
            config = LogConfig.from_env()
            assert config.default_level == logging.DEBUG
        finally:
            # Restore original value
            if original is not None:
                os.environ["LRGSGLIB_LOG_LEVEL"] = original
            else:
                os.environ.pop("LRGSGLIB_LOG_LEVEL", None)


class TestIntegration:
    """Integration tests for the full logging workflow."""

    def test_full_workflow(self):
        """Test complete logging workflow: enable, log, disable."""
        from lrgsglib.utils.tools.loglib import (
            get_logger,
            enable_logging,
            disable_logging,
            is_logging_enabled,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            # Initially disabled
            disable_logging()
            assert not is_logging_enabled()

            # Enable logging
            enable_logging(
                session_id="integration_test",
                log_dir=Path(tmpdir),
                level=logging.DEBUG,
            )
            assert is_logging_enabled()

            # Log a message
            logger = get_logger("integration_module")
            logger.debug("Debug message")
            logger.info("Info message")
            logger.warning("Warning message")

            # Disable and verify
            disable_logging()
            assert not is_logging_enabled()

            # Check log file contents
            log_file = Path(tmpdir) / "lrgsglib_integration_test.log"
            assert log_file.exists()

            contents = log_file.read_text()
            assert "Debug message" in contents
            assert "Info message" in contents
            assert "Warning message" in contents

    def test_silent_by_default(self):
        """Logging should be silent by default (no files created)."""
        from lrgsglib.utils.tools.loglib import get_logger, disable_logging

        with tempfile.TemporaryDirectory() as tmpdir:
            # Ensure logging is disabled
            disable_logging()

            # Get a logger and log some messages
            logger = get_logger("silent_test")
            logger.info("This should not create a file")
            logger.warning("Neither should this")

            # No log files should be created in tmpdir
            log_files = list(Path(tmpdir).glob("*.log"))
            assert len(log_files) == 0

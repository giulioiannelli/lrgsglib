"""Mixin for C backend integration in binary dynamics systems.

This module provides shared infrastructure for running C-based simulators
from Python. The CBackendMixin class extracts common patterns from
IsingDynamics, VoterModel, and ContactProcess to eliminate duplication
while maintaining backward compatibility.

Usage
-----
Subclasses should:
1. Inherit from both CBackendMixin and BinDynSys (mixin first)
2. Set `_c_program_name_template` class attribute
3. Set `_allowed_c_keys` tuple if validation is needed
4. Implement `_build_c_arglist()` for model-specific arguments
5. Optionally override `_get_cleanup_paths()` for model-specific cleanup

Example
-------
>>> class MyDynamics(CBackendMixin, BinDynSys):
...     _c_program_name_template = "MySimulator{}"
...     _allowed_c_keys = ("C0", "C1")
...
...     def _build_c_arglist(self) -> list[str]:
...         return [f"{self.N}", f"{self.param}"]
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from ..config.const import LOG, LRGSG_LOG
from ..utils.basic.strings import join_non_empty


def _get_ccore_bin() -> str:
    """Get C core binary directory path at runtime.

    Reads from environment variable to allow test isolation.
    Falls back to configured constant if not set.
    """
    env_val = os.environ.get("LRGSG_CCORE_BIN")
    if env_val:
        return env_val
    # Fall back to configured value from lrgsg_env
    from ..config.lrgsg_env import LRGSG_CCORE_BIN

    return LRGSG_CCORE_BIN

if TYPE_CHECKING:
    from .BinDynSys import BinDynSys


class CBackendMixin:
    """Mixin providing C backend integration for dynamics simulations.

    This mixin extracts common patterns for running C simulators:
    - stderr logging setup
    - subprocess execution with stdout parsing
    - file cleanup utilities

    Subclasses must implement:
    - `_c_program_name_template: str` - e.g., "IsingSimulator{}"
    - `_build_c_arglist(self) -> list[str]` - Build argument list

    Optional overrides:
    - `_allowed_c_keys: tuple[str, ...]` - Whitelist of runlang values
    - `_get_cleanup_paths(self) -> list[Path | None]` - Paths to clean up
    - `_c_program_suffix(self) -> str` - Custom suffix extraction
    """

    # Subclasses should override these
    _c_program_name_template: str = ""
    _allowed_c_keys: tuple[str, ...] = ()

    # These are set during operation
    CbaseName: str = ""
    cprogram: list[str | Path] = []
    stderr_path: Path | None = None
    stderr_fopen: Any = None

    # ------------------------------------------------------------------
    # C program key validation
    # ------------------------------------------------------------------
    def _c_program_key(self: "BinDynSys") -> str:
        """Validate and normalize C backend key from runlang.

        Converts runlang values like 'c1b', 'C1B', 'C0' to normalized
        uppercase format (e.g., 'C1B').

        Returns
        -------
        str
            Normalized key like 'C0', 'C1', 'C1B', etc.

        Raises
        ------
        ValueError
            If runlang doesn't start with 'C' or isn't in _allowed_c_keys.
        """
        if not self.runlang or not self.runlang.upper().startswith("C"):
            raise ValueError("C backend requested but runlang is not a C variant.")
        suffix = self.runlang[1:]
        key = "C0" if suffix == "" or suffix == "0" else f"C{suffix.upper()}"
        if self._allowed_c_keys and key not in self._allowed_c_keys:
            raise ValueError(
                f"runlang '{self.runlang}' not supported for {self.__class__.__name__}. "
                f"Allowed: {self._allowed_c_keys}"
            )
        return key

    def _c_program_suffix(self: "BinDynSys") -> str:
        """Extract program suffix from runlang (e.g., 'C1c' -> '1c').

        Override this method for custom suffix extraction logic.

        Returns
        -------
        str
            Suffix to append to program name template.
        """
        key = self._c_program_key()
        return "0" if key == "C0" else key[1:].lower()

    # ------------------------------------------------------------------
    # Command building
    # ------------------------------------------------------------------
    def build_cprogram_command(self: "BinDynSys") -> None:
        """Build the C program command list.

        Sets `self.CbaseName` and `self.cprogram` using the template
        and argument list from `_build_c_arglist()`.
        """
        suffix = self._c_program_suffix()
        self.CbaseName = self._c_program_name_template.format(suffix)
        arglist = self._build_c_arglist()
        self.cprogram = [Path(_get_ccore_bin()) / self.CbaseName] + arglist

    def _build_c_arglist(self: "BinDynSys") -> list[str]:
        """Build argument list for C program.

        Subclasses MUST override this method to provide model-specific
        command-line arguments.

        Returns
        -------
        list[str]
            Command-line arguments for the C program.
        """
        raise NotImplementedError("Subclasses must implement _build_c_arglist()")

    # ------------------------------------------------------------------
    # Stderr logging
    # ------------------------------------------------------------------
    def setup_stderr_logging(self: "BinDynSys") -> None:
        """Set up stderr logging file for C backend.

        Creates a log file in LRGSG_LOG directory with format:
        err{CbaseName}_{N}_{run_id}_{out_suffix}.log
        """
        fname = join_non_empty(
            '_',
            f"err{self.CbaseName}",
            f"{self.N}",
            getattr(self, 'run_id', ''),
            getattr(self, 'out_suffix', ''),
        ) + LOG
        self.stderr_path = Path(LRGSG_LOG) / fname
        self.stderr_path.parent.mkdir(parents=True, exist_ok=True)
        self.stderr_fopen = open(self.stderr_path, 'w')

    def _close_stderr_handle(self: "BinDynSys") -> None:
        """Close stderr file handle if open."""
        if self.stderr_fopen is not None and not self.stderr_fopen.closed:
            self.stderr_fopen.close()
        self.stderr_fopen = None

    # ------------------------------------------------------------------
    # Program execution
    # ------------------------------------------------------------------
    def run_cprogram(self: "BinDynSys", verbose: bool = False) -> None:
        """Execute the C backend and parse output state.

        Runs the C program via subprocess, captures stdout as the final
        state configuration (as int8 array), and stores it in `self.s`.

        Parameters
        ----------
        verbose : bool, optional
            If True, may print additional information (default: False).

        Raises
        ------
        RuntimeError
            If C program command hasn't been initialized.
        FileNotFoundError
            If the C binary doesn't exist.
        RuntimeError
            If the returned state size doesn't match N.
        """
        if not self.cprogram:
            raise RuntimeError("C program command has not been initialised.")

        binary_path = Path(self.cprogram[0])
        if not binary_path.exists():
            self._close_stderr_handle()
            raise FileNotFoundError(
                f"C backend executable '{binary_path}' was not found. "
                "Build the C components (e.g. via `make c-make`) before running."
            )

        result = subprocess.run(
            self.cprogram,
            stderr=self.stderr_fopen,
            stdout=subprocess.PIPE,
            check=False,
        )
        self._close_stderr_handle()

        state = np.frombuffer(result.stdout, dtype=np.int8)
        if state.size != self.N:
            raise RuntimeError(
                f"C backend returned invalid state: expected {self.N}, got {state.size}"
            )
        self.s = state.copy()

    # ------------------------------------------------------------------
    # File cleanup
    # ------------------------------------------------------------------
    def _remove_file_safe(self, path: Path | None) -> None:
        """Safely remove a file if it exists.

        Parameters
        ----------
        path : Path or None
            File path to remove. If None, does nothing.
        """
        if path is None:
            return
        try:
            if path.exists():
                path.unlink()
        except (FileNotFoundError, OSError):
            pass

    def _get_cleanup_paths(self: "BinDynSys") -> list[Path | None]:
        """Return list of paths to clean up after C run.

        Override in subclasses to add model-specific output files.
        Default implementation returns the state output file path.

        Returns
        -------
        list[Path | None]
            Paths to remove during cleanup.
        """
        return [getattr(self, 'sfout', None)]

    def _remove_stderr(self: "BinDynSys") -> None:
        """Remove the stderr log file."""
        self._remove_file_safe(self.stderr_path)

    def remove_run_c_files(self: "BinDynSys", remove_stderr: bool = True) -> None:
        """Remove exported and generated files from C run.

        Parameters
        ----------
        remove_stderr : bool, optional
            If True, also removes the stderr log file (default: True).
        """
        for path in self._get_cleanup_paths():
            self._remove_file_safe(path)
        if remove_stderr:
            self._remove_stderr()

    # ------------------------------------------------------------------
    # Utility methods
    # ------------------------------------------------------------------
    def _get_datdir_arg(self: "BinDynSys") -> str:
        """Get data directory path suitable for C program arguments.

        Attempts to make the path relative to cwd for cleaner output.

        Returns
        -------
        str
            Data directory path as string.
        """
        try:
            return str(self.sg.path_sgdata.relative_to(Path.cwd()))
        except (ValueError, AttributeError):
            # Fall back to absolute path or default
            path = getattr(self.sg, 'path_sgdata', getattr(self.sg, 'path_data', Path.cwd()))
            return str(path)

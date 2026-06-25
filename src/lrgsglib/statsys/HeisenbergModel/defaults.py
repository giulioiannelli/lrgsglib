"""Defaults and registry constants for :mod:`HeisenbergModel`."""

from __future__ import annotations

# Registry key under which HeisenbergModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``POTTS_SOLVER_NAME``).
HEISENBERG_SOLVER_NAME: str = "heisenberg"

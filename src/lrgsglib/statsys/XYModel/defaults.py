"""Defaults and registry constants for :mod:`XYModel`."""

from __future__ import annotations

# Registry key under which XYModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``POTTS_SOLVER_NAME``).
XY_SOLVER_NAME: str = "xy"

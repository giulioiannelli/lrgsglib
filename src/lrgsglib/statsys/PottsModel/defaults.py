"""Defaults and registry constants for :mod:`PottsModel`."""

from __future__ import annotations

# Registry key under which PottsModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``VOTER_SOLVER_NAME``).
POTTS_SOLVER_NAME: str = "potts"

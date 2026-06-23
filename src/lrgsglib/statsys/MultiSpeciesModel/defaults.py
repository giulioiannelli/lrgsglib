"""Defaults and registry constants for :mod:`MultiSpeciesModel`."""

from __future__ import annotations

# Registry key under which MultiSpeciesModel's solver backends are registered
# in ``statsys._solver_engine`` (mirrors ``POTTS_SOLVER_NAME``).
MULTISPEC_SOLVER_NAME: str = "multispec"

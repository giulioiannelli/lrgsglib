"""Defaults and registry constants for :mod:`ReactionDiffusionModel`."""

from __future__ import annotations

# Registry key under which ReactionDiffusionModel's solver backends are
# registered in ``statsys._solver_engine`` (mirrors ``KURAMOTO_SOLVER_NAME``).
RD_SOLVER_NAME: str = "rd"

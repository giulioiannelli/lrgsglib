"""Defaults and registry constants for :mod:`CoupledODEModel`."""

from __future__ import annotations

# Registry key under which CoupledODEModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``KURAMOTO_SOLVER_NAME``).
CODE_SOLVER_NAME: str = "code"

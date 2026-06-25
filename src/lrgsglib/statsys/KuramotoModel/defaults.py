"""Defaults and registry constants for :mod:`KuramotoModel`."""

from __future__ import annotations

# Registry key under which KuramotoModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``VOTER_SOLVER_NAME``).
KURAMOTO_SOLVER_NAME: str = "kuramoto"

"""Defaults and registry constants for :mod:`CoupledODEModel`."""

from __future__ import annotations

# Registry key under which CoupledODEModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``KURAMOTO_SOLVER_NAME``).
CODE_SOLVER_NAME: str = "code"

# Fallback dynamics subdir under ``<path_data>`` when the graph does not
# expose the graph-anchored ``path_coupled_ode``.
CODE_DYN_SUBDIR: str = "coupled_ode"

# Defaults elided from the run dirname when unchanged.
CODE_COUPLING_DEFAULT: str = "linear"
CODE_LOCAL_DEFAULT: str = "none"
CODE_INTEGRATOR_DEFAULT: str = "rk4"

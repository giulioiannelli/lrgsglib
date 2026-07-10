"""Defaults and registry constants for :mod:`ReactionDiffusionModel`."""

from __future__ import annotations

# Registry key under which ReactionDiffusionModel's solver backends are
# registered in ``statsys._solver_engine`` (mirrors ``KURAMOTO_SOLVER_NAME``).
RD_SOLVER_NAME: str = "rd"

# Fallback dynamics subdir under ``<path_data>`` when the graph does not
# expose the graph-anchored ``path_reaction_diffusion``.
RD_DYN_SUBDIR: str = "reaction_diffusion"

# Defaults elided from the run dirname when unchanged.
RD_REACTION_DEFAULT: str = "fisher_kpp"
RD_INTEGRATOR_DEFAULT: str = "rk4"

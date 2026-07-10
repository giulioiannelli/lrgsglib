"""Defaults and registry constants for :mod:`KuramotoModel`."""

from __future__ import annotations

# Registry key under which KuramotoModel's solver backends are registered in
# ``statsys._solver_engine`` (mirrors ``VOTER_SOLVER_NAME``).
KURAMOTO_SOLVER_NAME: str = "kuramoto"

# Fallback dynamics subdir under ``<path_data>`` when the graph does not
# expose the graph-anchored ``path_kuramoto`` (mirrors POTTS_DYN_SUBDIR).
KURAMOTO_DYN_SUBDIR: str = "kuramoto"

# Default integrator (elided from the run dirname when unchanged).
KURAMOTO_INTEGRATOR_DEFAULT: str = "rk4"

# Basename of the order-parameter series written into the per-run dir
# (``r.bin``, float64) when ``save_order_param=True``.
KURAMOTO_ORDER_FBASE: str = "r"

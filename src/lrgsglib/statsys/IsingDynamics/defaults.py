"""Model-level defaults for IsingDynamics.

Single-sources the registry key under which IsingDynamics' solver backends are
registered in ``statsys._solver_engine`` so the model and its solvers agree.
"""

from __future__ import annotations

# Registry key under which IsingDynamics' solver backends (py/pb/cu/C) are
# registered in ``statsys._solver_engine``.
ISING_SOLVER_NAME: str = "ising"

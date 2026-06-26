"""Model-level defaults for IsingDynamics.

Single-sources the registry key under which IsingDynamics' solver backends are
registered in ``statsys._solver_engine`` so the model and its solvers agree.
"""

from __future__ import annotations

# Registry key under which IsingDynamics' solver backends (py/pb/cu/C) are
# registered in ``statsys._solver_engine``.
ISING_SOLVER_NAME: str = "ising"

# Probability of accepting a tie (ΔE == 0) spin flip in the Metropolis update.
# 1.0 is standard Metropolis (a tie is always accepted, exp(-0/T)=1); 0.0 freezes
# ties, so at zero temperature a configuration with no energy-lowering flip is
# absorbing (frozen); 0.5 is the T=0 Glauber rule. Default 1.0 preserves the
# historical behaviour.
ISING_TIE_FLIP_DEFAULT: float = 1.0

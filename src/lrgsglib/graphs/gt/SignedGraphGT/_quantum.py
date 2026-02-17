"""Quantum propagator — re-exports from shared location.

The canonical implementation lives in ``graphs._shared._quantum``.
"""

from ..._shared._quantum import (  # noqa: F401
    compute_quantum_propagator,
    quantum_walk_probabilities,
    quantum_observables_time_series,
)

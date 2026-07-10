"""Defaults for :mod:`ContDynSys` — the continuous-flow base class.

Output-file basenames for the Phase-C per-run directory layout shared by
the flow trio (Kuramoto / ReactionDiffusion / CoupledODE)::

    <dynpath>/<run dirname>/{cfg.json, sfin.bin[, sout.bin, r.bin]}
"""

from __future__ import annotations

#: Trajectory rows (``(steps, N)`` float64), written when ``savedyn``.
CONTDYN_SOUT_FBASE: str = "sout"
#: Final state (``(N,)`` float64), written whenever ``savedisk``.
CONTDYN_SFIN_FBASE: str = "sfin"

"""Single front door for the Ising scheme classes.

``IsingModel(sg, scheme="pt", **params)`` mirrors the ``ContactProcess(...)``
factory: one callable entry point that dispatches on ``scheme`` to the
scheme leaves (``IsingMetropolis`` / ``IsingSimulatedAnnealing`` /
``IsingParallelTempering`` / ``IsingCEM``). The ``scheme`` vocabulary is
the same short form used by the ``sch=`` filename token ('met', 'sa',
'pt', 'cem'), with the long spellings accepted as aliases.

Unlike the ContactProcess factory this dispatch is EXPLICIT — the scheme
is named, never inferred from which kwargs happen to be present — so a
typo'd parameter can never silently select a different sampler.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from .defaults import ISING_SCHEME_ALIASES, ISING_SCHEME_DEFAULT
from .IsingCEM import IsingCEM
from .IsingMetropolis import IsingMetropolis
from .IsingParallelTempering import IsingParallelTempering
from .IsingSimulatedAnnealing import IsingSimulatedAnnealing

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph
    from .IsingBase import IsingBase

__all__ = ["IsingModel"]

_SCHEME_CLASSES = {
    "metropolis": IsingMetropolis,
    "sa": IsingSimulatedAnnealing,
    "pt": IsingParallelTempering,
    "cem": IsingCEM,
}


def IsingModel(
    sg: "SignedGraph",
    scheme: str = ISING_SCHEME_DEFAULT,
    **params: Any,
) -> "IsingBase":
    """Build the Ising scheme class selected by ``scheme``.

    Parameters
    ----------
    sg : SignedGraph
        The substrate graph (couplings are read from its signed weights).
    scheme : str
        Which sampler: ``'metropolis'``/``'met'`` (default), ``'sa'``,
        ``'pt'`` or ``'cem'`` (long spellings accepted, case-insensitive).
    **params
        Forwarded verbatim to the selected class — e.g. ``T=2.27`` for
        Metropolis, ``T_init=10, T_final=0.01`` for SA,
        ``n_replicas=8, T_min=1.0, T_max=4.0`` for PT.

    Returns
    -------
    IsingBase
        The configured scheme instance; call ``.run()`` on it.
    """
    key = ISING_SCHEME_ALIASES.get(str(scheme).lower())
    if key is None:
        raise ValueError(
            f"unknown Ising scheme {scheme!r}; expected one of "
            f"{sorted(set(ISING_SCHEME_ALIASES))}"
        )
    return _SCHEME_CLASSES[key](sg, **params)

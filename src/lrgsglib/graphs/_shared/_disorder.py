"""Engine-neutral disorder specification for signed graphs.

A *disorder realization* assigns weights to a subset of a graph's edges. It
factors into two orthogonal axes:

- **support** -- *which* edges are affected. The ``nwDict`` vocabulary
  (``rand``, ``randXERR``, ``randZERR``, ``single``, ``singleXERR``,
  ``singleZERR``) plus ``all`` (every edge) and ``none`` / ``None`` (defer).
  Structured supports drive ``nwDict`` construction.
- **coupling law** -- *what value* those edges receive. ``flip`` (the classic
  ``w -> -w`` sign flip, the default) or a distribution
  (``gaussian`` / ``uniform`` / ``powerlaw`` / a callable) for continuous
  couplings (SK / spin-glass-style).

The :class:`Disorder` dataclass bundles ``support`` x ``law`` x ``params`` and
travels on the ``SignedGraph`` object (``sg.disorder``). Both engine bases call
:meth:`Disorder.draw` for distributional laws and special-case ``flip``; the
edge-resolution and weight-writing live on the engine bases (NX ``weight``
float attribute / GT ``sign`` edge property) since they are engine-specific.

The shape mirrors the other ``_shared`` helpers (``_nw_container``, ``_spec``):
one engine-neutral module imported by both backends.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Optional, Union

import numpy as np

from ...config.const import SG_COUPLING, SG_DISORDER, SG_PFLIP

__all__ = [
    "Disorder",
    "as_disorder",
    "register_coupling",
    "structured_build_flags",
    "STRUCTURED_SUPPORTS",
    "FLIP",
]

#: literal name of the discrete sign-flip law (``w -> -w``); handled specially
#: by the engine bases (it is a negation, not a fresh draw).
FLIP = "flip"

#: structured supports route through ``nwDict``; the rest (``rand``/``all``)
#: are resolved directly on the graph.
STRUCTURED_SUPPORTS = frozenset(
    {"randXERR", "randZERR", "single", "singleXERR", "singleZERR"}
)

#: which ``NwContainer`` build flags a structured support requires. ``randXERR``
#: is always built (no flag); the central ``single*`` patterns need a central
#: edge (``build_single``), and any ``*ZERR`` cell pattern needs ``build_zerr``.
#: Used by the engine bases to force a graph that defaults to the geometry-free
#: container (ER/BA/…) to build the requested pattern on demand, so structured
#: supports work on *every* subclass, not only lattices.
_SUPPORT_BUILD_FLAGS: Dict[str, tuple] = {
    "single": (True, False),
    "singleXERR": (True, False),
    "singleZERR": (True, True),
    "randZERR": (False, True),
    "randXERR": (False, False),
}


def structured_build_flags(support: Optional[str]) -> tuple:
    """``(build_single, build_zerr)`` a structured ``support`` needs built.

    Returns ``(False, False)`` for non-structured supports (``rand``/``all``/
    ``none``), so callers can pass the result unconditionally.
    """
    return _SUPPORT_BUILD_FLAGS.get(support or "", (False, False))

#: registry of distributional coupling laws: ``name -> (n, rng, **params) -> ndarray``.
#: ``rng`` is duck-typed -- either ``np.random`` (NX, global seed) or a
#: ``np.random.Generator`` (GT); both expose ``normal``/``uniform``/``pareto``.
_COUPLING_LAWS: Dict[str, Callable[..., np.ndarray]] = {}


def register_coupling(name: str) -> Callable[[Callable], Callable]:
    """Register a named coupling law. The function signature is
    ``fn(n: int, rng, **params) -> array-like of length n``."""

    def _deco(fn: Callable[..., np.ndarray]) -> Callable[..., np.ndarray]:
        _COUPLING_LAWS[name] = fn
        return fn

    return _deco


@register_coupling("gaussian")
def _gaussian(n: int, rng, mu: float = 0.0, sigma: float = 1.0) -> np.ndarray:
    """i.i.d. normal couplings ``N(mu, sigma)`` (sign included)."""
    return rng.normal(mu, sigma, n)


@register_coupling("uniform")
def _uniform(n: int, rng, low: float = -1.0, high: float = 1.0) -> np.ndarray:
    """i.i.d. uniform couplings on ``[low, high]``."""
    return rng.uniform(low, high, n)


@register_coupling("powerlaw")
def _powerlaw(
    n: int, rng, exponent: float = 2.5, signed: bool = True, xmin: float = 1.0
) -> np.ndarray:
    """Heavy-tailed couplings: magnitude ~ Pareto(exponent) scaled by ``xmin``,
    with a random sign when ``signed`` (default). ``exponent`` is the Pareto
    shape ``a`` (tail index)."""
    mag = (rng.pareto(exponent, n) + 1.0) * xmin
    if signed:
        sign = np.where(rng.uniform(0.0, 1.0, n) < 0.5, -1.0, 1.0)
        return mag * sign
    return mag


@dataclass
class Disorder:
    """Engine-neutral disorder spec: ``support`` x coupling ``law`` x ``params``.

    Parameters
    ----------
    support : str
        Which edges are affected -- an ``nwDict`` key, ``'all'``, or
        ``'none'``. Default ``'rand'`` (a random ``pflip`` fraction).
    law : str or callable
        Coupling law. ``'flip'`` (default) negates the sign (``w -> -w``); a
        registered name (``'gaussian'``/``'uniform'``/``'powerlaw'``) or a
        callable ``(n, rng, **params) -> array`` draws continuous couplings.
    pflip : float
        Fraction / intensity for random supports (ignored by ``'all'``).
    params : dict
        Keyword arguments forwarded to the coupling law.
    seed : int, optional
        Reserved for per-disorder reseeding; by default the graph's RNG is used.
    """

    support: str = SG_DISORDER
    law: Union[str, Callable[..., np.ndarray]] = SG_COUPLING
    pflip: float = SG_PFLIP
    params: Dict[str, Any] = field(default_factory=dict)
    seed: Optional[int] = None

    def __post_init__(self) -> None:
        if isinstance(self.law, str) and self.law != FLIP and self.law not in _COUPLING_LAWS:
            raise ValueError(
                f"unknown coupling law {self.law!r}; known: "
                f"{[FLIP, *sorted(_COUPLING_LAWS)]} or a callable"
            )

    @property
    def is_flip(self) -> bool:
        """True for the discrete sign-flip law (negation, not a fresh draw)."""
        return self.law == FLIP

    @property
    def is_structured(self) -> bool:
        """True when the support routes through ``nwDict``."""
        return self.support in STRUCTURED_SUPPORTS

    @property
    def is_none(self) -> bool:
        """True when this spec applies nothing (defer)."""
        return self.support in (None, "none", "")

    def draw(self, n: int, rng) -> np.ndarray:
        """Return ``n`` continuous coupling values for distributional laws.

        Not valid for ``'flip'`` (a negation, handled by the engine base).
        """
        if self.is_flip:
            raise ValueError("Disorder.draw() is not defined for the 'flip' law")
        fn = self.law if callable(self.law) else _COUPLING_LAWS[self.law]
        return np.asarray(fn(n, rng, **self.params), dtype=float)


def as_disorder(
    disorder: Optional[Union[str, Disorder]], pflip: float = SG_PFLIP
) -> Optional[Disorder]:
    """Coerce the ``disorder=`` constructor argument to a :class:`Disorder` or ``None``.

    - ``None`` / ``'none'`` / ``''`` -> ``None`` (defer; realize signs manually).
    - ``str`` -> ``Disorder(support=str, law='flip', pflip=pflip)`` (top-level
      ``pflip`` feeds the string/default path).
    - :class:`Disorder` -> used as-is (its own ``pflip`` wins).
    """
    if disorder is None:
        return None
    if isinstance(disorder, Disorder):
        return None if disorder.is_none else disorder
    if isinstance(disorder, str):
        if disorder.lower() in ("none", ""):
            return None
        return Disorder(support=disorder, law=SG_COUPLING, pflip=pflip)
    raise TypeError(
        f"disorder must be None, a str, or a Disorder; got {type(disorder).__name__}"
    )

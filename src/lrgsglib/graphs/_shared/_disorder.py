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

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Union

import numpy as np

from ...config.const import SG_COUPLING, SG_DISORDER, SG_PFLIP

__all__ = [
    "Disorder",
    "CompositeDisorder",
    "as_disorder",
    "register_coupling",
    "register_support",
    "registered_supports",
    "structured_build_flags",
    "STRUCTURED_SUPPORTS",
    "COMPOSITION_MODES",
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
    #: extra kwargs forwarded to a registered-support builder (kept separate
    #: from ``params``, which feeds the coupling law).
    support_params: Dict[str, Any] = field(default_factory=dict)

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

    # --- registered-support + composition surface -------------------------
    @property
    def is_registered_support(self) -> bool:
        """True when ``support`` names a :func:`register_support` builder."""
        return isinstance(self.support, str) and self.support in _SUPPORT_BUILDERS

    def build_support(self, sg, rng, on_g) -> list:
        """Resolve a registered support to its ``(u, v)`` edge tuples."""
        fn = _SUPPORT_BUILDERS[self.support]
        return list(fn(sg, self.pflip, rng, on_g, **self.support_params))

    def build_flags(self) -> tuple:
        """``(build_single, build_zerr)`` this support needs from the container."""
        return structured_build_flags(self.support)

    # Operator sugar: each yields a CompositeDisorder (see as_disorder / below).
    def __or__(self, other):
        return _compose(self, other, "union")

    def __and__(self, other):
        return _compose(self, other, "intersection")

    def __xor__(self, other):
        return _compose(self, other, "symmetric_difference")

    def __add__(self, other):
        return _compose(self, other, "overlay")

    def __sub__(self, other):
        return CompositeDisorder(
            [self, _coerce_component(other, self.pflip)], mode="difference"
        )


def as_disorder(disorder: Any, pflip: float = SG_PFLIP) -> Any:
    """Coerce the ``disorder=`` constructor argument to a spec or ``None``.

    - ``None`` / ``'none'`` / ``''`` -> ``None`` (defer; realize signs manually).
    - ``str`` -> ``Disorder(support=str, law='flip', pflip=pflip)`` (top-level
      ``pflip`` feeds the string/default path).
    - :class:`Disorder` / :class:`CompositeDisorder` -> used as-is.
    - ``list`` -> a :class:`CompositeDisorder`: ``[A, B]`` overlays; a list of
      ``(component, weight)`` pairs is a ``mixture``.
    - ``dict`` ``{support: weight}`` -> a ``mixture`` :class:`CompositeDisorder`
      (the ``{'randXERR': 0.5, 'randZERR': 0.5}`` "half-and-half" shorthand).
    """
    if disorder is None:
        return None
    if isinstance(disorder, CompositeDisorder):
        return None if disorder.is_none else disorder
    if isinstance(disorder, Disorder):
        return None if disorder.is_none else disorder
    if isinstance(disorder, str):
        if disorder.lower() in ("none", ""):
            return None
        return Disorder(support=disorder, law=SG_COUPLING, pflip=pflip)
    if isinstance(disorder, dict):
        comps = [_coerce_component(k, pflip) for k in disorder]
        weights = [float(w) for w in disorder.values()]
        return _mixture_with_budget(comps, weights, pflip)
    if isinstance(disorder, (list, tuple)):
        pairs = len(disorder) > 0 and all(
            isinstance(x, (list, tuple)) and len(x) == 2 for x in disorder
        )
        if pairs:
            comps = [_coerce_component(c, pflip) for c, _ in disorder]
            weights = [float(w) for _, w in disorder]
            return _mixture_with_budget(comps, weights, pflip)
        comps = [_coerce_component(c, pflip) for c in disorder]
        return CompositeDisorder(comps, mode="overlay")
    raise TypeError(
        "disorder must be None, a str, Disorder, CompositeDisorder, list, or "
        f"dict; got {type(disorder).__name__}"
    )


# =====================================================================
# Registered supports (§6.4): a third support class that resolves straight
# to edges and bypasses ``nwDict`` (registered names are deliberately NOT in
# STRUCTURED_SUPPORTS, so ``is_structured`` stays False and the container is
# never built). Mirrors ``register_coupling``.
# =====================================================================

#: registry of support builders: name -> fn(sg, pflip, rng, on_g, **params) -> edges
_SUPPORT_BUILDERS: Dict[str, Callable[..., Any]] = {}


def register_support(name: str) -> Callable[[Callable], Callable]:
    """Register a named disorder *support* builder.

    The function signature is ``fn(sg, pflip, rng, on_g, **params) -> iterable of
    (u, v) tuples`` in representation ``on_g``. Registered supports resolve
    directly to edges (bypassing ``nwDict``), so they compose with any coupling
    law and need no ``NwContainer`` machinery. Draw all randomness from ``rng``
    for reproducibility, and sort node/edge inputs before sampling so the two
    engines agree for a fixed seed.
    """

    def _deco(fn: Callable) -> Callable:
        _SUPPORT_BUILDERS[name] = fn
        return fn

    return _deco


def registered_supports() -> frozenset:
    """The set of currently-registered support names."""
    return frozenset(_SUPPORT_BUILDERS)


@register_support("hubXERR")
def _hub_xerr(sg, pflip, rng, on_g, k: int = 1):
    """X-defect (star) on the ``k`` highest-degree hubs — every incident edge.

    ``pflip`` is ignored: the hubs are chosen deterministically by degree. Pass
    ``support_params={'k': K}`` to flip the top-``K`` hubs' stars.
    """
    from ._nw_geometry import star_edges

    def _deg(n):
        return len(list(sg.get_graph_neighbors(n, on_g)))

    nodes = sorted(sg.get_nodes_list(on_g), key=lambda n: (_deg(n), str(n)), reverse=True)
    keys = set()
    for hub in nodes[: max(1, int(k))]:
        for (u, v) in star_edges(sg, hub, on_g):
            iu, iv = int(u), int(v)
            keys.add((iu, iv) if iu <= iv else (iv, iu))
    return list(keys)


# =====================================================================
# Composition (§6.8): combine disorders — set-algebra, seed-mixtures, layering.
# =====================================================================

_SET_MODES = frozenset({"union", "intersection", "difference", "symmetric_difference"})
#: composition modes understood by :class:`CompositeDisorder`.
COMPOSITION_MODES = frozenset({"overlay", "mixture"} | set(_SET_MODES))


def _coerce_component(c: Any, pflip: float) -> "Disorder":
    """A composite component: a :class:`Disorder` as-is, or a support string."""
    if isinstance(c, Disorder):
        return c
    if isinstance(c, str):
        return Disorder(support=c, law=SG_COUPLING, pflip=pflip)
    raise TypeError(
        f"composite component must be a str or Disorder; got {type(c).__name__}"
    )


def _mixture_with_budget(comps, weights, pflip):
    """A ``mixture`` composite whose total budget is the top-level ``pflip``
    (so ``{'randXERR': .5, 'randZERR': .5}`` at ``pflip=0.1`` splits ~0.1·N seeds)."""
    comp = CompositeDisorder(comps, mode="mixture", weights=weights)
    if pflip and pflip > 0:
        comp._pflip = float(pflip)
    return comp


def _compose(a: Any, b: Any, mode: str) -> "CompositeDisorder":
    """Combine two operands into a :class:`CompositeDisorder`, flattening
    same-mode unweighted composites so ``A | B | C`` stays flat."""
    default_pflip = next(
        (x.pflip for x in (a, b) if isinstance(x, Disorder)), SG_PFLIP
    )
    parts: List[Disorder] = []
    for x in (a, b):
        if isinstance(x, CompositeDisorder):
            if x.mode == mode and x.weights is None:
                parts.extend(x.components)
            else:
                raise NotImplementedError(
                    "cannot combine composites of differing modes via operators; "
                    "build CompositeDisorder(...) explicitly"
                )
        else:
            parts.append(_coerce_component(x, default_pflip))
    return CompositeDisorder(parts, mode=mode)


def _largest_remainder(total: int, weights: List[float]) -> List[int]:
    """Split ``total`` into per-weight integer counts summing exactly to it."""
    raw = [total * w for w in weights]
    base = [int(math.floor(x)) for x in raw]
    rem = total - sum(base)
    order = sorted(range(len(weights)), key=lambda i: raw[i] - base[i], reverse=True)
    for i in order[: max(0, rem)]:
        base[i] += 1
    return base


@dataclass
class CompositeDisorder:
    """Combination of several :class:`Disorder` specs under one ``mode``.

    Modes: ``overlay`` (apply each component in turn — idempotent for ``flip``),
    ``mixture`` (partition the seed budget by ``weights`` — the "half-X/half-Z"
    case; v1 restricted to the ``randXERR``/``randZERR`` family), and the
    edge-set algebra ``union`` / ``intersection`` / ``difference`` /
    ``symmetric_difference``. Duck-types :class:`Disorder`'s ``is_*`` / ``pflip``
    surface so the constructor sites treat it uniformly.
    """

    components: List["Disorder"]
    mode: str = "overlay"
    weights: Optional[List[float]] = None
    seed: Optional[int] = None

    def __post_init__(self) -> None:
        self.components = [
            c for c in self.components if isinstance(c, Disorder) and not c.is_none
        ]
        if not self.components:
            raise ValueError("CompositeDisorder needs a non-empty component")
        if self.mode not in COMPOSITION_MODES:
            raise ValueError(
                f"unknown composition mode {self.mode!r}; known: "
                f"{sorted(COMPOSITION_MODES)}"
            )
        if self.mode == "difference" and len(self.components) != 2:
            raise ValueError("'difference' is binary; give exactly 2 components")
        if self.weights is not None:
            if len(self.weights) != len(self.components):
                raise ValueError("weights length must match components")
            s = float(sum(self.weights))
            if s <= 0:
                raise ValueError("weights must sum to a positive value")
            self.weights = [w / s for w in self.weights]

    # --- duck-typed Disorder surface (so constructor sites need no branching) ---
    @property
    def is_none(self) -> bool:
        return not self.components

    @property
    def is_structured(self) -> bool:
        return any(c.is_structured for c in self.components)

    @property
    def is_registered_support(self) -> bool:
        return any(c.is_registered_support for c in self.components)

    @property
    def is_flip(self) -> bool:
        return all(c.is_flip for c in self.components)

    @property
    def pflip(self) -> float:
        # A ``mixture`` budget set from the top-level ``pflip`` (via as_disorder)
        # wins; otherwise fall back to the largest component intensity.
        override = getattr(self, "_pflip", None)
        if override is not None:
            return float(override)
        return max((c.pflip for c in self.components), default=0.0)

    @property
    def support(self) -> str:
        return "composite"

    def build_flags(self) -> tuple:
        """OR the container build flags across all components."""
        bs = bz = False
        for c in self.components:
            s, z = structured_build_flags(c.support)
            bs, bz = bs or s, bz or z
        return bs, bz

    # Operator sugar so chains like ``A | B | C`` stay flat.
    def __or__(self, other):
        return _compose(self, other, "union")

    def __and__(self, other):
        return _compose(self, other, "intersection")

    def __xor__(self, other):
        return _compose(self, other, "symmetric_difference")

    def __add__(self, other):
        return _compose(self, other, "overlay")


def plan_composite_ops(comp: "CompositeDisorder", sg, on_g, resolve_keys) -> list:
    """Reduce a non-overlay composite to ``[(edge_key_set, law_disorder), ...]``.

    ``edge_key_set`` is a set of canonical ``(min, max)`` int edge tuples;
    ``law_disorder`` is the component whose coupling law writes them.
    ``resolve_keys(d)`` is the engine's support->key-set resolver. ``overlay`` is
    handled by the engine directly (apply each component) and never reaches here.
    """
    comps = comp.components
    if comp.mode in _SET_MODES:
        sets = [set(resolve_keys(c)) for c in comps]
        combined = set(sets[0])
        for s in sets[1:]:
            if comp.mode == "union":
                combined |= s
            elif comp.mode == "intersection":
                combined &= s
            elif comp.mode == "difference":
                combined -= s
            else:  # symmetric_difference
                combined ^= s
        return [(combined, comps[0])] if combined else []
    if comp.mode == "mixture":
        return _plan_mixture(comp, sg, on_g)
    return []


def _plan_mixture(comp: "CompositeDisorder", sg, on_g) -> list:
    """Partition the structured seed budget by ``weights`` (half-X / half-Z)."""
    from ._nw_geometry import star_edges

    comps = comp.components
    supports = [c.support for c in comps]
    randfam = {"randXERR", "randZERR"}
    if not all(s in randfam for s in supports):
        raise ValueError(
            "mixture (v1) supports only the structured rand-family "
            f"{sorted(randfam)} (e.g. half randXERR / half randZERR); got "
            f"{supports}. Use mode='overlay' to combine other supports."
        )
    weights = comp.weights or [1.0 / len(comps)] * len(comps)
    seed = comp.seed if comp.seed is not None else getattr(sg, "_seed", None)
    rng = np.random.default_rng(seed)
    nodes = sorted(sg.get_nodes_list(on_g), key=str)
    nflip = int(getattr(sg, "nflip", 0) or 0)
    if nflip <= 0 or not nodes:
        return []
    nflip = min(nflip, len(nodes))
    perm = rng.permutation(len(nodes))
    seeds = [nodes[int(i)] for i in perm[:nflip]]
    counts = _largest_remainder(len(seeds), weights)
    ops: list = []
    idx = 0
    for c, k in zip(comps, counts):
        chunk = seeds[idx : idx + k]
        idx += k
        keys = set()
        for s in chunk:
            pat = (
                star_edges(sg, s, on_g)
                if c.support == "randXERR"
                else sg.cell_edges(s, on_g)
            )
            for (u, v) in pat:
                iu, iv = int(u), int(v)
                keys.add((iu, iv) if iu <= iv else (iv, iu))
        if keys:
            ops.append((keys, c))
    return ops

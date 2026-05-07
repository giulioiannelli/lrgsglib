"""Parametric substrate construction for protein TMD experiments.

A *substrate* here is a network whose signed Laplacian eigenvectors form the
basis ``Φ`` against which a protein coordinate vector is decomposed in the
sense of Iannelli & Villegas 2025 (``z = Φᵀx``). This module collapses the
three near-duplicate basis builders that previously lived in ``spectral.py``
into a single dataclass-driven primitive.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Literal, Optional, Tuple

import numpy as np

from ....config import (
    DEFAULT_BA_M,
    DEFAULT_ER_DEGREE,
    DEFAULT_RGG_RADIUS_3D,
    DEFAULT_SBM_BLOCKS,
    DEFAULT_SBM_PIN,
    DEFAULT_SBM_POUT,
    DEFAULT_WS_K,
    DEFAULT_WS_P,
)
from ....graphs import (
    BarabasiAlbert,
    ErdosRenyi,
    Lattice2D,
    Lattice3D,
    RandomGeometric,
    StochasticBlockModel,
    WattsStrogatz,
)

__all__ = [
    "SubstrateSpec",
    "substrate_signature",
    "build_substrate_basis",
    "spec_label",
    "DEFAULT_SUBSTRATES",
]

SubstrateFamily = Literal[
    "lattice2d", "lattice3d", "er", "ba", "ws", "sbm", "rgg"
]


@dataclass(frozen=True)
class SubstrateSpec:
    """Frozen description of a substrate family + its identity-defining knobs.

    The same ``(family, geo, pflip, extra_kw)`` always produces the same
    ``substrate_signature`` regardless of the protein it is sized for.
    The ``seed`` only matters for stochastic families and is hashed into the
    signature alongside the other fields.
    """

    family: SubstrateFamily
    geo: Optional[str] = None  # 'sqr','tri','hex' / 'sc','bcc','fcc'
    pflip: float = 0.0
    engine: Literal["nx", "gt"] = "nx"
    extra_kw: Tuple[Tuple[str, Any], ...] = field(default_factory=tuple)
    seed: Optional[int] = None

    def kw(self) -> Dict[str, Any]:
        return dict(self.extra_kw)

    def to_dict(self) -> Dict[str, Any]:
        d = asdict(self)
        d["extra_kw"] = list(self.extra_kw)
        return d


# ---------------------------------------------------------------------------
# Defaults — the 10 substrates from the CHL-SF04 plan
# ---------------------------------------------------------------------------

DEFAULT_SUBSTRATES: Tuple[SubstrateSpec, ...] = (
    SubstrateSpec(family="lattice2d", geo="sqr"),
    SubstrateSpec(family="lattice2d", geo="tri", pflip=1.0),
    SubstrateSpec(family="lattice2d", geo="hex"),
    SubstrateSpec(family="lattice3d", geo="sc"),
    SubstrateSpec(family="lattice3d", geo="bcc"),
    SubstrateSpec(family="lattice3d", geo="fcc"),
    SubstrateSpec(family="er"),
    SubstrateSpec(family="rgg"),
    SubstrateSpec(family="ba"),
    SubstrateSpec(family="sbm"),
)


# ---------------------------------------------------------------------------
# Signatures
# ---------------------------------------------------------------------------

def substrate_signature(spec: SubstrateSpec, n_target: int) -> str:
    """Stable 12-char SHA1 of a (spec, n_target) pair.

    The signature pins identity for on-disk caches: identical ``(spec,
    n_target)`` always produce the same basis up to seed-driven randomness.
    """
    payload = json.dumps(
        {**spec.to_dict(), "n_target": int(n_target)},
        sort_keys=True,
        default=str,
    )
    return hashlib.sha1(payload.encode()).hexdigest()[:12]


def spec_label(spec: SubstrateSpec) -> str:
    """Short human-readable label for figures / file names."""
    base = spec.family if spec.geo is None else f"{spec.family}_{spec.geo}"
    if spec.pflip > 0:
        base += f"_pf{spec.pflip:g}"
    return base


# ---------------------------------------------------------------------------
# Builders — one helper per family
# ---------------------------------------------------------------------------


def _build_lattice2d(spec: SubstrateSpec, n_target: int):
    geo = spec.geo or "sqr"
    # Hex lattices return roughly N ≈ side·(side-c) nodes; oversize.
    factor = 2.0 if geo == "hex" else 1.0
    side = max(4, int(math.ceil(math.sqrt(factor * n_target))))
    # Hex with PBC requires an even side (see L2D_ERRMSG_GEO).
    if geo == "hex" and side % 2 == 1:
        side += 1
    kw = dict(
        side1=side,
        geo=geo,
        pflip=spec.pflip,
        engine=spec.engine,
        **spec.kw(),
    )
    if spec.seed is not None:
        kw.setdefault("seed", spec.seed)
    return Lattice2D(**kw)


def _build_lattice3d(spec: SubstrateSpec, n_target: int):
    side = max(2, int(math.ceil(n_target ** (1.0 / 3.0))))
    kw = dict(
        dim=(side, side, side),
        geo=spec.geo or "sc",
        pflip=spec.pflip,
        engine=spec.engine,
        **spec.kw(),
    )
    if spec.seed is not None:
        kw.setdefault("seed", spec.seed)
    return Lattice3D(**kw)


def _build_er(spec: SubstrateSpec, n_target: int):
    extra = spec.kw()
    # Oversize to absorb the ~5–15% lost to extract_giant_component.
    n = int(extra.pop("n", int(math.ceil(1.25 * n_target))))
    k_avg = float(extra.pop("k_avg", DEFAULT_ER_DEGREE))
    p = float(extra.pop("p", k_avg / max(n - 1, 1)))
    return ErdosRenyi(
        n=n,
        p=p,
        pflip=spec.pflip,
        seed=spec.seed,
        engine=spec.engine,
        **extra,
    )


def _build_ba(spec: SubstrateSpec, n_target: int):
    extra = spec.kw()
    n = int(extra.pop("n", n_target))
    m = int(extra.pop("m", DEFAULT_BA_M))
    return BarabasiAlbert(
        n=n,
        m=m,
        pflip=spec.pflip,
        seed=spec.seed,
        engine=spec.engine,
        **extra,
    )


def _build_ws(spec: SubstrateSpec, n_target: int):
    extra = spec.kw()
    n = int(extra.pop("n", n_target))
    k = int(extra.pop("k", DEFAULT_WS_K))
    p = float(extra.pop("p", DEFAULT_WS_P))
    return WattsStrogatz(
        n=n,
        k=k,
        p=p,
        pflip=spec.pflip,
        seed=spec.seed,
        engine=spec.engine,
        **extra,
    )


def _build_sbm(spec: SubstrateSpec, n_target: int):
    extra = spec.kw()
    n_blocks = int(extra.pop("n_blocks", DEFAULT_SBM_BLOCKS))
    p_in = float(extra.pop("p_in", DEFAULT_SBM_PIN))
    p_out = float(extra.pop("p_out", DEFAULT_SBM_POUT))
    inflated = int(math.ceil(1.25 * n_target))
    block = int(math.ceil(inflated / n_blocks))
    sizes = [block] * n_blocks
    p_matrix = [
        [p_in if i == j else p_out for j in range(n_blocks)]
        for i in range(n_blocks)
    ]
    return StochasticBlockModel(
        sizes=sizes,
        p_matrix=p_matrix,
        pflip=spec.pflip,
        seed=spec.seed,
        engine=spec.engine,
        **extra,
    )


def _build_rgg(spec: SubstrateSpec, n_target: int):
    extra = spec.kw()
    # RGG with extract_giant_component drops 30–50 % at the default radius;
    # oversample and let the giant component absorb the loss.
    n = int(extra.pop("n", int(math.ceil(2.0 * n_target))))
    radius = float(extra.pop("radius", DEFAULT_RGG_RADIUS_3D))
    dim = int(extra.pop("dim", 3))
    return RandomGeometric(
        n=n,
        radius=radius,
        dim=dim,
        pflip=spec.pflip,
        seed=spec.seed,
        engine=spec.engine,
        **extra,
    )


_BUILDERS = {
    "lattice2d": _build_lattice2d,
    "lattice3d": _build_lattice3d,
    "er": _build_er,
    "ba": _build_ba,
    "ws": _build_ws,
    "sbm": _build_sbm,
    "rgg": _build_rgg,
}


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def build_substrate_basis(
    spec: SubstrateSpec,
    n_target: int,
    spectrum_backend: Optional[str] = None,
) -> Tuple[np.ndarray, Any, Dict[str, Any]]:
    """Build a substrate sized to host ``>= n_target`` nodes and return its basis.

    Parameters
    ----------
    spec :
        Frozen substrate description.
    n_target :
        Desired lower bound on the number of nodes ``N``. The resulting
        substrate is **not** truncated; the full basis ``Φ ∈ ℝ^{N×N}`` is
        returned. Ensure ``N >= n_target`` so a coordinate signal of size
        ``n_target`` can be zero-padded into the substrate's coordinate space.
    spectrum_backend :
        Backend forwarded to ``compute_laplacian_spectrum_weigV``. ``None`` →
        the library default (typically ``'cupy'`` if available else
        ``'scipy'``).

    Returns
    -------
    basis : ndarray (N, N)
        Eigenvector matrix (columns = eigenvectors, ascending eigenvalue).
    sg : SignedGraph
        The constructed substrate object (eigenvalues stored on ``sg.eigv``).
    sizing_meta : dict
        ``{"N", "n_target", "signature", "spec", "label", "eigvals"}``.
    """
    builder = _BUILDERS[spec.family]
    sg = builder(spec, n_target)
    if sg.N < n_target:
        # One retry with a 2x inflated target — covers stochastic shortfalls.
        sg = builder(spec, int(math.ceil(2 * n_target)))
        if sg.N < n_target:
            raise RuntimeError(
                f"Substrate {spec_label(spec)} produced N={sg.N} < n_target={n_target}; "
                f"increase the family-specific oversize factor."
            )
    if spec.pflip > 0.0:
        sg.flip_random_fract_edges()
    backend_kw: Dict[str, Any] = {}
    if spectrum_backend is not None:
        backend_kw["backend"] = spectrum_backend
    sg.compute_laplacian_spectrum_weigV(**backend_kw)
    basis = sg.get_sgspect_basis(max_factor=1)
    sizing_meta = {
        "N": int(sg.N),
        "n_target": int(n_target),
        "signature": substrate_signature(spec, n_target),
        "spec": spec.to_dict(),
        "label": spec_label(spec),
        "eigvals": np.asarray(sg.eigv, dtype=np.float64),
    }
    return basis, sg, sizing_meta

"""
Unified RandomGeometric factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import RandomGeometric
>>>
>>> rgg = RandomGeometric(n=200, radius=0.2, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.RandomGeometricNX import RandomGeometricNX


def _get_nx_impl():
    from .nx.RandomGeometricNX import RandomGeometricNX

    return RandomGeometricNX


def _get_gt_impl():
    from .gt.RandomGeometricGT import RandomGeometricGT

    return RandomGeometricGT


register_implementation(
    "RandomGeometric", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "RandomGeometric", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class RandomGeometric:
    """
    Create a random geometric signed graph.

    Nodes are placed uniformly at random in a unit hypercube and
    connected if within the given radius.

    Parameters
    ----------
    n : int
        Number of nodes.
    radius : float
        Connection radius.
    dim : int, default 2
        Dimension of the embedding space.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine to use ('nx' or 'gt').
    **kwargs
        Additional engine-specific parameters.

    Returns
    -------
    SignedGraphProtocol
        A random geometric graph instance.
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX, the most complete
    # backend) so editors give full, precise method navigation -- including when
    # the call uses ``**dict`` unpacking, which defeats @overload-based typing.
    def __new__(
        cls,
        n: Any,
        radius: Any,
        dim: Any = 2,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "RandomGeometricNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("RandomGeometric", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "n": n,
                "radius": radius,
                "dim": dim,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            impl_kwargs = {
                "n": n,
                "radius": radius,
                "dim": dim,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "n": n,
                "radius": radius,
                "dim": dim,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)

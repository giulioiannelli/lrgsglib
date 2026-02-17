"""
Engine selection and registration system for multi-backend graph support.

This module provides the infrastructure for selecting between different graph
backends (NetworkX, graph-tool, igraph) and registering implementations.

Example
-------
>>> from lrgsglib.graphs import set_default_engine, get_default_engine
>>> set_default_engine('gt')  # Use graph-tool by default
>>> get_default_engine()
<GraphEngine.GRAPHTOOL: 'gt'>
"""

from __future__ import annotations

import os
import warnings
from enum import Enum
from typing import TYPE_CHECKING, Callable, Optional, Type, TypeVar, Union

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol

T = TypeVar("T")


class GraphEngine(Enum):
    """
    Enumeration of supported graph backends.

    Attributes
    ----------
    NETWORKX : str
        NetworkX backend ('nx'). Default, always available.
    GRAPHTOOL : str
        graph-tool backend ('gt'). High performance, requires installation.
    IGRAPH : str
        igraph backend ('ig'). Future support.
    """

    NETWORKX = "nx"
    GRAPHTOOL = "gt"
    IGRAPH = "ig"


# Global default engine (can be overridden)
_default_engine: GraphEngine = GraphEngine.NETWORKX

# Registry: graph_type -> {engine -> implementation_class_or_factory}
# Values can be either a class or a callable that returns the class (lazy import)
_implementations: dict[str, dict[GraphEngine, Union[Type, Callable[[], Type]]]] = {}

# Cache for resolved lazy imports
_resolved_implementations: dict[str, dict[GraphEngine, Type]] = {}


def set_default_engine(engine: Union[str, GraphEngine]) -> None:
    """
    Set the default graph engine globally.

    This affects all subsequent graph creation calls that don't explicitly
    specify an engine.

    Parameters
    ----------
    engine : str or GraphEngine
        Engine to use. Can be 'nx', 'gt', 'ig', or a GraphEngine enum value.

    Raises
    ------
    ValueError
        If the engine string is not recognized.

    Example
    -------
    >>> set_default_engine('gt')
    >>> set_default_engine(GraphEngine.NETWORKX)
    """
    global _default_engine
    if isinstance(engine, str):
        engine = GraphEngine(engine)
    _default_engine = engine


def get_default_engine() -> GraphEngine:
    """
    Get the current default graph engine.

    Returns
    -------
    GraphEngine
        The currently configured default engine.
    """
    return _default_engine


def get_default_engine_from_env() -> GraphEngine:
    """
    Get the default engine from environment variable if set.

    Checks LRGSG_GRAPH_ENGINE environment variable.

    Returns
    -------
    GraphEngine
        Engine from environment or current default.
    """
    env_engine = os.environ.get("LRGSG_GRAPH_ENGINE")
    if env_engine is not None:
        try:
            return GraphEngine(env_engine.lower())
        except ValueError:
            warnings.warn(
                f"Invalid LRGSG_GRAPH_ENGINE value: {env_engine}. "
                f"Valid options: nx, gt, ig. Using default.",
                stacklevel=2,
            )
    return _default_engine


def register_implementation(
    graph_type: str,
    engine: GraphEngine,
    cls_or_factory: Union[Type[T], Callable[[], Type[T]]],
) -> None:
    """
    Register an implementation class for a graph type and engine.

    Parameters
    ----------
    graph_type : str
        Name of the graph type (e.g., 'Lattice2D', 'ErdosRenyi').
    engine : GraphEngine
        The engine this implementation supports.
    cls_or_factory : Type or Callable
        Either the implementation class directly, or a callable that returns
        the class (for lazy imports to avoid circular dependencies).

    Example
    -------
    >>> # Direct registration
    >>> register_implementation('Lattice2D', GraphEngine.NETWORKX, Lattice2DNX)
    >>>
    >>> # Lazy registration (preferred for avoiding circular imports)
    >>> def _get_lattice2d_gt():
    ...     from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT
    ...     return Lattice2DGT
    >>> register_implementation('Lattice2D', GraphEngine.GRAPHTOOL, _get_lattice2d_gt)
    """
    if graph_type not in _implementations:
        _implementations[graph_type] = {}
    _implementations[graph_type][engine] = cls_or_factory


def _resolve_implementation(graph_type: str, engine: GraphEngine) -> Type:
    """
    Resolve a potentially lazy implementation to its actual class.

    Parameters
    ----------
    graph_type : str
        Name of the graph type.
    engine : GraphEngine
        The engine to resolve for.

    Returns
    -------
    Type
        The resolved implementation class.

    Raises
    ------
    KeyError
        If the graph type or engine is not registered.
    """
    # Check cache first
    if graph_type in _resolved_implementations:
        if engine in _resolved_implementations[graph_type]:
            return _resolved_implementations[graph_type][engine]

    # Get the registered value
    impl = _implementations[graph_type][engine]

    # Resolve if it's a factory
    if callable(impl) and not isinstance(impl, type):
        impl = impl()

    # Cache the resolved class
    if graph_type not in _resolved_implementations:
        _resolved_implementations[graph_type] = {}
    _resolved_implementations[graph_type][engine] = impl

    return impl


def get_implementation(
    graph_type: str,
    engine: Optional[Union[str, GraphEngine]] = None,
    fallback: bool = True,
) -> Type:
    """
    Get the implementation class for a graph type and engine.

    Parameters
    ----------
    graph_type : str
        Name of the graph type (e.g., 'Lattice2D', 'ErdosRenyi').
    engine : str, GraphEngine, or None
        The engine to use. If None, uses the default engine (which may
        be configured via environment variable or set_default_engine).
    fallback : bool, default True
        If True and the requested engine is not available, fall back to
        NetworkX. If False, raise an error for unavailable engines.

    Returns
    -------
    Type
        The implementation class for the specified graph type and engine.

    Raises
    ------
    ValueError
        If the graph type is unknown or no implementation is available.

    Example
    -------
    >>> cls = get_implementation('Lattice2D', engine='gt')
    >>> lattice = cls(side1=100, geo='sqr')
    """
    # Resolve engine
    if engine is None:
        engine = get_default_engine_from_env()
    elif isinstance(engine, str):
        engine = GraphEngine(engine)

    # Check if graph type is registered
    if graph_type not in _implementations:
        raise ValueError(f"Unknown graph type: {graph_type}")

    impls = _implementations[graph_type]

    # Try requested engine first
    if engine in impls:
        try:
            return _resolve_implementation(graph_type, engine)
        except ImportError as e:
            if not fallback:
                raise
            warnings.warn(
                f"Failed to import {engine.value} implementation for "
                f"{graph_type}: {e}. Falling back to NetworkX.",
                stacklevel=2,
            )

    # Fallback chain: try other engines
    if fallback:
        fallback_order = [
            GraphEngine.NETWORKX,
            GraphEngine.GRAPHTOOL,
            GraphEngine.IGRAPH,
        ]
        for fb_engine in fallback_order:
            if fb_engine in impls and fb_engine != engine:
                try:
                    return _resolve_implementation(graph_type, fb_engine)
                except ImportError:
                    continue

    raise ValueError(
        f"No implementation available for {graph_type} with engine {engine.value}"
    )


def list_registered_types() -> list[str]:
    """
    List all registered graph types.

    Returns
    -------
    list[str]
        Names of all registered graph types.
    """
    return list(_implementations.keys())


def list_engines_for_type(graph_type: str) -> list[GraphEngine]:
    """
    List available engines for a graph type.

    Parameters
    ----------
    graph_type : str
        Name of the graph type.

    Returns
    -------
    list[GraphEngine]
        Engines that have implementations registered for this type.
    """
    if graph_type not in _implementations:
        return []
    return list(_implementations[graph_type].keys())


def is_engine_available(engine: Union[str, GraphEngine]) -> bool:
    """
    Check if a graph engine is available (importable).

    Parameters
    ----------
    engine : str or GraphEngine
        Engine to check.

    Returns
    -------
    bool
        True if the engine's library can be imported.
    """
    if isinstance(engine, str):
        engine = GraphEngine(engine)

    if engine == GraphEngine.NETWORKX:
        try:
            import networkx  # noqa: F401

            return True
        except ImportError:
            return False
    elif engine == GraphEngine.GRAPHTOOL:
        try:
            import graph_tool  # noqa: F401

            return True
        except ImportError:
            return False
    elif engine == GraphEngine.IGRAPH:
        try:
            import igraph  # noqa: F401

            return True
        except ImportError:
            return False
    return False


def available_engines() -> list[GraphEngine]:
    """
    List all engines that are currently importable.

    Returns
    -------
    list[GraphEngine]
        Engines whose libraries are installed and importable.
    """
    return [e for e in GraphEngine if is_engine_available(e)]

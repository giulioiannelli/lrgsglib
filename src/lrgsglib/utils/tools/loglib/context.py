"""Session context utilities for lrgsglib logging.

Provides automatic session ID generation from graph parameters and timestamps.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    import networkx as nx


@dataclass
class SessionContext:
    """Context for generating session identifiers.

    Builds descriptive session IDs from graph parameters (N, pflip, p1, p2, etc.)
    combined with timestamps for unique, meaningful log filenames.

    Parameters
    ----------
    prefix : str
        Session prefix (e.g., "mcg", "l2d", "ising").
    params : dict[str, Any]
        Key-value pairs for parameters to include in session ID.
    include_timestamp : bool
        Whether to append timestamp to session ID (default: True).

    Examples
    --------
    >>> ctx = SessionContext(prefix="mcg", params={"N": 4096, "p1": 0.8})
    >>> ctx.generate_id()
    'mcg_N=4096_p1=0.8_20260202_143052'
    """

    prefix: str
    params: dict[str, Any] = field(default_factory=dict)
    include_timestamp: bool = True

    def generate_id(self) -> str:
        """Generate session identifier string.

        Returns
        -------
        str
            Session ID combining prefix, parameters, and optional timestamp.
        """
        parts = [self.prefix]

        # Add parameters
        for key, value in self.params.items():
            if isinstance(value, float):
                # Format floats compactly
                parts.append(f"{key}={value:.3g}")
            else:
                parts.append(f"{key}={value}")

        # Add timestamp if requested
        if self.include_timestamp:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            parts.append(timestamp)

        return "_".join(parts)

    def add_param(self, key: str, value: Any) -> "SessionContext":
        """Add a parameter to the context (returns self for chaining).

        Parameters
        ----------
        key : str
            Parameter name.
        value : Any
            Parameter value.

        Returns
        -------
        SessionContext
            Self for method chaining.
        """
        self.params[key] = value
        return self


def session_from_graph(
    graph: Any,
    prefix: str = "graph",
    include_timestamp: bool = True,
) -> SessionContext:
    """Create SessionContext from a graph object.

    Extracts relevant parameters from the graph (N, pflip, p1, p2, geo, etc.)
    based on the graph type.

    Parameters
    ----------
    graph : Any
        A SignedGraph instance or networkx graph.
    prefix : str
        Prefix for the session ID (default: "graph").
    include_timestamp : bool
        Whether to include timestamp (default: True).

    Returns
    -------
    SessionContext
        Context with extracted graph parameters.

    Examples
    --------
    >>> from lrgsglib import MultiplicativeCascadeGraph
    >>> mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6)
    >>> ctx = session_from_graph(mc, prefix="mcg")
    >>> ctx.generate_id()
    'mcg_N=4096_p1=0.8_p2=0.6_20260202_143052'
    """
    params: dict[str, Any] = {}

    # Try to get node count
    if hasattr(graph, "N"):
        params["N"] = graph.N
    elif hasattr(graph, "gr") and "G" in getattr(graph, "gr", {}):
        try:
            params["N"] = graph.gr["G"].number_of_nodes()
        except (AttributeError, TypeError):
            pass
    elif hasattr(graph, "number_of_nodes"):
        params["N"] = graph.number_of_nodes()

    # Extract common signed graph parameters
    param_attrs = [
        "pflip",  # SignedGraph frustration parameter
        "p1", "p2", "p3", "p4",  # MultispectralGraph parameters
        "fraction",  # MCG fraction
        "iterations",  # MCG iterations
        "geo",  # Lattice geometry
        "side1", "side2",  # Lattice dimensions
        "dim",  # 3D lattice dimensions
        "p", "k",  # ErdosRenyi parameters
    ]

    for attr in param_attrs:
        if hasattr(graph, attr):
            value = getattr(graph, attr)
            if value is not None:
                # Skip default/uninteresting values
                if attr == "pflip" and value == 0.0:
                    continue
                if attr == "fraction" and value == 0.4:
                    continue
                params[attr] = value

    return SessionContext(
        prefix=prefix,
        params=params,
        include_timestamp=include_timestamp,
    )


def session_from_dynamics(
    dynamics: Any,
    prefix: str = "dynamics",
    include_timestamp: bool = True,
) -> SessionContext:
    """Create SessionContext from a dynamics object.

    Extracts parameters from dynamics objects (IsingDynamics, ContactProcess, etc.)

    Parameters
    ----------
    dynamics : Any
        A dynamics instance (IsingDynamics, ContactProcessSIR, etc.).
    prefix : str
        Prefix for the session ID (default: "dynamics").
    include_timestamp : bool
        Whether to include timestamp (default: True).

    Returns
    -------
    SessionContext
        Context with extracted dynamics parameters.
    """
    params: dict[str, Any] = {}

    # Get graph info first
    if hasattr(dynamics, "sg") and dynamics.sg is not None:
        if hasattr(dynamics.sg, "N"):
            params["N"] = dynamics.sg.N

    # Dynamics-specific parameters
    dynamics_attrs = [
        "T",  # Temperature (Ising)
        "steps",  # Simulation steps
        "gamma", "mu",  # Contact process rates
        "p0",  # Initial density
        "runlang",  # Backend
    ]

    for attr in dynamics_attrs:
        if hasattr(dynamics, attr):
            value = getattr(dynamics, attr)
            if value is not None:
                params[attr] = value

    return SessionContext(
        prefix=prefix,
        params=params,
        include_timestamp=include_timestamp,
    )

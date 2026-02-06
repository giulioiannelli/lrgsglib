"""
Abstract protocol definitions for signed graph implementations.

This module defines the interface that all graph implementations must satisfy
to be interoperable across different backends (NetworkX, graph-tool, igraph).

The SignedGraphProtocol uses Python's Protocol for structural subtyping,
meaning any class that implements the required methods is considered compatible
without explicit inheritance.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    import networkx as nx


@runtime_checkable
class SignedGraphProtocol(Protocol):
    """
    Protocol that all signed graph implementations must satisfy.

    This protocol defines the minimal interface required for a graph class
    to work with the unified graphs API. Implementations can provide
    additional methods beyond this protocol.

    Required Properties
    -------------------
    N : int
        Number of nodes in the graph.
    num_edges : int
        Number of edges in the graph.
    pflip : float
        Fraction of edges that are negative (0.0 to 1.0).
    seed : int or None
        Random seed used for graph generation.

    Required Methods
    ----------------
    flip_random_fract_edges(fraction=None)
        Flip a fraction of edges to negative.
    get_signed_laplacian() -> ndarray
        Get the signed Laplacian matrix.
    get_signed_adjacency() -> ndarray
        Get the signed adjacency matrix.
    count_negative_edges() -> int
        Count edges with negative sign.
    count_positive_edges() -> int
        Count edges with positive sign.
    to_networkx() -> nx.Graph
        Convert to NetworkX graph.

    Example
    -------
    >>> def analyze_graph(sg: SignedGraphProtocol) -> float:
    ...     '''Works with any graph implementation.'''
    ...     L = sg.get_signed_laplacian()
    ...     eigenvalues = np.linalg.eigvalsh(L)
    ...     return eigenvalues[1]  # First non-zero eigenvalue
    """

    # === Required Properties ===

    @property
    def N(self) -> int:
        """Number of nodes in the graph."""
        ...

    @property
    def Ne(self) -> int:
        """Number of edges in the graph (NX naming convention)."""
        ...

    @property
    def num_edges(self) -> int:
        """Number of edges in the graph (GT naming convention, alias for Ne)."""
        ...

    @property
    def pflip(self) -> float:
        """Fraction of edges that are negative."""
        ...

    @property
    def seed(self) -> Optional[int]:
        """Random seed used for generation."""
        ...

    # === Required Methods ===

    def flip_random_fract_edges(self, fraction: Optional[float] = None) -> None:
        """
        Flip a fraction of edges to negative.

        Parameters
        ----------
        fraction : float or None
            Fraction of edges to flip. If None, uses self.pflip.
        """
        ...

    def get_signed_laplacian(self) -> NDArray[np.floating]:
        """
        Get the signed Laplacian matrix L_s = D_s - A.

        Where D_s is the diagonal matrix of absolute degrees and A is
        the signed adjacency matrix.

        Returns
        -------
        ndarray
            NxN signed Laplacian matrix.
        """
        ...

    def get_signed_adjacency(self) -> NDArray[np.floating]:
        """
        Get the signed adjacency matrix.

        Returns
        -------
        ndarray
            NxN signed adjacency matrix with +1/-1 edge weights.
        """
        ...

    def count_negative_edges(self) -> int:
        """
        Count edges with negative sign.

        Returns
        -------
        int
            Number of negative edges.
        """
        ...

    def count_positive_edges(self) -> int:
        """
        Count edges with positive sign.

        Returns
        -------
        int
            Number of positive edges.
        """
        ...

    def to_networkx(self) -> "nx.Graph":
        """
        Convert to NetworkX graph.

        Returns
        -------
        nx.Graph
            NetworkX graph with 'weight' edge attribute for signs.
        """
        ...


@runtime_checkable
class SpectralGraphProtocol(SignedGraphProtocol, Protocol):
    """
    Extended protocol for graphs with spectral computation capabilities.

    This protocol extends SignedGraphProtocol with methods for eigenvalue
    and eigenvector computations.
    """

    def compute_laplacian_spectrum(
        self,
        k: Optional[int] = None,
        backend: str = "numpy",
    ) -> NDArray[np.floating]:
        """
        Compute eigenvalues of the signed Laplacian.

        Parameters
        ----------
        k : int or None
            Number of eigenvalues to compute. If None, compute all.
        backend : str
            Computation backend: 'numpy', 'scipy', 'cupy'.

        Returns
        -------
        ndarray
            Array of eigenvalues, sorted ascending.
        """
        ...

    def compute_laplacian_spectrum_weigV(
        self,
        backend: str = "numpy",
    ) -> tuple[NDArray[np.floating], NDArray[np.floating]]:
        """
        Compute eigenvalues and eigenvectors of the signed Laplacian.

        Parameters
        ----------
        backend : str
            Computation backend: 'numpy', 'scipy', 'cupy'.

        Returns
        -------
        tuple[ndarray, ndarray]
            (eigenvalues, eigenvectors) where eigenvectors are column-major.
        """
        ...


@runtime_checkable
class LatticeGraphProtocol(SignedGraphProtocol, Protocol):
    """
    Protocol for lattice-based graphs with geometry and dimension info.

    This protocol extends SignedGraphProtocol with properties specific
    to lattice structures.
    """

    @property
    def side1(self) -> int:
        """Size in first dimension."""
        ...

    @property
    def side2(self) -> int:
        """Size in second dimension."""
        ...

    @property
    def geo(self) -> str:
        """Geometry type: 'sqr', 'tri', 'hex', etc."""
        ...

    @property
    def periodic(self) -> bool:
        """Whether periodic boundary conditions are used."""
        ...


def is_signed_graph(obj: Any) -> bool:
    """
    Check if an object satisfies the SignedGraphProtocol.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if the object implements SignedGraphProtocol.
    """
    return isinstance(obj, SignedGraphProtocol)


def is_lattice_graph(obj: Any) -> bool:
    """
    Check if an object satisfies the LatticeGraphProtocol.

    Parameters
    ----------
    obj : Any
        Object to check.

    Returns
    -------
    bool
        True if the object implements LatticeGraphProtocol.
    """
    return isinstance(obj, LatticeGraphProtocol)

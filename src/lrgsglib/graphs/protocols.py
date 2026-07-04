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
    from pathlib import Path

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

    # === Dynamics-Critical Methods ===

    def get_neighbors_with_weights(self, node: int) -> list[tuple[int, float]]:
        """Return neighbors of *node* with signed edge weights.

        Parameters
        ----------
        node : int
            Node index.

        Returns
        -------
        list[tuple[int, float]]
            List of ``(neighbor_index, sign_weight)`` pairs.
        """
        ...

    def get_edges_with_weights(self) -> list[tuple[int, int, float]]:
        """Return all edges with their signed weights.

        Returns
        -------
        list[tuple[int, int, float]]
            List of ``(u, v, weight)`` triples.
        """
        ...

    def get_neighbor_indices(self, node: int) -> list[int]:
        """Return neighbor indices without weights.

        Parameters
        ----------
        node : int
            Node index.

        Returns
        -------
        list[int]
            Neighbor node indices.
        """
        ...

    @property
    def gr(self) -> dict:
        """Engine-agnostic graph-representation dict (``gr['G']`` / ``gr['H']``)."""
        ...

    def get_node_attributes(self, attr: str = "pos", on_g: Optional[str] = None) -> dict:
        """Return a ``{node: value}`` dict for vertex attribute *attr*.

        Parameters
        ----------
        attr : str
            Node attribute / vertex-property name (default ``'pos'``).
        on_g : str, optional
            Graph representation; ignored on single-representation engines.
        """
        ...


@runtime_checkable
class DynamicsGraphProtocol(SignedGraphProtocol, Protocol):
    """Protocol for graphs used with dynamics simulations.

    Extends :class:`SignedGraphProtocol` with file I/O and path methods
    required by the C backend pipeline.
    """

    @property
    def syshapePth(self) -> str:
        """Parameter string for directory/file naming."""
        ...

    @property
    def path_sgdata(self) -> "Path":
        """Root data directory for this graph instance."""
        ...

    @property
    def path_data(self) -> "Path":
        """Base data directory the graph's whole output tree hangs off."""
        ...

    def get_p_fname(self, who: str, out_suffix: str = "", ext: str = ".bin") -> "str | Path":
        """Build a parametric filename (``<who>_p=<pflip>[_suffix]<ext>``).

        Parameters
        ----------
        who : str
            Filename prefix (e.g. ``'s'``, ``'h'``, ``'m'``).
        out_suffix : str
            Additional suffix.
        ext : str
            File extension (default ``'.bin'``).

        Returns
        -------
        str or Path
            Constructed filename.
        """
        ...

    def get_eigV_bin_check(
        self,
        which: int = 0,
        reshaped: bool = False,
        backend: Optional[str] = None,
        typf: type = np.float64,
        transpose: bool = True,
        flip_to_pos: bool = True,
    ) -> NDArray:
        """Binarized eigenvector (computed on demand); used by spectral ICs."""
        ...

    def _export_edgel_bin(self, exName: str = "") -> None:
        """Write edge list to binary file for C backend."""
        ...

    def remove_exported_files(self) -> None:
        """Remove all exported binary files."""
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

    def draw(self, ax: Optional[Any] = None, **kwargs: Any) -> Any:
        """Draw the lattice using its stored node positions.

        Engine-agnostic: identical signature and result on the NetworkX and
        graph-tool backends. Requires the lattice to have been built with
        ``with_positions=True``. See
        :func:`lrgsglib.graphs._shared._draw.draw` for the full signature.
        """
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

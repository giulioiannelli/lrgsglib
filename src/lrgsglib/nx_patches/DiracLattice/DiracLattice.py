"""Abstract base class for Dirac lattice structures."""

from __future__ import annotations

import networkx as nx

from ..MultispectralGraph import MultispectralGraph
from . import _spectral


class DiracLatticeGraph(MultispectralGraph):
    """Abstract base class for Dirac structure graphs (comb, brush).

    Dirac structures are hierarchical graphs with a base graph and fiber
    graphs attached to each base node. This structure allows efficient
    spectral computation by separating base and fiber spectra.
    """

    def __init__(self, periodic: bool = True, **kwargs):
        """Initialize Dirac lattice graph.

        Parameters
        ----------
        periodic : bool
            Use periodic boundary conditions
        **kwargs
            Additional arguments passed to MultispectralGraph
        """
        self.periodic = periodic

        # Dirac-specific attributes
        self.dirac_structure = None
        self.base_graph: nx.Graph | None = None
        self.fiber_graph: nx.Graph | None = None
        self._base_eigenvalues = None
        self._fiber_laplacian = None

        super().__init__(**kwargs)

    @property
    def base_eigenvalues(self):
        """Cached base eigenvalues."""
        return self._base_eigenvalues

    @base_eigenvalues.setter
    def base_eigenvalues(self, value):
        self._base_eigenvalues = value

    @property
    def fiber_laplacian(self):
        """Cached fiber Laplacian."""
        return self._fiber_laplacian

    @fiber_laplacian.setter
    def fiber_laplacian(self, value):
        self._fiber_laplacian = value

    def is_dirac_structure(self) -> bool:
        """Always True for DiracLatticeGraph subclasses."""
        return True


# Attach spectral methods to class
DiracLatticeGraph.compute_base_spectrum = _spectral.compute_base_spectrum
DiracLatticeGraph.compute_fiber_laplacian = _spectral.compute_fiber_laplacian
DiracLatticeGraph.compute_dirac_spectrum_separated = (
    _spectral.compute_dirac_spectrum_separated
)
DiracLatticeGraph.get_base_fiber_dimensions = _spectral.get_base_fiber_dimensions

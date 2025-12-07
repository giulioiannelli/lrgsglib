"""
Quantum extension of SignedLaplacianAnalysis for quantum RG flows.

This module provides the QuantumSignedLaplacianAnalysis class which extends
the classical SignedLaplacianAnalysis to support quantum propagators e^(-itL)
and quantum information-theoretic observables.
"""

import numpy as np
from typing import Optional
from numpy.typing import NDArray

from .quantum import (
    compute_quantum_observables_from_eigenvalues,
    quantum_classical_divergence,
)

# Import base class - will be imported from core module
# from ...core import SignedLaplacianAnalysis
# For now, we'll import the necessary components
try:
    from ...config.const import LRSG_ENTROPY_STEP, DEFAULT_MAX_THRESHOLD
    from ...nx_patches.SignedGraph import SignedGraph
except ImportError:
    # Fallback for testing
    LRSG_ENTROPY_STEP = 600
    DEFAULT_MAX_THRESHOLD = 0.1

__all__ = ["QuantumSignedLaplacianAnalysis"]


class QuantumSignedLaplacianAnalysis:
    """
    Quantum extension of SignedLaplacianAnalysis.

    Analyzes quantum propagators e^(-itL) and quantum information-theoretic
    properties including von Neumann entropy, quantum coherence, and
    interference effects.

    This class can be used standalone or in conjunction with the classical
    SignedLaplacianAnalysis for comparative studies.

    Attributes
    ----------
    sg : SignedGraph
        The signed graph to analyze.
    tTs_quantum : NDArray
        Logarithmic time grid for quantum evolution.
    init_node : int
        Initial node for localized quantum states.
    von_neumann_entropy : NDArray or None
        von Neumann entropy profile S_vN(t).
    quantum_coherence : NDArray or None
        Quantum coherence profile C(t).
    quantum_purity : NDArray or None
        Purity profile Tr(ρ²).
    participation_ratio : NDArray or None
        Inverse participation ratio profile.
    prob_distributions : NDArray or None
        Probability distribution P_j(t) for all nodes and times.
    quantum_classical_div : NDArray or None
        Divergence between quantum and classical evolution.

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> from lrgsglib.utils.lrg.quantum_core import QuantumSignedLaplacianAnalysis
    >>> l = Lattice2D(side=20, pflip=0.1)
    >>> l.flip_random_fract_edges()
    >>> qrg = QuantumSignedLaplacianAnalysis(l, init_node=200)
    >>> qrg.compute_all_quantum_observables()
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(qrg.tTs_quantum, qrg.von_neumann_entropy)
    >>> plt.xlabel("Quantum time")
    >>> plt.ylabel("von Neumann entropy")
    >>> plt.show()
    """

    # Quantum observables storage
    von_neumann_entropy: Optional[NDArray] = None
    quantum_coherence: Optional[NDArray] = None
    quantum_purity: Optional[NDArray] = None
    participation_ratio: Optional[NDArray] = None
    prob_distributions: Optional[NDArray] = None
    quantum_classical_div: Optional[NDArray] = None

    def __init__(
        self,
        sg: SignedGraph,
        tstep: int = LRSG_ENTROPY_STEP,
        tlex: float = -2,
        thex: float = 5,
        init_node: int = 0,
        init_type: str = "localized",
        initspect: bool = True,
    ) -> None:
        """
        Initialize quantum LRG analysis.

        Parameters
        ----------
        sg : SignedGraph
            Signed graph to analyze.
        tstep : int, optional
            Number of logarithmic time steps (default: LRSG_ENTROPY_STEP).
        tlex : float, optional
            Lower exponent for time grid (log scale) (default: -2).
        thex : float, optional
            Upper exponent for time grid (log scale) (default: 5).
        init_node : int, optional
            Initial node for localized quantum states (default: 0).
        init_type : str, optional
            Initial state type: "localized", "uniform", "thermal"
            (default: "localized").
        initspect : bool, optional
            Compute eigendecomposition on initialization (default: True).
        """
        self.sg = sg
        self.tstep = tstep
        self.tlex = tlex
        self.thex = thex
        self.init_node = init_node
        self.init_type = init_type

        # Quantum time grid
        self.tTs_quantum = np.logspace(self.tlex, self.thex, self.tstep)

        # Compute eigendecomposition if requested
        if initspect and not hasattr(self.sg, "eigv"):
            self.__initSpectrum__()

    def __initSpectrum__(self) -> None:
        """Compute eigendecomposition if not already done."""
        self.sg.compute_k_eigvV(backend="numpy")

    def compute_quantum_entropy(self) -> None:
        """
        Compute von Neumann entropy S_vN(t) = -Tr(ρ(t) log ρ(t)).

        Also computes coherence, purity, and participation ratio in a single
        pass for efficiency.

        Stores results in:
        - self.von_neumann_entropy
        - self.quantum_coherence
        - self.quantum_purity
        - self.participation_ratio
        - self.prob_distributions

        Raises
        ------
        ValueError
            If eigendecomposition has not been computed.

        Examples
        --------
        >>> from lrgsglib import Lattice2D
        >>> from lrgsglib.utils.lrg.quantum_core import QuantumSignedLaplacianAnalysis
        >>> l = Lattice2D(side=10)
        >>> qrg = QuantumSignedLaplacianAnalysis(l, init_node=50)
        >>> qrg.compute_quantum_entropy()
        >>> print(qrg.von_neumann_entropy.shape)
        (600,)
        """
        if not hasattr(self.sg, "eigv") or not hasattr(self.sg, "eigV"):
            raise ValueError(
                "Eigendecomposition required. "
                "Call compute_k_eigvV() or set initspect=True."
            )

        obs = compute_quantum_observables_from_eigenvalues(
            self.sg.eigv,
            self.sg.eigV,
            self.sg.N,
            self.tTs_quantum,
            init_type=self.init_type,
            init_node=self.init_node,
        )

        self.von_neumann_entropy = obs["von_neumann_entropy"]
        self.quantum_coherence = obs["coherence"]
        self.quantum_purity = obs["purity"]
        self.participation_ratio = obs["participation_ratio"]
        self.prob_distributions = obs["probability_distribution"]

    def compute_quantum_classical_divergence(self) -> None:
        """
        Compute divergence between quantum and classical evolution.

        D(t) = ||ρ_Q(t) - ρ_C(τ)||_F

        Uses Frobenius norm by default. The classical comparison uses thermal
        state at the same time parameter.

        Stores result in self.quantum_classical_div.

        Raises
        ------
        ValueError
            If eigendecomposition has not been computed.

        Examples
        --------
        >>> from lrgsglib import Lattice2D
        >>> from lrgsglib.utils.lrg.quantum_core import QuantumSignedLaplacianAnalysis
        >>> l = Lattice2D(side=10)
        >>> qrg = QuantumSignedLaplacianAnalysis(l, init_node=50)
        >>> qrg.compute_quantum_classical_divergence()
        >>> print(qrg.quantum_classical_div.shape)
        (600,)
        """
        if not hasattr(self.sg, "eigv") or not hasattr(self.sg, "eigV"):
            raise ValueError(
                "Eigendecomposition required. "
                "Call compute_k_eigvV() or set initspect=True."
            )

        div_profile = np.zeros(len(self.tTs_quantum))

        for i, t in enumerate(self.tTs_quantum):
            div_profile[i] = quantum_classical_divergence(
                self.sg.eigv,
                self.sg.eigV,
                t,
                tau=t,  # Use same time for comparison
                init_node=self.init_node,
                metric="frobenius",
            )

        self.quantum_classical_div = div_profile

    def compute_all_quantum_observables(self) -> None:
        """
        Compute all quantum observables in one pass.

        This is the main entry point for quantum analysis. It computes:
        - von Neumann entropy
        - Quantum coherence
        - Purity
        - Participation ratio
        - Probability distributions
        - Quantum-classical divergence

        Examples
        --------
        >>> from lrgsglib import Lattice2D
        >>> from lrgsglib.utils.lrg.quantum_core import QuantumSignedLaplacianAnalysis
        >>> l = Lattice2D(side=10)
        >>> qrg = QuantumSignedLaplacianAnalysis(l, init_node=50)
        >>> qrg.compute_all_quantum_observables()
        >>> # All observables are now available
        >>> print(qrg.von_neumann_entropy is not None)
        True
        >>> print(qrg.quantum_classical_div is not None)
        True
        """
        self.compute_quantum_entropy()  # Computes entropy, coherence, purity, etc.
        self.compute_quantum_classical_divergence()

    def get_quantum_entropy_derivative(self) -> NDArray:
        """
        Compute derivative of von Neumann entropy with respect to log(t).

        Returns
        -------
        entropy_derivative : NDArray
            d(S_vN) / d(log(t))

        Raises
        ------
        ValueError
            If von_neumann_entropy has not been computed.

        Notes
        -----
        Similar to classical capacity computation, this derivative can reveal
        quantum RG flow features.

        Examples
        --------
        >>> qrg = QuantumSignedLaplacianAnalysis(l, init_node=50)
        >>> qrg.compute_quantum_entropy()
        >>> dS_dt = qrg.get_quantum_entropy_derivative()
        """
        if self.von_neumann_entropy is None:
            raise ValueError(
                "von Neumann entropy not computed. "
                "Call compute_quantum_entropy() first."
            )

        return np.gradient(self.von_neumann_entropy, np.log(self.tTs_quantum))

    def get_coherence_derivative(self) -> NDArray:
        """
        Compute derivative of quantum coherence with respect to log(t).

        Returns
        -------
        coherence_derivative : NDArray
            d(C) / d(log(t))

        Raises
        ------
        ValueError
            If quantum_coherence has not been computed.

        Notes
        -----
        The coherence derivative can identify time scales where quantum
        interference effects are strongest.

        Examples
        --------
        >>> qrg = QuantumSignedLaplacianAnalysis(l, init_node=50)
        >>> qrg.compute_quantum_entropy()
        >>> dC_dt = qrg.get_coherence_derivative()
        """
        if self.quantum_coherence is None:
            raise ValueError(
                "Quantum coherence not computed. "
                "Call compute_quantum_entropy() first."
            )

        return np.gradient(self.quantum_coherence, np.log(self.tTs_quantum))

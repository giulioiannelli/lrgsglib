from __future__ import annotations

import numpy as np
from ..FullyConnected import FullyConnected
from .generators_scs import generate_scs_generalized_graph


class SCSGeneralizedNN(FullyConnected):
    """
    Generalized Sompolinsky-Crisanti-Sommers fully connected network.

    The class extends :class:`~lrgsglib.nx_patches.FullyConnected.FullyConnected`
    by drawing weighted edges from the correlated Gaussian ensemble described
    by Martorell et al. The correlation parameter ``gamma`` (in ``[-1, 1]``)
    tunes the reciprocity between pairs of directed edges, ranging from fully
    symmetric (``gamma = 1``) to fully anti-correlated (``gamma = -1``).
    """

    def __init__(
        self,
        N: int,
        gamma: float,
        J0: float = 0.0,
        J: float = 1.0,
        g: float = 1.0,
        diagonal=None,
        sgpathn: str = None,
        stdFnameSFFX: str = None,
        only_const_mode: bool = False,
        **kwargs,
    ):
        self.gamma = gamma
        self.J0 = J0
        self.J = J
        self.g = g
        self.diagonal = diagonal
        self.only_const_mode = only_const_mode
        
        # Extract seed from kwargs to create RNG before parent initialization
        seed = kwargs.get('seed', None)
        self._rng = np.random.default_rng(seed)
        
        # Set sgpathn before calling parent __init__
        # Use parent's FC_SGPATH if sgpathn not provided, or SCS-specific path
        from ...config.const import SCS_SGPATH, SCS_PHTABB, FC_SGPATH
        if sgpathn is None:
            sgpathn = SCS_SGPATH if SCS_SGPATH else FC_SGPATH
        
        # Store sgpathn before parent init, but don't pass it (parent will set it)
        _sgpathn = SCS_PHTABB if not sgpathn else sgpathn
        
        super().__init__(
            N=N, 
            sgpathn=_sgpathn,
            stdFnameSFFX=stdFnameSFFX if stdFnameSFFX is not None else "",
            only_const_mode=only_const_mode, 
            **kwargs
        )
        base_fname = f"scs_nn_N={self._N}_gamma={self.gamma:.3g}"
        suffix = getattr(self, "peq_str", "")
        self.std_fname = f"{base_fname}_{suffix}" if suffix else base_fname

    def __init_fullyconnected__(self) -> None:
        self.G = generate_scs_generalized_graph(
            self._N,
            self.gamma,
            J0=self.J0,
            J=self.J,
            g=self.g,
            diagonal=self.diagonal,
            rng=self._rng,
        )

    def __repr__(self) -> str:
        return (
            "SCSGeneralizedNN("
            f"N={self._N}, gamma={self.gamma:.3g}, J0={self.J0:.3g}, J={self.J:.3g}, g={self.g:.3g}"
            ")"
        )

    def __str__(self) -> str:
        return self.__repr__()

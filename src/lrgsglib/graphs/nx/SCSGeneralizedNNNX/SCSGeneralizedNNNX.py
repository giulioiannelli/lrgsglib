from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
from ..FullyConnectedNX import FullyConnectedNX
from .generators_scs import generate_scs_generalized_graph


class SCSGeneralizedNNNX(FullyConnectedNX):
    """
    Generalized Sompolinsky-Crisanti-Sommers fully connected network.

    The class extends :class:`~lrgsglib.graphs.nx.FullyConnectedNX.FullyConnectedNX`
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

        # Set sgpathn with intermediate folder (like Lattice2D does with L2D_PATH_DICT)
        # This ensures output goes to data/<workdir>/scs_nn/ instead of data/<workdir>/
        from ....config.const import SCS_SGPATH, SCS_PHTABB, FC_SGPATH
        if sgpathn is None:
            sgpathn = SCS_SGPATH if SCS_SGPATH else FC_SGPATH

        # Add intermediate folder 'scs_nn' to match Lattice2D structure
        if sgpathn:
            _sgpathn = str(Path(sgpathn) / SCS_PHTABB)
        else:
            _sgpathn = SCS_PHTABB

        super().__init__(
            N=N,
            sgpathn=_sgpathn,
            stdFnameSFFX=stdFnameSFFX if stdFnameSFFX is not None else "",
            only_const_mode=only_const_mode,
            **kwargs
        )
        # std_fname should only contain parameters not already in directory path
        # (scs_nn and N are already in the path: data/.../scs_nn/.../N=<N>/)
        base_fname = f"gamma={self.gamma:.3g}"
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

    def is_symmetric(self) -> bool:
        """Whether the weight matrix is symmetric (gamma=1, undirected)."""
        return not self.G.is_directed()

    def get_weight_matrix(self) -> np.ndarray:
        """Return the weight matrix as a numpy array of shape (N, N)."""
        return nx.to_numpy_array(self.G, weight="weight")

    def __repr__(self) -> str:
        return (
            "SCSGeneralizedNNNX("
            f"N={self._N}, gamma={self.gamma:.3g}, J0={self.J0:.3g}, J={self.J:.3g}, g={self.g:.3g}"
            ")"
        )

    def __str__(self) -> str:
        return self.__repr__()

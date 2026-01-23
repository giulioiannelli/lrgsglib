"""
ConfigurationModel: Random graph with prescribed degree sequence.

The configuration model generates a random graph with a specific
degree sequence, allowing study of degree distribution effects
decoupled from other structural properties.

Example
-------
>>> from lrgsglib.nx_patches.random import ConfigurationModel
>>> degrees = [3, 3, 3, 3, 4, 4, 4, 4, 5, 5]  # sum must be even
>>> cm = ConfigurationModel(degree_sequence=degrees, pflip=0.1)
>>> cm.flip_random_fract_edges()
"""

from typing import Sequence

import networkx as nx
import numpy as np

from .random_graph import RandomGraph


# Constants
CM_PHTABB = "cm"
CM_SGPATH = ""
CM_STDFN = ""


class ConfigurationModel(RandomGraph):
    """
    Signed random graph with prescribed degree sequence.

    Uses the configuration model to generate a simple graph (no self-loops
    or multi-edges) with the specified degree sequence.

    Parameters
    ----------
    degree_sequence : Sequence[int]
        Degree sequence. Sum must be even. Length determines n.
    sgpathn : str, default CM_SGPATH
        Subpath for storing graph data.
    stdFnameSFFX : str, default CM_STDFN
        Filename suffix for exports.
    only_const_mode : bool, default False
        If True, only populate metadata without building graph.
    **kwargs : Any
        Forwarded to SignedGraph (pflip, seed, etc.).

    Attributes
    ----------
    degree_sequence : list[int]
        The prescribed degree sequence.
    mean_degree : float
        Mean of the degree sequence.

    Notes
    -----
    The configuration model is useful for:
    - Studying effects of degree heterogeneity
    - Creating null models with matched degree sequences
    - Testing community detection under degree constraints

    The generated graph may have self-loops and multi-edges removed,
    so the actual degrees may be slightly less than prescribed.

    Examples
    --------
    >>> from lrgsglib.nx_patches.random import ConfigurationModel
    >>> # Power-law degree sequence
    >>> import numpy as np
    >>> degrees = np.random.zipf(2.5, 100)
    >>> degrees = degrees[degrees < 50]  # cap max degree
    >>> if sum(degrees) % 2: degrees[0] += 1  # ensure even sum
    >>> cm = ConfigurationModel(degree_sequence=degrees, seed=42)
    """

    def __init__(
        self,
        degree_sequence: Sequence[int],
        sgpathn: str = CM_SGPATH,
        stdFnameSFFX: str = CM_STDFN,
        only_const_mode: bool = False,
        **kwargs,
    ) -> None:
        self.degree_sequence = list(degree_sequence)

        if sum(self.degree_sequence) % 2 != 0:
            raise ValueError("Sum of degree sequence must be even")

        self.mean_degree = np.mean(self.degree_sequence)
        self._init_std_fname(stdFnameSFFX)
        sgpath = CM_PHTABB if not sgpathn else sgpathn

        super().__init__(
            n=len(degree_sequence),
            sgpathn=sgpath,
            std_fname=self.std_fname,
            only_const_mode=only_const_mode,
            extract_giant_component=True,
            **kwargs,
        )

    def _init_std_fname(self, suffix: str = "") -> None:
        self.std_fname = CM_PHTABB + suffix

    def _generate_graph(self) -> nx.Graph:
        """Generate configuration model graph, removing self-loops and multi-edges."""
        G = nx.configuration_model(self.degree_sequence)
        # Remove self-loops and multi-edges
        G = nx.Graph(G)
        G.remove_edges_from(nx.selfloop_edges(G))
        return G

    def _compute_syshapePth(self) -> str:
        return f"N={self.n}_kavg={self.mean_degree:.2f}"

    def get_degree_statistics(self) -> dict:
        """
        Return statistics of the prescribed degree sequence.

        Returns
        -------
        dict
            Dictionary with mean, std, min, max of degree sequence.
        """
        return {
            'mean': float(np.mean(self.degree_sequence)),
            'std': float(np.std(self.degree_sequence)),
            'min': int(np.min(self.degree_sequence)),
            'max': int(np.max(self.degree_sequence)),
        }

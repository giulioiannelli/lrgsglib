"""Hierarchical Modular Network graph class."""

from __future__ import annotations

import networkx as nx

from ..common import *
from ..MultispectralGraph import MultispectralGraph
from .generators_hmn import hierarchical_modular_network


# Default constants for HMN
HMN_LEVELS = 3
HMN_BRANCHING = 2
HMN_LEAF_NODES = 8
HMN_P_INTRA = 0.8
HMN_P_RATIO = 0.1
HMN_STDFN = ""
HMN_SGPATH = "msg_hmn"


class HierarchicalModularNetwork(MultispectralGraph):
    """Hierarchical Modular Network with nested community structure.

    Creates a signed graph with hierarchical modularity where nodes are
    organized into nested communities. Connection probability decreases
    exponentially with hierarchical distance.

    Parameters
    ----------
    levels : int
        Number of hierarchical levels (tree depth). Default: 3
    branching : int
        Number of submodules per module. Default: 2
    leaf_nodes : int
        Number of nodes in each leaf (lowest level) module. Default: 8
    p_intra : float
        Edge probability within leaf modules. Default: 0.8
    p_ratio : float
        Decay ratio for edge probability between levels.
        p(level k) = p_intra * p_ratio^k. Default: 0.1
    pflip : float
        Fraction of edges to flip negative (frustration parameter). Default: 0.0
    seed : int, optional
        Random seed for reproducibility.
    **kwargs
        Additional arguments passed to MultispectralGraph/SignedGraph.

    Attributes
    ----------
    levels : int
        Number of hierarchical levels.
    branching : int
        Branching factor.
    leaf_nodes : int
        Nodes per leaf module.
    p_intra : float
        Intra-module connection probability.
    p_ratio : float
        Inter-level probability ratio.
    n_modules : int
        Total number of leaf modules = branching^levels.

    Examples
    --------
    >>> from lrgsglib.nx_patches.HierarchicalModular import HierarchicalModularNetwork
    >>> hmn = HierarchicalModularNetwork(levels=2, branching=2, leaf_nodes=4)
    >>> hmn.N  # Total nodes = 2^2 * 4 = 16
    16
    >>> hmn.n_modules
    4

    >>> # With frustration
    >>> hmn = HierarchicalModularNetwork(levels=2, branching=2, leaf_nodes=4, pflip=0.3)
    >>> hmn.flip_random_fract_edges()
    >>> hmn.Ne_n > 0  # Some negative edges
    True

    Notes
    -----
    The hierarchical structure creates a natural multiscale organization:
    - Level 0: Tight communities (leaf modules) with p_intra
    - Level 1: Module pairs with p_intra * p_ratio
    - Level k: Larger modules with p_intra * p_ratio^k

    This model is useful for studying:
    - Community detection algorithms
    - Diffusion on modular networks
    - Spectral properties of hierarchical structures
    """

    def __init__(
        self,
        levels: int = HMN_LEVELS,
        branching: int = HMN_BRANCHING,
        leaf_nodes: int = HMN_LEAF_NODES,
        p_intra: float = HMN_P_INTRA,
        p_ratio: float = HMN_P_RATIO,
        stdFnameSFFX: str = HMN_STDFN,
        sgpathn: str = HMN_SGPATH,
        out_suffix: str = "",
        **kwargs,
    ):
        # Store parameters
        self.levels = levels
        self.branching = branching
        self.leaf_nodes = leaf_nodes
        self.p_intra = p_intra
        self.p_ratio = p_ratio
        self.out_suffix = out_suffix

        # Computed attributes
        self.n_modules = branching ** levels

        # Build syshapePth before calling parent
        from ..common import DEFAULT_P_FSTR_FMT
        fmtfunc = lambda x: f"{x:{DEFAULT_P_FSTR_FMT}}"

        # Format: L={levels}_B={branching}_N={leaf_nodes}_pi={p_intra}_pr={p_ratio}
        self.syshapePth = (
            f"L={levels}_"
            f"B={branching}_"
            f"Nl={leaf_nodes}_"
            f"pi={fmtfunc(p_intra)}_"
            f"pr={fmtfunc(p_ratio)}"
        )

        if out_suffix:
            self.syshapePth += f"_{out_suffix}"

        super().__init__(
            stdFnameSFFX=stdFnameSFFX,
            sgpathn=sgpathn,
            **kwargs
        )

    def _generate(self) -> nx.Graph:
        """Generate hierarchical modular network."""
        seed = getattr(self, 'seed', None)
        H = hierarchical_modular_network(
            levels=self.levels,
            branching=self.branching,
            leaf_nodes=self.leaf_nodes,
            p_intra=self.p_intra,
            p_ratio=self.p_ratio,
            seed=seed,
        )
        self.H = H
        return H

    def get_expected_num_nodes(self) -> int:
        """Return expected number of nodes."""
        return self.n_modules * self.leaf_nodes

    def get_module_for_node(self, node: int) -> int:
        """Return the leaf module index for a given node."""
        return node // self.leaf_nodes

    def get_module_nodes(self, module_id: int) -> list[int]:
        """Return list of nodes in a given leaf module."""
        start = module_id * self.leaf_nodes
        return list(range(start, start + self.leaf_nodes))

    def __repr__(self) -> str:
        return (
            f"HierarchicalModularNetwork("
            f"N={self.N}, levels={self.levels}, "
            f"branching={self.branching}, leaf_nodes={self.leaf_nodes})"
        )

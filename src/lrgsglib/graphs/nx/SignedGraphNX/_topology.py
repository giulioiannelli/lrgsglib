import networkx as nx
import numpy as np
from scipy.sparse import spdiags, identity as scsp_identity, csr_array
from typing import Any, TYPE_CHECKING

from ....config.const import (
    SG_REPR,
    SG_LAPL_SIGNED,
    SG_LAPL_RW,
    SG_LAPL_SYM,
    SG_LAPL_TYPES,
    SG_LAPL_DEFAULT_TYPE,
)

if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX




def get_adjacency_matrix_for(self, on_g: str = SG_REPR):
    """Get the adjacency matrix for a specific graph representation."""
    return self.adjacency_matrices.get(on_g)

def get_degree_matrix_for(self, on_g: str = SG_REPR):
    """Get the degree matrix for a specific graph representation."""
    return self.degree_matrices.get(on_g)

def get_laplacian_matrix_for(self, on_g: str = SG_REPR):
    """Get the Laplacian matrix for a specific graph representation."""
    return self.laplacian_matrices.get(on_g)

def get_signed_degree_matrix_for(self, on_g: str = SG_REPR):
    """Get the signed degree matrix for a specific graph representation."""
    return self.signed_degree_matrices.get(on_g)

def get_signed_laplacian_matrix_for(self, on_g: str = SG_REPR):
    """Get the signed Laplacian matrix for a specific graph representation."""
    return self.signed_laplacian_matrices.get(on_g)

def get_laplacian(self: "SignedGraphNX"):
    return self.degm - self.adj


def get_signed_laplacian(self: "SignedGraphNX"):
    return self.sdeg - self.adj


def get_signed_rw_laplacian(
    self: "SignedGraphNX", sym: bool = True, format: str = 'csr'
):
    """Random-walk normalized signed Laplacian (Kunegis 2010, §3.3/§3.4).

    Parameters
    ----------
    sym : bool, default True
        ``True``  -> ``L_sym = I - D_s^{-1/2} A D_s^{-1/2}`` (symmetric, PSD).
        ``False`` -> ``L_rw  = I - D_s^{-1} A`` (random walk, NOT symmetric).
    format : str, default 'csr'
        Sparse format of the returned matrix.

    Notes
    -----
    Built from the cached signed adjacency ``self.adj``; ``D_s`` is the signed
    (absolute) degree matrix. The two variants are isospectral
    (``L_rw = D_s^{-1/2} L_sym D_s^{1/2}``), so their eigenvalues coincide and
    are real in ``[0, 2]``. Isolated nodes (zero signed degree) get a zero
    inverse, leaving that row equal to the identity row (eigenvalue 1).
    """
    if self.adj is None:
        self.upd_graph_matrices()
    A = self.adj
    n = A.shape[0]
    deg = np.asarray(abs(A).sum(axis=1)).ravel()
    nz = deg > 0
    inv = np.zeros_like(deg, dtype=float)
    if sym:
        inv[nz] = 1.0 / np.sqrt(deg[nz])
        D = csr_array(spdiags(inv, 0, n, n, format=format))
        return csr_array(scsp_identity(n, format=format)) - D @ A @ D
    inv[nz] = 1.0 / deg[nz]
    D = csr_array(spdiags(inv, 0, n, n, format=format))
    return csr_array(scsp_identity(n, format=format)) - D @ A


def _laplacian_operator(
    self: "SignedGraphNX", laplacian_type: str = SG_LAPL_DEFAULT_TYPE
):
    """Resolve ``laplacian_type`` to ``(matrix, is_symmetric)``.

    ``'signed'`` -> combinatorial ``self.slp`` (symmetric, default behavior).
    ``'sym'``    -> symmetric normalized ``L_sym`` (symmetric).
    ``'rw'``     -> random-walk ``L_rw`` (NOT symmetric; needs a general solver).
    """
    if laplacian_type == SG_LAPL_SIGNED:
        if self.slp is None:
            self.upd_graph_matrices()
        return self.slp, True
    if laplacian_type == SG_LAPL_SYM:
        return self.get_signed_rw_laplacian(sym=True), True
    if laplacian_type == SG_LAPL_RW:
        return self.get_signed_rw_laplacian(sym=False), False
    raise ValueError(
        f"laplacian_type must be one of {SG_LAPL_TYPES}, got {laplacian_type!r}"
    )


def get_signed_laplacian_embedding(self: "SignedGraphNX", k: int = 2):
    return self.eigV[:k]


def make_rescaled_signed_laplacian(
    self: "SignedGraphNX",
    MODE: str = 'field'
):
    if MODE == 'field':
        self.resLp = self.slp - self.eigv[0] * scsp_identity(self.N)
    elif MODE == 'double':
        from scipy.linalg import eigvalsh

        self.resLp = self.slp - np.array([self.eigv[0]])
        new_eigv0 = eigvalsh(
            self.resLp.astype(np.float64), 
            subset_by_index=[0, 0]
        )
        self.resLp = self.resLp - new_eigv0 * np.identity(self.N)


def nodes_in(self: "SignedGraphNX", on_g: str = SG_REPR) -> list:
    """
    Get list of nodes in specified graph representation.

    Parameters
    ----------
    on_g : str, default SG_REPR
        Graph representation identifier.

    Returns
    -------
    list
        List of node identifiers.
    """
    return list(self.gr[on_g].nodes())


def get_nodes_list(self: "SignedGraphNX", on_g: str = None) -> list:
    """
    Get list of nodes for a specific graph representation.

    Parameters
    ----------
    on_g : str, optional
        Graph representation identifier. If None, uses self.on_g.

    Returns
    -------
    list
        List of nodes in the specified graph representation.
    """
    on_g = on_g or self.on_g
    return list(self.gr[on_g].nodes())


def get_node_attributes(
        self: "SignedGraphNX",
        attr: str = 'pos', 
        on_g: str = SG_REPR
):
    return nx.get_node_attributes(self.gr[on_g], attr)


def get_edge_data(
        self: "SignedGraphNX",
        u: Any, 
        v: Any, 
        thedata: str = 'weight',
        on_g: str = SG_REPR
):
    edge_data = self.gr[on_g].get_edge_data(u, v)
    if edge_data is None:
        raise KeyError(f"Edge ({u}, {v}) not found in graph {on_g}")
    return edge_data[thedata]


def get_edge_mapping(
    self: "SignedGraphNX",
    edge,
    target_g,
    on_g: str = SG_REPR
):
    try:
        return self.map_edge[target_g][on_g][edge]
    except KeyError:
        rev_edge = edge[::-1]
        return self.map_edge[target_g][on_g].get(rev_edge, None)


def get_edge_color(
        self: "SignedGraphNX",
        pec = "blue", 
        nec = "red",
        thedata: str = 'weight', 
        on_g: str = SG_REPR,
        continuous: bool = False, 
        cmap: str = 'coolwarm'
):
    from matplotlib.pyplot import get_cmap
    
    def map_values(value):
        if continuous:
            # Normalize value to [0, 1]
            if max_val != min_val:
                norm_value = (value - min_val) / (max_val - min_val)
            else:
                norm_value = 0.5  # If all values are the same, use middle of colormap
            return norm_value  # Return normalized value, not RGBA
        return nec if value == -1 else pec if value == 1 else value
    
    arr = nx.get_edge_attributes(self.gr[on_g], thedata)

    if continuous:
        values = np.array(list(arr.values()))
        min_val, max_val = values.min(), values.max()

    return list(map(map_values, arr.values()))


def get_graph_neighbors(
        self: "SignedGraphNX",
        node: Any, 
        on_g: str = SG_REPR
):
    return list(self.gr[on_g].neighbors(node))


def get_adjacency_matrix(
        self: "SignedGraphNX",
        on_g: str = SG_REPR, 
        weight: str = 'weight', 
        format: str = 'csr'
):
    return nx.to_scipy_sparse_array(
        self.gr[on_g], 
        weight=weight, 
        format=format
    )


def get_degree_matrix(self: "SignedGraphNX", format: str = 'csr'):
    diag_values = self.adj.sum(axis=1)
    return spdiags(diag_values, 0, *self.adj.shape, format=format)


def get_abs_degree_matrix(self: "SignedGraphNX", format: str = 'csr'):
    diag_values = abs(self.adj).sum(axis=1)
    return spdiags(
        diag_values, 0, *self.adj.shape, 
        format=format
    )


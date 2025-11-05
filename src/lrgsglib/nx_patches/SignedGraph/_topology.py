import networkx as nx
import numpy as np
from scipy.sparse import spdiags, identity as scsp_identity
from typing import Any, TYPE_CHECKING

from ...config.const import SG_REPR

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph




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

def get_laplacian(self: "SignedGraph"):
    return self.degm - self.adj


def get_signed_laplacian(self: "SignedGraph"):
    return self.sdeg - self.adj


def get_signed_laplacian_embedding(self: "SignedGraph", k: int = 2):
    return self.eigV[:k]


def make_rescaled_signed_laplacian(
    self: "SignedGraph",
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


def get_nodes_list(self: "SignedGraph", on_g: str = None) -> list:
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
        self: "SignedGraph",
        attr: str = 'pos', 
        on_g: str = SG_REPR
):
    return nx.get_node_attributes(self.gr[on_g], attr)


def get_edge_data(
        self: "SignedGraph",
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
    self: "SignedGraph",
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
        self: "SignedGraph",
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
            norm_value = (value - min_val) / (max_val - min_val)
            return get_cmap(cmap)(norm_value)  # Get color from cmap
        return nec if value == -1 else pec if value == 1 else value
    
    arr = nx.get_edge_attributes(self.gr[on_g], thedata)

    if continuous:
        values = np.array(list(arr.values()))
        min_val, max_val = values.min(), values.max()

    return list(map(map_values, arr.values()))


def get_graph_neighbors(
        self: "SignedGraph",
        node: Any, 
        on_g: str = SG_REPR
):
    return list(self.gr[on_g].neighbors(node))


def get_adjacency_matrix(
        self: "SignedGraph",
        on_g: str = SG_REPR, 
        weight: str = 'weight', 
        format: str = 'csr'
):
    return nx.to_scipy_sparse_array(
        self.gr[on_g], 
        weight=weight, 
        format=format
    )


def get_degree_matrix(self: "SignedGraph", format: str = 'csr'):
    diag_values = self.adj.sum(axis=1)
    return spdiags(diag_values, 0, *self.adj.shape, format=format)


def get_abs_degree_matrix(self: "SignedGraph", format: str = 'csr'):
    diag_values = abs(self.adj).sum(axis=1)
    return spdiags(
        diag_values, 0, *self.adj.shape, 
        format=format
    )


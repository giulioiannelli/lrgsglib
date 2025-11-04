"""
Test to verify that type hints allow VSCode/Pylance to resolve cross-file definitions.

When you open this file in VSCode with Pylance enabled:
1. Hover over `self.get_adjacency_matrix` - should show definition from _topology.py
2. Ctrl+Click (Go to Definition) on any self.method should jump to the submodule
3. Autocomplete on `self.` should show all available methods from all submodules
"""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def example_function_using_self_methods(self: "SignedGraph"):
    """
    Example showing that with type hints, Pylance can now resolve:
    - self.get_adjacency_matrix (from _topology.py)
    - self.upd_graph_matrices (from _representations.py)  
    - self.flip_sel_edges (from _ongraph.py)
    - self.get_eigV_cluster_sizes (from _partitioning.py)
    
    Try: Hover over any of these method calls to see definitions!
    Try: Ctrl+Click to jump to definition in the submodule file!
    """
    # From _topology.py
    adj = self.get_adjacency_matrix(on_g=self.on_g)
    
    # From _representations.py
    self.upd_graph_matrices()
    
    # From _ongraph.py
    self.flip_sel_edges([(0, 1)])
    
    # From _partitioning.py  
    cluster_sizes = self.get_eigV_cluster_sizes(which=0)
    
    # Properties defined in SignedGraph should also work
    num_nodes = self.N
    num_edges = self.Ne
    
    return adj, cluster_sizes, num_nodes, num_edges

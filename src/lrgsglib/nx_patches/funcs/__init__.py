"""Utility functions used across :mod:`lrgsglib.nx_patches`."""

from .base import (
    fast_set_weights_from_matrix,
    remove_edges,
    rewire_edges_optimized,
    get_giant_component,
    get_giant_component_leftoff,
)
from .neighbors import (
    get_kth_order_neighbours,
    get_neighbors_at_distance,
    get_smallest_cycle_graph_node,
)
from .spectral import (
    signed_laplacian_matrix,
    signed_spectral_layout,
)
from .lattice import (
    LatticeND_graph_FastPatch,
    LatticeND_graph_with_dilution,
)
from .thresholding import (
    compute_threshold_stats,
    compute_threshold_stats_fast,
    find_exact_detachment_threshold,
    select_threshold_and_graph,
    threshold_graph,
)

__all__ = [
    'fast_set_weights_from_matrix',
    'remove_edges',
    'find_exact_detachment_threshold',
    'get_kth_order_neighbours',
    'get_neighbors_at_distance',
    'get_smallest_cycle_graph_node',
    'signed_laplacian_matrix',
    'signed_spectral_layout',
    'LatticeND_graph_FastPatch',
    'rewire_edges_optimized',
    'get_giant_component',
    'get_giant_component_leftoff',
    'compute_threshold_stats',
    'compute_threshold_stats_fast',
    'select_threshold_and_graph',
    'threshold_graph',
]

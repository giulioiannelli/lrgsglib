"""Laplacian Renormalization Group utilities.

This package provides tools for spectral analysis, entropy computation,
clustering, and quantum propagation on signed graphs.
"""

# Re-export public APIs from submodules
from .spectral import (
    get_graph_lspectrum,
    get_graph_lspectrum_rw,
    compute_laplacian_properties,
)
# Import from infocomm subpackage (split for modularity)
from .infocomm import (
    extract_ultrametric_matrix,
    lapl_dists,
    entropy,
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_slq,
    compute_entropy_observables_expm_multiply,
    compute_renyi_observables_from_eigenvalues,
)
from .ising import (
    compose_product_ising_states,
    compose_weighted_ising_state,
    compose_xor_ising_state,
    compute_ising_pairwise_energy,
    spin_overlap,
    spin_matching_fraction,
    spin_matching_fraction_fromovp,
    compute_spin_overlap_series,
    compute_spin_match_series,
    ising_spinglass_pmJ_2D_Tcrit,
)
from .clustering import (
    MakeLinkageMatrix,
    compute_normalized_linkage,
    compute_optimal_threshold,
    circular_layout_by_cluster,
)
from .quantum import (
    compute_quantum_propagator_spectral,
    compute_quantum_propagator_matrix,
    quantum_density_matrix_evolution,
    quantum_probability_distribution,
    compute_quantum_coherence,
    von_neumann_entropy,
    quantum_classical_divergence,
    compute_interference_visibility,
    compute_quantum_observables_from_eigenvalues,
)
from .quantum_core import QuantumSignedLaplacianAnalysis

__all__ = [
    # spectral
    "get_graph_lspectrum",
    "get_graph_lspectrum_rw",
    "compute_laplacian_properties",
    # infocomm
    "extract_ultrametric_matrix",
    "lapl_dists",
    "entropy",
    "compute_entropy_observables_from_eigenvalues",
    "compute_entropy_observables_slq",
    "compute_entropy_observables_expm_multiply",
    "compute_renyi_observables_from_eigenvalues",
    # ising
    "compose_product_ising_states",
    "compose_weighted_ising_state",
    "compose_xor_ising_state",
    "compute_ising_pairwise_energy",
    "spin_overlap",
    "spin_matching_fraction",
    "spin_matching_fraction_fromovp",
    "compute_spin_overlap_series",
    "compute_spin_match_series",
    "ising_spinglass_pmJ_2D_Tcrit",
    # clustering
    "MakeLinkageMatrix",
    "compute_normalized_linkage",
    "compute_optimal_threshold",
    "circular_layout_by_cluster",
    # quantum
    "compute_quantum_propagator_spectral",
    "compute_quantum_propagator_matrix",
    "quantum_density_matrix_evolution",
    "quantum_probability_distribution",
    "compute_quantum_coherence",
    "von_neumann_entropy",
    "quantum_classical_divergence",
    "compute_interference_visibility",
    "compute_quantum_observables_from_eigenvalues",
    "QuantumSignedLaplacianAnalysis",
]
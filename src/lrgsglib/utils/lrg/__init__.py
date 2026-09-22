"""Laplacian Renormalization Group utilities.

This package provides tools for spectral analysis, entropy computation,
clustering, and quantum propagation on signed graphs.
"""

from .clustering import (
    MakeLinkageMatrix,
    circular_layout_by_cluster,
    compute_normalized_linkage,
    compute_optimal_threshold,
    dendrogram_leaf_node_colors,
    log_dendrogram,
)

# Import from infocomm subpackage (split for modularity)
from .infocomm import (
    compute_entropy_observables_expm_multiply,
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_slq,
    compute_renyi_observables_from_eigenvalues,
    entropy,
    extract_ultrametric_matrix,
    lapl_dists,
    specific_heat_tau_window,
)
from .ising import (
    compose_product_ising_states,
    compose_weighted_ising_state,
    compose_xor_ising_state,
    compute_ising_pairwise_energy,
    compute_spin_match_series,
    compute_spin_overlap_series,
    ising_spinglass_pmJ_2D_Tcrit,
    spin_matching_fraction,
    spin_matching_fraction_fromovp,
    spin_overlap,
)
from .percolation import (
    build_cycle_dual,
    build_face_syndrome,
    build_signed_cycle_dual,
    build_vertex_syndrome,
    compute_pareto_point,
    compute_pareto_point_direct,
    compute_pareto_point_general,
    default_p_range,
    defect_percolation_sweep,
    dual_percolation_sweep,
    dual_percolation_sweep_general,
    find_threshold,
    giant_component_fraction,
    identify_frustrated_faces,
    identify_frustrated_vertices,
    percolation_sweep,
)
from .quantum import (
    compute_interference_visibility,
    compute_quantum_coherence,
    compute_quantum_observables_from_eigenvalues,
    compute_quantum_propagator_matrix,
    compute_quantum_propagator_spectral,
    quantum_classical_divergence,
    quantum_density_matrix_evolution,
    quantum_probability_distribution,
    von_neumann_entropy,
)
from .quantum_core import QuantumSignedLaplacianAnalysis
from .screening import (
    compute_face_statistics,
    compute_graph_properties,
    distance_above_pareto_front,
    fujii_reference_data,
    generate_gog_configs,
    generate_lattice_configs,
    generate_random_configs,
    plot_face_distributions,
    plot_pareto_front,
    plot_pc_landscape,
    plot_pc_vs_dimension,
    run_general_screening_batch,
    run_screening_batch,
    screen_topology,
    screen_topology_general,
)

# Re-export public APIs from submodules
from .spectral import (
    compute_laplacian_properties,
    get_graph_lspectrum,
    get_graph_lspectrum_rw,
)
from .spectral_rg import compute_frustration_fraction  # deprecated alias
from .spectral_rg import (
    agmon_geodesic_distance,
    build_reduced_graph,
    compute_eigenmode_sign_distance,
    compute_reduced_spectrum,
    compute_signed_diffusion_distance,
    find_rg_scales,
    frustration_index,
    lattice_plaquettes,
    negative_edge_fraction,
    partition_at_scale,
    rg_flow_observables,
    spectral_frustration,
    spectral_rg_flow,
    spectral_rg_step,
)

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
    "specific_heat_tau_window",
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
    "log_dendrogram",
    "dendrogram_leaf_node_colors",
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
    # percolation (eigenmode method)
    "percolation_sweep",
    "find_threshold",
    "dual_percolation_sweep",
    "compute_pareto_point",
    "default_p_range",
    # percolation (direct defect method)
    "identify_frustrated_faces",
    "identify_frustrated_vertices",
    "giant_component_fraction",
    "build_face_syndrome",
    "build_vertex_syndrome",
    "defect_percolation_sweep",
    "compute_pareto_point_direct",
    # generalized dual
    "build_cycle_dual",
    "build_signed_cycle_dual",
    "dual_percolation_sweep_general",
    "compute_pareto_point_general",
    # screening
    "compute_face_statistics",
    "screen_topology",
    "run_screening_batch",
    "generate_gog_configs",
    "fujii_reference_data",
    "distance_above_pareto_front",
    "plot_pareto_front",
    "plot_face_distributions",
    # spectral RG
    "negative_edge_fraction",
    "frustration_index",
    "spectral_frustration",
    "lattice_plaquettes",
    "compute_frustration_fraction",  # deprecated alias
    "compute_signed_diffusion_distance",
    "compute_eigenmode_sign_distance",
    "find_rg_scales",
    "partition_at_scale",
    "build_reduced_graph",
    "compute_reduced_spectrum",
    "spectral_rg_step",
    "spectral_rg_flow",
    "rg_flow_observables",
    # general topology screening
    "compute_graph_properties",
    "screen_topology_general",
    "run_general_screening_batch",
    "generate_lattice_configs",
    "generate_random_configs",
    "plot_pc_landscape",
    "plot_pc_vs_dimension",
]

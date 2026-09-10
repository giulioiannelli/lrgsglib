"""LRG / spectral analysis surface, reconstruction kernels, protein-TMD."""
# LRG / spectral utilities most commonly needed alongside SignedGraph methods.
from ..utils.lrg import (
    get_graph_lspectrum,
    compute_entropy_observables_from_eigenvalues,
    specific_heat_tau_window,
    lapl_dists,
    extract_ultrametric_matrix,
    MakeLinkageMatrix,
    compute_normalized_linkage,
    compute_optimal_threshold,
    circular_layout_by_cluster,
    log_dendrogram,
    dendrogram_leaf_node_colors,
    compute_signed_diffusion_distance,
    compute_eigenmode_sign_distance,
    agmon_geodesic_distance,
)

# Spectral-reconstruction kernel (used by CHL Chladni-state notebooks).
from ..utils.basic.linalg import compute_recon_ultra, compute_mse_from_recon

# Protein-TMD primitives (CHL-SF04 family — see utils/recon/protein/).
from ..utils.recon.protein import (
    SubstrateSpec,
    build_substrate_basis,
    substrate_signature,
    spec_label,
    DEFAULT_SUBSTRATES,
    download_pdb,
    extract_ca_coordinates,
    extract_atoms_coordinates,
    coords_to_pdb_string_with_structure,
    assign_secondary_structure_from_coords,
    pad_protein_coordinates,
    safe_reconstruct_coordinates_from_features,
    protein_to_coordinate_feature_vector,
    kabsch_align,
    kabsch_rmsd,
    q3_score,
    # nb_helpers — used by the CHL-SF04 notebook
    load_protein_corpus,
    load_tmd_bundles,
    interp_mse_curve,
    mse_per_substrate,
    pflip_curves,
    auc_per_substrate,
    q3_curve_per_substrate,
    q3_emergence_per_substrate,
    ss_segments,
    layout_positions_2d,
    render_substrate_panel,
    reconstruct_protein,
)

LINKAGE_METHOD = 'average'   # UPGMA — classical LRG choice
CMAP_CLUSTERS  = 'tab20'

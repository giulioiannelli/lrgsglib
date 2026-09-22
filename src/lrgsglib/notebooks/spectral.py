"""LRG / spectral analysis surface, reconstruction kernels, protein-TMD."""

# LRG / spectral utilities most commonly needed alongside SignedGraph methods.
# Spectral-reconstruction kernel (used by CHL Chladni-state notebooks).
from ..utils.basic.linalg import compute_mse_from_recon, compute_recon_ultra
from ..utils.lrg import (
    MakeLinkageMatrix,
    agmon_geodesic_distance,
    circular_layout_by_cluster,
    compute_eigenmode_sign_distance,
    compute_entropy_observables_from_eigenvalues,
    compute_normalized_linkage,
    compute_optimal_threshold,
    compute_signed_diffusion_distance,
    dendrogram_leaf_node_colors,
    extract_ultrametric_matrix,
    get_graph_lspectrum,
    lapl_dists,
    log_dendrogram,
    specific_heat_tau_window,
)

# Protein-TMD primitives (CHL-SF04 family — see utils/recon/protein/).
from ..utils.recon.protein import (  # nb_helpers — used by the CHL-SF04 notebook
    DEFAULT_SUBSTRATES,
    SubstrateSpec,
    assign_secondary_structure_from_coords,
    auc_per_substrate,
    build_substrate_basis,
    coords_to_pdb_string_with_structure,
    download_pdb,
    extract_atoms_coordinates,
    extract_ca_coordinates,
    interp_mse_curve,
    kabsch_align,
    kabsch_rmsd,
    layout_positions_2d,
    load_protein_corpus,
    load_tmd_bundles,
    mse_per_substrate,
    pad_protein_coordinates,
    pflip_curves,
    protein_to_coordinate_feature_vector,
    q3_curve_per_substrate,
    q3_emergence_per_substrate,
    q3_score,
    reconstruct_protein,
    render_substrate_panel,
    safe_reconstruct_coordinates_from_features,
    spec_label,
    ss_segments,
    substrate_signature,
)

LINKAGE_METHOD = "average"  # UPGMA — classical LRG choice
CMAP_CLUSTERS = "tab20"

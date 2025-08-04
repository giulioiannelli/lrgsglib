"""Spectral reconstruction utilities."""

from __future__ import annotations

from typing import Optional, Tuple, Dict, Any, List
from pathlib import Path

import numpy as np

from ....nx_patches import Lattice2D, Lattice3D, ErdosRenyi
from .feature_extraction import create_enhanced_feature_vector
from .reconstruction import safe_reconstruct_coordinates_from_features

__all__ = [
    "create_spectral_basis",
    "spectral_reconstruction_with_topology",
    "create_adaptive_spectral_basis",
    "create_fixed_dimension_spectral_basis",
    "fixed_dimension_protein_reconstruction",
    "analyze_proteins_with_residue_range",
]


def create_spectral_basis(
    topology: str,
    n_features: int = 1024,
    k_eigenvectors: Optional[int] = None,
    path_data=None,
    kw_args: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, Any]:
    """Create a spectral basis matrix for a given network topology."""
    if kw_args is None:
        kw_args = {}
    if path_data is not None and "path_data" not in kw_args:
        kw_args["path_data"] = path_data

    match topology:
        case "l2d_square":
            kw_args.setdefault("side1", 32)
            kw_args.setdefault("geo", "sqr")
            sg = Lattice2D(**kw_args)
        case "l2d_hex":
            kw_args.setdefault("side1", 32)
            kw_args.setdefault("side2", 32)
            kw_args.setdefault("geo", "hex")
            sg = Lattice2D(**kw_args)
        case "l2d_anti_tri":
            kw_args.setdefault("side1", 32)
            kw_args.setdefault("geo", "tri")
            kw_args.setdefault("pflip", 1)
            sg = Lattice2D(**kw_args)
            sg.flip_random_fract_edges()
        case "l3d":
            kw_args.setdefault("dim", (16, 8, 8))
            sg = Lattice3D(**kw_args)
        case "l3d_sg":
            kw_args.setdefault("dim", (16, 8, 8))
            kw_args.setdefault("pflip", 0.35)
            sg = Lattice3D(**kw_args)
            sg.flip_random_fract_edges()
        case "er":
            kw_args.setdefault("n", n_features)
            kw_args.setdefault("p", 0.01)
            kw_args.setdefault("seed", 42)
            sg = ErdosRenyi(**kw_args)
        case "er_sg":
            kw_args.setdefault("n", n_features)
            kw_args.setdefault("p", 0.01)
            kw_args.setdefault("seed", 42)
            kw_args.setdefault("pflip", 0.35)
            sg = ErdosRenyi(**kw_args)
            sg.flip_random_fract_edges()
        case "er_dense":
            kw_args.setdefault("n", n_features)
            kw_args.setdefault("p", 0.05)
            kw_args.setdefault("seed", 42)
            sg = ErdosRenyi(**kw_args)
        case "l2d_sqr_rew":
            kw_args.setdefault("side1", 32)
            kw_args.setdefault("geo", "sqr")
            kw_args.setdefault("prew", 0.1)
            kw_args.setdefault("seed", 42)
            sg = Lattice2D(**kw_args)
        case _:
            raise ValueError(f"Unknown topology: {topology}")

    sg.compute_laplacian_spectrum_weigV()
    basis = sg.get_sgspect_basis(max_factor=1)

    original_k = basis.shape[1]
    if k_eigenvectors is not None:
        if k_eigenvectors < original_k:
            print(f"  \U0001F4CA Truncating basis: {original_k} → {k_eigenvectors} eigenvectors")
            basis = basis[:, :k_eigenvectors]
        else:
            print(
                f"  \u26A0\uFE0F Requested {k_eigenvectors} eigenvectors, but only {original_k} available - using all"
            )
    else:
        print(f"  \U0001F4CA Using full spectral basis: {original_k} eigenvectors")

    return basis, sg


def spectral_reconstruction_with_topology(
    features: np.ndarray,
    topology: str,
    k_eigenvectors: Optional[int] = None,
    path_data=None,
    kw_args: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, np.ndarray, Any]:
    """Reconstruct a feature vector using a network topology."""
    basis, sg = create_spectral_basis(topology, len(features), k_eigenvectors, path_data, kw_args)
    basis_size = basis.shape[0]
    feature_size = len(features)

    if feature_size != basis_size:
        print(f"  \U0001F4CF Dimension mismatch: features={feature_size}, basis={basis_size}")
        if feature_size < basis_size:
            padded = np.zeros(basis_size)
            padded[:feature_size] = features
            features = padded
        else:
            features = features[:basis_size]

    spectral_coeffs = basis.T @ features
    reconstructed = basis @ spectral_coeffs

    if len(reconstructed) != feature_size:
        if len(reconstructed) < feature_size:
            padded = np.zeros(feature_size)
            padded[: len(reconstructed)] = reconstructed
            reconstructed = padded
        else:
            reconstructed = reconstructed[:feature_size]

    return reconstructed, basis, sg


def create_adaptive_spectral_basis(
    topology: str,
    n_residues: int,
    k_eigenvectors: Optional[int] = None,
    path_data=None,
    kw_args: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, Any]:
    """Create a spectral basis adapting the network size to the number of residues."""
    dof = 3 * n_residues
    if kw_args is None:
        kw_args = {}
    if path_data is not None and "path_data" not in kw_args:
        kw_args["path_data"] = path_data
    print(f"  \U0001F527 Adapting {topology} topology for {n_residues} residues ({dof} DOF)")
    match topology:
        case "l2d_square":
            side = int(np.ceil(np.sqrt(dof)))
            kw_args_adapted = {"side1": side, "geo": "sqr", **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            print(f"    \U0001F4D0 2D Square: {side}×{side} = {side**2} nodes (target: {dof})")
        case "l2d_hex":
            side = int(np.ceil(np.sqrt(dof)))
            kw_args_adapted = {"side1": side, "side2": side, "geo": "hex", **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            print(f"    \U0001F4D0 2D Hexagonal: {side}×{side} = {side**2} nodes (target: {dof})")
        case "l2d_anti_tri":
            side = int(np.ceil(np.sqrt(dof)))
            kw_args_adapted = {"side1": side, "geo": "tri", "pflip": 1, **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            sg.flip_random_fract_edges()
            print(f"    \U0001F4D0 2D Anti-triangular: {side}×{side} = {side**2} nodes (target: {dof})")
        case "l3d":
            side = int(np.ceil(dof ** (1 / 3)))
            dim = (side + 1, side, side) if side**3 < dof * 0.8 else (side, side, side)
            kw_args_adapted = {"dim": dim, **kw_args}
            sg = Lattice3D(**kw_args_adapted)
            actual_nodes = dim[0] * dim[1] * dim[2]
            print(f"    \U0001F4D0 3D Lattice: {dim[0]}×{dim[1]}×{dim[2]} = {actual_nodes} nodes (target: {dof})")
        case "l3d_sg":
            side = int(np.ceil(dof ** (1 / 3)))
            dim = (side + 1, side, side) if side**3 < dof * 0.8 else (side, side, side)
            kw_args_adapted = {"dim": dim, "pflip": 0.35, **kw_args}
            sg = Lattice3D(**kw_args_adapted)
            sg.flip_random_fract_edges()
            actual_nodes = dim[0] * dim[1] * dim[2]
            print(f"    \U0001F4D0 3D Small-world: {dim[0]}×{dim[1]}×{dim[2]} = {actual_nodes} nodes (target: {dof})")
        case "er":
            p = max(0.005, min(0.02, 10.0 / dof))
            kw_args_adapted = {"n": dof, "p": p, "seed": 42, **kw_args}
            sg = ErdosRenyi(**kw_args_adapted)
            print(f"    \U0001F4D0 Erd\u0151s-R\u00e9nyi: {dof} nodes, p={p:.4f}")
        case "er_sg":
            p = max(0.005, min(0.02, 10.0 / dof))
            kw_args_adapted = {"n": dof, "p": p, "seed": 42, "pflip": 0.35, **kw_args}
            sg = ErdosRenyi(**kw_args_adapted)
            sg.flip_random_fract_edges()
            print(f"    \U0001F4D0 Erd\u0151s-R\u00e9nyi SG: {dof} nodes, p={p:.4f}")
        case "er_dense":
            p = max(0.02, min(0.1, 50.0 / dof))
            kw_args_adapted = {"n": dof, "p": p, "seed": 42, **kw_args}
            sg = ErdosRenyi(**kw_args_adapted)
            print(f"    \U0001F4D0 Erd\u0151s-R\u00e9nyi Dense: {dof} nodes, p={p:.4f}")
        case "l2d_sqr_rew":
            side = int(np.ceil(np.sqrt(dof)))
            kw_args_adapted = {"side1": side, "geo": "sqr", "prew": 0.1, "seed": 42, **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            print(f"    \U0001F4D0 2D Square Rewired: {side}×{side} = {side**2} nodes (target: {dof})")
        case _:
            raise ValueError(f"Unknown topology: {topology}")

    sg.compute_laplacian_spectrum_weigV()
    basis = sg.get_sgspect_basis(max_factor=1)

    original_k = basis.shape[1]
    if k_eigenvectors is not None:
        if k_eigenvectors < original_k:
            print(f"    \U0001F4CA Truncating basis: {original_k} → {k_eigenvectors} eigenvectors")
            basis = basis[:, :k_eigenvectors]
        else:
            print(
                f"    \u26A0\uFE0F Requested {k_eigenvectors} eigenvectors, but only {original_k} available - using all"
            )
    else:
        print(f"    \U0001F4CA Using full spectral basis: {original_k} eigenvectors")

    return basis, sg


def create_fixed_dimension_spectral_basis(
    topology: str,
    target_residues: int,
    k_eigenvectors: Optional[int] = None,
    path_data=None,
    kw_args: Optional[Dict[str, Any]] = None,
) -> Tuple[np.ndarray, Any]:
    """Create a spectral basis sized for ``target_residues`` residues."""
    target_dof = 3 * target_residues
    if kw_args is None:
        kw_args = {}
    if path_data is not None and "path_data" not in kw_args:
        kw_args["path_data"] = path_data
    print(f"  \U0001F527 Creating {topology} topology for {target_residues} residues ({target_dof} DOF)")
    match topology:
        case "l2d_square":
            side = int(np.ceil(np.sqrt(target_dof)))
            kw_args_adapted = {"side1": side, "geo": "sqr", **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            print(f"    \U0001F4D0 2D Square: {side}×{side} = {side**2} nodes")
        case "l2d_hex":
            side = int(np.ceil(np.sqrt(target_dof)))
            kw_args_adapted = {"side1": side, "side2": side, "geo": "hex", **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            print(f"    \U0001F4D0 2D Hexagonal: {side}×{side} = {side**2} nodes")
        case "l2d_anti_tri":
            side = int(np.ceil(np.sqrt(target_dof)))
            kw_args_adapted = {"side1": side, "geo": "tri", "pflip": 1, **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            sg.flip_random_fract_edges()
            print(f"    \U0001F4D0 2D Anti-triangular: {side}×{side} = {side**2} nodes")
        case "l3d":
            side = int(np.ceil(target_dof ** (1 / 3)))
            dim = (side + 1, side, side) if side**3 < target_dof * 0.8 else (side, side, side)
            kw_args_adapted = {"dim": dim, **kw_args}
            sg = Lattice3D(**kw_args_adapted)
            actual_nodes = dim[0] * dim[1] * dim[2]
            print(f"    \U0001F4D0 3D Lattice: {dim[0]}×{dim[1]}×{dim[2]} = {actual_nodes} nodes")
        case "l3d_sg":
            side = int(np.ceil(target_dof ** (1 / 3)))
            dim = (side + 1, side, side) if side**3 < target_dof * 0.8 else (side, side, side)
            kw_args_adapted = {"dim": dim, "pflip": 0.35, **kw_args}
            sg = Lattice3D(**kw_args_adapted)
            sg.flip_random_fract_edges()
            actual_nodes = dim[0] * dim[1] * dim[2]
            print(f"    \U0001F4D0 3D Small-world: {dim[0]}×{dim[1]}×{dim[2]} = {actual_nodes} nodes")
        case "er":
            p = max(0.005, min(0.02, 10.0 / target_dof))
            kw_args_adapted = {"n": target_dof, "p": p, "seed": 42, **kw_args}
            sg = ErdosRenyi(**kw_args_adapted)
            print(f"    \U0001F4D0 Erd\u0151s-R\u00e9nyi: {target_dof} nodes, p={p:.4f}")
        case "er_sg":
            p = max(0.005, min(0.02, 10.0 / target_dof))
            kw_args_adapted = {"n": target_dof, "p": p, "seed": 42, "pflip": 0.35, **kw_args}
            sg = ErdosRenyi(**kw_args_adapted)
            sg.flip_random_fract_edges()
            print(f"    \U0001F4D0 Erd\u0151s-R\u00e9nyi SG: {target_dof} nodes, p={p:.4f}")
        case "er_dense":
            p = max(0.02, min(0.1, 50.0 / target_dof))
            kw_args_adapted = {"n": target_dof, "p": p, "seed": 42, **kw_args}
            sg = ErdosRenyi(**kw_args_adapted)
            print(f"    \U0001F4D0 Erd\u0151s-R\u00e9nyi Dense: {target_dof} nodes, p={p:.4f}")
        case "l2d_sqr_rew":
            side = int(np.ceil(np.sqrt(target_dof)))
            kw_args_adapted = {"side1": side, "geo": "sqr", "prew": 0.1, "seed": 42, **kw_args}
            sg = Lattice2D(**kw_args_adapted)
            print(f"    \U0001F4D0 2D Square Rewired: {side}×{side} = {side**2} nodes")
        case _:
            raise ValueError(f"Unknown topology: {topology}")

    sg.compute_laplacian_spectrum_weigV()
    basis = sg.get_sgspect_basis(max_factor=1)

    original_k = basis.shape[1]
    if k_eigenvectors is not None:
        if k_eigenvectors < original_k:
            print(f"    \U0001F4CA Truncating basis: {original_k} → {k_eigenvectors} eigenvectors")
            basis = basis[:, :k_eigenvectors]
        else:
            print(
                f"    \u26A0\uFE0F Requested {k_eigenvectors} eigenvectors, but only {original_k} available - using all"
            )
    else:
        print(f"    \U0001F4CA Using full spectral basis: {original_k} eigenvectors")

    return basis, sg


def pad_protein_coordinates(coords: np.ndarray, target_residues: int) -> np.ndarray:
    """Pad or truncate coordinates to a target number of residues."""
    n_residues = coords.shape[0]
    if n_residues == target_residues:
        return coords
    if n_residues < target_residues:
        padded_coords = np.zeros((target_residues, 3))
        padded_coords[:n_residues] = coords
        if n_residues > 0:
            centroid = coords.mean(axis=0)
            padded_coords[n_residues:] = centroid
        return padded_coords
    return coords[:target_residues]


def fixed_dimension_protein_reconstruction(
    coords: np.ndarray,
    topology: Optional[str] = None,
    target_residues: Optional[int] = None,
    k_eigenvectors: Optional[int] = None,
    path_data=None,
    kw_args: Optional[Dict[str, Any]] = None,
    basis: Optional[np.ndarray] = None,
    sg: Optional[Any] = None,
) -> Tuple[np.ndarray, np.ndarray, Any]:
    """
    Reconstruct coordinates using a spectral basis with fixed dimension.
    
    Args:
        coords: Input protein coordinates
        topology: Network topology (required if basis/sg not provided)
        target_residues: Target number of residues (required if basis/sg not provided)
        k_eigenvectors: Number of eigenvectors to use
        path_data: Path to data directory
        kw_args: Additional keyword arguments
        basis: Precomputed spectral basis (optional, computed if not provided)
        sg: Precomputed network structure (optional, computed if not provided)
    
    Returns:
        Tuple of (reconstructed_coords, basis, sg)
    """
    # Validate arguments
    if basis is None or sg is None:
        if topology is None or target_residues is None:
            raise ValueError("Either (basis, sg) or (topology, target_residues) must be provided")
        
        # Compute basis and sg if not provided
        if target_residues is None:
            target_residues = len(coords)
        basis, sg = create_fixed_dimension_spectral_basis(
            topology, target_residues, k_eigenvectors, path_data, kw_args
        )
    else:
        # When basis is provided, infer target_residues from basis size unless explicitly provided
        # This ensures each topology uses its optimal dimension
        if target_residues is None:
            target_residues = basis.shape[0] // 3  # Assuming 3D coordinates
    
    padded_coords = pad_protein_coordinates(coords, target_residues)
    enhanced_features, metadata = create_enhanced_feature_vector(
        padded_coords, n_features=3 * target_residues
    )
    
    basis_size = basis.shape[0]
    feature_size = len(enhanced_features)
    
    # Align feature size with basis size
    if feature_size != basis_size:
        if feature_size < basis_size:
            padded_features = np.zeros(basis_size)
            padded_features[:feature_size] = enhanced_features
            enhanced_features = padded_features
        else:
            enhanced_features = enhanced_features[:basis_size]
    
    # Perform spectral reconstruction
    spectral_coeffs = basis.T @ enhanced_features
    reconstructed_features = basis @ spectral_coeffs
    
    # Align reconstructed features back to original feature size
    if len(reconstructed_features) != feature_size:
        if len(reconstructed_features) < feature_size:
            padded = np.zeros(feature_size)
            padded[: len(reconstructed_features)] = reconstructed_features
            reconstructed_features = padded
        else:
            reconstructed_features = reconstructed_features[:feature_size]
    
    # Reconstruct coordinates
    reconstructed_coords = safe_reconstruct_coordinates_from_features(
        reconstructed_features, metadata
    )
    
    # Trim to original protein size
    original_size = coords.shape[0]
    reconstructed_coords = reconstructed_coords[:original_size]
    
    return reconstructed_coords, basis, sg


def analyze_proteins_with_residue_range(
    pdb_ids: List[str], 
    path_pdb_storage: Path,
    topologies: List[str],
    topology_names: Dict[str, str],
    min_residues: Optional[int] = None, 
    max_residues: Optional[int] = None,
    fixed_basis_size: Optional[int] = None,
    path_data=None
) -> Tuple[List[Dict], Dict[str, np.ndarray]]:
    """
    Analyze proteins with spectral reconstruction for different topologies.
    
    Args:
        pdb_ids: List of PDB IDs to analyze
        path_pdb_storage: Path to PDB storage directory
        topologies: List of topology names to test
        topology_names: Dictionary mapping topology codes to display names
        min_residues: Minimum number of residues (optional)
        max_residues: Maximum number of residues (optional)
        fixed_basis_size: Use fixed basis size instead of adaptive sizing (optional)
        path_data: Path for storing data
        
    Returns:
        Tuple of (selected_proteins_info, mse_results_dict)
    """
    from .io import download_pdb, extract_ca_coordinates
    
    print(f"🧬 PROTEIN SPECTRAL RECONSTRUCTION ANALYSIS")
    print("=" * 70)
    
    # Determine filter criteria
    if min_residues is None and max_residues is None:
        print("📊 No residue filters specified - using all proteins")
        filter_description = "all proteins"
    elif min_residues is not None and max_residues is not None:
        print(f"📊 Extracting proteins with {min_residues}-{max_residues} residues...")
        filter_description = f"{min_residues}-{max_residues} residues"
    elif min_residues is not None:
        print(f"📊 Extracting proteins with ≥{min_residues} residues...")
        filter_description = f"≥{min_residues} residues"
    else:  # max_residues is not None
        print(f"📊 Extracting proteins with ≤{max_residues} residues...")
        filter_description = f"≤{max_residues} residues"
    
    # Basis sizing mode
    if fixed_basis_size is not None:
        print(f"🔧 Using fixed basis size: {fixed_basis_size} features")
        sizing_mode = "fixed"
    else:
        print("🔧 Using adaptive basis sizing based on largest protein")
        sizing_mode = "adaptive"
    
    # Step 1: Extract proteins based on residue criteria
    selected_proteins = []
    
    for pdb_id in pdb_ids:
        try:
            pdb_content = download_pdb(pdb_id, path_pdb_storage)
            coords = extract_ca_coordinates(pdb_content)
            n_residues = len(coords)
            
            # Apply filters
            include_protein = True
            if min_residues is not None and n_residues < min_residues:
                include_protein = False
            if max_residues is not None and n_residues > max_residues:
                include_protein = False
            
            if include_protein:
                selected_proteins.append({
                    'pdb_id': pdb_id, 
                    'coords': coords, 
                    'n_residues': n_residues
                })
                print(f"  ✅ {pdb_id}: {n_residues} residues")
            else:
                reason = "too small" if (min_residues and n_residues < min_residues) else "too large"
                print(f"  ❌ {pdb_id}: {n_residues} residues ({reason})")
                
        except Exception as e:
            print(f"  ❌ Error with {pdb_id}: {e}")
    
    print(f"\n📈 Found {len(selected_proteins)} proteins matching criteria ({filter_description})")
    
    # Fallback if no proteins found
    if len(selected_proteins) == 0:
        print("❌ No proteins found matching criteria. Using first 5 available proteins instead.")
        fallback_min = 50  # Reasonable fallback
        for pdb_id in pdb_ids:
            try:
                pdb_content = download_pdb(pdb_id, path_pdb_storage)
                coords = extract_ca_coordinates(pdb_content)
                n_residues = len(coords)
                
                if n_residues >= fallback_min:
                    selected_proteins.append({
                        'pdb_id': pdb_id, 
                        'coords': coords, 
                        'n_residues': n_residues
                    })
                    if len(selected_proteins) >= 5:
                        break
            except:
                continue
    
    # Determine network sizing
    if selected_proteins:
        if fixed_basis_size is not None:
            min_features = fixed_basis_size
            print(f"🎯 Using fixed network sizing: {fixed_basis_size} features")
        else:
            protein_sizes = [p['n_residues'] for p in selected_proteins]
            target_residues = max(protein_sizes)  # Use largest protein for network sizing
            min_features = target_residues * 3
            print(f"🎯 Network sizing based on largest protein: {target_residues} residues ({min_features} features)")
    else:
        min_features = fixed_basis_size if fixed_basis_size is not None else 600
        print(f"⚠️ Using default network sizing: {min_features} features")
    
    # Step 2: Analyze each topology
    mse_protein = {}
    
    print(f"\n🌐 Computing spectral bases for each topology ({min_features} features)...")
    
    for topology in topologies:
        print(f'\n🔧 Processing topology: {topology_names.get(topology, topology)}')
        
        if fixed_basis_size is not None:
            # Use fixed dimension topology creation
            basis, sg = _create_topology_with_fixed_features(topology, fixed_basis_size, path_data)
        else:
            # Use adaptive topology creation
            basis, sg = _create_topology_with_min_features(topology, min_features, path_data)
            
        print(f"  ✅ Network: {sg.N} nodes, Basis: {basis.shape}")
        
        # Compute MSE curves for each protein
        protein_mse_curves = _compute_protein_mse_curves(
            selected_proteins[:10], basis, topology_names.get(topology, topology)
        )
        
        if protein_mse_curves:
            mse_protein[topology] = np.array(protein_mse_curves)
            print(f"  📊 Stored MSE curves for {len(protein_mse_curves)} proteins")
    
    print(f"\n✅ SPECTRAL RECONSTRUCTION ANALYSIS COMPLETED!")
    print(f"📊 Processed topologies: {list(mse_protein.keys())}")
    print(f"🧬 Analyzed {len(selected_proteins)} proteins ({filter_description})")
    print(f"🔧 Basis sizing: {sizing_mode} ({min_features} features)")
    
    return selected_proteins, mse_protein


def _create_topology_with_min_features(topology: str, min_features: int, path_data=None):
    """Create a topology with at least min_features nodes using only_const_mode."""
    kw_args = {'path_data': path_data, 'only_const_mode': True}
    
    match topology:
        case "l2d_square":
            side = int(np.ceil(np.sqrt(min_features)))
            kw_args.update({'side1': side, 'geo': 'sqr'})
            sg = Lattice2D(**kw_args)
        case "l2d_hex":
            side = int(np.ceil(np.sqrt(min_features)))
            kw_args.update({'side1': side, 'side2': side, 'geo': 'hex'})
            sg = Lattice2D(**kw_args)
        case "l2d_anti_tri":
            side = int(np.ceil(np.sqrt(min_features)))
            kw_args.update({'side1': side, 'geo': 'tri', 'pflip': 1})
            sg = Lattice2D(**kw_args)
        case "l3d":
            side = int(np.ceil(min_features ** (1/3)))
            kw_args.update({'dim': (side, side, side)})
            sg = Lattice3D(**kw_args)
        case "l3d_sg":
            side = int(np.ceil(min_features ** (1/3)))
            kw_args.update({'dim': (side, side, side), 'pflip': 0.35})
            sg = Lattice3D(**kw_args)
        case "er":
            kw_args.update({'n': min_features, 'p': 0.01, 'seed': 42})
            sg = ErdosRenyi(**kw_args)
        case "er_sg":
            kw_args.update({'n': min_features, 'p': 0.01, 'seed': 42, 'pflip': 0.35})
            sg = ErdosRenyi(**kw_args)
        case "er_dense":
            kw_args.update({'n': min_features, 'p': 0.05, 'seed': 42})
            sg = ErdosRenyi(**kw_args)
        case _:
            raise ValueError(f"Unknown topology: {topology}")
    
    # Get expected nodes
    expected_nodes = sg.get_expected_num_nodes()
    print(f"  📊 Expected nodes: {expected_nodes} (target: ≥{min_features})")
    
    # Reinitialize without only_const_mode
    kw_args['only_const_mode'] = False
    
    match topology:
        case "l2d_square":
            sg = Lattice2D(**kw_args)
        case "l2d_hex":
            sg = Lattice2D(**kw_args)
        case "l2d_anti_tri":
            sg = Lattice2D(**kw_args)
            sg.flip_random_fract_edges()
        case "l3d":
            sg = Lattice3D(**kw_args)
        case "l3d_sg":
            sg = Lattice3D(**kw_args)
            sg.flip_random_fract_edges()
        case "er":
            sg = ErdosRenyi(**kw_args)
        case "er_sg":
            sg = ErdosRenyi(**kw_args)
            sg.flip_random_fract_edges()
        case "er_dense":
            sg = ErdosRenyi(**kw_args)
    
    # Compute spectral basis
    sg.compute_laplacian_spectrum_weigV()
    basis = sg.get_sgspect_basis(max_factor=1)
    
    return basis, sg


def _create_topology_with_fixed_features(topology: str, fixed_features: int, path_data=None):
    """Create a topology with exactly the specified number of features."""
    kw_args = {'path_data': path_data, 'only_const_mode': True}
    
    print(f"  🔧 Creating {topology} with exactly {fixed_features} features...")
    
    match topology:
        case "l2d_square":
            side = int(np.ceil(np.sqrt(fixed_features)))
            kw_args.update({'side1': side, 'geo': 'sqr'})
            sg = Lattice2D(**kw_args)
        case "l2d_hex":
            side = int(np.ceil(np.sqrt(fixed_features)))
            kw_args.update({'side1': side, 'side2': side, 'geo': 'hex'})
            sg = Lattice2D(**kw_args)
        case "l2d_anti_tri":
            side = int(np.ceil(np.sqrt(fixed_features)))
            kw_args.update({'side1': side, 'geo': 'tri', 'pflip': 1})
            sg = Lattice2D(**kw_args)
        case "l3d":
            side = int(np.ceil(fixed_features ** (1/3)))
            kw_args.update({'dim': (side, side, side)})
            sg = Lattice3D(**kw_args)
        case "l3d_sg":
            side = int(np.ceil(fixed_features ** (1/3)))
            kw_args.update({'dim': (side, side, side), 'pflip': 0.35})
            sg = Lattice3D(**kw_args)
        case "er":
            kw_args.update({'n': fixed_features, 'p': 0.01, 'seed': 42})
            sg = ErdosRenyi(**kw_args)
        case "er_sg":
            kw_args.update({'n': fixed_features, 'p': 0.01, 'seed': 42, 'pflip': 0.35})
            sg = ErdosRenyi(**kw_args)
        case "er_dense":
            kw_args.update({'n': fixed_features, 'p': 0.05, 'seed': 42})
            sg = ErdosRenyi(**kw_args)
        case "l2d_sqr_rew":
            side = int(np.ceil(np.sqrt(fixed_features)))
            kw_args.update({'side1': side, 'geo': 'sqr', 'prew': 0.1, 'seed': 42})
            sg = Lattice2D(**kw_args)
        case _:
            raise ValueError(f"Unknown topology: {topology}")
    
    # Get expected nodes and show comparison
    expected_nodes = sg.get_expected_num_nodes()
    print(f"    📊 Expected: {expected_nodes} nodes (target: {fixed_features})")
    
    # Reinitialize without only_const_mode to actually create the network
    kw_args['only_const_mode'] = False
    
    match topology:
        case "l2d_square":
            sg = Lattice2D(**kw_args)
        case "l2d_hex":
            sg = Lattice2D(**kw_args)
        case "l2d_anti_tri":
            sg = Lattice2D(**kw_args)
            sg.flip_random_fract_edges()
        case "l3d":
            sg = Lattice3D(**kw_args)
        case "l3d_sg":
            sg = Lattice3D(**kw_args)
            sg.flip_random_fract_edges()
        case "er":
            sg = ErdosRenyi(**kw_args)
        case "er_sg":
            sg = ErdosRenyi(**kw_args)
            sg.flip_random_fract_edges()
        case "er_dense":
            sg = ErdosRenyi(**kw_args)
        case "l2d_sqr_rew":
            sg = Lattice2D(**kw_args)
    
    # Compute spectral basis
    sg.compute_laplacian_spectrum_weigV()
    basis = sg.get_sgspect_basis(max_factor=1)
    
    print(f"    ✅ Actual network: {sg.N} nodes, Basis: {basis.shape}")
    
    return basis, sg


def _compute_protein_mse_curves(proteins: List[Dict], basis: np.ndarray, topology_name: str) -> List[List[float]]:
    """Compute MSE curves for proteins using the given basis."""
    protein_mse_curves = []
    
    for protein in proteins:
        coords = protein['coords']
        pdb_id = protein['pdb_id']
        
        # Flatten coordinates to feature vector
        features = coords.flatten()
        
        # Pad or truncate to match basis size
        if len(features) < basis.shape[0]:
            padded_features = np.zeros(basis.shape[0])
            padded_features[:len(features)] = features
            features = padded_features
        else:
            features = features[:basis.shape[0]]
        
        # Compute MSE curve using i/N eigenvectors
        mse_curve = []
        N = basis.shape[1]
        
        for i in range(0, N, max(1, N//100)):  # Sample 100 points
            # Reconstruct using first i eigenvectors
            if i == 0:
                reconstructed = np.zeros_like(features)
            else:
                basis_i = basis[:, :i]
                coeffs = basis_i.T @ features
                reconstructed = basis_i @ coeffs
            
            # Compute MSE
            mse = np.mean((features - reconstructed)**2)
            mse_curve.append(mse)
        
        protein_mse_curves.append(mse_curve)
        print(f"    ✅ {pdb_id}: MSE curve computed")
    
    return protein_mse_curves


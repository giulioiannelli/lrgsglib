"""Feature extraction utilities for protein coordinates."""

from __future__ import annotations

import os
import tempfile
from typing import Optional, Tuple, List, Dict, Any

import numpy as np
from scipy.spatial.distance import pdist, squareform

# Bio.PDB imports for PDB processing
try:
    from Bio.PDB import PDBParser
    BIO_PDB_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    BIO_PDB_AVAILABLE = False

__all__ = [
    "protein_to_distance_feature_vector",
    "create_enhanced_feature_vector",
    "assign_secondary_structure_from_coords",
    "extract_comprehensive_protein_features",
    "create_comprehensive_feature_vector",
    "protein_to_coordinate_feature_vector",
]


def protein_to_distance_feature_vector(coords: np.ndarray, n_features: int = 784) -> np.ndarray:
    """Convert coordinates to a fixed-size distance feature vector."""
    if len(coords) == 0:
        return np.zeros(n_features)
    coords_centered = coords - coords.mean(axis=0)
    dist_matrix = squareform(pdist(coords_centered))
    upper_tri_indices = np.triu_indices(len(coords), k=1)
    distances = dist_matrix[upper_tri_indices]
    if len(distances) >= n_features:
        feature_vector = distances[:n_features]
    else:
        feature_vector = np.zeros(n_features)
        feature_vector[: len(distances)] = distances
    max_val = feature_vector.max()
    if max_val > 0:
        feature_vector /= max_val
    return feature_vector


def create_enhanced_feature_vector(
    coords: np.ndarray,
    residue_types: Optional[List[str]] = None,
    secondary_structure: Optional[List[str]] = None,
    n_features: int = 1024,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Create an enhanced feature vector including coordinates and biology."""
    coords_centered = coords - coords.mean(axis=0)
    coord_mean = coords.mean(axis=0)
    coord_feature_size = int(n_features * 0.75)
    coords_flat = coords_centered.flatten()
    coord_features = np.zeros(coord_feature_size)
    if len(coords_flat) <= coord_feature_size:
        coord_features[: len(coords_flat)] = coords_flat
        subsampling_indices = None
    else:
        indices = np.linspace(0, len(coords_flat) - 1, coord_feature_size, dtype=int)
        coord_features = coords_flat[indices]
        subsampling_indices = indices
    bio_feature_size = n_features - coord_feature_size
    bio_features = np.zeros(bio_feature_size)
    if secondary_structure:
        ss_dict = {"H": 0, "E": 1, "C": 2}
        for ss in secondary_structure:
            if ss in ss_dict:
                bio_features[ss_dict[ss]] += 1
        bio_features[:3] /= len(secondary_structure)
    if residue_types:
        aa_dict = {
            "ALA": 0,
            "ARG": 1,
            "ASN": 2,
            "ASP": 3,
            "CYS": 4,
            "GLN": 5,
            "GLU": 6,
            "GLY": 7,
            "HIS": 8,
            "ILE": 9,
            "LEU": 10,
            "LYS": 11,
            "MET": 12,
            "PHE": 13,
            "PRO": 14,
            "SER": 15,
            "THR": 16,
            "TRP": 17,
            "TYR": 18,
            "VAL": 19,
        }
        for residue in residue_types:
            if residue in aa_dict:
                bio_features[3 + aa_dict[residue]] += 1
        bio_features[3:23] /= len(residue_types)
    if len(coords_centered) > 0:
        rog = np.sqrt(np.mean(np.sum(coords_centered ** 2, axis=1)))
        bio_features[23] = rog / 50.0
        if len(coords_centered) > 2:
            cov_matrix = np.cov(coords_centered.T)
            eigenvals = np.linalg.eigvals(cov_matrix)
            eigenvals = np.sort(eigenvals)[::-1]
            bio_features[24:27] = eigenvals / (eigenvals[0] + 1e-8)
    full_features = np.concatenate([coord_features, bio_features])
    feature_scale = np.max(np.abs(full_features))
    if feature_scale > 0:
        full_features = full_features / feature_scale
    else:
        feature_scale = 1.0
    metadata = {
        "coord_mean": coord_mean,
        "coords_shape": coords.shape,
        "coord_feature_size": coord_feature_size,
        "feature_scale": feature_scale,
        "subsampling_indices": subsampling_indices,
        "original_coords": coords,
        "residue_types": residue_types,
        "secondary_structure": secondary_structure,
    }
    return full_features, metadata




def assign_secondary_structure_from_coords(coords: np.ndarray) -> List[str]:
    """Assign a rudimentary secondary structure from geometry."""
    n_residues = len(coords)
    ss_assignments = ["C"] * n_residues
    if n_residues < 4:
        return ss_assignments
    for i in range(2, n_residues - 2):
        v1 = coords[i] - coords[i - 1]
        v2 = coords[i + 1] - coords[i]
        v3 = coords[i + 2] - coords[i + 1]
        v1_norm = v1 / (np.linalg.norm(v1) + 1e-8)
        v2_norm = v2 / (np.linalg.norm(v2) + 1e-8)
        v3_norm = v3 / (np.linalg.norm(v3) + 1e-8)
        angle1 = np.arccos(np.clip(np.dot(v1_norm, v2_norm), -1, 1))
        angle2 = np.arccos(np.clip(np.dot(v2_norm, v3_norm), -1, 1))
        avg_angle = (angle1 + angle2) / 2
        if avg_angle < np.pi / 3:
            ss_assignments[i] = "E"
        elif avg_angle > 2 * np.pi / 3:
            ss_assignments[i] = "H"
    smoothed = ss_assignments.copy()
    for i in range(1, n_residues - 1):
        if ss_assignments[i - 1] == ss_assignments[i + 1] and ss_assignments[i - 1] != "C":
            smoothed[i] = ss_assignments[i - 1]
    return smoothed


def extract_comprehensive_protein_features(
    pdb_content: str, pdb_id: str, n_features: int = 1024
) -> Tuple[np.ndarray, Dict[str, Any], np.ndarray]:
    """Parse a PDB string and produce a comprehensive feature vector."""
    if not BIO_PDB_AVAILABLE:
        raise ImportError("Bio.PDB required for comprehensive PDB processing. Install with: pip install biopython")

    parser = PDBParser(QUIET=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
        tmp.write(pdb_content)
        tmp_path = tmp.name
    try:
        structure = parser.get_structure("protein", tmp_path)
        coordinates, residue_types, secondary_structure = [], [], []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if "CA" in residue:
                        coordinates.append(residue["CA"].get_coord())
                        residue_types.append(residue.get_resname())
                        secondary_structure.append("C")
        if not coordinates:
            return np.zeros(n_features), {}, np.empty((0, 3))
        coordinates = np.array(coordinates)
        feature_vector = create_comprehensive_feature_vector(
            coordinates, residue_types, secondary_structure, n_features
        )
        struct_info = {
            "coordinates": coordinates,
            "residue_types": residue_types,
            "secondary_structure": secondary_structure,
            "pdb_id": pdb_id,
            "n_residues": len(coordinates),
        }
        return feature_vector, struct_info, coordinates
    except Exception as e:  # pragma: no cover - best effort
        print(f"Error processing {pdb_id}: {e}")
        return np.zeros(n_features), {}, np.empty((0, 3))
    finally:
        os.unlink(tmp_path)


def create_comprehensive_feature_vector(
    coords: np.ndarray,
    residue_types: List[str],
    secondary_structure: List[str],
    n_features: int,
) -> np.ndarray:
    """Create a comprehensive feature vector from coordinates and metadata."""
    if len(coords) == 0:
        return np.zeros(n_features)
    coords_centered = coords - coords.mean(axis=0)
    dist_matrix = squareform(pdist(coords_centered))
    upper_tri_indices = np.triu_indices(len(coords), k=1)
    distances = dist_matrix[upper_tri_indices]
    AA_DICT = {
        "ALA": 0,
        "ARG": 1,
        "ASN": 2,
        "ASP": 3,
        "CYS": 4,
        "GLN": 5,
        "GLU": 6,
        "GLY": 7,
        "HIS": 8,
        "ILE": 9,
        "LEU": 10,
        "LYS": 11,
        "MET": 12,
        "PHE": 13,
        "PRO": 14,
        "SER": 15,
        "THR": 16,
        "TRP": 17,
        "TYR": 18,
        "VAL": 19,
    }
    aa_composition = np.zeros(20)
    for residue in residue_types:
        if residue in AA_DICT:
            aa_composition[AA_DICT[residue]] += 1
    if len(residue_types) > 0:
        aa_composition /= len(residue_types)
    ss_composition = np.zeros(3)
    SS_DICT = {"H": 0, "E": 1, "C": 2}
    for ss in secondary_structure:
        if ss in SS_DICT:
            ss_composition[SS_DICT[ss]] += 1
    if len(secondary_structure) > 0:
        ss_composition /= len(secondary_structure)
    structural_features: List[float] = []
    if len(coords_centered) > 0:
        rog = np.sqrt(np.mean(np.sum(coords_centered ** 2, axis=1)))
        structural_features.append(rog)
        if len(coords_centered) > 2:
            cov_matrix = np.cov(coords_centered.T)
            eigenvals = np.sort(np.linalg.eigvals(cov_matrix))[::-1]
            asphericity = eigenvals[0] - 0.5 * (eigenvals[1] + eigenvals[2])
            structural_features.extend([asphericity, eigenvals[0], eigenvals[1], eigenvals[2]])
        contact_threshold = 8.0
        contacts = (dist_matrix < contact_threshold).sum()
        contact_density = contacts / (len(coords) * (len(coords) - 1)) if len(coords) > 1 else 0
        structural_features.append(contact_density)
    max_distance_features = n_features // 2
    distance_features = distances[: min(len(distances), max_distance_features)]
    all_features = np.concatenate([distance_features, aa_composition, ss_composition, structural_features])
    if len(all_features) >= n_features:
        feature_vector = all_features[:n_features]
    else:
        feature_vector = np.zeros(n_features)
        feature_vector[: len(all_features)] = all_features
    max_val = feature_vector.max()
    if max_val > 0:
        feature_vector /= max_val
    return feature_vector


def protein_to_coordinate_feature_vector(coords: np.ndarray, n_features: int = 1024) -> np.ndarray:
    """
    PROPER coordinate feature extraction from protein structure.
    
    This function extracts actual x,y,z coordinates as features, which can be 
    meaningfully reconstructed back to protein structures.
    
    Args:
        coords: Protein coordinates array of shape (n_atoms, 3)
        n_features: Target number of features
        
    Returns:
        Feature vector of length n_features containing coordinate information
    """
    if len(coords) == 0:
        return np.zeros(n_features)
    
    # Center the coordinates
    coords_centered = coords - coords.mean(axis=0)
    
    # Flatten to 1D: [x1, y1, z1, x2, y2, z2, ...]
    coords_flat = coords_centered.flatten()
    
    # Create feature vector
    if len(coords_flat) <= n_features:
        # Pad with zeros if we have fewer coordinates than features
        feature_vector = np.zeros(n_features)
        feature_vector[:len(coords_flat)] = coords_flat
    else:
        # Subsample coordinates if we have more than n_features
        # Take every k-th coordinate to maintain spatial structure
        step = len(coords_flat) // n_features
        indices = np.arange(0, len(coords_flat), step)[:n_features]
        feature_vector = coords_flat[indices]
    
    # Normalize to reasonable range
    max_coord = np.abs(feature_vector).max()
    if max_coord > 0:
        feature_vector = feature_vector / max_coord
    
    return feature_vector



"""Low level utilities for reconstructing coordinates from feature vectors."""

from __future__ import annotations

from typing import Dict, Any, Optional

import numpy as np

__all__ = [
    "reconstruct_coordinates_from_features",
    "safe_reconstruct_coordinates_from_features",
]


def reconstruct_coordinates_from_features(
    reconstructed_features: np.ndarray, metadata: Dict[str, Any]
) -> Optional[np.ndarray]:
    """Reconstruct coordinates from an enhanced feature vector."""
    if metadata is None:
        return None

    feature_scale = metadata["feature_scale"]
    denormalized = reconstructed_features * feature_scale
    coord_feature_size = metadata["coord_feature_size"]
    coord_features = denormalized[:coord_feature_size]
    coords_shape = metadata["coords_shape"]
    n_atoms, n_dims = coords_shape

    if metadata["subsampling_indices"] is not None:
        full_coords_flat = np.zeros(n_atoms * n_dims)
        indices = metadata["subsampling_indices"]
        full_coords_flat[indices] = coord_features
        for i in range(len(full_coords_flat)):
            if i not in indices:
                distances = np.abs(indices - i)
                nearest_idx = indices[np.argmin(distances)]
                nearest_val = coord_features[np.argmin(distances)]
                full_coords_flat[i] = nearest_val
        coords_reconstructed = full_coords_flat.reshape(coords_shape)
    else:
        coords_flat_size = n_atoms * n_dims
        if len(coord_features) >= coords_flat_size:
            coords_reconstructed = coord_features[:coords_flat_size].reshape(coords_shape)
        else:
            padded_coords = np.zeros(coords_flat_size)
            padded_coords[: len(coord_features)] = coord_features
            coords_reconstructed = padded_coords.reshape(coords_shape)

    coords_reconstructed += metadata["coord_mean"]
    return coords_reconstructed


def safe_reconstruct_coordinates_from_features(
    reconstructed_features: np.ndarray, metadata: Dict[str, Any]
) -> np.ndarray:
    """
    FIXED: Proper coordinate reconstruction from feature vectors.
    This function now correctly handles coordinate features, not distances.
    """
    if metadata is None:
        return np.zeros((10, 3))
    
    try:
        # Handle the case where metadata is from coordinate-based feature extraction
        if 'normalization_factor' in metadata:
            # This is from our new coordinate-based approach
            coords_flat = reconstructed_features * metadata['normalization_factor']
            
            if 'subsampling_indices' in metadata:
                # Reconstruct from subsampled coordinates
                full_coords_flat = np.zeros(metadata['coords_flat_size'])
                indices = metadata['subsampling_indices']
                n_available = min(len(coords_flat), len(indices))
                full_coords_flat[indices[:n_available]] = coords_flat[:n_available]
                
                # Interpolate missing coordinates
                for i in range(metadata['coords_flat_size']):
                    if i not in indices[:n_available]:
                        distances = np.abs(indices[:n_available] - i)
                        if len(distances) > 0:
                            nearest_idx = np.argmin(distances)
                            full_coords_flat[i] = coords_flat[nearest_idx]
                coords_flat = full_coords_flat
            else:
                # Remove padding
                coords_flat = coords_flat[:metadata['coords_flat_size']]
            
            # Reshape to coordinates
            coords_reconstructed = coords_flat.reshape(metadata['n_atoms'], 3)
            coords_reconstructed += metadata['coord_mean']
            return coords_reconstructed
        
        # Legacy handling for old metadata format (enhanced features)
        else:
            feature_scale = metadata.get("feature_scale", 1.0)
            denormalized = reconstructed_features * feature_scale
            coord_feature_size = metadata.get("coord_feature_size", len(reconstructed_features))
            coord_features = denormalized[:coord_feature_size]
            coords_shape = metadata.get("coords_shape", (len(coord_features)//3, 3))
            n_atoms, n_dims = coords_shape
            
            if metadata.get("subsampling_indices") is not None:
                full_coords_flat = np.zeros(n_atoms * n_dims)
                indices = metadata["subsampling_indices"]
                max_coords_to_use = min(len(coord_features), len(indices))
                full_coords_flat[indices[:max_coords_to_use]] = coord_features[:max_coords_to_use]
                for i in range(len(full_coords_flat)):
                    if i not in indices[:max_coords_to_use]:
                        if len(indices[:max_coords_to_use]) > 0:
                            distances = np.abs(indices[:max_coords_to_use] - i)
                            nearest_idx_pos = np.argmin(distances)
                            if nearest_idx_pos < len(coord_features):
                                full_coords_flat[i] = coord_features[nearest_idx_pos]
                coords_reconstructed = full_coords_flat.reshape(coords_shape)
            else:
                coords_flat_size = n_atoms * n_dims
                if len(coord_features) >= coords_flat_size:
                    coords_reconstructed = coord_features[:coords_flat_size].reshape(coords_shape)
                else:
                    padded_coords = np.zeros(coords_flat_size)
                    padded_coords[: len(coord_features)] = coord_features
                    coords_reconstructed = padded_coords.reshape(coords_shape)
            
            coord_mean = metadata.get("coord_mean", np.zeros(3))
            coords_reconstructed += coord_mean
            return coords_reconstructed
            
    except Exception as e:  # pragma: no cover - best effort
        print(f"  Warning: Coordinate reconstruction failed ({e})")
        print("  Returning original coordinates for visualization")
        return metadata.get("original_coords", np.zeros((10, 3)))

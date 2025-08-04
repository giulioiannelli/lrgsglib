import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde
import plotly.graph_objects as go


def create_volume_visualization(data_sets, X, Y, Z, colors=None, colors_solid=None,
                               # Volume appearance parameters
                               opacity=0.4,
                               surface_count=2,
                               # Density filtering parameters
                               significance_threshold=0.05,
                               contrast_power=1.2,
                               min_threshold_percentile=60,
                               max_threshold_percentile=90,
                               # Point visualization parameters
                               sample_size=300,
                               min_point_size=4,
                               max_point_size=25):
    """
    Create 3D volume visualization with full parameter control.
    
    Parameters:
    -----------
    data_sets : list of tuples
        List of (x_data, y_data, z_data, color_name, label) tuples for each class
    X, Y, Z : ndarray
        3D coordinate grids for density evaluation
    colors : list, optional
        List of color strings for volume rendering (with alpha). 
        Default: ['rgba(255,0,0,0.3)', 'rgba(0,255,0,0.3)', 'rgba(0,0,255,0.3)']
    colors_solid : list, optional
        List of solid color strings for scatter points.
        Default: ['red', 'green', 'blue']
    
    KEY PARAMETERS FOR TUNING:
    =========================
    
    🚫 TO AVOID BIG SURFACES AROUND ISOLATED POINTS:
    ------------------------------------------------
    significance_threshold=0.1-0.2  # Higher = stricter filtering
    min_threshold_percentile=70-80  # Higher = only show dense cores
    
    🎨 VOLUME APPEARANCE:
    --------------------
    opacity=0.3-0.6                # Lower = more transparent
    surface_count=1-3               # Fewer = cleaner appearance
    significance_threshold=0.01-0.3 # Lower = more volume, Higher = less volume
    
    ✨ CONTRAST & SHARPNESS:
    -----------------------
    contrast_power=0.8-2.0          # Higher = sharper contrast
    min_threshold_percentile=40-80  # Higher = fewer surfaces
    max_threshold_percentile=85-95  # Controls surface extent
    
    📊 POINT VISUALIZATION:
    ----------------------
    sample_size=200-500             # More points = more detail
    min_point_size=3-6              # Smallest point size
    max_point_size=15-30            # Largest point size
    
    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The 3D volume visualization figure
    """
    
    # Set default colors if not provided
    if colors is None:
        colors = ['rgba(255,0,0,0.3)', 'rgba(0,255,0,0.3)', 'rgba(0,0,255,0.3)']
    if colors_solid is None:
        colors_solid = ['red', 'green', 'blue']
    
    fig = go.Figure()
    
    for idx, (x_data, y_data, z_data, color_name, label) in enumerate(data_sets):
        
        # Step 1: Create 3D density volume
        density = create_3d_density_volume(x_data, y_data, z_data, X, Y, Z)
        
        # Check if we got meaningful density values
        if density.max() < 1e-10:
            continue
        
        # Step 2: Normalize density with custom parameters
        density_norm, iso_min, iso_max = process_density_3d_visualization(
            density,
            significance_threshold=significance_threshold,
            contrast_power=contrast_power,
            min_threshold_percentile=min_threshold_percentile,
            max_threshold_percentile=max_threshold_percentile
        )
        
        # Step 3: Create volume plot with custom settings
        color_idx = idx % len(colors)  # Handle cases with more data_sets than colors
        fig.add_trace(go.Volume(
            x=X.flatten(),
            y=Y.flatten(),
            z=Z.flatten(),
            value=density_norm.flatten(),
            isomin=iso_min,
            isomax=iso_max,
            opacity=opacity,
            surface_count=surface_count,
            colorscale=[[0, colors[color_idx]], [1, colors[color_idx]]],
            name=label,
            showscale=False,
            caps=dict(x_show=False, y_show=False, z_show=False),
        ))
        
        # Step 4: Add density-weighted scatter points
        sample_idx, point_sizes, normalized_densities = calculate_point_densities(
            x_data, y_data, z_data,
            sample_size=sample_size,
            min_point_size=min_point_size,
            max_point_size=max_point_size
        )
        
        color_solid_idx = idx % len(colors_solid)  # Handle cases with more data_sets than colors
        fig.add_trace(go.Scatter3d(
            x=x_data[sample_idx],
            y=y_data[sample_idx],
            z=z_data[sample_idx],
            mode='markers',
            marker=dict(
                size=point_sizes,
                color=colors_solid[color_solid_idx],
                opacity=0.8,
                line=dict(width=0.5, color='white')
            ),
            name=f'{label} (density-weighted)',
            showlegend=True,
            hovertemplate=f'<b>{label}</b><br>' +
                         'X: %{x:.3f}<br>' +
                         'Y: %{y:.3f}<br>' +
                         'Z: %{z:.3f}<br>' +
                         'Density: %{text:.3f}<extra></extra>',
            text=normalized_densities
        ))
    
    return fig


def create_3d_density_volume(x_data, y_data, z_data, X, Y, Z):
    """
    Create 3D density volume using kernel density estimation.
    
    Parameters:
    -----------
    x_data, y_data, z_data : array-like
        Coordinates of data points
    X, Y, Z : ndarray
        3D coordinate grids for density evaluation
        
    Returns:
    --------
    density : ndarray
        3D density volume on the coordinate grid
    """
    # Combine coordinates and ensure they are finite
    points = np.vstack([x_data, y_data, z_data])
    
    # Remove any infinite or NaN values
    valid_mask = np.isfinite(points).all(axis=0)
    if not np.any(valid_mask):
        return np.zeros(X.shape)
    
    points = points[:, valid_mask]
    # print(f"    Using {points.shape[1]} valid points out of {len(x_data)}")
    
    # Calculate data scale to adjust bandwidth appropriately
    data_scale = np.std(points, axis=1)
    mean_scale = np.mean(data_scale)
    # print(f"    Data scale (std): {data_scale}, mean scale: {mean_scale:.4f}")
    
    # Handle case with too few points
    if points.shape[1] < 10:
        # print("    Too few points, using simple Gaussian")
        mean_point = np.mean(points, axis=1)
        grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
        distances = np.linalg.norm(grid_points - mean_point[:, np.newaxis], axis=0)
        density = np.exp(-distances**2 / (2 * (mean_scale * 0.5)**2))
        return density.reshape(X.shape)
    
    try:
        # Use automatic bandwidth selection with adjustment
        kde = gaussian_kde(points)
        
        # Adjust bandwidth if too small
        if hasattr(kde, 'factor') and kde.factor < 0.1:
            # print(f"    Auto bandwidth too small ({kde.factor:.6f}), adjusting...")
            kde.set_bandwidth(0.15)
        
        # print(f"    Using bandwidth: {kde.factor:.6f}")
        
        # Evaluate KDE on grid
        grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
        density = kde(grid_points).reshape(X.shape)
        
        # Apply light smoothing
        density = gaussian_filter(density, sigma=0.3)
        density = np.nan_to_num(density, nan=0.0, posinf=0.0, neginf=0.0)
        
        # print(f"    KDE successful, density range: {density.min():.8f} to {density.max():.8f}")
        
        # Check for zero densities
        if density.max() < 1e-10:
            # print("    KDE produced zero densities, using fallback...")
            raise ValueError("Zero density from KDE")
        
    except Exception as e:
        # print(f"    KDE failed, using fallback method: {e}")
        
        # Fallback: gaussian mixture around data points
        grid_points = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
        density = np.zeros(grid_points.shape[1])
        
        # Use data-scale adjusted bandwidth
        bandwidth = max(0.3, mean_scale * 0.4)
        # print(f"    Using fallback bandwidth: {bandwidth:.4f}")
        
        # Sample points to avoid slowdown
        sample_points = points[:, ::max(1, points.shape[1]//200)]
        # print(f"    Using {sample_points.shape[1]} sample points for fallback")
        
        for point in sample_points.T:
            distances = np.linalg.norm(grid_points - point[:, np.newaxis], axis=0)
            density += np.exp(-distances**2 / (2 * bandwidth**2))
        
        density = density.reshape(X.shape)
        # print(f"    Fallback successful, density range: {density.min():.8f} to {density.max():.8f}")
    
    return density


def process_density_3d_visualization(density, 
                                     significance_threshold=0.05,  # Filter isolated points
                                     contrast_power=1.2,           # Enhance contrast
                                     min_threshold_percentile=60,  # Conservative thresholds
                                     max_threshold_percentile=90,
                                     default_iso_min=0.3,         # Fallback values
                                     default_iso_max=0.9,
                                     min_separation=0.2):          # Minimum iso separation
    """
    Normalize density values for visualization with full parameter control.
    
    Parameters:
    -----------
    density : ndarray
        Raw density values
    significance_threshold : float, default=0.05
        Percentage of max density below which points are filtered out (0.0-1.0)
        - Higher values (0.1-0.2): Very strict, only core regions shown
        - Medium values (0.05-0.1): Balanced filtering
        - Lower values (0.01-0.05): More permissive, shows more outliers
        - 0.0: No filtering (shows all points including isolated ones)
    
    contrast_power : float, default=1.2
        Power transform for contrast enhancement (>0)
        - Higher values (1.5-2.0): Strong contrast, emphasizes cores
        - Medium values (1.0-1.5): Moderate contrast enhancement
        - Lower values (0.5-1.0): Softer contrast
        - 1.0: No contrast enhancement (linear)
    
    min_threshold_percentile : float, default=60
        Percentile for isomin threshold (0-100)
        - Higher values (70-90): Only show highest density regions
        - Medium values (40-70): Balanced coverage
        - Lower values (10-40): Show more low-density regions
    
    max_threshold_percentile : float, default=90
        Percentile for isomax threshold (0-100)
        - Should be > min_threshold_percentile
        - 95-99: Include almost all significant density
        - 85-95: Exclude highest peaks for cleaner view
    
    default_iso_min/max : float, default=0.3/0.9
        Fallback threshold values when percentile calculation fails
        
    min_separation : float, default=0.2
        Minimum separation between iso_min and iso_max
        - Higher values (0.3-0.5): Ensure good visual separation
        - Lower values (0.1-0.2): Allow closer thresholds
        
    Returns:
    --------
    density_norm : ndarray
        Normalized density values [0, 1]
    iso_min, iso_max : float
        Adaptive threshold values for isosurfaces
    """
    density_min = density.min()
    density_max = density.max()
    density_range = density_max - density_min
    
    # print(f"  Raw density range: {density_min:.8f} to {density_max:.8f} (range: {density_range:.8f})")
    
    # Apply significance threshold to filter out isolated points
    threshold_value = density_max * significance_threshold
    # print(f"  Significance threshold ({significance_threshold*100:.1f}% of max): {threshold_value:.8f}")
    
    # Zero out low-significance regions
    density_filtered = np.where(density > threshold_value, density, 0)
    
    # Count filtering results
    kept_points = np.sum(density_filtered > 0)
    total_points = density.size
    # print(f"  Kept {kept_points}/{total_points} points ({100*kept_points/total_points:.1f}%)")
    # print(f"  Filtered out {total_points-kept_points} low-density points")
    
    # Recalculate range after filtering
    density_min_filt = density_filtered.min()
    density_max_filt = density_filtered.max()
    density_range_filt = density_max_filt - density_min_filt
    
    # Smart normalization on filtered data
    if density_range_filt > 1e-10:
        density_norm = (density_filtered - density_min_filt) / density_range_filt
    else:
        # print(f"  No significant density variation after filtering")
        return np.zeros_like(density), default_iso_min, default_iso_max
    
    # Apply contrast enhancement
    density_norm = np.power(density_norm, contrast_power)
    
    # print(f"  Filtered density range: {density_norm.min():.6f} to {density_norm.max():.6f}")
    # print(f"  Non-zero values after filtering: {np.sum(density_norm > 0)} / {density_norm.size}")
    
    # Calculate adaptive thresholds
    non_zero_density = density_norm[density_norm > 0.01]
    if len(non_zero_density) > 50:
        iso_min = np.percentile(non_zero_density, min_threshold_percentile)
        iso_max = np.percentile(non_zero_density, max_threshold_percentile)
        # print(f"  Using percentile thresholds: {iso_min:.6f} to {iso_max:.6f}")
    else:
        iso_min = default_iso_min
        iso_max = default_iso_max
        # print(f"  Using default thresholds: {iso_min} to {iso_max}")
    
    # Ensure minimum separation
    if iso_max - iso_min < min_separation:
        iso_max = min(1.0, iso_min + min_separation)
        # print(f"  Adjusted thresholds for separation: {iso_min:.6f} to {iso_max:.6f}")
    
    return density_norm, iso_min, iso_max


def calculate_point_densities(x_data, y_data, z_data, 
                             sample_size=300,
                             min_point_size=4,
                             max_point_size=25,
                             density_radius=0.3):
    """
    Calculate local density for sampled points with full parameter control.
    
    Parameters:
    -----------
    x_data, y_data, z_data : array-like
        Point coordinates
    sample_size : int, default=300
        Number of points to sample for visualization
        - Higher values (500-1000): More detailed representation
        - Lower values (100-300): Faster rendering, less detail
        
    min_point_size : float, default=4
        Minimum point size for low-density points
        
    max_point_size : float, default=25  
        Maximum point size for high-density points
        
    density_radius : float, default=0.3
        Radius for fallback distance-based density calculation
        - Larger values: Smoother density estimates
        - Smaller values: More local density variations
        
    Returns:
    --------
    sample_idx : ndarray
        Indices of sampled points
    point_sizes : ndarray
        Point sizes proportional to density
    normalized_densities : ndarray
        Normalized density values for hover text
    """
    sample_size = min(sample_size, len(x_data))
    sample_idx = np.random.choice(len(x_data), sample_size, replace=False)
    
    sample_points = np.column_stack([x_data[sample_idx], y_data[sample_idx], z_data[sample_idx]])
    points_for_kde = np.vstack([x_data, y_data, z_data])
    
    try:
        kde_for_points = gaussian_kde(points_for_kde)
        
        # Adjust bandwidth if too small
        if hasattr(kde_for_points, 'factor') and kde_for_points.factor < 0.1:
            kde_for_points.set_bandwidth(0.15)
            
        local_densities = kde_for_points(sample_points.T)
        
        if local_densities.max() < 1e-10:
            raise ValueError("Zero densities from KDE")
            
    except Exception as e:
        # print(f"  Using distance-based density for points: {e}")
        # Fallback: count nearby points
        local_densities = np.zeros(len(sample_points))
        for i, point in enumerate(sample_points):
            distances = np.linalg.norm(points_for_kde.T - point, axis=1)
            nearby_count = np.sum(distances < density_radius)
            local_densities[i] = nearby_count
    
    # Normalize densities
    min_density, max_density = local_densities.min(), local_densities.max()
    
    if max_density - min_density > 1e-10:
        normalized_densities = (local_densities - min_density) / (max_density - min_density)
    else:
        # Random variation if all densities are the same
        normalized_densities = 0.5 + 0.2 * (np.random.random(len(local_densities)) - 0.5)
    
    normalized_densities = np.clip(normalized_densities, 0.0, 1.0)
    point_sizes = min_point_size + (max_point_size - min_point_size) * normalized_densities
    
    return sample_idx, point_sizes, normalized_densities


def create_volume_visualization_simple(data_sets, 
                                     grid_resolution=30,
                                     grid_padding=0.2,
                                     colors=None,
                                     colors_solid=None,
                                     **kwargs):
    """
    Simplified wrapper for create_volume_visualization that automatically creates the grid.
    
    Parameters:
    -----------
    data_sets : list of tuples
        List of (x_data, y_data, z_data, color_name, label) tuples for each class
    grid_resolution : int, default=30
        Resolution of the 3D grid for density estimation
    grid_padding : float, default=0.2
        Padding around the data bounds for the grid
    colors : list, optional
        List of color strings for volume rendering (with alpha)
    colors_solid : list, optional
        List of solid color strings for scatter points
    **kwargs : dict
        Additional parameters passed to create_volume_visualization
        
    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The 3D volume visualization figure
    """
    
    # Extract all coordinates
    all_x = np.concatenate([data[0] for data in data_sets])
    all_y = np.concatenate([data[1] for data in data_sets])
    all_z = np.concatenate([data[2] for data in data_sets])
    
    # Define grid bounds
    x_min, x_max = all_x.min() - grid_padding, all_x.max() + grid_padding
    y_min, y_max = all_y.min() - grid_padding, all_y.max() + grid_padding
    z_min, z_max = all_z.min() - grid_padding, all_z.max() + grid_padding
    
    # Create coordinate grids
    x_grid = np.linspace(x_min, x_max, grid_resolution)
    y_grid = np.linspace(y_min, y_max, grid_resolution)
    z_grid = np.linspace(z_min, z_max, grid_resolution)
    X, Y, Z = np.meshgrid(x_grid, y_grid, z_grid, indexing='ij')
    
    # Set default colors if not provided
    if colors is None:
        colors = ['rgba(255,0,0,0.3)', 'rgba(0,255,0,0.3)', 'rgba(0,0,255,0.3)']
    if colors_solid is None:
        colors_solid = ['red', 'green', 'blue']
    
    # Call the main function
    return create_volume_visualization(
        data_sets, X, Y, Z, colors, colors_solid, **kwargs
    )


def update_camera_view(fig, zoom_level=1.0, rotation_x=0, rotation_y=0, rotation_z=0):
    """
    Update the camera view with zoom and rotation controls.
    
    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
        The figure to update
    zoom_level : float, default=1.0
        Zoom level (0.5 = very close, 1.0 = normal, 2.0 = far away)
    rotation_x : float, default=0
        Rotation angle in degrees around X-axis (pitch - tilt up/down)
    rotation_y : float, default=0
        Rotation angle in degrees around Y-axis (yaw - turn left/right)
    rotation_z : float, default=0
        Rotation angle in degrees around Z-axis (roll - rotate clockwise/counter-clockwise)
        
    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The updated figure with new camera position
        
    Examples:
    ---------
    # Zoom in close
    fig_zoomed = update_camera_view(fig, zoom_level=2.0)
    
    # Rotate to front angle
    fig_front = update_camera_view(fig, zoom_level=1.2, rotation_y=45)
    
    # Top-down perspective
    fig_top = update_camera_view(fig, zoom_level=1.5, rotation_x=70)
    
    # Custom view (zoomed + rotated)
    fig_custom = update_camera_view(fig, zoom_level=1.8, rotation_x=20, rotation_y=30)
    """
    
    # Convert degrees to radians
    rx, ry, rz = np.radians([rotation_x, rotation_y, rotation_z])
    
    # Base camera position (adjust this for different base views)
    base_eye = np.array([1.5, 1.5, 1.5])
    
    # Apply zoom (smaller values = closer)
    zoomed_eye = base_eye / zoom_level
    
    # Apply rotations (simplified rotation around origin)
    # X rotation (pitch)
    cos_x, sin_x = np.cos(rx), np.sin(rx)
    y_rot = zoomed_eye[1] * cos_x - zoomed_eye[2] * sin_x
    z_rot = zoomed_eye[1] * sin_x + zoomed_eye[2] * cos_x
    zoomed_eye[1], zoomed_eye[2] = y_rot, z_rot
    
    # Y rotation (yaw)
    cos_y, sin_y = np.cos(ry), np.sin(ry)
    x_rot = zoomed_eye[0] * cos_y + zoomed_eye[2] * sin_y
    z_rot = -zoomed_eye[0] * sin_y + zoomed_eye[2] * cos_y
    zoomed_eye[0], zoomed_eye[2] = x_rot, z_rot
    
    # Z rotation (roll)
    cos_z, sin_z = np.cos(rz), np.sin(rz)
    x_rot = zoomed_eye[0] * cos_z - zoomed_eye[1] * sin_z
    y_rot = zoomed_eye[0] * sin_z + zoomed_eye[1] * cos_z
    zoomed_eye[0], zoomed_eye[1] = x_rot, y_rot
    
    # Update the figure
    fig.update_layout(
        scene_camera=dict(
            eye=dict(x=zoomed_eye[0], y=zoomed_eye[1], z=zoomed_eye[2]),
            center=dict(x=0, y=0, z=0),
            up=dict(x=0, y=0, z=1)
        )
    )
    
    return fig


# Additional camera presets for convenience
CAMERA_PRESETS = {
    'default': dict(eye=dict(x=1.5, y=1.5, z=1.5)),
    'front_view': dict(eye=dict(x=0, y=-2.5, z=0)),
    'side_view': dict(eye=dict(x=2.5, y=0, z=0)),
    'top_view': dict(eye=dict(x=0, y=0, z=2.5)),
    'diagonal_high': dict(eye=dict(x=2, y=2, z=2)),
    'paper_view': dict(eye=dict(x=1.8, y=1.2, z=1.5), center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1)),
    'presentation': dict(eye=dict(x=2.2, y=1.8, z=1.2), center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1)),
    'close_up': dict(eye=dict(x=1.0, y=1.0, z=1.0), center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1)),
    'far_out': dict(eye=dict(x=3.0, y=3.0, z=3.0), center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1))
}


def apply_camera_preset(fig, preset_name='default'):
    """
    Apply a predefined camera preset to the figure.
    
    Parameters:
    -----------
    fig : plotly.graph_objects.Figure
        The figure to update
    preset_name : str, default='default'
        Name of the preset to apply. Available presets:
        'default', 'front_view', 'side_view', 'top_view', 'diagonal_high',
        'paper_view', 'presentation', 'close_up', 'far_out'
        
    Returns:
    --------
    fig : plotly.graph_objects.Figure
        The updated figure with new camera position
    """
    if preset_name not in CAMERA_PRESETS:
        print(f"Warning: Preset '{preset_name}' not found. Available presets: {list(CAMERA_PRESETS.keys())}")
        preset_name = 'default'
    
    fig.update_layout(scene_camera=CAMERA_PRESETS[preset_name])
    return fig

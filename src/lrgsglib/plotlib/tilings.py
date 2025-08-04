import numpy as np
#
from typing import Any, Union
from numpy.typing import NDArray
#
from .const_plotlib import *
#
__all__ = [
    'plot_honeycomb_grid',
    'plot_honeycomb_grid_fast',
    'plot_hex_tiling_from_nodes',
    'plot_hex_tiling_from_pos',
    'plot_octagonal_square_tiling_from_nodes',
    'plot_octagonal_square_tiling_from_pos',
]
#
def plot_honeycomb_grid(
    data: NDArray,
    ax: Any,
    triangle_size: float = 1.0,
    cmap: Union[Colormap, str] = 'viridis',
    rotation: float = 0.0
) -> None:
    """
    Plot a rotated honeycomb (triangular) grid with equal side lengths.

    Parameters
    ----------
    data : NDArray
        2D array of values.
    triangle_size : float
        Side length of each equilateral triangle.
    cmap : Colormap or str
    rotation : float
        Angle in degrees to rotate the entire grid (counter-clockwise).
    """
    n_rows, n_cols = data.shape
    cmap_fun = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    h = (np.sqrt(3) / 2) * triangle_size

    # prepare extents
    min_x, max_x = np.inf, -np.inf
    min_y, max_y = np.inf, -np.inf
    theta = np.deg2rad(rotation)
    c, s = np.cos(theta), np.sin(theta)

    for i in range(n_rows):
        for j in range(n_cols):
            # base position
            x0 = j * triangle_size / 2
            y0 = i * h

            # pick orientation
            if (i + j) % 2 == 0:
                verts = [
                    ( x0,       y0 + h/2),
                    ( x0 - triangle_size/2, y0 - h/2),
                    ( x0 + triangle_size/2, y0 - h/2),
                ]
            else:
                verts = [
                    ( x0,       y0 - h/2),
                    ( x0 - triangle_size/2, y0 + h/2),
                    ( x0 + triangle_size/2, y0 + h/2),
                ]

            # rotate verts
            verts_r = [ (x*c - y*s, x*s + y*c) for x, y in verts ]

            # update extents
            xs, ys = zip(*verts_r)
            min_x, max_x = min(min_x, *xs), max(max_x, *xs)
            min_y, max_y = min(min_y, *ys), max(max_y, *ys)

            # draw
            color = cmap_fun(data[i, j])
            ax.add_patch(Polygon(verts_r, facecolor=color, edgecolor='none'))

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_aspect('equal')


def plot_hex_tiling_from_nodes(x, y, vals, ax=None, cmap='viridis'):
    """
    Creates a tiling of pointy-topped hexagons (one per node) on a rectangular grid.
    
    The tiling is computed so that the overall pattern exactly fits the bounding box
    defined by the input x and y coordinates.
    
    For pointy-topped hexagons:
      - Width = √3 * R, Height = 2 * R.
      - Horizontal spacing = √3 * R.
      - Vertical spacing = 1.5 * R.
      - Even rows: center_x = x_min + (√3/2)*R + col*(√3*R)
      - Odd rows:  center_x = x_min + √3*R + col*(√3*R)
      - For all rows: center_y = y_min + R + row*(1.5*R)
    
    Parameters:
      x, y : array-like
          Node coordinates.
      vals : array-like
          A value per node used for coloring.
      ax : matplotlib.axes.Axes, optional
          Axis to plot on; a new figure is created if None.
      cmap : str or Colormap, optional
          Colormap for the hexagons.
    
    Returns:
      ax : matplotlib.axes.Axes
          The axis with the hexagon tiling.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
    # Use actual node positions as hexagon centers directly
    hex_centers = np.column_stack([x, y])
    
    # Bounding box of input data.
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    
    # Calculate radius based on nearest neighbor distance
    # For triangular lattice, find minimum distance between points
    min_dist = np.inf
    n_sample = min(100, len(hex_centers))  # Sample for efficiency
    for i in range(n_sample):
        for j in range(i+1, n_sample):
            dist = np.sqrt((hex_centers[i, 0] - hex_centers[j, 0])**2 + 
                          (hex_centers[i, 1] - hex_centers[j, 1])**2)
            if dist < min_dist:
                min_dist = dist
    R = min_dist * 0.6  # Scale factor for good tiling
    
    # Sort centers and node data lexicographically (by x then y).
    hex_sort_idx = np.lexsort((hex_centers[:,1], hex_centers[:,0]))
    node_sort_idx = np.lexsort((np.array(y), np.array(x)))
    hex_centers_sorted = hex_centers[hex_sort_idx]
    vals_sorted = np.array(vals)[node_sort_idx]
    
    norm = Normalize(vmin=np.min(vals), vmax=np.max(vals))
    colormap = plt.get_cmap(cmap)
    
    # Build one hexagon template (pointy-topped) - much faster approach
    angles = np.linspace(0, 2*np.pi, 7)[:-1] + np.pi/6
    base_hex = np.column_stack([np.cos(angles)*R, np.sin(angles)*R])

    # Generate all hexagon polygons using the same positions
    polys = [
        base_hex + np.array([cx, cy])
        for cx, cy in hex_centers_sorted
    ]

    # Map values to colors
    facecolors = colormap(norm(vals_sorted))

    # Create & add a single PolyCollection - much faster than individual patches
    coll = PolyCollection(
        polys,
        facecolors=facecolors,
        edgecolors='none',
        antialiased=False
    )
    ax.add_collection(coll)
    
    # Compute exact extents to fill the axis without padding
    left = np.min(hex_centers_sorted[:,0] - (np.sqrt(3)/2)*R)
    right = np.max(hex_centers_sorted[:,0] + (np.sqrt(3)/2)*R)
    bottom = np.min(hex_centers_sorted[:,1] - R)
    top = np.max(hex_centers_sorted[:,1] + R)
    
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    ax.set_aspect('equal')

    return ax

def plot_hex_tiling_from_pos(pos, vals, ax=None, cmap='viridis'):
    """
    Extracts x and y coordinates from a position dictionary and calls
    plot_hex_tiling_from_nodes to create a hexagon tiling.

    Parameters:
      pos : dict
          Dictionary of node positions (e.g. {node: (x, y), ...}).
      vals : array-like
          A value per node used for coloring.
      ax : matplotlib.axes.Axes, optional
          Axis to plot on; if None, a new figure is created.
      cmap : str or Colormap, optional
          Colormap for the hexagons.

    Returns:
      ax : matplotlib.axes.Axes
          The axis with the hexagon tiling.
    """
    points = np.array(list(pos.values()))
    x, y = points.T
    return plot_hex_tiling_from_nodes(x, y, vals, ax=ax, cmap=cmap)


def plot_honeycomb_grid_fast(
    data: NDArray,
    ax: Any,
    triangle_size: float = 1.0,
    cmap: Union[Colormap, str] = 'viridis',
    rotation: float = 0.0
) -> None:
    """
    Fast version of plot_honeycomb_grid using PolyCollection.
    
    Plot a rotated honeycomb (triangular) grid with equal side lengths.
    Much faster than the original as it creates all triangles at once.

    Parameters
    ----------
    data : NDArray
        2D array of values. MUST be 2D with shape matching the desired grid dimensions.
        Common mistake: using flat (1D) arrays will cause errors or incorrect output.
    ax : matplotlib.axes.Axes
        The axes to plot on.
    triangle_size : float
        Side length of each equilateral triangle.
    cmap : Colormap or str
        Colormap for the triangles.
    rotation : float
        Angle in degrees to rotate the entire grid (counter-clockwise).
        
    Raises
    ------
    ValueError
        If data is not a 2D array.
        
    Note
    ----
    This function requires a 2D array matching the grid shape. If you have flat data,
    reshape it first: data.reshape(n_rows, n_cols).
    """
    # Input validation
    if data.ndim != 2:
        raise ValueError(
            f"Input data must be a 2D array, got {data.ndim}D array with shape {data.shape}. "
            f"If you have flat data, reshape it first: data.reshape(n_rows, n_cols)"
        )
    
    n_rows, n_cols = data.shape
    cmap_fun = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    h = (np.sqrt(3) / 2) * triangle_size

    # Rotation setup
    theta = np.deg2rad(rotation)
    c, s = np.cos(theta), np.sin(theta)

    # Pre-create triangle templates
    # Upward pointing triangle template
    up_tri_template = np.array([
        [0,       h/2],
        [-triangle_size/2, -h/2],
        [triangle_size/2, -h/2],
    ])
    
    # Downward pointing triangle template
    down_tri_template = np.array([
        [0,       -h/2],
        [-triangle_size/2, h/2],
        [triangle_size/2, h/2],
    ])

    # Rotate templates if needed
    if rotation != 0.0:
        rot_matrix = np.array([[c, -s], [s, c]])
        up_tri_template = up_tri_template @ rot_matrix.T
        down_tri_template = down_tri_template @ rot_matrix.T

    # Collect all polygons and colors
    polys = []
    colors = []
    
    for i in range(n_rows):
        for j in range(n_cols):
            # Base position
            x0 = j * triangle_size / 2
            y0 = i * h

            # Apply rotation to position
            if rotation != 0.0:
                x0_rot = x0 * c - y0 * s
                y0_rot = x0 * s + y0 * c
                x0, y0 = x0_rot, y0_rot

            # Pick orientation
            if (i + j) % 2 == 0:
                triangle = up_tri_template + np.array([x0, y0])
            else:
                triangle = down_tri_template + np.array([x0, y0])

            polys.append(triangle)
            colors.append(data[i, j])

    # Normalize colors
    norm = Normalize(vmin=np.min(data), vmax=np.max(data))
    facecolors = cmap_fun(norm(colors))

    # Create PolyCollection - much faster than individual patches
    coll = PolyCollection(
        polys,
        facecolors=facecolors,
        edgecolors='none',
        antialiased=False
    )
    ax.add_collection(coll)

    # Calculate bounds from all vertices
    all_verts = np.vstack(polys)
    min_x, max_x = np.min(all_verts[:, 0]), np.max(all_verts[:, 0])
    min_y, max_y = np.min(all_verts[:, 1]), np.max(all_verts[:, 1])

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_aspect('equal')


def plot_octagonal_square_tiling_from_nodes(x, y, vals, ax=None, cmap='viridis', 
                                           edgecolors='black', linewidths=0.1, padding=0.0):
    """
    Creates a triangular tiling for the octagonal-square (rhomb-octagonal) lattice.
    
    For each node, creates a triangle with one vertex at the rhomb center 
    and the other two vertices at the centers of the two octagons that the node belongs to.
    
    Parameters:
    -----------
    x, y : array-like
        x and y coordinates of the nodes.
    vals : array-like
        A value per node used for coloring the triangles.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on; if None, a new figure is created.
    cmap : str or Colormap, optional
        Colormap for the triangles.
    edgecolors : str or array-like, optional
        Colors for triangle edges. Use 'none' or 'face' to remove borders.
    linewidths : float, optional
        Width of triangle edges.
    padding : float, optional
        Padding fraction to add around the plot (default: 0.0 for no padding).

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axis with the triangle tiling.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    
    # Convert to numpy arrays for easier manipulation
    x = np.array(x)
    y = np.array(y)
    vals = np.array(vals)
    n_nodes = len(x)
    
    # Lattice parameters (should match the generator)
    rhomb_size = 0.4
    rhomb_spacing = 4.0 * rhomb_size  # 1.6
    
    # For each node, determine:
    # 1. Which rhomb it belongs to (rhomb center)
    # 2. Which two octagons it belongs to (octagon centers)
    
    triangles = []
    triangle_colors = []
    
    for i in range(n_nodes):
        node_x, node_y = x[i], y[i]
        
        # Find rhomb center by rounding node position to nearest rhomb grid point
        rhomb_col = round(node_x / rhomb_spacing)
        rhomb_row = round(node_y / rhomb_spacing)
        rhomb_center_x = rhomb_col * rhomb_spacing
        rhomb_center_y = rhomb_row * rhomb_spacing
        
        # Determine which node this is within the rhomb (0=bottom, 1=left, 2=right, 3=top)
        # based on its offset from the rhomb center
        offset_x = node_x - rhomb_center_x
        offset_y = node_y - rhomb_center_y
        
        # Tolerance for floating point comparison
        tol = 1e-6
        
        if abs(offset_y + rhomb_size) < tol and abs(offset_x) < tol:
            # Node 0: bottom
            # Belongs to octagons at: (rhomb_col-0.5, rhomb_row-0.5) and (rhomb_col+0.5, rhomb_row-0.5)
            oct1_x = (rhomb_col - 0.5) * rhomb_spacing
            oct1_y = (rhomb_row - 0.5) * rhomb_spacing
            oct2_x = (rhomb_col + 0.5) * rhomb_spacing
            oct2_y = (rhomb_row - 0.5) * rhomb_spacing
        elif abs(offset_x + rhomb_size) < tol and abs(offset_y) < tol:
            # Node 1: left
            # Belongs to octagons at: (rhomb_col-0.5, rhomb_row-0.5) and (rhomb_col-0.5, rhomb_row+0.5)
            oct1_x = (rhomb_col - 0.5) * rhomb_spacing
            oct1_y = (rhomb_row - 0.5) * rhomb_spacing
            oct2_x = (rhomb_col - 0.5) * rhomb_spacing
            oct2_y = (rhomb_row + 0.5) * rhomb_spacing
        elif abs(offset_x - rhomb_size) < tol and abs(offset_y) < tol:
            # Node 2: right
            # Belongs to octagons at: (rhomb_col+0.5, rhomb_row-0.5) and (rhomb_col+0.5, rhomb_row+0.5)
            oct1_x = (rhomb_col + 0.5) * rhomb_spacing
            oct1_y = (rhomb_row - 0.5) * rhomb_spacing
            oct2_x = (rhomb_col + 0.5) * rhomb_spacing
            oct2_y = (rhomb_row + 0.5) * rhomb_spacing
        elif abs(offset_y - rhomb_size) < tol and abs(offset_x) < tol:
            # Node 3: top
            # Belongs to octagons at: (rhomb_col-0.5, rhomb_row+0.5) and (rhomb_col+0.5, rhomb_row+0.5)
            oct1_x = (rhomb_col - 0.5) * rhomb_spacing
            oct1_y = (rhomb_row + 0.5) * rhomb_spacing
            oct2_x = (rhomb_col + 0.5) * rhomb_spacing
            oct2_y = (rhomb_row + 0.5) * rhomb_spacing
        else:
            # Fallback: shouldn't happen with correct rhomb-octagonal lattice
            print(f"Warning: Node {i} at ({node_x:.3f}, {node_y:.3f}) doesn't match expected rhomb positions")
            # Use rhomb center as octagon centers (degenerate triangle)
            oct1_x = oct2_x = rhomb_center_x
            oct1_y = oct2_y = rhomb_center_y
        
        # Create triangle with vertices: rhomb center, octagon 1 center, octagon 2 center
        triangle = np.array([
            [rhomb_center_x, rhomb_center_y],  # Rhomb center
            [oct1_x, oct1_y],                  # First octagon center
            [oct2_x, oct2_y]                   # Second octagon center
        ])
        
        triangles.append(triangle)
        triangle_colors.append(vals[i])
    
    # Normalize colors
    norm = Normalize(vmin=np.min(vals), vmax=np.max(vals))
    colormap = plt.get_cmap(cmap)
    facecolors = colormap(norm(triangle_colors))
    
    # Create PolyCollection
    coll = PolyCollection(
        triangles,
        facecolors=facecolors,
        edgecolors=edgecolors,
        linewidths=linewidths,
        antialiased=True
    )
    ax.add_collection(coll)
    
    # Set axis limits based on triangle extents
    all_verts = np.vstack(triangles)
    min_x, max_x = np.min(all_verts[:, 0]), np.max(all_verts[:, 0])
    min_y, max_y = np.min(all_verts[:, 1]), np.max(all_verts[:, 1])
    
    # Add padding if specified
    if padding > 0:
        padding_x = (max_x - min_x) * padding
        padding_y = (max_y - min_y) * padding
        ax.set_xlim(min_x - padding_x, max_x + padding_x)
        ax.set_ylim(min_y - padding_y, max_y + padding_y)
    else:
        ax.set_xlim(min_x, max_x)
        ax.set_ylim(min_y, max_y)
    ax.set_aspect('equal')
    
    return ax


def plot_octagonal_square_tiling_from_pos(pos, vals, ax=None, cmap='viridis', 
                                         edgecolors='black', linewidths=0.1):
    """
    Extracts x and y coordinates from a position dictionary and calls
    plot_octagonal_square_tiling_from_nodes to create the octagonal-square tiling.

    Parameters:
    -----------
    pos : dict
        Dictionary of node positions (e.g. {node: (x, y), ...}).
    vals : array-like
        A value per node used for coloring.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on; if None, a new figure is created.
    cmap : str or Colormap, optional
        Colormap for the triangles.
    edgecolors : str or array-like, optional
        Colors for triangle edges. Use 'none' or 'face' to remove borders.
    linewidths : float, optional
        Width of triangle edges.

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axis with the triangle tiling.
    """
    # Sort positions by node ID to ensure consistent ordering
    sorted_items = sorted(pos.items())
    points = np.array([point for _, point in sorted_items])
    x, y = points.T
    return plot_octagonal_square_tiling_from_nodes(x, y, vals, ax=ax, cmap=cmap, 
                                                  edgecolors=edgecolors, linewidths=linewidths)
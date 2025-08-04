#
#                           ,((((((((((((((((((((((((.
#                      @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                   @@@@@@                              @@@@@@
#                  @@@@@@                               *@@@@@/
#    /@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                                #    TRASH CODE
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#   *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@,
#    @@@@@@@@@@@@@@    @@@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@
#    &@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@&
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    .@@@@@@@@@@@@@
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    (@@@@@@@@@@@@@
#     @@@@@@@@@@@@@     @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    @@@@@@@@@@@@@@
#     @@@@@@@@@@@@@.    @@@@@@@@@@@@@@     @@@@@@@@@@@@@@    @@@@@@@@@@@@@@
#     @@@@@@@@@@@@@/    @@@@@@@@@@@@@@     @@@@@@@@@@@@@&    @@@@@@@@@@@@@@
#     *@@@@@@@@@@@@@    @@@@@@@@@@@@@@     @@@@@@@@@@@@@(    @@@@@@@@@@@@@.
#      @@@@@@@@@@@@@    &@@@@@@@@@@@@@     @@@@@@@@@@@@@.    @@@@@@@@@@@@@
#      @@@@@@@@@@@@@    *@@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@
#      &@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@#
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@
#       @@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#       *@@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@     @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@@    @@@@@@@@@@@@@     @@@@@@@@@@@@@    @@@@@@@@@@@@
#        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#         @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&
#          @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&
#
# def generate_maxpercdiff_colormap(
#     number_of_distinct_colors: int = 80, number_of_shades: int = 7
# ):
#     """
#     Generates a perceptually distinct colormap.
# 
#     Parameters:
#     -----------
#     number_of_distinct_colors : int, optional
#         The number of distinct colors in the colormap. Default is 80.
# 
#     Returns:
#     --------
#     ListedColormap
#         A ListedColormap object representing the generated colormap.
# 
#     Notes:
#     ------
#     This function generates a perceptually distinct colormap using a saw-tooth 
#         pattern in the HSV color space. The number of distinct colors can be ]
#         customized by providing the `number_of_distinct_colors` parameter.
# 
#     Reference:
#     -----------
#     - Based on the "saw-like" pattern technique for generating distinct colors.
#     - HSV colormap is used to create cyclic color variations.
# 
#     Example:
#     --------
#     >>> colormap = generate_colormap(100)
#     >>> import matplotlib.pyplot as plt
#     >>> import numpy as np
#     >>> data = np.random.rand(10, 10)
#     >>> plt.imshow(data, cmap=colormap)
#     >>> plt.colorbar()
#     >>> plt.show()
#     """
# 
#     no_distinct_colors_wmultshades = int(
#         np.ceil(number_of_distinct_colors / number_of_shades) * number_of_shades
#     )
# 
#     # Create an array with uniformly drawn floats taken from <0, 1) partition
#     linearly_distributed_nums = (
#         np.arange(no_distinct_colors_wmultshades)
#         / no_distinct_colors_wmultshades
#     )
# 
#     # We are going to reorganize monotonically growing numbers in such a way that there will be a single array with a saw-like pattern
#     #     but each sawtooth is slightly higher than the one before
#     # First divide linearly_distributed_nums into number_of_shades sub-arrays containing linearly distributed numbers
#     arr_by_shade_rows = linearly_distributed_nums.reshape(
#         number_of_shades, no_distinct_colors_wmultshades // number_of_shades
#     )
# 
#     # Transpose the above matrix (columns become rows) - as a result, each row contains a sawtooth with values slightly higher than the row above
#     arr_by_shade_columns = arr_by_shade_rows.T
# 
#     # Keep the number of sawtooths for later
#     number_of_partitions = arr_by_shade_columns.shape[0]
# 
#     # Flatten the above matrix - join each row into a single array
#     nums_distributed_like_rising_saw = arr_by_shade_columns.reshape(-1)
# 
#     # HSV colormap is cyclic (https://matplotlib.org/tutorials/colors/colormaps.html#cyclic), we'll use this property
#     initial_cm = hsv(nums_distributed_like_rising_saw)
# 
#     lower_partitions_half = number_of_partitions // 2
#     upper_partitions_half = number_of_partitions - lower_partitions_half
# 
#     # Modify the lower half in such a way that colors towards the beginning of the partition are darker
#     # First colors are affected more, colors closer to the middle are affected less
#     lower_half = lower_partitions_half * number_of_shades
#     for i in range(3):
#         initial_cm[0:lower_half, i] *= np.arange(0.2, 1, 0.8 / lower_half)
# 
#     # Modify the second half in such a way that colors towards the end of the partition are less intense and brighter
#     # Colors closer to the middle are affected less, colors closer to the end are affected more
#     for i in range(3):
#         for j in range(upper_partitions_half):
#             modifier = (
#                 np.ones(number_of_shades)
#                 - initial_cm[
#                     lower_half
#                     + j * number_of_shades : lower_half
#                     + (j + 1) * number_of_shades,
#                     i,
#                 ]
#             )
#             modifier = j * modifier / upper_partitions_half
#             initial_cm[
#                 lower_half
#                 + j * number_of_shades : lower_half
#                 + (j + 1) * number_of_shades,
#                 i,
#             ] += modifier
#
#     return ListedColormap(initial_cm)
#
# def plot_line_with_style(ax, x_coords, y_coords, color, kwargs_lines, cpec):
#     """
#     Helper function to plot a line with a specific color and optional style.
#     """
#     kwargs_lines["color"] = color
#     if color == cpec:
#         with rc_context({'path.sketch': (5, 15, 1)}):
#             ax.plot(x_coords, y_coords, **kwargs_lines)
#     else:
#         ax.plot(x_coords, y_coords, **kwargs_lines)
#
# def pos_evolving_grid_rw(n_steps: int, p: float, init: str = 'random'):
#     pos_state = np.zeros((n_steps, 3))  # Now also includes state in the last dimension]
#     if init == 'random':
#         pos_state[0, 2] = np.random.choice([-1, +1])  # Initial state
#     elif init == 'fixed':
#         pos_state[0, 2] = +1
#     global edge_signs 
#     edge_signs = {}
# 
#     # Function to update position and state
#     def update_position_state(current_pos_state, direction, p):
#         x, y, state = current_pos_state
#         key = (x, y, direction)
#         if key not in edge_signs:
#             edge_signs[key] = -1 if np.random.rand() < p else 1
#         sign = edge_signs[key]
# 
#         # Update position based on direction
#         if direction == 0:  # Move up
#             y += 1
#         elif direction == 1:  # Move down
#             y -= 1
#         elif direction == 2:  # Move left
#             x -= 1
#         elif direction == 3:  # Move right
#             x += 1
# 
#         # Flip state if necessary
#         if sign == -1:
#             state *= -1
# 
#         return np.array([x, y, state])
# 
#     # Generate the walk
#     for i in range(1, n_steps):
#         direction = np.random.randint(0, 4)  # Choose direction
#         pos_state[i] = update_position_state(pos_state[i-1], direction, p)
#     return pos_state
# 
# def average_evolving_rw(replica: int = 10**2, n_steps: int = 10**4, p=0.3, init: str = 'random'):
#     cumulative_walk = {}
#     for _ in range(replica):
#         walk = pos_evolving_grid_rw(n_steps, p,  init = 'fixed')
#         for x, y, value in walk:
#             key = (int(x), int(y))
#             if key in cumulative_walk:
#                 cumulative_walk[key] += value
#             else:
#                 cumulative_walk[key] = value

#     # Find extents for the 2D array
#     min_x = min(key[0] for key in cumulative_walk)
#     max_x = max(key[0] for key in cumulative_walk)
#     min_y = min(key[1] for key in cumulative_walk)
#     max_y = max(key[1] for key in cumulative_walk)

#     # Create and fill the array
#     array_shape = (int(max_y - min_y + 1), int(max_x - min_x + 1))
#     walk_array = np.zeros(array_shape)
#     for (x, y), value in cumulative_walk.items():
#         walk_array[y - min_y, x - min_x] = value / replica  # Average the states
#     return walk_array
#
# def plot_evolving_grid_rw(n_steps: int = 10**4, p: float = 0.103, ax: Axes = plt.Axes, col1: Any = 'red', col2: Any = 'blue', init: str = 'random'):
#     # Initialize a 3D array: (n_steps, 2 positions, 1 state)
#     colors = {1: col1, -1: col2}
#     pos_state = pos_evolving_grid_rw(n_steps, p, init = init)
#     # Plot the steps with appropriate color based on the state
#     for i in range(1, n_steps):
#         ax.plot(pos_state[i-1:i+1, 0], pos_state[i-1:i+1, 1], color=colors[pos_state[i, 2]], linewidth=2)
#
# def perform_random_walks(G, steps, N):
#     """
#     Perform N random walks on graph G, each with a given number of steps,
#     and update node states based on edge weights.
#     """
#     node_states = {node: 0 for node in G.nodes()}  # Initial node states
#     for _ in range(N):
#         current_state = np.random.choice([-1, 1])
#         current_node = list(G.nodes())[np.random.randint(len(G))]  # Start at the origin for each walk
#         node_states[current_node] += current_state
#         for _ in range(steps):
#             neighbors = list(G.neighbors(current_node))
#             next_node = neighbors[np.random.randint(len(neighbors))]
#             edge_weight = G.edges[current_node, next_node]['weight']
#             # Update the state of the next node based on edge weight
#             current_state = current_state*edge_weight
#             node_states[next_node] += current_state
#             current_node = next_node
#     return node_states
#
# def visualize_node_states(node_states, n, m, ax: Axes = plt.Axes):
#     """
#     Visualize the final states of nodes as a heatmap.
#     """
#     # Convert node states to a 2D array
#     state_array = np.zeros((n, m))
#     for (x, y), state in node_states.items():
#         state_array[x, y] = state
#
#     ax.imshow(state_array)
#     ax.axis('off')  # Hide the axes
#
# def plot_honeycomb_grid(data: NDArray, fig: Any, ax: Any, triangle_size: float = 1.0, cmap: Union[Colormap, str] = credcblu) -> None:
#     """
#     Plot a triangular grid representing the honeycomb structure based on the provided 2D array data.

#     Parameters
#     ----------
#     data : np.ndarray
#         A 2D numpy array representing values at each point in the honeycomb structure.
#     fig : matplotlib.figure.Figure
#         A matplotlib figure object to plot on.
#     ax : matplotlib.axes._axes.Axes
#         A matplotlib axes object to plot on.
#     triangle_size : float, optional
#         The side length of the equilateral triangles in the grid. Default is 1.0.
#
#     Returns
#     -------
#     None
#     """
#     import matplotlib.patches as patches
#     # Determine number of rows and columns from data
#     n_rows, n_cols = data.shape
#     if isinstance(cmap, str):
#         cmap_fun = plt.get_cmap(cmap)
#     else:
#         cmap_fun = cmap
#     # Height of an equilateral triangle
#     triangle_height = (np.sqrt(3) / 2) * triangle_size
#
#     # Loop over each row and column to create triangular patches
#     for row in range(n_rows):
#         for col in range(n_cols):
#             # Determine the x and y positions for the center of each triangle
#             x_base = col * triangle_size / 2
#             y_base = row * triangle_height
#
#             # Define the vertices of the triangles (upward and downward)
#             if (row + col) % 2 == 0:
#                 # Upward-pointing triangle
#                 vertices = [
#                     (x_base, y_base + (triangle_height / 2)),
#                     (x_base - (triangle_size / 2), y_base - (triangle_height / 2)),
#                     (x_base + (triangle_size / 2), y_base - (triangle_height / 2)),
#                 ]
#             else:
#                 # Downward-pointing triangle
#                 vertices = [
#                     (x_base, y_base - (triangle_height / 2)),
#                     (x_base - (triangle_size / 2), y_base + (triangle_height / 2)),
#                     (x_base + (triangle_size / 2), y_base + (triangle_height / 2)),
#                 ]
#
#             # Get color based on data
#             color = cmap_fun(data[row, col])

#             # Create and add the triangle patch to the axis
#             triangle = patches.Polygon(vertices, facecolor=color)
#             ax.add_patch(triangle)
#
#     # Set the limits and aspect ratio to ensure proper display
#     ax.set_xlim(-triangle_size / 2, (n_cols / 2 + 0.5) * triangle_size)
#     ax.set_ylim(-triangle_height / 2, (n_rows + 0.5) * triangle_height)
#     ax.set_aspect('equal')
#
#     # Hide the axes for better visualization
#     ax.axis('off')
# def plot_tiangular_grid(data, fig, ax, triangle_size: float = 1.0) -> None:
#     """
#     Plot a planar triangular grid by placing a hexagon at each node.
#     Each hexagon is colored based on the corresponding value in the provided 2D array.
    
#     The centers of the hexagons are arranged on a triangular (offset) grid. For a pointy-topped
#     hexagon tiling the horizontal distance between adjacent centers is defined by:
#         distance = √3 * R,  with  R = triangle_size / √3.
#     This yields centers computed as:
#         x_center = col * triangle_size + (row % 2) * (triangle_size / 2)
#         y_center = row * ((√3)/2 * triangle_size)
    
#     Parameters
#     ----------
#     data : np.ndarray
#         A 2D numpy array representing values at each node of the triangular grid.
#     fig : matplotlib.figure.Figure
#         A matplotlib figure object to plot on.
#     ax : matplotlib.axes._axes.Axes
#         A matplotlib axes object to plot on.
#     triangle_size : float, optional
#         The grid spacing parameter. It is used both to determine the positions of nodes and
#         to set the size of each hexagon (with the hexagon circumradius being triangle_size/√3).
#         Default is 1.0.
    
#     Returns
#     -------
#     None
#     """
#     import numpy as np
#     import matplotlib.patches as patches

#     # Determine number of rows and columns in the data array.
#     n_rows, n_cols = data.shape

#     # For a pointy-topped hexagon, the distance between centers horizontally equals √3 * R.
#     # We choose R so that √3 * R = triangle_size.
#     hex_radius = triangle_size / np.sqrt(3)

#     # Loop over the grid and add a hexagon at each node.
#     for row in range(n_rows):
#         for col in range(n_cols):
#             # Calculate the center of the hexagon.
#             x_center = col * triangle_size + (row % 2) * (triangle_size / 2)
#             y_center = row * ((np.sqrt(3) / 2) * triangle_size)
            
#             # Get color based on the data value (assuming credcblu is defined elsewhere).
#             color = credcblu(data[row, col])
            
#             # Create a pointy-topped hexagon (using orientation=π/2 so that one vertex is at the top).
#             hexagon = patches.RegularPolygon(
#                 (x_center, y_center),
#                 numVertices=6,
#                 radius=hex_radius,
#                 orientation=np.pi/2,
#                 facecolor=color
#             )
#             ax.add_patch(hexagon)

#     # Compute plot limits to ensure all hexagons are visible.
#     # For x: the leftmost hexagon (from an even row) is at 0 - hex_radius;
#     # for rows with an offset (odd rows) the rightmost center is at (n_cols-1)*triangle_size + triangle_size/2.
#     x_min = -hex_radius
#     if n_rows > 1:
#         x_max = (n_cols - 1) * triangle_size + (triangle_size / 2) + hex_radius
#     else:
#         x_max = (n_cols - 1) * triangle_size + hex_radius

#     # For y: the lowest center is at 0, so subtract hex_radius;
#     # the highest center is at (n_rows-1) * ((√3)/2 * triangle_size) plus hex_radius.
#     y_min = -hex_radius
#     y_max = (n_rows - 1) * ((np.sqrt(3) / 2) * triangle_size) + hex_radius

#     ax.set_xlim(x_min, x_max)
#     ax.set_ylim(y_min, y_max)
#     ax.set_aspect('equal')
#     ax.axis('off')

# def plot_hex_tiling_from_nodes(x, y, vals, ax=None, cmap='viridis'):
#     if ax is None:
#         fig, ax = plt.subplots(figsize=(8, 8))
#     x_min, x_max = np.min(x), np.max(x)
#     y_min, y_max = np.min(y), np.max(y)
#     unique_x = np.unique(x); unique_y = np.unique(y)
#     cols, rows = len(unique_x), len(unique_y)
#     Rw = (x_max - x_min) / ((cols + 0.5)*np.sqrt(3))
#     Rh = (y_max - y_min) / (1.5*(rows-1) + 2)
#     R = min(Rw, Rh)
#     # precompute single hexagon template (centered at 0,0)
#     angles = np.linspace(0, 2*np.pi, 7)[:-1] + np.pi/6
#     base_hex = np.column_stack([np.cos(angles)*R, np.sin(angles)*R])
#     polys = []
#     for r in range(rows):
#         for c in range(cols):
#             cx = x_min + ((np.sqrt(3)/2) + c*np.sqrt(3) + (r%2)* (np.sqrt(3)/2))*R
#             cy = y_min + R + r*(1.5*R)
#             polys.append(base_hex + [cx, cy])
#     # sort lexicographically if you need to match vals order
#     hex_centers = np.array([p.mean(axis=0) for p in polys])
#     idx = np.lexsort((hex_centers[:,1], hex_centers[:,0]))
#     polys = [polys[i] for i in idx]
#     vals_sorted = np.array(vals)[ np.lexsort((y, x)) ][idx]
#     norm = Normalize(np.min(vals), np.max(vals))
#     cmap = plt.get_cmap(cmap)
#     facecolors = cmap(norm(vals_sorted))
#     coll = PolyCollection(
#         polys,
#         facecolors=facecolors,
#         edgecolors='none',
#         antialiased=False
#     )
#     ax.add_collection(coll)
#     ax.autoscale()
#     ax.set_aspect('equal')



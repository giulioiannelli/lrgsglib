import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

# Import the specific module directly
from lrgsglib.plotlib.tilings import plot_tri_tiling_from_nodes, plot_honeycomb_grid

# Create a simple test grid
n_rows, n_cols = 6, 8
data = np.random.rand(n_rows, n_cols)

# Create coordinate arrays for the new function
x_coords = []
y_coords = []
vals = []

triangle_size = 1.0
h = (np.sqrt(3) / 2) * triangle_size

for i in range(n_rows):
    for j in range(n_cols):
        x_coords.append(j * triangle_size / 2)
        y_coords.append(i * h)
        vals.append(data[i, j])

# Test the new optimized function
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

# Original slow function
plot_honeycomb_grid(data, fig, ax1, triangle_size=triangle_size)
ax1.set_title('Original plot_honeycomb_grid')

# New fast function
plot_tri_tiling_from_nodes(x_coords, y_coords, vals, ax=ax2, triangle_size=triangle_size)
ax2.set_title('New plot_tri_tiling_from_nodes')

plt.tight_layout()
plt.savefig('triangle_tiling_comparison.png', dpi=150, bbox_inches='tight')
plt.show()

print("Test completed successfully!")

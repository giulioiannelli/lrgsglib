Visualization
=============

This guide covers the visualization capabilities of lrgsglib for plotting
signed graphs, lattice configurations, spectral data, and 3D structures.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

lrgsglib provides visualization tools through the ``plotlib`` module, which
includes:

- **Lattice plots**: 2D and 3D lattice visualizations with spin configurations
- **Color utilities**: Custom colormaps and color schemes for signed graphs
- **3D volume plots**: Density visualizations using Plotly
- **Spectral plots**: Eigenvalue distributions and entropy profiles
- **Axis utilities**: Formatting helpers for publication-quality figures

Basic Graph Visualization
-------------------------

Using NetworkX
~~~~~~~~~~~~~~

For basic graph drawing, use NetworkX's built-in functions with lrgsglib graphs:

.. code-block:: python

   import matplotlib.pyplot as plt
   import networkx as nx
   from lrgsglib.graphs.ErdosRenyi import ErdosRenyi

   # Create signed graph
   er = ErdosRenyi(nnodes=50, prob=0.15, pflip=0.3, seed=42)
   er.flip_random_fract_edges()

   # Get the NetworkX graph
   G = er.gr['G']

   # Color edges by sign
   edge_colors = ['blue' if G[u][v].get('weight', 1) > 0 else 'red'
                  for u, v in G.edges()]

   # Draw with spring layout
   plt.figure(figsize=(10, 8))
   pos = nx.spring_layout(G, seed=42)
   nx.draw(G, pos,
           node_size=100,
           node_color='lightgray',
           edge_color=edge_colors,
           width=1.5,
           with_labels=False)
   plt.title(f'Signed Erdos-Renyi (N={er.N}, pflip={er.pflip})')
   plt.savefig('signed_graph.png', dpi=150)
   plt.show()

Lattice Visualizations
----------------------

2D Lattice with Positions
~~~~~~~~~~~~~~~~~~~~~~~~~

Lattice2D objects have built-in position information:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from lrgsglib.graphs.Lattice2D import Lattice2D

   # Create lattice
   lattice = Lattice2D(side1=16, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   # Get coordinate graph (H has positions)
   H = lattice.gr['H']

   # Extract positions from node labels (they are (x, y) tuples)
   pos = {node: node for node in H.nodes()}

   # Color edges by sign
   edge_colors = ['blue' if H[u][v].get('weight', 1) > 0 else 'red'
                  for u, v in H.edges()]

   # Plot
   fig, ax = plt.subplots(figsize=(10, 10))
   nx.draw(H, pos,
           ax=ax,
           node_size=50,
           node_color='black',
           edge_color=edge_colors,
           width=1.0)
   ax.set_aspect('equal')
   ax.set_title(f'2D Lattice ({lattice.side1}x{lattice.side1})')
   plt.savefig('lattice2d.png', dpi=150)
   plt.show()

Spin Configuration Heatmaps
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Visualize spin states as a 2D heatmap:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from lrgsglib.graphs.Lattice2D import Lattice2D
   from lrgsglib.statsys.IsingDynamics import IsingDynamics

   # Create lattice and run Ising dynamics
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   ising = IsingDynamics(sg=lattice, T=2.0, steps=2000, runlang='py')
   ising.init_ising_dynamics()
   ising.run(verbose=False)

   # Reshape spin configuration to 2D
   spins_2d = ising.s.reshape(lattice.side1, lattice.side1)

   # Plot heatmap
   fig, ax = plt.subplots(figsize=(8, 8))
   im = ax.imshow(spins_2d, cmap='RdBu', vmin=-1, vmax=1)
   ax.set_title(f'Ising Configuration at T={ising.T}')
   plt.colorbar(im, ax=ax, label='Spin', shrink=0.8)
   plt.savefig('ising_config.png', dpi=150)
   plt.show()

Eigenvector Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~

Visualize eigenvector components on the lattice:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from lrgsglib.graphs.Lattice2D import Lattice2D

   # Create lattice and compute eigenvectors
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()
   lattice.compute_laplacian_spectrum_weigV()

   # Reshape Fiedler vector (second eigenvector) to 2D
   fiedler = lattice.eigV[1]
   fiedler_2d = fiedler.reshape(lattice.side1, lattice.side1)

   # Plot
   fig, axes = plt.subplots(1, 2, figsize=(14, 6))

   # Eigenvector heatmap
   im1 = axes[0].imshow(fiedler_2d, cmap='RdBu')
   axes[0].set_title(f'Fiedler Vector (eigenvalue={lattice.eigv[1]:.4f})')
   plt.colorbar(im1, ax=axes[0], shrink=0.8)

   # Ground state
   ground = lattice.eigV[0].reshape(lattice.side1, lattice.side1)
   im2 = axes[1].imshow(ground, cmap='viridis')
   axes[1].set_title(f'Ground State (eigenvalue={lattice.eigv[0]:.4f})')
   plt.colorbar(im2, ax=axes[1], shrink=0.8)

   plt.tight_layout()
   plt.savefig('eigenvectors.png', dpi=150)
   plt.show()

Spectral Visualization
----------------------

Eigenvalue Distribution
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from lrgsglib.graphs.Lattice2D import Lattice2D

   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.3, seed=42)
   lattice.flip_random_fract_edges()
   lattice.compute_laplacian_spectrum()

   fig, axes = plt.subplots(1, 3, figsize=(15, 4))

   # Histogram
   axes[0].hist(lattice.eigv, bins=50, edgecolor='black', alpha=0.7)
   axes[0].axvline(x=0, color='r', linestyle='--', label='Zero')
   axes[0].set_xlabel('Eigenvalue')
   axes[0].set_ylabel('Count')
   axes[0].set_title('Eigenvalue Distribution')
   axes[0].legend()

   # Sorted spectrum (staircase)
   axes[1].plot(np.sort(lattice.eigv), 'b-', linewidth=0.5)
   axes[1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
   axes[1].set_xlabel('Index')
   axes[1].set_ylabel('Eigenvalue')
   axes[1].set_title('Spectral Staircase')

   # Density of states
   density, bins = np.histogram(lattice.eigv, bins=100, density=True)
   bin_centers = (bins[:-1] + bins[1:]) / 2
   axes[2].fill_between(bin_centers, density, alpha=0.5)
   axes[2].plot(bin_centers, density, 'b-', linewidth=1)
   axes[2].set_xlabel('Eigenvalue')
   axes[2].set_ylabel('Density')
   axes[2].set_title('Density of States')

   plt.tight_layout()
   plt.savefig('spectrum.png', dpi=150)
   plt.show()

Entropy Profiles
~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from lrgsglib.graphs.MultiplicativeCascade import MultiplicativeCascadeGraph
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create graph and compute spectrum
   mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, fraction=0.4, seed=42)
   L, eigenvalues = get_graph_lspectrum(mc.gr['G'], library='numpy')

   # Compute entropy observables
   entropy, heat, var, tau = compute_entropy_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=mc.N,
       steps=300,
   )

   # Publication-quality figure
   fig, axes = plt.subplots(1, 3, figsize=(15, 4))

   # Entropy
   axes[0].semilogx(tau, entropy, 'b-', linewidth=2)
   axes[0].set_xlabel(r'$\tau$', fontsize=12)
   axes[0].set_ylabel(r'$1 - S(\tau)/\log N$', fontsize=12)
   axes[0].set_title('Von Neumann Entropy')
   axes[0].grid(True, alpha=0.3)

   # Specific heat
   ds = 2.0 * np.mean(heat[-30:])
   axes[1].semilogx(tau, heat, 'r-', linewidth=2)
   axes[1].axhline(y=ds/2, color='k', linestyle='--',
                   label=f'$d_s/2 = {ds/2:.2f}$')
   axes[1].set_xlabel(r'$\tau$', fontsize=12)
   axes[1].set_ylabel(r'$C(\tau)$', fontsize=12)
   axes[1].set_title('Specific Heat')
   axes[1].legend()
   axes[1].grid(True, alpha=0.3)

   # Variance
   axes[2].loglog(tau, var, 'g-', linewidth=2)
   axes[2].set_xlabel(r'$\tau$', fontsize=12)
   axes[2].set_ylabel(r'$\sigma^2(\tau)$', fontsize=12)
   axes[2].set_title('Eigenvalue Variance')
   axes[2].grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('entropy_profile.png', dpi=300, bbox_inches='tight')
   plt.show()

3D Visualization
----------------

3D Lattice Structure
~~~~~~~~~~~~~~~~~~~~

For 3D lattices, use Matplotlib's 3D plotting:

.. code-block:: python

   import matplotlib.pyplot as plt
   from mpl_toolkits.mplot3d import Axes3D
   import numpy as np
   from lrgsglib.graphs.Lattice3D import Lattice3D

   # Create 3D lattice
   lattice3d = Lattice3D(side1=6, geo='cub', pflip=0.2, seed=42)
   lattice3d.flip_random_fract_edges()

   # Extract node positions and edge signs
   H = lattice3d.gr['H']
   nodes = list(H.nodes())
   coords = np.array(nodes)

   fig = plt.figure(figsize=(10, 8))
   ax = fig.add_subplot(111, projection='3d')

   # Plot nodes
   ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
              c='black', s=50, alpha=0.8)

   # Plot edges colored by sign
   for u, v in H.edges():
       weight = H[u][v].get('weight', 1)
       color = 'blue' if weight > 0 else 'red'
       alpha = 0.5 if weight > 0 else 0.8
       ax.plot([u[0], v[0]], [u[1], v[1]], [u[2], v[2]],
               color=color, alpha=alpha, linewidth=1)

   ax.set_xlabel('X')
   ax.set_ylabel('Y')
   ax.set_zlabel('Z')
   ax.set_title(f'3D Cubic Lattice (pflip={lattice3d.pflip})')
   plt.savefig('lattice3d.png', dpi=150)
   plt.show()

Plotly 3D Volume Visualization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For advanced 3D visualization, use the plotlib.plot3d module:

.. code-block:: python

   import numpy as np
   from lrgsglib.plotlib.plot3d import create_volume_visualization

   # Create sample 3D density data (e.g., from eigenvector analysis)
   # This example creates synthetic data for demonstration
   np.random.seed(42)
   n_points = 500

   # Class 1: cluster centered at (2, 2, 2)
   x1 = np.random.randn(n_points) + 2
   y1 = np.random.randn(n_points) + 2
   z1 = np.random.randn(n_points) + 2

   # Class 2: cluster centered at (-2, -2, -2)
   x2 = np.random.randn(n_points) - 2
   y2 = np.random.randn(n_points) - 2
   z2 = np.random.randn(n_points) - 2

   # Create coordinate grids
   grid_resolution = 30
   x_range = np.linspace(-5, 5, grid_resolution)
   y_range = np.linspace(-5, 5, grid_resolution)
   z_range = np.linspace(-5, 5, grid_resolution)
   X, Y, Z = np.meshgrid(x_range, y_range, z_range)

   # Prepare data sets
   data_sets = [
       (x1, y1, z1, 'red', 'Cluster 1'),
       (x2, y2, z2, 'blue', 'Cluster 2'),
   ]

   # Create visualization
   fig = create_volume_visualization(
       data_sets, X, Y, Z,
       opacity=0.4,
       surface_count=2,
       sample_size=200,
   )

   fig.update_layout(
       title='3D Density Visualization',
       scene=dict(
           xaxis_title='X',
           yaxis_title='Y',
           zaxis_title='Z',
       )
   )

   fig.write_html('volume_3d.html')
   fig.show()

Custom Colormaps
----------------

lrgsglib provides custom colormaps for signed graph visualization:

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np
   from lrgsglib.plotlib.colormaps import restr_twilight
   from lrgsglib.plotlib.color import set_alpha_torgb

   # Get color palette for multiple datasets
   n_datasets = 5
   colors = restr_twilight(np.linspace(0, 1, n_datasets))

   # Plot with custom colors
   fig, ax = plt.subplots(figsize=(10, 6))
   x = np.linspace(0, 10, 100)

   for i, c in enumerate(colors):
       y = np.sin(x + i * 0.5) + i * 0.5
       # Add transparency
       c_alpha = set_alpha_torgb(c, 0.8)
       ax.plot(x, y, color=c, linewidth=2, label=f'Series {i+1}')
       ax.fill_between(x, y - 0.2, y + 0.2, color=c, alpha=0.2)

   ax.legend()
   ax.set_xlabel('x')
   ax.set_ylabel('y')
   ax.set_title('Custom Colormap Example')
   plt.savefig('custom_colors.png', dpi=150)
   plt.show()

Dynamics Animation
------------------

For time-evolving visualizations, create animations:

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib.animation as animation
   import numpy as np
   from lrgsglib.graphs.Lattice2D import Lattice2D
   from lrgsglib.statsys.IsingDynamics import IsingDynamics

   # Create lattice and run dynamics with snapshots
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   ising = IsingDynamics(
       sg=lattice,
       T=2.27,  # Critical temperature
       steps=200,
       runlang='py',
   )
   ising.init_ising_dynamics()

   # Collect snapshots
   snapshots = [ising.s.copy().reshape(32, 32)]
   for _ in range(199):
       ising.metropolis(np.random.randint(lattice.N))
       if _ % 10 == 0:  # Save every 10 steps
           snapshots.append(ising.s.copy().reshape(32, 32))

   # Create animation
   fig, ax = plt.subplots(figsize=(8, 8))
   im = ax.imshow(snapshots[0], cmap='RdBu', vmin=-1, vmax=1)
   ax.set_title('Ising Dynamics at Critical Temperature')

   def update(frame):
       im.set_array(snapshots[frame])
       ax.set_title(f'Step {frame * 10}')
       return [im]

   ani = animation.FuncAnimation(fig, update, frames=len(snapshots),
                                  interval=100, blit=True)
   ani.save('ising_animation.gif', writer='pillow', fps=10)
   plt.show()

High-Level Animation API
^^^^^^^^^^^^^^^^^^^^^^^^

For convenience, lrgsglib provides a high-level animation API. The simplest
entry points are bound directly on the lattice (engine-agnostic) --
``lat.animate_states(states)`` and ``lat.animate_largest_cluster(states)``. The
underlying builders live in the structural, graph-type module
:mod:`lrgsglib.graphs._shared.animation` and the generic codec primitives in
:mod:`lrgsglib.plotlib.animation`:

.. code-block:: python

   import matplotlib.pyplot as plt
   from lrgsglib.graphs import Lattice2D
   from lrgsglib.graphs._shared.animation import make_lattice2d_animation
   from lrgsglib.plotlib.animation import save_animation

   # Create lattice
   lattice = Lattice2D(side1=32, geo='sqr', seed=42)

   # Generate example frames (e.g., from simulation)
   import numpy as np
   frames = [np.random.randn(lattice.N) for _ in range(50)]

   # Create animation with one call
   fig, ax = plt.subplots(figsize=(8, 8))
   result = make_lattice2d_animation(
       lattice=lattice,
       fig=fig,
       ax=ax,
       frames=frames,
       interval_ms=100,
       cmap='viridis',
       add_colorbar=True,
       autoscale=True,  # Adjust color limits per frame
   )

   # Save to file
   save_animation(result.animation, 'lattice_anim.gif', fps=10, dpi=150)
   plt.show()

The ``make_lattice2d_animation`` function returns a ``LatticeAnimationResult``
containing the animation, image, and optional colorbar objects.

**Parameters:**

- ``lattice``: A Lattice2D instance (needs ``syshape`` and ``N`` attributes)
- ``frames``: Sequence of 1D vectors (length N) or 2D arrays (shape syshape)
- ``interval_ms``: Delay between frames in milliseconds
- ``cmap``: Matplotlib colormap name
- ``autoscale``: If True, update color limits per frame
- ``vmin``, ``vmax``: Fixed color limits when autoscale=False

Publication-Quality Figures
---------------------------

Tips for creating publication-ready figures:

.. code-block:: python

   import matplotlib.pyplot as plt
   import matplotlib as mpl

   # Set publication style
   plt.rcParams.update({
       'font.size': 12,
       'font.family': 'serif',
       'axes.labelsize': 14,
       'axes.titlesize': 14,
       'xtick.labelsize': 12,
       'ytick.labelsize': 12,
       'legend.fontsize': 11,
       'figure.dpi': 150,
       'savefig.dpi': 300,
       'savefig.bbox': 'tight',
   })

   # Example figure with LaTeX labels
   fig, ax = plt.subplots(figsize=(6, 4))

   x = np.linspace(0, 10, 100)
   ax.plot(x, np.exp(-x/3), 'b-', linewidth=2, label=r'$e^{-\tau/3}$')
   ax.set_xlabel(r'$\tau$')
   ax.set_ylabel(r'$\rho(\tau)$')
   ax.set_title('Diffusion Decay')
   ax.legend()
   ax.grid(True, alpha=0.3)

   plt.savefig('publication_figure.pdf', format='pdf')
   plt.show()

Best Practices
--------------

1. **Use appropriate colormaps**: ``RdBu`` for bipolar data (spins), ``viridis``
   for continuous data

2. **Set figure size early**: Define ``figsize`` based on publication requirements

3. **Use vector formats**: Save as PDF or SVG for publications

4. **Include colorbars**: Always show the scale for heatmaps

5. **Set seeds for layouts**: Use ``seed`` parameter in ``spring_layout`` for
   reproducible graph visualizations

6. **Adjust DPI**: Use 150+ DPI for screen, 300+ for print

API Reference
-------------

.. seealso::

   - :mod:`lrgsglib.plotlib.lattices` - Lattice visualization utilities
   - :mod:`lrgsglib.plotlib.plot3d` - 3D volume visualization
   - :mod:`lrgsglib.plotlib.colormaps` - Custom colormaps
   - :mod:`lrgsglib.plotlib.color` - Color manipulation utilities

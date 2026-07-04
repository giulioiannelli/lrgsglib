Quick Start
===========

This guide provides a quick introduction to the basic functionality of lrgsglib. For more detailed information, see the :doc:`user_guide/index`.

Creating Your First Signed Graph
---------------------------------

Let's start by creating a simple signed Erdős-Rényi random graph:

.. code-block:: python

   import numpy as np
   from lrgsglib.graphs import ErdosRenyi

   # Set random seed for reproducibility
   np.random.seed(42)

   # Create a signed Erdős-Rényi graph
   N = 100  # Number of nodes
   p = 0.1  # Edge probability
   q = 0.3  # Fraction of negative edges

   # Any SignedGraph subclass auto-flips a ``pflip`` fraction of edges at
   # construction; ``flip_random_fract_edges()`` (re)applies the disorder.
   er = ErdosRenyi(n=N, p=p, pflip=q)
   er.flip_random_fract_edges()
   G = er.gr['G']  # underlying networkx graph

   print(f"Created graph with {G.number_of_nodes()} nodes")
   print(f"Number of edges: {G.number_of_edges()}")

   # Count positive and negative edges (sign is stored under the 'weight' key)
   positive_edges = sum(1 for u, v in G.edges() if G[u][v]['weight'] > 0)
   negative_edges = sum(1 for u, v in G.edges() if G[u][v]['weight'] < 0)

   print(f"Positive edges: {positive_edges}")
   print(f"Negative edges: {negative_edges}")

Working with 2D Lattices
-------------------------

Create and manipulate 2D lattices with signed edges:

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   import numpy as np

   # Create a square lattice with some negative edges
   np.random.seed(12345)

   lattice = Lattice2D(
       side1=20,           # 20x20 lattice
       geo='sqr',          # Square geometry
       pflip=0.2,          # 20% of edges will be negative
       with_positions=True # Include node positions
   )

   # Randomly flip edges to negative
   lattice.flip_random_fract_edges()

   print(f"Lattice size: {lattice.side1}x{lattice.side2}")
   print(f"Total nodes: {lattice.N}")
   print(f"Total edges: {lattice.gr['G'].number_of_edges()}")

Computing Spectral Properties
------------------------------

Analyze the Laplacian spectrum of signed graphs:

.. code-block:: python

   from lrgsglib.utils.lrg import spectral
   import numpy as np
   from lrgsglib.graphs import Lattice2D

   # Create a small lattice
   np.random.seed(42)
   lattice = Lattice2D(side1=10, geo='sqr', pflip=0.3)
   lattice.flip_random_fract_edges()

   # Compute Laplacian spectrum
   L, eigenvalues = spectral.get_graph_lspectrum(
       lattice.gr['G'],
       library='numpy'
   )

   print(f"Spectrum shape: {eigenvalues.shape}")
   print(f"Smallest eigenvalue: {eigenvalues.min():.4f}")
   print(f"Largest eigenvalue: {eigenvalues.max():.4f}")

   # Compute Laplacian properties
   props = spectral.compute_laplacian_properties(lattice.gr['G'])
   print(f"Laplacian properties computed")

Running Contact Process Simulations
------------------------------------

Simulate the contact process on a signed lattice:

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   from lrgsglib.statsys import ContactProcess
   import numpy as np

   # Create a lattice
   np.random.seed(42)
   lattice = Lattice2D(side1=50, geo='sqr', pflip=0.2)
   lattice.flip_random_fract_edges()

   # Initialize contact process
   cp = ContactProcess(
       sg=lattice,
       gamma=0.5,              # Recovery rate
       activation='relu',       # Activation function
       state_type='binary',     # Binary states (0 or 1)
       runlang='C1c'           # Use C implementation
   )

   # Set initial condition (5% infected)
   np.random.seed(99)
   initial_state = np.random.choice(
       [0, 1],
       size=lattice.N,
       p=[0.95, 0.05]
   )
   cp.init_contact_dynamics(custom=initial_state)

   # Run simulation
   steps = 10000
   cp.run(steps=steps, verbose=True, tqdm_on=True)

   # Check final state
   final_density = np.mean(cp.s)
   print(f"Final infected density: {final_density:.4f}")

Visualizing Results
-------------------

Plot signed graphs and lattices:

.. code-block:: python

   from lrgsglib.plotlib import lattices
   import matplotlib.pyplot as plt

   # Draw a schematic signed square lattice with a random negative-edge pattern
   fig, ax = plt.subplots(figsize=(8, 8))
   lattices.scheme_Lattice2DSquared(
       ax,
       side1=10,
       side2=10,
       mode='rand',
       pflip=0.3,
   )
   plt.title("Signed 2D Lattice")
   plt.show()

Working with 3D Lattices
-------------------------

Create and analyze 3D cubic lattices:

.. code-block:: python

   from lrgsglib.graphs import Lattice3D
   import numpy as np

   # Create a 3D cubic lattice
   np.random.seed(42)

   lattice_3d = Lattice3D(
       dim=(10, 10, 10),   # 10x10x10 lattice
       geo='sc',           # Simple-cubic geometry
       pflip=0.2           # 20% of edges will be negative
   )

   # Flip edges
   lattice_3d.flip_random_fract_edges()

   sx, sy, sz = lattice_3d.dim
   print(f"3D Lattice: {sx}x{sy}x{sz}")
   print(f"Total nodes: {lattice_3d.N}")
   print(f"Total edges: {lattice_3d.gr['G'].number_of_edges()}")

Using Different Graph Types
----------------------------

lrgsglib supports various graph topologies:

.. code-block:: python

   from lrgsglib.graphs import (
       ErdosRenyi,
       Lattice2D,
       Lattice3D,
   )
   import numpy as np

   np.random.seed(42)

   # Erdős-Rényi random graph
   er = ErdosRenyi(n=100, p=0.1, pflip=0.3)
   er.flip_random_fract_edges()
   er_graph = er.gr['G']

   # 2D triangular lattice
   tri_lattice = Lattice2D(side1=15, geo='tri', pflip=0.2)
   tri_lattice.flip_random_fract_edges()

   # 2D hexagonal lattice
   hex_lattice = Lattice2D(side1=15, geo='hex', pflip=0.2)
   hex_lattice.flip_random_fract_edges()

   # 3D simple-cubic lattice
   cubic_3d = Lattice3D(dim=(8, 8, 8), geo='sc', pflip=0.2)
   cubic_3d.flip_random_fract_edges()

   print(f"ER graph: {er_graph.number_of_nodes()} nodes")
   print(f"Triangular lattice: {tri_lattice.N} nodes")
   print(f"Hexagonal lattice: {hex_lattice.N} nodes")
   print(f"3D cubic: {cubic_3d.N} nodes")

Choosing Simulation Backends
-----------------------------

lrgsglib provides multiple backends for simulations:

.. code-block:: python

   from lrgsglib.statsys import ContactProcess
   from lrgsglib.graphs import Lattice2D
   import numpy as np

   np.random.seed(42)
   lattice = Lattice2D(side1=30, geo='sqr', pflip=0.2)
   lattice.flip_random_fract_edges()

   # Python backend (slower, more flexible)
   cp_python = ContactProcess(
       sg=lattice,
       gamma=0.5,
       runlang='py'
   )

   # C backend (faster, optimized)
   cp_c = ContactProcess(
       sg=lattice,
       gamma=0.5,
       runlang='C1c'
   )

   print(f"Python backend initialized: {cp_python.runlang}")
   print(f"C backend initialized: {cp_c.runlang}")

Next Steps
----------

Now that you've seen the basics, explore more advanced features:

- :doc:`user_guide/graphs` - Detailed guide on signed graph creation
- :doc:`user_guide/lattices` - Working with 2D and 3D lattices
- :doc:`user_guide/spectral` - Spectral analysis and renormalization
- :doc:`user_guide/dynamics` - Statistical physics simulations
- :doc:`user_guide/plotting` - Visualization techniques
- :doc:`user_guide/examples` - Complete end-to-end examples

Or dive into the API documentation:

- :doc:`api/index` - Complete API reference

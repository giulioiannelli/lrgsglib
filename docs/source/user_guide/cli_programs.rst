CLI Program Reference
=====================

lrgsglib includes standalone programs for batch computation on clusters
or workstations. Each program combines a graph type with a computation
type and writes results to the ``data/`` directory tree.

.. contents:: Contents
   :local:
   :depth: 2

Common Flags
------------

All programs share a set of common flags:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Flag
     - Short
     - Description
   * - ``--graph_engine``
     - ``-ge``
     - Graph backend: ``nx`` (NetworkX, default) or ``gt`` (graph-tool)
   * - ``--number_of_averages``
     - ``-na``
     - Number of disorder realizations to average over
   * - ``--verbose``
     - ``-v``
     - Enable verbose logging output
   * - ``--remove_files``
     - ``-rf``
     - Remove intermediate binary files after processing

Spectral Programs
-----------------

L2D_SlaplSpect.py
~~~~~~~~~~~~~~~~~

Compute Laplacian eigenvalues/eigenvectors for 2D lattices.

.. code-block:: bash

   # Full spectrum (all eigenvalues)
   python src/L2D_SlaplSpect.py 16 0.3 --mode eigvals --howmany 0 -na 5 -g sqr

   # With GT backend
   python src/L2D_SlaplSpect.py 16 0.3 --mode eigvals --howmany 0 -na 2 -ge gt

   # Entropy mode
   python src/L2D_SlaplSpect.py 32 0.0 --mode entropy -na 3 --backend scipy

**Positional args:** ``side`` (lattice size), ``pflip`` (frustration fraction)

MCG_SlaplSpect.py
~~~~~~~~~~~~~~~~~

Compute spectra for multiplicative cascade graphs.

.. code-block:: bash

   # Eigenvalue computation
   python src/MCG_SlaplSpect.py 0.8 0.6 0 0 0.4 3 --mode eigvals --howmany 0 -na 2

   # Entropy computation
   python src/MCG_SlaplSpect.py 0.8 0.6 0 0 0.4 3 --mode entropy -na 2

   # With GT backend
   python src/MCG_SlaplSpect.py 0.8 0.6 0 0 0.4 3 --mode eigvals -na 1 -ge gt

**Positional args:** ``p1``, ``p2``, ``p3``, ``p4`` (cascade probabilities),
``fraction``, ``iterations``

Dynamics Programs
-----------------

L2D_IsingDynamics.py
~~~~~~~~~~~~~~~~~~~~

Run Ising model simulations on 2D lattices.

.. code-block:: bash

   # Python backend
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.27 -rl py -na 1 -ic random

   # C backend with energy/magnetization output
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.27 -rl C1b -na 1 -ic random -rf

   # Pybind11 backend
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.0 -rl pb_met -na 1

   # GT engine with pybind11
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.0 -rl pb_met -na 1 -ge gt

   # Cluster algorithms
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.27 -rl wolff -na 1

**Key flags:** ``-T`` (temperature), ``-rl`` (runlang backend),
``-ic`` (initial condition: ``random``, ``ground_state_0``, etc.)

ER_IsingDynamics.py
~~~~~~~~~~~~~~~~~~~

Run Ising dynamics on Erdos-Renyi random graphs.

.. code-block:: bash

   python src/ER_IsingDynamics.py 100 5 0.0 -T 2.0 -rl C1 -na 2 -ic ground_state_0

   # GT engine
   python src/ER_IsingDynamics.py 100 5 0.0 -T 2.0 -rl py -na 1 -ge gt

**Positional args:** ``N`` (nodes), ``k`` (average degree), ``pflip``

L2D_ContactProcess.py
~~~~~~~~~~~~~~~~~~~~~

Contact process dynamics (excitation-inhibition or SIR) on 2D lattices.

.. code-block:: bash

   # C backend (EI dynamics)
   python src/L2D_ContactProcess.py 16 0.0 -dy EI -rl C1c -ga 1.5 -na 3 -rf

   # Python backend
   python src/L2D_ContactProcess.py 16 0.0 -dy EI -rl py -ga 1.5 -na 2

   # GT engine
   python src/L2D_ContactProcess.py 16 0.0 -dy EI -rl py -ga 1.5 -na 2 -ge gt

**Key flags:** ``-dy`` (dynamics type: ``EI`` or ``SIR``),
``-ga`` (coupling gamma), ``-rl`` (runlang)

TransCluster Programs
---------------------

L2D_TransCluster.py
~~~~~~~~~~~~~~~~~~~

Compute transient cluster distributions and order parameters on 2D lattices.

.. code-block:: bash

   # Cluster distribution mode
   python src/L2D_TransCluster.py 16 0.3 --mode pCluster -na 10 -g sqr -c rand

   # Order parameter mode
   python src/L2D_TransCluster.py 16 0.3 --mode ordParam -na 10 -g sqr -c rand

   # GT engine
   python src/L2D_TransCluster.py 16 0.3 --mode pCluster -na 5 -ge gt

**Modes:** ``pCluster`` (cluster size distribution), ``ordParam`` (Pinf, gap)

L3D_TransCluster.py
~~~~~~~~~~~~~~~~~~~

Same as 2D but for 3D lattice geometries. Additional flags for edge weight
models and bond dilution.

SCS_TransCluster.py
~~~~~~~~~~~~~~~~~~~

Transient clustering on SCS generalized networks (scale-free, non-lattice).

Reconstruction Programs
-----------------------

L2D_Recon.py / L3D_Recon.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Reconstruct Laplacian eigenmodes from Ising dynamics spin configurations.

.. code-block:: bash

   # 2D lattice reconstruction
   python src/L2D_Recon.py 16 0.0 -T 0.1 -rl C1b -na 5

   # 3D lattice reconstruction
   python src/L3D_Recon.py 8 0.0 -T 0.1 -rl C1b -na 3

   # GT engine (uses pybind11 automatically)
   python src/L2D_Recon.py 16 0.0 -T 0.1 -rl pb_met -na 3 -ge gt

Graph Engine Compatibility
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 10 10 60

   * - Program
     - NX
     - GT
     - Notes
   * - L2D_SlaplSpect
     - Yes
     - Yes
     - Full spectral support
   * - MCG_SlaplSpect
     - Yes
     - Yes
     - Via MultispectralGraphGT
   * - L2D_IsingDynamics
     - Yes
     - Yes
     - GT: py/pb/cu/cluster backends only
   * - ER_IsingDynamics
     - Yes
     - Yes
     - GT: py/pb/cu backends only
   * - L2D_ContactProcess
     - Yes
     - Yes
     - GT: Python backend only
   * - L2D_TransCluster
     - Yes
     - Yes
     - GT: ``cell=rand`` only
   * - L3D_TransCluster
     - Yes
     - No
     - NX only (complex edge weight models)
   * - SCS_TransCluster
     - Yes
     - No
     - SCS networks are NX-only
   * - L2D_Recon
     - Yes
     - Yes
     - GT: requires pb/py backend
   * - L3D_Recon
     - Yes
     - Yes
     - GT: requires pb/py backend

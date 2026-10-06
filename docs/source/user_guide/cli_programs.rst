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
     - Graph backend: ``nx`` (NetworkX, default) or ``gt`` (graph-tool).
       Where the graph-tool path is not ported, ``gt`` exits with code 2
       and says why (see :ref:`cli-engine-compat`)
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
   python src/L2D_SlaplSpect.py 16 -p 0.3 --mode eigvals --howmany 0 -na 5 -g sqr

   # With GT backend
   python src/L2D_SlaplSpect.py 16 -p 0.3 --mode eigvals --howmany 0 -na 2 -ge gt

   # Entropy mode
   python src/L2D_SlaplSpect.py 32 -p 0.0 --mode entropy -na 3 --backend scipy

**Positional args:** ``side``/``L`` (lattice size).
**Key flags:** ``-p``/``--pflip`` (frustration fraction, default 0.0)

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
   python src/L2D_IsingDynamics.py 16 -p 0.0 -T 2.27 -rl py -na 1 -ic random

   # C backend with energy/magnetization output
   python src/L2D_IsingDynamics.py 16 -p 0.0 -T 2.27 -rl C1b -na 1 -ic random -rf

   # Pybind11 backend
   python src/L2D_IsingDynamics.py 16 -p 0.0 -T 2.0 -rl pb_met -na 1

   # GT engine with pybind11
   python src/L2D_IsingDynamics.py 16 -p 0.0 -T 2.0 -rl pb_met -na 1 -ge gt

   # Cluster algorithms
   python src/L2D_IsingDynamics.py 16 -p 0.0 -T 2.27 -rl wolff -na 1

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
   python src/L2D_ContactProcess.py 16 -p 0.0 -dy EI -rl C1c -ga 1.5 -na 3 -rf

   # Python backend
   python src/L2D_ContactProcess.py 16 -p 0.0 -dy EI -rl py -ga 1.5 -na 2

   # GT engine
   python src/L2D_ContactProcess.py 16 -p 0.0 -dy EI -rl py -ga 1.5 -na 2 -ge gt

**Key flags:** ``-dy`` (dynamics type: ``EI`` or ``SIR``),
``-ga`` (coupling gamma), ``-rl`` (runlang)

TransCluster Programs
---------------------

L2D_TransCluster.py
~~~~~~~~~~~~~~~~~~~

Compute transient cluster distributions and order parameters on 2D lattices.

.. code-block:: bash

   # Cluster distribution mode
   python src/L2D_TransCluster.py 16 -p 0.3 --mode pCluster -na 10 -g sqr -c rand

   # Order parameter mode
   python src/L2D_TransCluster.py 16 -p 0.3 --mode ordParam -na 10 -g sqr -c rand

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

.. warning::

   ``L3D_Recon.py`` currently stops with a ``TypeError``: the kernel passes
   ``compute`` to the 3D lattice constructor. ``L2D_Recon.py`` runs as shown.

.. code-block:: bash

   # 2D lattice reconstruction
   python src/L2D_Recon.py 16 -p 0.0 -T 0.1 -rl C1b -na 5

   # 3D lattice reconstruction
   python src/L3D_Recon.py 8 -p 0.0 -T 0.1 -rl C1b -na 3

.. _cli-engine-compat:

Graph Engine Compatibility
--------------------------

``-ge gt`` exits with code 2 and a one-line reason wherever a graph-tool run
would drop an option or is not ported yet. The SignedRW bubble and
time-to-touch programs leave ``-ge`` unset by default and then follow
``LRGSG_GRAPH_ENGINE`` (else ``nx``), as the graph factories do. The other
programs default to ``nx``, and ``RRG_TransCluster`` to ``gt``.

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
     - GT refuses ``-m eigvec_dist``, the default mode, so pass ``-m``.
       GT file names carry ``_ge=gt``
   * - MCG_SlaplSpect
     - Yes
     - Yes
     - GT refuses ``-m eigvec_dist``. GT file names carry ``_ge=gt``
   * - L2D_IsingDynamics
     - Yes
     - Yes
     - GT: py/pb/cu/cluster backends only; refuses a cell other than
       ``rand``, ``-wd`` and ``--save-results``
   * - L3D_IsingDynamics
     - Yes
     - Yes
     - GT: as L2D_IsingDynamics, and refuses ``--pdil``
   * - ER_IsingDynamics
     - Yes
     - Yes
     - GT: py/pb/cu backends only (C runlangs are refused); refuses a cell
       other than ``rand`` and ``-wd``
   * - L2D_ContactProcess
     - Yes
     - Yes
     - GT: Python backend only; refuses a cell other than ``rand`` and ``-wd``
   * - L2D_VoterModel, L2D_SignedRW
     - Yes
     - Yes
     - GT refuses a cell other than ``rand`` and ``-wd``
   * - L2D_SignedRW_Bubbles, L2D_SignedRW_BubblesRand,
       L2D_SignedRW_TimeToTouch
     - Yes
     - Yes
     - GT refuses ``-wd``, so pass ``-wd ''`` (BubblesRand and
       TimeToTouch set one by default)
   * - RRG_TransCluster
     - Yes
     - Yes
     - Defaults to ``gt``
   * - L2D_TransCluster
     - Yes
     - No
     - Refused: graph-tool eigenvectors lack the majority flip
   * - L3D_TransCluster
     - Yes
     - No
     - NX only (complex edge weight models)
   * - SCS_TransCluster, SCS_TransCluster_sweep
     - Yes
     - No
     - SCS networks are NX-only
   * - L2D_Recon, L3D_Recon
     - Yes
     - No
     - Refused: the kernel passes NetworkX-only options to the lattice
   * - L2D_SignedRW_Bubbles_Scan
     - Yes
     - No
     - Refused: ``--skip-existing`` looks for the NetworkX output layout

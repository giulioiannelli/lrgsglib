Theoretical Background
======================

This section provides theoretical background on the mathematical and physical
concepts underlying lrgsglib.

.. toctree::
   :maxdepth: 2

   signed_graphs
   lrg_theory
   stat_phys

Overview
--------

lrgsglib implements the **Laplacian Renormalization Group for Signed Graphs**,
a framework that connects network structure to dynamical processes through
spectral analysis.

Core Concepts
~~~~~~~~~~~~~

**Signed Graphs**
   Networks where edges carry signs (+1 or -1) representing positive or
   negative interactions. Unlike standard graphs, signed graphs can exhibit
   frustration when constraints cannot all be satisfied.

**Spectral Analysis**
   The eigenvalues and eigenvectors of the signed Laplacian encode structural
   information about the network, including frustration levels, community
   structure, and diffusion properties.

**Laplacian Renormalization**
   A coarse-graining procedure that reveals how network properties change
   across scales, connecting microscopic structure to macroscopic behavior.

**Statistical Physics**
   Models like the Ising model, contact process, and voter model study how
   local interactions on networks lead to emergent collective behavior.

Mathematical Framework
----------------------

The Signed Laplacian
~~~~~~~~~~~~~~~~~~~~

For a signed graph G with N nodes and signed adjacency matrix A:

.. math::

   L = D - A

where D is the absolute degree matrix. Key properties:

- Symmetric for undirected graphs
- Can have negative eigenvalues (frustration signature)
- Spectral gap controls mixing/diffusion times

See :doc:`signed_graphs` for complete theory.

Information from Spectra
~~~~~~~~~~~~~~~~~~~~~~~~

The spectral density matrix:

.. math::

   \rho(\tau) = \frac{e^{-\tau L}}{\text{Tr}[e^{-\tau L}]}

yields the Von Neumann entropy:

.. math::

   S(\tau) = -\text{Tr}[\rho \log \rho]

The spectral dimension :math:`d_s` characterizes scaling:

.. math::

   S(\tau) \sim \frac{d_s}{2} \log \tau

See :doc:`lrg_theory` for the renormalization framework.

Dynamical Models
~~~~~~~~~~~~~~~~

**Ising Model:**

.. math::

   H = -\sum_{(i,j)} J_{ij} s_i s_j, \quad s_i \in \{-1, +1\}

**Contact Process:**

- Infection: S → I at rate λ × (infected neighbors)
- Recovery: I → S at rate μ

**Voter Model:**

- Agent adopts neighbor's opinion (modified by edge sign)

See :doc:`stat_phys` for detailed formulations.

Quick Reference
---------------

Key Equations
~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Quantity
     - Definition
     - Equation
   * - Signed Laplacian
     - L = D - A
     - :math:`L_{ij} = \delta_{ij} d_i - A_{ij}`
   * - Spectral entropy
     - Von Neumann entropy
     - :math:`S = -\sum_i p_i \log p_i`
   * - Spectral dimension
     - Scaling exponent
     - :math:`d_s = 2 \lim_{\tau \to \infty} C(\tau)`
   * - Frustration index
     - Min edges to balance
     - :math:`\phi(G) = \min_{\sigma'} |E_{\text{wrong}}|`

Key Results
~~~~~~~~~~~

- **Balanced graphs**: All eigenvalues ≥ 0, λ₁ = 0
- **Frustrated graphs**: Negative eigenvalues, ``|λ₁|`` measures frustration
- **Spectral gap**: Controls relaxation times: τ ∼ 1/Δλ
- **Localization**: Frustration can localize eigenvectors

Connections to Physics
----------------------

Spin Glasses
~~~~~~~~~~~~

Signed graphs naturally model spin glass systems:

- Positive edges = ferromagnetic (align)
- Negative edges = antiferromagnetic (anti-align)
- Frustration = competing interactions

The ground state problem maps to balance optimization.

Network Science
~~~~~~~~~~~~~~~

Spectral methods connect to:

- Community detection (spectral clustering)
- Diffusion processes (random walks)
- Epidemic spreading (contact process)
- Opinion dynamics (voter model)

Renormalization
~~~~~~~~~~~~~~~

The LRG framework draws from statistical physics:

- Coarse-graining preserves essential structure
- Fixed points identify universal behavior
- Spectral dimension generalizes embedding dimension

Implementation Notes
--------------------

lrgsglib provides:

1. **Graph types** - SignedGraph, Lattice2D/3D, ErdosRenyi, etc.
2. **Spectral analysis** - Eigenvalue/eigenvector computation with multiple backends
3. **Entropy calculation** - Von Neumann entropy and specific heat
4. **Dynamics simulation** - Ising, contact process, voter model with C backends

See the :doc:`/user_guide/index` for practical usage.

References
----------

Key papers on signed graphs and spectral methods:

.. note::

   A comprehensive bibliography is planned for future releases.
   Key topics include:

   - Structural balance theory (Heider, Cartwright-Harary)
   - Spectral graph theory (Chung, Mohar)
   - Spin glasses (Edwards-Anderson, Sherrington-Kirkpatrick)
   - Contact process (Harris, Liggett)
   - Renormalization group (Wilson, Kadanoff)

Further Reading
---------------

- :doc:`/user_guide/spectral` - Practical spectral analysis
- :doc:`/user_guide/dynamics` - Running simulations
- :doc:`/user_guide/examples` - End-to-end examples
- :doc:`/api/index` - API reference

Theoretical Background
======================

This section provides theoretical background on the mathematical and physical concepts underlying lrgsglib.

.. note::
   This section is currently under development. For detailed theoretical background,
   please refer to the academic papers listed in the references below.

Topics
------

The theory section will cover:

- Laplacian Renormalization Group framework
- Signed graph theory and balance
- Spectral properties of signed graphs
- Statistical physics models on networks
- Frustration and ground states

Laplacian Renormalization Group
--------------------------------

The Laplacian Renormalization Group (LRG) for signed graphs is a theoretical framework for analyzing the spectral properties of networks with both positive and negative edges.

Key concepts:

- **Signed Laplacian**: Generalization of the graph Laplacian to signed graphs
- **Frustration**: Measure of incompatibility with balanced partitions
- **Spectral flow**: Evolution of eigenvalues under renormalization
- **Coarse-graining**: Reduction of graph complexity while preserving spectral properties

Signed Graph Theory
-------------------

Signed graphs extend ordinary graphs by allowing edge weights of +1 (positive) or -1 (negative).

- **Balance**: A graph is balanced if it can be partitioned into two sets with all positive edges within sets and negative edges between sets
- **Frustration index**: Minimum number of edges to remove to achieve balance
- **Structural balance theory**: Social and physical interpretations

Statistical Physics on Networks
--------------------------------

lrgsglib implements several statistical physics models:

- **Ising model**: Spin dynamics with ferromagnetic/antiferromagnetic interactions
- **Contact process**: Epidemic spreading with recovery
- **Voter model**: Opinion dynamics and consensus formation

References
----------

Key papers on the Laplacian Renormalization Group and signed graphs:

.. note::
   References will be added in the theory section development phase.

Coming Soon
-----------

Detailed sections on:

- Mathematical formulation of the LRG
- Spectral properties and eigenvalue bounds
- Statistical physics models: formulation and analysis
- Numerical methods and algorithms
- Connection to renormalization group theory in physics

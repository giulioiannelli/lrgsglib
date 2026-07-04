statsys - Statistical Physics Simulations
==========================================

The ``statsys`` module provides Python wrappers and interfaces for statistical physics simulations on networks.

.. currentmodule:: lrgsglib.statsys

Overview
--------

This module includes simulation systems organized by state type:

**Binary state** (spins :math:`\pm 1`):

- **Ising Dynamics**: Spin dynamics with Metropolis, simulated annealing, and parallel tempering
- **Contact Process**: Epidemic spreading (SIR) and excitation-inhibition (EI) dynamics
- **Voter Model**: Opinion dynamics and consensus formation

**Continuous state** (real-valued):

- **Kuramoto Model**: Coupled oscillator synchronization
- **Reaction-Diffusion Model**: Pattern formation on networks
- **Coupled ODE Model**: General coupled ODE systems on graphs

**Vector state** (multi-component):

- **Potts Model**: :math:`q`-state spin model
- **XY Model**: Planar spin model
- **Heisenberg Model**: 3D spin model
- **MultiSpecies Model**: Multi-component population dynamics

**Other**:

- **Signed Random Walks**: Dynamics on signed networks

Base Classes
------------

DynSys (Abstract Base)
^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.DynSys
   :members:
   :undoc-members:
   :show-inheritance:

BinDynSys (Binary State)
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.BinDynSys
   :members:
   :undoc-members:
   :show-inheritance:

ContDynSys (Continuous State)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.ContDynSys
   :members:
   :undoc-members:
   :show-inheritance:

VecDynSys (Vector State)
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.VecDynSys
   :members:
   :undoc-members:
   :show-inheritance:

Solver Registry
---------------

Every dynamics model dispatches its ``run()`` through a solver backend
selected via ``runlang=``. The registry maps a ``(model, backend)`` pair to a
:class:`~lrgsglib.statsys.Solver` implementation; the available backends are
enumerated by :class:`~lrgsglib.statsys.SolverBackend` (``py``, ``pb``,
``np``, ``cu``, ``c``). New backends can be registered without modifying the
model classes.

.. autoclass:: lrgsglib.statsys.Solver
   :members:

.. autoclass:: lrgsglib.statsys.SolverBackend
   :members:

.. autofunction:: lrgsglib.statsys.get_solver

.. autofunction:: lrgsglib.statsys.register_solver

.. autofunction:: lrgsglib.statsys.list_solvers_for_model

.. autofunction:: lrgsglib.statsys.list_solver_models

.. autofunction:: lrgsglib.statsys.is_backend_available

Binary State Models
-------------------

Ising Dynamics
^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.IsingDynamics
   :members:
   :undoc-members:
   :show-inheritance:

Contact Process
^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.ContactProcess
   :members:
   :undoc-members:
   :show-inheritance:

Voter Model
^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.VoterModel
   :members:
   :undoc-members:
   :show-inheritance:

Continuous State Models
-----------------------

Kuramoto Model
^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.KuramotoModel
   :members:
   :undoc-members:
   :show-inheritance:

Reaction-Diffusion Model
^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.ReactionDiffusionModel
   :members:
   :undoc-members:
   :show-inheritance:

Coupled ODE Model
^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.CoupledODEModel
   :members:
   :undoc-members:
   :show-inheritance:

Vector State Models
-------------------

Potts Model
^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.PottsModel
   :members:
   :undoc-members:
   :show-inheritance:

XY Model
^^^^^^^^

.. automodule:: lrgsglib.statsys.XYModel
   :members:
   :undoc-members:
   :show-inheritance:

Heisenberg Model
^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.HeisenbergModel
   :members:
   :undoc-members:
   :show-inheritance:

MultiSpecies Model
^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.MultiSpeciesModel
   :members:
   :undoc-members:
   :show-inheritance:

Other Models
------------

Signed Random Walks
^^^^^^^^^^^^^^^^^^^

.. automodule:: lrgsglib.statsys.SignedRW
   :members:
   :undoc-members:
   :show-inheritance:

"""Statistical physics systems for dynamics on signed graphs.

This module contains all dynamics models organized in a folder-per-class
pattern.  Each model lives in its own subdirectory with optional ``ccore/``
subfolder for C backend source and binaries.
"""

# Solver-backend registry (engine-style dispatch for dynamics)
from ._solver import Solver, SolverBackend
from ._solver_engine import (
    get_solver,
    is_backend_available,
    list_solver_models,
    list_solvers_for_model,
    register_solver,
)

# Base classes
from .DynSys import DynSys
from .BinDynSys import BinDynSys
from .ContDynSys import ContDynSys
from .VecDynSys import VecDynSys

# Binary-state dynamics
from .IsingDynamics import IsingDynamics, IsingModel
from .ContactProcess import (
    ContactProcess,
    ContactProcessBase,
    ContactProcessEI,
    ContactProcessSIR,
)
from .SignedRW import SignedRW
from .VoterModel import VoterModel

# Continuous-state dynamics
from .KuramotoModel import KuramotoModel
from .ReactionDiffusionModel import ReactionDiffusionModel
from .CoupledODEModel import CoupledODEModel

# Vector-state dynamics
from .PottsModel import PottsModel
from .XYModel import XYModel
from .HeisenbergModel import HeisenbergModel
from .MultiSpeciesModel import MultiSpeciesModel

__all__ = [
    # Solver registry
    "Solver",
    "SolverBackend",
    "get_solver",
    "register_solver",
    "list_solvers_for_model",
    "list_solver_models",
    "is_backend_available",
    # Base classes
    "DynSys",
    "BinDynSys",
    "ContDynSys",
    "VecDynSys",
    # Binary
    "IsingModel",
    "IsingDynamics",
    "ContactProcess",
    "ContactProcessBase",
    "ContactProcessEI",
    "ContactProcessSIR",
    "SignedRW",
    "VoterModel",
    # Continuous
    "KuramotoModel",
    "ReactionDiffusionModel",
    "CoupledODEModel",
    # Vector
    "PottsModel",
    "XYModel",
    "HeisenbergModel",
    "MultiSpeciesModel",
]

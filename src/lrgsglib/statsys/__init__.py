"""Statistical physics systems for dynamics on signed graphs."""

# Base classes
from .DynSys import DynSys
from .BinDynSys import BinDynSys
from .ContDynSys import ContDynSys
from .VecDynSys import VecDynSys

# Binary-state dynamics
from .IsingDynamics import IsingDynamics
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
    # Base classes
    "DynSys",
    "BinDynSys",
    "ContDynSys",
    "VecDynSys",
    # Binary
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

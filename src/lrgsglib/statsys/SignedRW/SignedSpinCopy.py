"""Signed spin-copy dynamics on a signed graph.

At each step a node copies a random neighbour's spin multiplied by the
connecting edge sign. This is **not** a position-tracked random walker;
for the walker family see :mod:`SignedWalker`.
"""

from __future__ import annotations

import random
import subprocess
from typing import TYPE_CHECKING

import numpy as np

from ..BinDynSys import BinDynSys
from .._solver import SolverBackend
from .._solver_engine import get_solver
from .defaults import SRW_COPY_SOLVER_NAME

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph


class SignedSpinCopy(BinDynSys):
    """Spin-copy-with-sign dynamics.

    At each step, node ``nd`` sets its spin to
    ``s[nd] = w[nd, nn] * s[nn]`` where ``nn`` is a uniformly random
    neighbour of ``nd`` and ``w`` is the signed edge weight.
    """

    dyn_UVclass = "signed_spin_copy"
    id_string_signedrw = ""

    def __init__(self, sg: "SignedGraph", **kwargs) -> None:
        dynpath = getattr(sg, "path_srw", None)
        super().__init__(sg, dynpath=dynpath, **kwargs)

    def ds1step(self, nd: int) -> None:
        """One spin-copy-with-sign update at node ``nd``."""
        nodedict = dict(self.sg.H[nd])
        neighs = list(nodedict.keys())
        nn = np.random.choice(neighs)
        self.s[nd] = nodedict[nn]["weight"] * self.s[nn]

    def run_py(self) -> None:
        """Run the dynamics with the Python backend."""
        dsNstep = self.dsNstep()
        nodes = list(self.sg.H.nodes())
        for _ in range(self.steps):
            smp = random.sample(nodes, self.sg.N)
            self.s_t.append(self.s.copy())
            dsNstep(smp)

    def run_c(
        self,
        adjfname: str = "",
        out_suffix: str = "",
        eqSTEP: int = 0,
    ) -> None:
        """Run using the C backend (not yet implemented)."""
        if eqSTEP:
            self.eqSTEP_def = eqSTEP
        if adjfname == "":
            adjfname = self.sg.std_fname
        out_suffix = out_suffix + self.id_string_signedrw
        if out_suffix == "":
            out_suffix = "''"
        self.cprogram = [
            f"{self.sg.N}",
            f"{self.sg.pflip}",
            f"{self.eqSTEP_def}",
            str(self.sg.path_data),
            adjfname,
            out_suffix,
        ]
        subprocess.call(self.cprogram)

    def _resolve_backend(self) -> SolverBackend:
        """Map the runlang code to a solver family (``py`` / ``C``)."""
        if self.runlang.startswith("py"):
            return SolverBackend.PY
        if self.runlang.startswith("C"):
            return SolverBackend.C
        raise ValueError(
            f"SignedSpinCopy supports 'py' or 'C' backends, got {self.runlang!r}."
        )

    def run(self, **kwargs) -> None:
        """Dispatch to the Python or C backend through the shared solver registry."""
        backend = self._resolve_backend()
        solver = get_solver(SRW_COPY_SOLVER_NAME, backend)
        solver.supports(self)
        solver.execute(self, **kwargs)

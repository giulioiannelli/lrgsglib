"""Signed random walk dynamics on graphs."""

from __future__ import annotations

import random
import subprocess
from typing import TYPE_CHECKING

import numpy as np

from .BinDynSys import BinDynSys

if TYPE_CHECKING:
    from ..nx_patches import SignedGraph


class SignedRW(BinDynSys):
    """Signed random walk dynamics.

    At each step, a node copies the spin of a neighbor multiplied by the
    edge sign (weight). This creates correlated spin patterns based on
    the underlying signed graph structure.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph to run dynamics on.
    **kwargs
        Additional arguments passed to BinDynSys.
    """

    dyn_UVclass = "signed_rw"
    id_string_signedrw = ""

    def __init__(self, sg: "SignedGraph", **kwargs) -> None:
        dynpath = getattr(sg, "path_signedrw", None)
        super().__init__(sg, dynpath=dynpath, **kwargs)

    def ds1step(self, nd: int) -> None:
        """Perform one signed random walk step for node nd.

        The node copies the spin of a random neighbor, multiplied by
        the edge sign between them.
        """
        nodedict = dict(self.sg.H[nd])
        neighs = list(nodedict.keys())
        nn = np.random.choice(neighs)
        self.s[nd] = nodedict[nn]["weight"] * self.s[nn]

    def run_py(self) -> None:
        """Run signed random walk using Python backend."""
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
        """Run signed random walk using C backend.

        Note: C backend is incomplete and may not work.
        """
        if eqSTEP:
            self.eqSTEP_def = eqSTEP
        if adjfname == "":
            adjfname = self.sg.std_fname
        out_suffix = out_suffix + self.id_string_signedrw
        if out_suffix == "":
            out_suffix = "''"
        self.cprogram = [
            # TODO: Add path to C binary
            f"{self.sg.N}",
            f"{self.sg.pflip}",
            f"{self.eqSTEP_def}",
            str(self.sg.path_data),
            adjfname,
            out_suffix,
        ]
        subprocess.call(self.cprogram)

    def run(self, **kwargs) -> None:
        """Run signed random walk dynamics.

        Uses Python or C backend based on self.runlang setting.
        """
        if self.runlang.startswith("py"):
            self.run_py()
        elif self.runlang.startswith("C"):
            self.run_c(**kwargs)

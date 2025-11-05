"""Voter-model dynamics on signed graphs."""

from __future__ import annotations

import random
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import tqdm

from .BinDynSys import BinDynSys
from ..config.const import LOG, LRGSG_CCORE_BIN, LRGSG_LOG
from ..utils.basic.strings import join_non_empty
from ..utils.tools.chronometer import time_function_accumulate


class VoterModel(BinDynSys):
    """Binary voter dynamics with optional C backend."""

    dyn_UVclass = "voter_model"
    magn: list[float] = []
    s_t: list[np.ndarray] = []

    def __init__(
        self,
        sg,
        *,
        eqSTEP: int = 20,
        save_magnetization: bool = False,
        upd_mode: str = "asynchronous",
        freq: int = 10,
        nSampleLog: int = 100,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, 'path_voter', None)
        super().__init__(sg, dynpath=dynpath, **kwargs)
        self.eqSTEP = int(eqSTEP)
        self.save_magnetization = save_magnetization
        self.upd_mode = upd_mode
        self.freq = freq
        self.nSampleLog = nSampleLog
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.stderr_path: Path | None = None
        self.stderr_fopen = None
        self.cprogram: list[str | Path] = []
        self.out_id: str = self.out_suffix
        self.magn_path: Path | None = None

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset cached observables collected during a run."""

        self.magn = []
        self.s_t = []

    def init_voter_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the spin configuration and export data if required."""

        self.reset_observables()
        self.init_s(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            self.sg._export_edgel_bin(exName=self.run_id)
            self.sg.export_adj_bin(exName=self.run_id)
        self.sini = self.s.copy()

    def check_attribute(self) -> None:
        try:
            getattr(self, "CbaseName")
        except AttributeError:
            self.init_voter_dynamics()

    def initialize_run_parameters(self, eqSTEP: int | None = None) -> None:
        if eqSTEP is not None:
            self.eqSTEP = int(eqSTEP)
        self.simtime = self.eqSTEP

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, nd: int) -> None:
        neigh = self.sg.gr[self.sg.on_g][nd]
        if not neigh:
            return
        neighbours = list(neigh)
        choice = neighbours[np.random.randint(len(neighbours))]
        weight = neigh[choice].get("weight", 1.0)
        edge_sign = -1 if weight < 0 else 1
        self.s[nd] = np.int8(edge_sign * self.s[choice])

    def voter_sampling(self, tqdm_on: bool) -> None:
        dsNstep = self.dsNstep()
        nodes = list(self.sg.gr[self.sg.on_g].nodes())
        iterator = tqdm.tqdm(range(self.eqSTEP)) if tqdm_on else range(self.eqSTEP)
        for _ in iterator:
            if self.save_magnetization:
                self.magn.append(float(np.sum(self.s)) / float(self.N))
            if self.savedyn:
                self.s_t.append(self.s.copy())
            dsNstep(random.sample(nodes, len(nodes)))

    def run_py(self, tqdm_on: bool = False) -> None:
        self.voter_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def build_cprogram_command(self) -> None:
        self.CbaseName = f"VoterSimulator{self.runlang[-1]}"
        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except ValueError:
            datdir = self.sg.path_sgdata
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self.out_suffix
        self.magn_path = self.dynpath / self.sg.get_p_fname('m', self.out_id)
        
        # Base arguments common to all simulators
        arglist = [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.eqSTEP}",
            str(datdir),
            syshape,
            self.run_id,
            self.out_id,
        ]
        
        # Add extra arguments for VoterSimulator1 and above
        if self.runlang[-1] != "0":
            arglist.append(f"{self.nSampleLog}")
        
        self.cprogram = [LRGSG_CCORE_BIN / self.CbaseName] + arglist

    def setup_stderr_logging(self) -> None:
        fname = join_non_empty(
            '_',
            f"err{self.CbaseName}",
            f"{self.N}",
            self.run_id,
            self.out_suffix,
        ) + LOG
        self.stderr_path = LRGSG_LOG / fname
        self.stderr_path.parent.mkdir(parents=True, exist_ok=True)
        self.stderr_fopen = open(self.stderr_path, 'w')

    def run_cprogram(self, verbose: bool = False) -> None:
        if not self.cprogram:
            raise RuntimeError("C program command has not been initialised.")
        binary_path = Path(self.cprogram[0])
        if not binary_path.exists():
            if self.stderr_fopen and not self.stderr_fopen.closed:
                self.stderr_fopen.close()
            raise FileNotFoundError(
                f"C backend executable '{binary_path}' was not found. Build the "
                "C components (e.g. via `make c-make`) before running the C backend."
            )
        result = subprocess.run(
            self.cprogram,
            stderr=self.stderr_fopen,
            stdout=subprocess.PIPE,
            check=False,
        )
        if self.stderr_fopen and not self.stderr_fopen.closed:
            self.stderr_fopen.close()
        state = np.frombuffer(result.stdout, dtype=np.int8)
        if state.size != self.N:
            raise RuntimeError("C backend returned an invalid spin configuration.")
        self.s = state.copy()
        if self.magn_path and self.magn_path.exists():
            self.magn = np.fromfile(self.magn_path, dtype=np.float64).tolist()
        try:
            self.sfout.unlink()
        except FileNotFoundError:
            pass

    def _remove_magn(self) -> None:
        if not self.magn_path:
            return
        try:
            self.magn_path.unlink()
        except FileNotFoundError:
            pass

    def _remove_stderr(self) -> None:
        if not self.stderr_path:
            return
        try:
            self.stderr_path.unlink()
        except FileNotFoundError:
            pass

    def remove_run_c_files(self, remove_stderr: bool = True) -> None:
        self._remove_sfout()
        self._remove_magn()
        if remove_stderr:
            self._remove_stderr()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        eqSTEP: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        self.check_attribute()
        self.initialize_run_parameters(eqSTEP)
        if self.runlang.startswith("C"):
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.voter_sampling(tqdm_on)

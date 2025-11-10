"""Contact process dynamics for signed graphs."""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import tqdm

from .BinDynSys import BinDynSys
from ..config.const import LOG, LRGSG_CCORE_BIN, LRGSG_LOG
from ..nx_patches import SignedGraph
from ..utils.basic.strings import join_non_empty
from ..utils.tools.chronometer import time_function_accumulate


class ContactProcess(BinDynSys):
    """Simple contact process dynamics on a weighted signed graph."""

    dyn_UVclass = "contact_process"
    density: list[float] = []
    s_t: list[np.ndarray] = []

    def __init__(
        self,
        sg: SignedGraph,
        mu: float = 1.0,
        *,
        save_density: bool = False,
        **kwargs: Any,
    ) -> None:
        dynpath = getattr(sg, "path_contact", None)
        if dynpath is None:
            dynpath = getattr(sg, "path_lrgsg", getattr(sg, "path_data", Path.cwd()))
        super().__init__(sg, dynpath=dynpath, **kwargs)
        self.mu = float(mu)
        self.save_density = save_density
        self.reset_observables()
        self.sini: np.ndarray | None = None
        self.stderr_path: Path | None = None
        self.stderr_fopen = None
        self.cprogram: list[str | Path] = []
        self.out_id: str = self.out_suffix

    # ------------------------------------------------------------------
    # Utility helpers
    # ------------------------------------------------------------------
    def reset_observables(self) -> None:
        """Reset cached observables collected during a run."""

        self.density = []
        self.s_t = []

    def _iter_neighbour_data(self, node: int) -> Iterable[tuple[int, float]]:
        graph = self.sg.gr[self.sg.on_g]
        if hasattr(graph, "is_directed") and graph.is_directed():
            iterator = graph.predecessors(node)
            for neighbour in iterator:
                data = graph[neighbour][node]
                weight = data.get("weight", 1.0)
                yield neighbour, weight
        else:
            iterator = graph.neighbors(node)
            for neighbour in iterator:
                data = graph[node][neighbour]
                weight = data.get("weight", 1.0)
                yield neighbour, weight

    def init_contact_dynamics(self, custom: Any = None, exName: str = "") -> None:
        """Initialise the state and export data if required."""

        self.reset_observables()
        self.init_state(custom)
        self.s = self.s.astype(np.int8, copy=False)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            if self.rndStr and not exName:
                exName = self.run_id
            self.sg._export_edgel_bin(exName=self.run_id)
        self.sini = self.s.copy()

    def init_state(self, custom: Any = None) -> None:
        match self.ic:
            case "uniform" | "random" | "rand":
                self.s = np.random.choice([0, 1], size=self.sg.N).astype(np.int8)
            case "homogeneous" | "homo":
                self.s = np.ones(self.sg.N, dtype=np.int8)
            case "delta":
                self.s = np.zeros(self.sg.N, dtype=np.int8)
                self.s[np.random.randint(self.sg.N)] = 1
            case "custom":
                if custom is None:
                    raise ValueError("A custom state must be provided for 'custom' initial condition.")
                self.s = np.asarray(custom, dtype=np.int8)
            case _:
                raise ValueError("Invalid initial condition for ContactProcess.")

    def check_attribute(self) -> None:
        try:
            getattr(self, "CbaseName")
        except AttributeError:
            self.init_contact_dynamics()

    def initialize_run_parameters(self, steps: int | None = None) -> None:
        if steps is not None:
            self.simtime = int(steps)

    # ------------------------------------------------------------------
    # Python dynamics
    # ------------------------------------------------------------------
    def ds1step(self, node: int) -> None:
        neighbours = list(self._iter_neighbour_data(node))
        if self.s[node]:
            rate = self.mu
            for neighbour, weight in neighbours:
                if weight < 0.0 and self.s[neighbour]:
                    rate += -weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(0)
        else:
            rate = 0.0
            for neighbour, weight in neighbours:
                if weight > 0.0 and self.s[neighbour]:
                    rate += weight
            prob = 1.0 - np.exp(-rate)
            if prob > 0.0 and np.random.random() < prob:
                self.s[node] = np.int8(1)

    def contact_sampling(self, tqdm_on: bool) -> None:
        dsNstep = self.dsNstep()
        nodes = np.arange(self.N)
        iterator = tqdm.tqdm(range(self.simtime)) if tqdm_on else range(self.simtime)
        for _ in iterator:
            if self.save_density:
                self.density.append(float(np.mean(self.s)))
            if self.savedyn:
                self.s_t.append(self.s.copy())
            np.random.shuffle(nodes)
            dsNstep(nodes)

    def run_py(self, tqdm_on: bool = False) -> None:
        self.contact_sampling(tqdm_on)

    # ------------------------------------------------------------------
    # C backend integration
    # ------------------------------------------------------------------
    def build_cprogram_command(self) -> None:
        self.CbaseName = f"ContactSimulator{self.runlang[-1]}"
        try:
            datdir = self.sg.path_sgdata.relative_to(Path.cwd())
        except ValueError:
            datdir = self.sg.path_sgdata
        syshape = getattr(self.sg, "syshapePth", f"N={self.N}")
        self.out_id = self.out_suffix
        arglist = [
            f"{self.N}",
            f"{self.sg.pflip:.12g}",
            f"{self.mu:.12g}",
            f"{self.simtime}",
            str(datdir),
            syshape,
            self.run_id,
            self.out_id,
        ]
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
            self.stderr_fopen = None
            raise FileNotFoundError(
                f"C backend executable '{binary_path}' was not found. "
                "Build the C components (e.g. via `make c-make`) before running the C backend."
            )
        result = subprocess.run(
            self.cprogram,
            stderr=self.stderr_fopen,
            stdout=subprocess.PIPE,
            check=False,
        )
        if self.stderr_fopen and not self.stderr_fopen.closed:
            self.stderr_fopen.close()
        self.stderr_fopen = None
        state = np.frombuffer(result.stdout, dtype=np.int8)
        if state.size != self.N:
            raise RuntimeError("C backend returned an invalid state configuration.")
        self.s = state.copy()

    def _remove_sfout(self) -> None:
        try:
            self.sfout.unlink()
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
        if remove_stderr:
            self._remove_stderr()

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------
    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        steps: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
    ) -> None:
        self.check_attribute()
        self.initialize_run_parameters(steps)
        if self.runlang.startswith("C"):
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            self.contact_sampling(tqdm_on)

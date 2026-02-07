"""Ising model dynamics with Metropolis, Simulated Annealing, and Parallel Tempering."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import tqdm

from ._c_backend import CBackendMixin
from .BinDynSys import BinDynSys
from ..config.const import BIN, SG_REPR
from ..config.funcs import build_pT_fname
from ..graphs.nx.funcs import get_kth_order_neighbours
from ..utils.lrg.ising import compute_ising_pairwise_energy
from ..utils.statsys import boltzmann_factor
from ..utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ..graphs.nx import SignedGraphNX as SignedGraph

# Type aliases for cooling/ladder schedules
CoolingSchedule = Literal["linear", "exponential", "logarithmic", "custom"]
LadderType = Literal["geometric", "linear", "custom"]


class IsingDynamics(CBackendMixin, BinDynSys):
    """Ising model dynamics with Metropolis, Simulated Annealing, and Parallel Tempering.

    Naming convention for C backends:
        - Number = Algorithm: 1=Metropolis, 3=SA, 4=PT
        - Letter = Output: a=final, b=E/M, c=snapshots, d=clusters, e=eigvec, f=exchange

    Parameters
    ----------
    sg : SignedGraph
        The signed graph to run dynamics on.
    T : float, optional
        Temperature for equilibrium simulations.
    steps : int, optional
        Number of Monte Carlo sweeps.
    simref : float, optional
        Size-normalized time (steps = simref * N).
    sa_enabled : bool, optional
        Enable simulated annealing mode.
    pt_enabled : bool, optional
        Enable parallel tempering mode.
    **kwargs
        Additional arguments passed to BinDynSys.

    Examples
    --------
    >>> ising = IsingDynamics(lattice, T=2.0, steps=1000, runlang='C1b')
    >>> ising.init_ising_dynamics()
    >>> ising.run(verbose=False)

    Notes
    -----
    C backend variants:
        - ``runlang="C1b"``: Metropolis with E/M output
        - ``runlang="C3b"``: Simulated Annealing with E/M output
        - ``runlang="C4b"``: Parallel Tempering with E/M output
    """

    dyn_UVclass = "ising_dynamics"

    # CBackendMixin configuration
    _c_program_name_template = "IsingSimulator{}"
    _allowed_c_keys: tuple[str, ...] = (
        # New Metropolis variants
        "C1B",
        # Simulated Annealing variants
        "C3B",
        # Parallel Tempering variants
        "C4B",
        # Legacy aliases (deprecated but supported)
        "C0", "C1", "C2", "C3", "C4", "C5",
    )
    magn = []
    ene = []
    s_t = []
    Ising_clusters = []
    k_B = 1

    def __init__(
        self,
        sg: "SignedGraph",
        T: float = 0.,
        *,
        NoClust: int = 1,
        thrmSTEP: int = 20,
        eqSTEP: int = 20,
        freq: int = 10,
        steps: int | None = None,
        simref: float | None = None,
        nstepsIsing: int | None = None,
        save_magnetization: bool = False,
        upd_mode: str = "asynchronous",
        # Simulated Annealing parameters
        sa_enabled: bool = False,
        T_init: float = 10.0,
        T_final: float = 0.01,
        cooling_schedule: CoolingSchedule = "exponential",
        cooling_rate: float = 0.95,
        steps_per_T: int = 100,
        n_temperatures: int = 100,
        custom_T_schedule: np.ndarray | None = None,
        # Parallel Tempering parameters
        pt_enabled: bool = False,
        n_replicas: int = 8,
        T_min: float = 0.5,
        T_max: float = 5.0,
        T_ladder_type: LadderType = "geometric",
        custom_T_ladder: np.ndarray | None = None,
        steps_per_exchange: int = 10,
        n_exchanges: int = 1000,
        **kwargs
    ) -> None:
        dynpath = getattr(sg, 'path_ising', None)
        resolved_steps = steps
        if resolved_steps is None:
            resolved_steps = nstepsIsing
        if resolved_steps is None:
            resolved_steps = eqSTEP
        super(IsingDynamics, self).__init__(sg, dynpath=dynpath, steps=resolved_steps, simref=simref, **kwargs)
        self.T = T
        self.thrmSTEP = thrmSTEP
        self.freq = freq
        self.save_magnetization = save_magnetization
        self.NoClust = NoClust
        self.NeigV = -1
        self.upd_mode = upd_mode

        # Simulated Annealing parameters
        self.sa_enabled = sa_enabled
        self.T_init = T_init
        self.T_final = T_final
        self.cooling_schedule = cooling_schedule
        self.cooling_rate = cooling_rate
        self.steps_per_T = steps_per_T
        self.n_temperatures = n_temperatures
        self.custom_T_schedule = custom_T_schedule

        # Parallel Tempering parameters
        self.pt_enabled = pt_enabled
        self.n_replicas = n_replicas
        self.T_min = T_min
        self.T_max = T_max
        self.T_ladder_type = T_ladder_type
        self.custom_T_ladder = custom_T_ladder
        self.steps_per_exchange = steps_per_exchange
        self.n_exchanges = n_exchanges

    @property
    def eqSTEP(self) -> int:
        """Compatibility alias for equilibration sweeps."""

        return self.steps

    @eqSTEP.setter
    def eqSTEP(self, value: int) -> None:
        self._set_time_controls(steps=value)

    @property
    def nstepsIsing(self) -> int:
        """Compatibility alias mapping legacy field to the sweep count."""

        return self.steps

    @nstepsIsing.setter
    def nstepsIsing(self, value: int) -> None:
        self._set_time_controls(steps=value)
    #
    def neigh_ene(self, neigh: list) -> float:
        return np.sum(neigh) / len(neigh)
    #
    def neigh_wghtmagn(self, node: int,on_g: str = SG_REPR) -> list:
        nd = dict(self.sg.gr[on_g][node])
        return [w["weight"] * self.s[nn] for nn, w in nd.items()]
    #
    def metropolis(self, node):
        neigh = self.neigh_wghtmagn(node)
        neighene = self.neigh_ene(neigh)
        DeltaE = 2 * self.s[node] * neighene
        if DeltaE < 0:
            self.flip_spin(node)
        elif np.random.uniform() < boltzmann_factor(DeltaE, self.T, self.k_B):
            self.flip_spin(node)
    #
    def calc_full_energy(self) -> float:
        neigh_energies = np.array([
            self.neigh_ene(self.neigh_wghtmagn(node)) for node in range(self.sg.N)
        ])
        return -np.dot(self.s, neigh_energies)
    #
    def init_ising_dynamics(self, custom: Any = None, exName: str = ""):
        self._check_c_backend_or_fallback()
        self.init_s(custom)
        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.setup_stderr_logging()
            self.export_s_init()
            self.export_hfield()
            if self.rndStr and not exName:
                exName = self.run_id
            exName_arg = exName if exName else self.run_id
            # Export uses build_p_fname which handles underscore via join_non_empty
            self.sg._export_edgel_bin(exName=exName_arg)
        self.sini = self.s.copy()
    #
    def check_attribute(self) -> None:
        """Initialize dynamics if not already done."""
        # Check sini (set at end of init_ising_dynamics) rather than CbaseName
        # because CbaseName has a class-level default from CBackendMixin
        if not hasattr(self, 'sini') or self.sini is None:
            self.init_ising_dynamics()

    def initialize_run_parameters(
        self,
        T_ising,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
    ):
        if T_ising:
            self.T = T_ising
        chosen_steps = steps if steps is not None else eqSTEP
        self._set_time_controls(steps=chosen_steps, simref=simref)

    # ------------------------------------------------------------------
    # C backend configuration (overrides CBackendMixin)
    # ------------------------------------------------------------------
    def _c_program_key(self) -> str:
        """Convert runlang to internal key (e.g., 'c1b' -> 'C1B').

        Overrides mixin to default to C1B when no suffix is provided.
        """
        if not self.runlang or not self.runlang.upper().startswith("C"):
            raise ValueError("C backend requires runlang starting with 'C'")
        suffix = self.runlang[1:]
        key = f"C{suffix.upper()}" if suffix else "C1B"  # Default to C1B
        if self._allowed_c_keys and key not in self._allowed_c_keys:
            raise ValueError(
                f"runlang '{self.runlang}' not supported for {self.__class__.__name__}. "
                f"Allowed: {self._allowed_c_keys}"
            )
        return key

    def _c_program_suffix(self) -> str:
        """Extract program suffix for IsingSimulator (handles legacy variants).

        Legacy single-digit variants (C0-C5) map to simple digits.
        New variants (C1B, C3B, etc.) use lowercase suffix.
        """
        key = self._c_program_key()
        legacy_map = {"C0": "0", "C1": "1", "C2": "2", "C3": "3", "C4": "4", "C5": "5"}
        if key in legacy_map:
            return legacy_map[key]
        return key[1:].lower()  # "C1B" -> "1b", "C3B" -> "3b"

    def _is_sa_variant(self) -> bool:
        """Check if current runlang is simulated annealing (new variant only).

        Note: C3 is a legacy variant that uses equilibrium arguments.
        Only C3B uses the new simulated annealing argument format.
        """
        key = self._c_program_key()
        return key.startswith("C3") and key != "C3"  # C3 is legacy

    def _is_pt_variant(self) -> bool:
        """Check if current runlang is parallel tempering."""
        key = self._c_program_key()
        return key.startswith("C4") and key != "C4"  # C4 is legacy

    # ------------------------------------------------------------------
    # Argument building (model-specific)
    # ------------------------------------------------------------------
    def _build_c_arglist(self) -> list[str]:
        """Build argument list based on variant type."""
        if self._is_pt_variant():
            return self._build_pt_arglist()
        elif self._is_sa_variant():
            return self._build_sa_arglist()
        else:
            return self._build_equilibrium_arglist()

    def _build_equilibrium_arglist(self) -> list[str]:
        """Arguments for standard Metropolis equilibrium variants."""
        return [
            f"{self.N}",
            f"{self.T:.3g}",
            f"{self.sg.pflip:.3g}",
            f"{self.NoClust}",
            f"{self.thrmSTEP:.3g}",
            f"{self.steps}",
            self._get_datdir_arg(),
            self.sg.syshapePth,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
            self.upd_mode,
            f"{self.freq}",
            f"{self.NeigV}"
        ]

    def _build_sa_arglist(self) -> list[str]:
        """Arguments for simulated annealing variants."""
        return [
            f"{self.N}",
            f"{self.sg.pflip:.3g}",
            f"{self.T_init:.6g}",
            f"{self.T_final:.6g}",
            self.cooling_schedule,
            f"{self.cooling_rate:.6g}",
            f"{self.steps_per_T}",
            f"{self.n_temperatures}",
            self._get_datdir_arg(),
            self.sg.syshapePth,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
            self.upd_mode,
        ]

    def _build_pt_arglist(self) -> list[str]:
        """Arguments for parallel tempering variants."""
        ladder_type_map = {"geometric": "0", "linear": "1", "custom": "2"}
        return [
            f"{self.N}",
            f"{self.sg.pflip:.3g}",
            f"{self.n_replicas}",
            f"{self.T_min:.6g}",
            f"{self.T_max:.6g}",
            ladder_type_map.get(self.T_ladder_type, "0"),
            f"{self.steps_per_exchange}",
            f"{self.n_exchanges}",
            self._get_datdir_arg(),
            self.sg.syshapePth,
            self._c_suffix_arg(self.run_id),
            self._c_suffix_arg(self.out_suffix),
            self.upd_mode,
        ]
    # ------------------------------------------------------------------
    # C backend execution (overrides CBackendMixin)
    # ------------------------------------------------------------------
    def run_cprogram(self, verbose: bool = False) -> None:
        """Execute C backend and read model-specific output files.

        Calls parent implementation for subprocess execution and state parsing,
        then reads energy/magnetization files for variants that output them.
        """
        super().run_cprogram(verbose)
        # Read energy/magnetization files for 'b' variants
        if self._output_includes_em():
            self._read_c_em_output()

    def _get_cleanup_paths(self) -> list[Path | None]:
        """Return paths to clean up after C run.

        Includes model-specific output files (energy, magnetization, etc.)
        in addition to standard state output.
        """
        paths: list[Path | None] = [
            getattr(self, 'sfout', None),
            getattr(self, 'hfout', None),
        ]
        # Add E/M output files if they exist
        if hasattr(self, '_c_output_paths'):
            paths.extend(self._c_output_paths)
        return paths
    #
    def _output_includes_em(self) -> bool:
        """Check if current variant outputs energy/magnetization time series."""
        key = self._c_program_key()
        return key.endswith('B') or key in ('C0',)
    #
    def _read_c_em_output(self) -> None:
        """Read energy and magnetization binary files from C backend output."""
        # Determine the temperature parameter for filename
        key = self._c_program_key()
        if self._is_sa_variant():
            T_for_fname = self.T_init
        elif self._is_pt_variant():
            T_for_fname = self.T_min
        else:
            T_for_fname = self.T

        # Build filenames matching C backend output format
        ene_fname = build_pT_fname(
            'ene', self.sg.pflip, T_for_fname,
            out_suffix=self.out_suffix, ext=BIN
        )
        magn_fname = build_pT_fname(
            'm', self.sg.pflip, T_for_fname,
            out_suffix=self.out_suffix, ext=BIN
        )

        ene_path = self.dynpath / ene_fname
        magn_path = self.dynpath / magn_fname

        # Read energy file
        if ene_path.exists():
            self.ene = np.fromfile(ene_path, dtype=np.float64)
            if self._is_pt_variant():
                # Reshape to (n_replicas, n_exchanges)
                self.ene = self.ene.reshape(self.n_replicas, -1)
        else:
            self.ene = np.array([])

        # Read magnetization file
        if magn_path.exists():
            self.magn = np.fromfile(magn_path, dtype=np.float64)
            if self._is_pt_variant():
                # Reshape to (n_replicas, n_exchanges)
                self.magn = self.magn.reshape(self.n_replicas, -1)
        else:
            self.magn = np.array([])

        # Store paths for cleanup
        self._c_output_paths = [ene_path, magn_path]

        # For PT, also read temperature ladder
        if self._is_pt_variant():
            tladder_fname = build_pT_fname(
                'Tladder', self.sg.pflip, T_for_fname,
                out_suffix=self.out_suffix, ext=BIN
            )
            tladder_path = self.dynpath / tladder_fname
            if tladder_path.exists():
                self.T_ladder = np.fromfile(tladder_path, dtype=np.float64)
                self._c_output_paths.append(tladder_path)

        # For SA, read temperature schedule
        if self._is_sa_variant():
            tsched_fname = build_pT_fname(
                'Tsched', self.sg.pflip, T_for_fname,
                out_suffix=self.out_suffix, ext=BIN
            )
            tsched_path = self.dynpath / tsched_fname
            if tsched_path.exists():
                self.sa_temps = np.fromfile(tsched_path, dtype=np.float64)
                self._c_output_paths.append(tsched_path)
    #
    def metropolis_sampling(self, tqdm_on):
        metropolis_1step = np.vectorize(self.metropolis, excluded="self")
        if self.save_magnetization:
            def save_magn_array():
                self.s_t.append(self.s)
        else:
            def save_magn_array():
                pass

        sample = list(range(self.sg.N))
        iterator = tqdm.tqdm(range(self.steps), desc="Metropolis") if tqdm_on \
            else range(self.steps)
        self.ene = []
        self.magn = []
        for _ in iterator:
            self.magn.append(np.sum(self.s))
            self.ene.append(self.calc_full_energy())
            metropolis_1step(sample)
            save_magn_array()

    # =========================================================================
    # Simulated Annealing Methods
    # =========================================================================

    def _generate_temperature_schedule(self) -> np.ndarray:
        """Generate temperature array based on cooling_schedule parameter."""
        if self.custom_T_schedule is not None:
            return self.custom_T_schedule

        n = self.n_temperatures

        if self.cooling_schedule == "linear":
            return np.linspace(self.T_init, self.T_final, n)
        elif self.cooling_schedule == "exponential":
            return self.T_init * (self.cooling_rate ** np.arange(n))
        elif self.cooling_schedule == "logarithmic":
            return self.T_init / np.log(2 + np.arange(n))
        else:
            raise ValueError(f"Unknown cooling schedule: {self.cooling_schedule}")

    def simulated_annealing_sampling(self, tqdm_on: bool = True) -> None:
        """Pure Python simulated annealing implementation."""
        T_schedule = self._generate_temperature_schedule()

        self.sa_energy = []
        self.sa_magn = []
        self.sa_temps = T_schedule.copy()

        outer_iter = tqdm.tqdm(T_schedule, desc="SA") if tqdm_on else T_schedule

        for T in outer_iter:
            self.T = T
            for _ in range(self.steps_per_T):
                # Single MC sweep
                for node in range(self.N):
                    self.metropolis(node)

            # Record observables at this temperature
            self.sa_energy.append(self.calc_full_energy())
            self.sa_magn.append(np.sum(self.s) / self.N)

        self.sa_energy = np.array(self.sa_energy)
        self.sa_magn = np.array(self.sa_magn)

    def get_sa_observables_by_temperature(self) -> dict:
        """Reshape SA observables by temperature."""
        return {
            'ene': self.sa_energy,
            'magn': self.sa_magn,
            'T': self.sa_temps
        }

    # =========================================================================
    # Parallel Tempering Methods
    # =========================================================================

    def _generate_T_ladder(self) -> np.ndarray:
        """Generate temperature ladder for PT."""
        if self.custom_T_ladder is not None:
            return self.custom_T_ladder

        n = self.n_replicas
        if self.T_ladder_type == "geometric":
            return self.T_min * (self.T_max / self.T_min) ** (np.arange(n) / (n - 1))
        elif self.T_ladder_type == "linear":
            return np.linspace(self.T_min, self.T_max, n)
        else:
            raise ValueError(f"Unknown T_ladder_type: {self.T_ladder_type}")

    def _attempt_exchange(self, replicas: list, i: int, j: int, T_ladder: np.ndarray) -> bool:
        """Attempt replica exchange using Metropolis criterion."""
        E_i = self._calc_replica_energy(replicas[i])
        E_j = self._calc_replica_energy(replicas[j])
        delta_beta = 1.0 / T_ladder[i] - 1.0 / T_ladder[j]
        delta_E = E_i - E_j

        if delta_beta * delta_E <= 0 or np.random.random() < np.exp(delta_beta * delta_E):
            # Swap spin configurations
            replicas[i], replicas[j] = replicas[j].copy(), replicas[i].copy()
            return True
        return False

    def _calc_replica_energy(self, s: np.ndarray) -> float:
        """Calculate energy for a given spin configuration."""
        # Temporarily swap s and calculate energy
        old_s = self.s
        self.s = s
        energy = self.calc_full_energy()
        self.s = old_s
        return energy

    def parallel_tempering_sampling(self, tqdm_on: bool = True) -> None:
        """Pure Python parallel tempering implementation."""
        T_ladder = self._generate_T_ladder()
        n_rep = self.n_replicas

        # Initialize replicas (each gets a copy of initial state)
        replicas = [self.s.copy() for _ in range(n_rep)]

        # Storage for observables
        self.pt_energy = np.zeros((n_rep, self.n_exchanges))
        self.pt_magn = np.zeros((n_rep, self.n_exchanges))
        self.pt_exchanges = []  # List of (round, i, j, accepted)

        outer_iter = tqdm.tqdm(range(self.n_exchanges), desc="PT") if tqdm_on else range(self.n_exchanges)

        for ex_round in outer_iter:
            # Run MC sweeps on all replicas
            for r in range(n_rep):
                self.T = T_ladder[r]
                self.s = replicas[r]
                for _ in range(self.steps_per_exchange):
                    for node in range(self.N):
                        self.metropolis(node)
                replicas[r] = self.s.copy()

                # Record observables
                self.pt_energy[r, ex_round] = self._calc_replica_energy(replicas[r])
                self.pt_magn[r, ex_round] = np.mean(replicas[r])

            # Attempt exchanges (even-odd alternation)
            start = ex_round % 2
            for i in range(start, n_rep - 1, 2):
                accepted = self._attempt_exchange(replicas, i, i + 1, T_ladder)
                self.pt_exchanges.append((ex_round, i, i + 1, accepted))

        # Store final states
        self.pt_final_states = replicas
        self.T_ladder = T_ladder

    def get_pt_exchange_rate(self) -> np.ndarray:
        """Compute exchange acceptance rate between each pair of adjacent replicas."""
        n_rep = self.n_replicas
        rates = np.zeros(n_rep - 1)
        counts = np.zeros(n_rep - 1)

        for (_, i, j, accepted) in self.pt_exchanges:
            idx = min(i, j)
            counts[idx] += 1
            if accepted:
                rates[idx] += 1

        return rates / np.maximum(counts, 1)

    def compute_energy(self, spins: np.ndarray | None = None) -> float:
        """Compute Ising energy from spin configuration.

        Parameters
        ----------
        spins : np.ndarray, optional
            Spin configuration to evaluate. If None, uses current state `self.s`.

        Returns
        -------
        float
            Pairwise Ising energy E = -∑_{(u,v)} w_uv * s_u * s_v
        """
        if spins is None:
            spins = self.s
        edges = list(self.sg.gr[self.sg.on_g].edges(data='weight', default=1.0))
        return compute_ising_pairwise_energy(spins, edges)

    # =========================================================================
    # Run Method
    # =========================================================================

    @time_function_accumulate(auto_log=False)
    def run(
        self,
        tqdm_on: bool = True,
        T_ising: float = None,
        steps: int | None = None,
        simref: float | None = None,
        eqSTEP: int | None = None,
        verbose: bool = False,
        clean_export: bool = True,
        sa_mode: bool = False,
        pt_mode: bool = False,
    ):
        """
        Run Ising dynamics simulation.

        Parameters
        ----------
        tqdm_on : bool
            Show progress bar.
        T_ising : float, optional
            Override temperature.
        steps : int, optional
            Override number of steps.
        simref : float, optional
            Reference simulation time.
        eqSTEP : int, optional
            Legacy parameter for equilibration steps.
        verbose : bool
            Verbose output.
        clean_export : bool
            Remove exported files after run.
        sa_mode : bool
            Run simulated annealing.
        pt_mode : bool
            Run parallel tempering.
        """
        self.check_attribute()
        self.initialize_run_parameters(T_ising, steps=steps, simref=simref, eqSTEP=eqSTEP)

        # Determine mode
        use_sa = sa_mode or self.sa_enabled
        use_pt = pt_mode or self.pt_enabled

        if self.runlang.startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()
        else:
            # Python backends
            if use_pt:
                self.parallel_tempering_sampling(tqdm_on)
            elif use_sa:
                self.simulated_annealing_sampling(tqdm_on)
            else:
                self.metropolis_sampling(tqdm_on)
        

    #
    def find_ising_clusters(self, import_cl: bool = False):
        #can be easily reworked
        if import_cl:
            for i in range(self.NoClust):
                self.Ising_clusters.append(
                    np.fromfile(
                        f"{self.sg.path_ising}cl{i}_{self.sg.std_fname}.bin",
                        dtype=int,
                    )
                )
            self.numIsing_cl = len(self.Ising_clusters)
        if self.Ising_clusters:
            return
        #
        self.sg.compute_k_eigvV(k=self.NoClust)
        eigVbin = self.sg.get_eigV_bin_check_list()
        #
        self.Ising_clusters = []
        for j in range(self.NoClust):
            lnodes = list(self.sg.H.nodes())
            lnodes_tmp = lnodes[:]

            def recursive_search(seed, magn_i, clustertmp):
                neighs = get_kth_order_neighbours(self.sg.H, seed, 1)
                neighs = np.array(
                    [e for e in neighs if e not in set(clustertmp)]
                )
                if not neighs.size:
                    return
                samecluster = np.array(eigVbin[j][neighs] == magn_i)
                if not samecluster.any():
                    return
                neighs_samecluster = list(neighs[samecluster])
                clustertmp.extend(neighs_samecluster)
                for ss in neighs_samecluster:
                    recursive_search(ss, magn_i, clustertmp)

            allclusters = []
            for i in lnodes:
                if i not in lnodes_tmp:
                    continue
                if not lnodes_tmp:
                    break
                #
                clustertmp = []
                clustertmp.extend([i])
                #
                recursive_search(i, eigVbin[j][i], clustertmp)
                lnodes_tmp = [e for e in lnodes_tmp if e not in set(clustertmp)]
                allclusters.append(clustertmp)
            allclusters.sort(key=len, reverse=True)
            self.Ising_clusters.append(allclusters[0])
        self.numIsing_cl = len(self.Ising_clusters)
        if self.runlang.startswith("C"):
            self.sg.export_ising_clust()

    #
    def mapping_nodes_to_clusters(self):
        if not self.Ising_clusters:
            self.find_ising_clusters()
        loc = [x for x in range(len(self.Ising_clusters))]
        self.loc = loc
        node_with_inherclust = [
            [[j, loc[i]] for j in clus]
            for i, clus in enumerate(self.Ising_clusters)
        ]
        self.node_with_inherclust = node_with_inherclust
        node_inherclust_flat = [i for j in node_with_inherclust for i in j]
        self.node_inherclust_flat = node_inherclust_flat
        sorted_list = sorted(node_inherclust_flat, key=lambda x: x[0])
        self.sorted_list = sorted_list
        result_array = np.empty((self.sg.side1, self.sg.side2), dtype=object)
        self.result_array = result_array

        # Fill the result_array with tuples from sorted_list
        for i, sublist in enumerate(sorted_list):
            row, col = divmod(
                i, self.sg.side1
            )  # Calculate the row and column index
            result_array[row, col] = sublist[1]
        self.mapping = result_array


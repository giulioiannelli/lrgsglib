"""Ising model dynamics with Metropolis, Simulated Annealing, and Parallel Tempering."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np
import tqdm

from .._c_backend import CBackendMixin
from ..BinDynSys import BinDynSys
from ...config.const import BIN, SG_REPR
from ...config.funcs import build_pT_fname
from ...utils.lrg.ising import compute_ising_pairwise_energy
from ...utils.statsys import boltzmann_factor
from ...utils.tools.chronometer import time_function_accumulate

if TYPE_CHECKING:
    from ...graphs.protocols import SignedGraphProtocol as SignedGraph

# Type aliases for cooling/ladder schedules
CoolingSchedule = Literal["linear", "exponential", "logarithmic", "custom"]
LadderType = Literal["geometric", "linear", "custom"]

# ---------------------------------------------------------------------------
# runlang naming convention
# ---------------------------------------------------------------------------
# New systematic format: C<digit><letters>
#   digit:   0 = Metropolis, 1 = SA, 2 = PT
#   letters: E = energy+magnetization, S = snapshots, K = cluster,
#            V = eigvec overlap, H = external field
#
# Bare digit (e.g. C0) defaults to E (energy+magnetization).
# Compound: C0ES = Metropolis + energy/magnetization + snapshots.
#
# Non-C backends: py, pb_met, cu_met, etc.

import warnings as _warnings

# Algorithm digit → backend key
_ISING_ALGORITHMS: dict[str, str] = {
    "0": "c_met",
    "1": "c_sa",
    "2": "c_pt",
}

# Save-mode letter → C binary flag
_ISING_SAVE_LETTERS: dict[str, str] = {
    "E": "em",
    "S": "snap",
    "K": "clust",
    "V": "eigvec",
    "H": "hfield",
}

_ISING_DEFAULT_SAVE = "em"

# Deprecated legacy codes → (backend, save_flags, new_equivalent)
_ISING_DEPRECATED: dict[str, tuple[str, str, str]] = {
    "C0":  ("c_met", "em",              "C0E"),
    "C1":  ("c_met", "clust",           "C0K"),
    "C1B": ("c_met", "em",              "C0E"),
    "C2":  ("c_met", "snap",            "C0S"),
    "C3":  ("c_met", "em,snap",         "C0ES"),
    "C3B": ("c_sa",  "em",              "C1E"),
    "C4":  ("c_met", "clust,eigvec",    "C0KV"),
    "C4B": ("c_pt",  "em",              "C2E"),
    "C5":  ("c_met", "em,snap,hfield",  "C0ESH"),
}

# Convenience aliases (non-C backends)
_RUNLANG_ALIASES: dict[str, str] = {
    "py":       "py_met",
    "python":   "py_met",
    "Python":   "py_met",
    "pb":       "pb_met",
    "pybind":   "pb_met",
    "cu":       "cu_met",
    "cuda":     "cu_met",
    "gpu":      "cu_met",
    "wolff":    "py_wolff",
    "sw":       "py_sw",
    "topo_met": "py_topo_met",
    "topo_fca": "py_topo_fca",
}

# Canonical algorithm → unified binary name
_C_UNIFIED_BINARY: dict[str, str] = {
    "c_met": "IsingMetropolis",
    "c_sa":  "IsingSimulatedAnnealing",
    "c_pt":  "IsingParallelTempering",
}


def _parse_ising_c_runlang(code: str) -> tuple[str, str]:
    """Parse C<digit><letters> → (backend, save_flags).

    Parameters
    ----------
    code : str
        Runlang code starting with 'C' (e.g., 'C0E', 'C0ES', 'C1E').

    Returns
    -------
    tuple[str, str]
        (backend_key, comma_separated_save_flags)

    Raises
    ------
    ValueError
        If the code is invalid.
    """
    upper = code.strip().upper()
    if not upper.startswith("C"):
        raise ValueError(f"Expected C-backend code, got '{code}'")

    rest = upper[1:]
    if not rest:
        raise ValueError(f"Empty runlang code: '{code}'")

    digit = rest[0]
    letters = rest[1:]

    # New-format code: digit is a valid algorithm, letters are all valid save letters
    if digit in _ISING_ALGORITHMS:
        if letters == "" or all(c in _ISING_SAVE_LETTERS for c in letters):
            backend = _ISING_ALGORITHMS[digit]
            if letters:
                flags = ",".join(_ISING_SAVE_LETTERS[c] for c in letters)
            else:
                flags = _ISING_DEFAULT_SAVE
            return backend, flags

    # Check deprecated legacy codes
    full_code = f"C{rest}"
    if full_code in _ISING_DEPRECATED:
        backend, flags, new_equiv = _ISING_DEPRECATED[full_code]
        _warnings.warn(
            f"runlang='{code}' is deprecated. Use '{new_equiv}' instead.",
            DeprecationWarning,
            stacklevel=4,
        )
        return backend, flags

    raise ValueError(
        f"Unknown Ising runlang code: '{code}'. "
        f"Valid: C<digit><letters> where digit={{0,1,2}} and "
        f"letters⊆{{E,S,K,V,H}}, or legacy codes."
    )


def _resolve_runlang(raw: str) -> str:
    """Normalise a runlang string to the canonical form.

    Accepts new systematic codes (``'C0E'``), legacy codes (``'C1b'``),
    and named backends (``'pb_met'``).
    """
    key = raw.strip()
    # Already canonical (contains underscore)?
    if "_" in key:
        return key
    # Non-C convenience alias?
    upper = key.upper()
    if upper in _RUNLANG_ALIASES:
        return _RUNLANG_ALIASES[upper]
    if key in _RUNLANG_ALIASES:
        return _RUNLANG_ALIASES[key]
    # C-backend code: parse via new scheme (handles both new and deprecated)
    if upper.startswith("C"):
        backend, _ = _parse_ising_c_runlang(key)
        return backend
    # Pass through
    return key


class IsingDynamics(CBackendMixin, BinDynSys):
    """Ising model dynamics with Metropolis, Simulated Annealing, and Parallel Tempering.

    Naming convention for ``runlang``
    ---------------------------------
    **C subprocess format:** ``C<digit><letters>``

    ===== ============================================
    digit algorithm
    ===== ============================================
    ``0`` Glauber-Metropolis at fixed *T*
    ``1`` Simulated Annealing
    ``2`` Parallel Tempering (Replica Exchange)
    ===== ============================================

    ====== ============================================
    letter save mode
    ====== ============================================
    ``E``  energy + magnetization time series
    ``S``  spin snapshots
    ``K``  cluster magnetization
    ``V``  eigenvector overlap
    ``H``  external field (read from file)
    ====== ============================================

    Bare digit defaults to ``E``. Combine letters for multiple outputs:
    ``C0ES`` = Metropolis + energy/magnetization + snapshots.

    **Other backends:**

    ========= ========================================
    backend   description
    ========= ========================================
    ``pb_met``  pybind11 Metropolis (fastest, any graph)
    ``cu_met``  CUDA GPU Metropolis (lattices only)
    ``py``      Pure-Python Metropolis
    ========= ========================================

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
    >>> # Pybind11 Metropolis (fastest, engine-agnostic):
    >>> ising = IsingDynamics(lattice, T=2.0, steps=1000, runlang='pb_met')
    >>> ising.init_ising_dynamics()
    >>> ising.run(verbose=False)

    >>> # C subprocess Metropolis with energy+magnetization output:
    >>> ising = IsingDynamics(lattice, T=2.0, steps=1000, runlang='C0E')
    """

    dyn_UVclass = "ising_dynamics"

    # CBackendMixin configuration
    _c_bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    _c_program_name_template = ""  # Overridden by build_cprogram_command
    _allowed_c_keys: tuple[str, ...] = ()
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
        save: str | None = None,
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
        # Topological Metropolis parameters
        topo_n_modes: int = 40,
        topo_sigma_init: float = 0.15,
        topo_proposal_mode: str = "weighted",
        topo_chunk_size: int = 5000,
        topo_polish: bool = True,
        topo_polish_sweeps: int = 50,
        # Topological Field Cooling Annealing parameters
        topo_tau: float = 1.0,
        topo_field_strength: float = 1.0,
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
        self.save = save  # output mode for unified C binaries
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

        # Topological Metropolis parameters
        self.topo_n_modes = topo_n_modes
        self.topo_sigma_init = topo_sigma_init
        self.topo_proposal_mode = topo_proposal_mode
        self.topo_chunk_size = topo_chunk_size
        self.topo_polish = topo_polish
        self.topo_polish_sweeps = topo_polish_sweeps

        # Topological Field Cooling Annealing parameters
        self.topo_tau = topo_tau
        self.topo_field_strength = topo_field_strength

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
    def neigh_wghtmagn(self, node: int) -> list:
        neighbors = self.sg.get_neighbors_with_weights(node)
        return [w * self.s[nn] for nn, w in neighbors]
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
        rl = _resolve_runlang(self.runlang)
        if rl.startswith("c_") or self.runlang.startswith("C"):
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

    def _resolve_c_backend(self) -> str:
        """Determine canonical backend+algorithm from runlang.

        Returns one of: ``'c_met'``, ``'c_sa'``, ``'c_pt'``.
        """
        upper = self.runlang.strip().upper()
        if upper.startswith("C") and "_" not in self.runlang:
            backend, _ = _parse_ising_c_runlang(self.runlang)
            return backend
        # Canonical name: extract backend+algorithm
        rl = _resolve_runlang(self.runlang)
        parts = rl.split("_")
        if len(parts) >= 2 and parts[0] == "c":
            return f"c_{parts[1]}"
        return "c_met"

    def _resolve_save_flags(self) -> str:
        """Determine output mode flags.

        Uses explicit ``self.save`` if set, otherwise derives from
        runlang code (new systematic or deprecated legacy).
        Defaults to ``"em"``.
        """
        if self.save:
            return self.save
        upper = self.runlang.strip().upper()
        if upper.startswith("C") and "_" not in self.runlang:
            _, flags = _parse_ising_c_runlang(self.runlang)
            return flags
        return _ISING_DEFAULT_SAVE

    def _c_program_key(self) -> str:
        """Convert runlang to internal key for unified binaries."""
        if not self.runlang:
            raise ValueError("C backend requires a runlang value")
        backend = self._resolve_c_backend()
        key_map = {
            "c_met": "UNIFIED_MET",
            "c_sa": "UNIFIED_SA",
            "c_pt": "UNIFIED_PT",
        }
        return key_map.get(backend, "UNIFIED_MET")

    def build_cprogram_command(self) -> None:
        """Build command for unified C binary."""
        backend = self._resolve_c_backend()
        binary_name = _C_UNIFIED_BINARY.get(backend, "IsingMetropolis")
        self.CbaseName = binary_name
        arglist = self._build_c_arglist()
        self.cprogram = [self._get_bin_dir() / binary_name] + arglist

    def _is_sa_variant(self) -> bool:
        """Check if current runlang is simulated annealing."""
        return self._resolve_c_backend() == "c_sa"

    def _is_pt_variant(self) -> bool:
        """Check if current runlang is parallel tempering."""
        return self._resolve_c_backend() == "c_pt"

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
        """Arguments for unified IsingMetropolis binary."""
        save_flags = self._resolve_save_flags()
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
            f"{self.NeigV}",
            save_flags,
        ]

    def _build_sa_arglist(self) -> list[str]:
        """Arguments for unified IsingSimulatedAnnealing binary."""
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
        """Arguments for unified IsingParallelTempering binary."""
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
    # Cluster Algorithm Methods (Python fallback)
    # =========================================================================

    def wolff_sampling(self, tqdm_on: bool = True) -> None:
        """Pure Python Wolff single-cluster algorithm."""
        self.ene = []
        self.magn = []
        self.cluster_sizes = []

        iterator = tqdm.tqdm(range(self.steps), desc="Wolff") if tqdm_on \
            else range(self.steps)

        for _ in iterator:
            self.magn.append(np.sum(self.s))
            self.ene.append(self.calc_full_energy())
            cluster_size = self._wolff_sweep()
            self.cluster_sizes.append(cluster_size)

    def _wolff_sweep(self) -> int:
        """Execute Wolff steps until ~N spins have been flipped."""
        total_flipped = 0
        while total_flipped < self.N:
            total_flipped += self._wolff_step()
        return total_flipped

    def _wolff_step(self) -> int:
        """Single Wolff cluster step: grow from random seed, flip."""
        seed = np.random.randint(0, self.N)
        seed_spin = self.s[seed]
        inv_T = 1.0 / self.T if self.T > 0 else 1e12

        in_cluster = np.zeros(self.N, dtype=np.bool_)
        in_cluster[seed] = True
        queue = [seed]
        head = 0

        while head < len(queue):
            node = queue[head]
            head += 1
            for nb, w in self.sg.get_neighbors_with_weights(node):
                if in_cluster[nb]:
                    continue
                # Bond activation: only when aligned relative to coupling
                if self.s[node] * self.s[nb] * (1 if w > 0 else -1) <= 0:
                    continue
                p_add = 1.0 - np.exp(-2.0 * abs(w) * inv_T)
                if np.random.random() < p_add:
                    in_cluster[nb] = True
                    queue.append(nb)

        # Flip entire cluster
        cluster_nodes = np.array(queue, dtype=np.intp)
        self.s[cluster_nodes] *= -1
        return len(queue)

    def swendsen_wang_sampling(self, tqdm_on: bool = True) -> None:
        """Pure Python Swendsen-Wang multi-cluster algorithm."""
        self.ene = []
        self.magn = []
        self.cluster_sizes = []

        iterator = tqdm.tqdm(range(self.steps), desc="SW") if tqdm_on \
            else range(self.steps)

        for _ in iterator:
            self.magn.append(np.sum(self.s))
            self.ene.append(self.calc_full_energy())
            n_clusters = self._sw_step()
            self.cluster_sizes.append(n_clusters)

    def _sw_step(self) -> int:
        """Single Swendsen-Wang step with union-find."""
        inv_T = 1.0 / self.T if self.T > 0 else 1e12
        N = self.N

        # Union-find (array-based)
        parent = np.arange(N)
        rank = np.zeros(N, dtype=np.int32)

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]  # path halving
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            rx, ry = find(x), find(y)
            if rx == ry:
                return
            if rank[rx] < rank[ry]:
                parent[rx] = ry
            elif rank[rx] > rank[ry]:
                parent[ry] = rx
            else:
                parent[ry] = rx
                rank[rx] += 1

        # Phase 1: Bond activation
        for u, v, w in self.sg.get_edges_with_weights():
            # Only activate when aligned relative to coupling
            if self.s[u] * self.s[v] * (1 if w > 0 else -1) <= 0:
                continue
            p_bond = 1.0 - np.exp(-2.0 * abs(w) * inv_T)
            if np.random.random() < p_bond:
                union(u, v)

        # Phase 2: Assign random flip bits per cluster root
        flip_bit = np.random.randint(0, 2, size=N, dtype=np.int8)

        # Phase 3: Flip clusters
        roots = set()
        for i in range(N):
            root = find(i)
            roots.add(root)
            if flip_bit[root]:
                self.s[i] = np.int8(-self.s[i])

        return len(roots)

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
        edges = self.sg.get_edges_with_weights()
        return compute_ising_pairwise_energy(spins, edges)

    # =========================================================================
    # Graph CSR builder (for pybind11 backend)
    # =========================================================================

    @staticmethod
    def _build_graph_csr_from_arrays(
        src: np.ndarray,
        dst: np.ndarray,
        wts: np.ndarray,
        N: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build CSR arrays from edge arrays (vectorized, no Python loops).

        Parameters
        ----------
        src, dst : ndarray[int64]
            Source and target vertex indices (directed edges, one per GT
            edge — this method symmetrises them).
        wts : ndarray[float64]
            Edge weights corresponding to ``(src, dst)`` pairs.
        N : int
            Number of nodes.

        Returns
        -------
        neigh_indices, neigh_weights, neigh_ptr
            CSR representation of the undirected graph.
        """
        # Symmetrise: add both directions
        all_src = np.concatenate([src, dst])
        all_dst = np.concatenate([dst, src])
        all_wts = np.concatenate([wts, wts])

        # Sort by source node
        order = np.argsort(all_src, kind="stable")
        all_dst = all_dst[order]
        all_wts = all_wts[order]
        all_src = all_src[order]

        # Build row pointers via searchsorted
        nodes = np.arange(N + 1, dtype=np.int64)
        ptr = np.searchsorted(all_src, nodes, side="left")

        return (
            all_dst.astype(np.int64),
            all_wts.astype(np.float64),
            ptr,
        )

    def _build_graph_csr(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Build CSR-like arrays from the graph for the pybind11 backend.

        If the graph provides ``get_edge_arrays()`` (GT graphs), uses a
        vectorized numpy path that avoids per-node Python loops.
        Otherwise falls back to the per-node loop via
        ``get_neighbors_with_weights()``.

        Returns
        -------
        neigh_indices : ndarray[int64]
            Concatenated neighbor indices.
        neigh_weights : ndarray[float64]
            Concatenated edge weights (signs).
        neigh_ptr : ndarray[int64]
            Row pointers (N+1,).
        """
        if hasattr(self.sg, "get_edge_arrays"):
            src, dst, wts = self.sg.get_edge_arrays()
            return self._build_graph_csr_from_arrays(
                src, dst, wts, self.N
            )

        # Fallback: per-node loop (fast for NX dict-of-dicts)
        indices_list: list[int] = []
        weights_list: list[float] = []
        ptr = [0]
        for node in range(self.N):
            nw = self.sg.get_neighbors_with_weights(node)
            for nn, w in nw:
                indices_list.append(nn)
                weights_list.append(w)
            ptr.append(len(indices_list))
        return (
            np.array(indices_list, dtype=np.int64),
            np.array(weights_list, dtype=np.float64),
            np.array(ptr, dtype=np.int64),
        )

    # =========================================================================
    # Pybind11 backend
    # =========================================================================

    @staticmethod
    def _load_native_module():
        """Import the compiled _ising_native pybind11 module."""
        try:
            from ..IsingDynamics.ccore import _ising_native  # type: ignore[import-untyped]
            return _ising_native
        except ImportError as exc:
            raise RuntimeError(
                "Pybind11 Ising backend not available. Build the C extensions "
                "with `pip install -e .` or `make all` first."
            ) from exc

    def _run_pybind_met(self) -> None:
        """Run Metropolis via pybind11 native kernels."""
        mod = self._load_native_module()
        ni, nw, nptr = self._build_graph_csr()
        s_out, ene, magn = mod.metropolis_sampling(
            self.s.astype(np.int8),
            ni, nw, nptr,
            self.field.astype(np.float64),
            float(self.T),
            int(self.steps),
            int(self.seed),
            self.upd_mode,
        )
        self.s = s_out
        self.ene = ene.tolist()
        self.magn = magn.tolist()

    def _run_pybind_sa(self) -> None:
        """Run Simulated Annealing via pybind11 native kernels."""
        mod = self._load_native_module()
        ni, nw, nptr = self._build_graph_csr()
        s_out, ene, magn, T_sched = mod.sa_sampling(
            self.s.astype(np.int8),
            ni, nw, nptr,
            self.field.astype(np.float64),
            float(self.T_init),
            float(self.T_final),
            self.cooling_schedule,
            float(self.cooling_rate),
            int(self.steps_per_T),
            int(self.n_temperatures),
            int(self.seed),
            self.upd_mode,
        )
        self.s = s_out
        self.sa_energy = np.asarray(ene)
        self.sa_magn = np.asarray(magn)
        self.sa_temps = np.asarray(T_sched)
        # Also set ene/magn for consistent API
        self.ene = ene.tolist()
        self.magn = magn.tolist()

    def _run_pybind_pt(self) -> None:
        """Run Parallel Tempering via pybind11 native kernels."""
        mod = self._load_native_module()
        ni, nw, nptr = self._build_graph_csr()
        T_ladder = self._generate_T_ladder()
        s_out, ene_2d, magn_2d, exchanges = mod.pt_sampling(
            self.s.astype(np.int8),
            ni, nw, nptr,
            self.field.astype(np.float64),
            T_ladder.astype(np.float64),
            int(self.steps_per_exchange),
            int(self.n_exchanges),
            int(self.seed),
            self.upd_mode,
        )
        self.s = s_out
        self.T_ladder = T_ladder
        self.pt_energy = np.asarray(ene_2d)
        self.pt_magn = np.asarray(magn_2d)
        self.ene = np.asarray(ene_2d)
        self.magn = np.asarray(magn_2d)
        # Rebuild exchange list for API compatibility
        exch = np.asarray(exchanges)
        self.pt_exchanges = []
        n_rep = len(T_ladder)
        for ex_round in range(exch.shape[0]):
            start = ex_round % 2
            for i in range(start, n_rep - 1, 2):
                self.pt_exchanges.append(
                    (ex_round, i, i + 1, bool(exch[ex_round, i]))
                )

    # =========================================================================
    # CuPy (GPU) backend
    # =========================================================================

    def _run_cupy_met(self) -> None:
        """Run Metropolis via CuPy GPU kernels."""
        from ._cupy_ising import cupy_metropolis

        ni, nw, nptr = self._build_graph_csr()
        final_spins, ene_trace, magn_trace = cupy_metropolis(
            self.s.astype(np.int8),
            ni, nw, nptr,
            float(self.T),
            int(self.steps),
            int(self.seed),
        )
        self.s = final_spins
        self.ene = ene_trace
        self.magn = magn_trace

    def _run_cupy_sa(self) -> None:
        """Run Simulated Annealing via CuPy GPU kernels."""
        from ._cupy_ising import cupy_sa

        ni, nw, nptr = self._build_graph_csr()
        final_spins, ene_trace, magn_trace, temps = cupy_sa(
            self.s.astype(np.int8),
            ni, nw, nptr,
            float(self.T_init),
            float(self.T_final),
            self.cooling_schedule,
            float(self.cooling_rate),
            int(self.steps_per_T),
            int(self.n_temperatures),
            int(self.seed),
        )
        self.s = final_spins
        self.sa_energy = np.array(ene_trace)
        self.sa_magn = np.array(magn_trace)
        self.sa_temps = np.array(temps)
        self.ene = ene_trace
        self.magn = magn_trace

    def _run_cupy_pt(self) -> None:
        """Run Parallel Tempering via CuPy GPU kernels."""
        from ._cupy_ising import cupy_pt

        ni, nw, nptr = self._build_graph_csr()
        T_ladder = self._generate_T_ladder()
        final_spins, ene_2d, magn_2d, exchanges = cupy_pt(
            self.s.astype(np.int8),
            ni, nw, nptr,
            T_ladder.astype(np.float64),
            int(self.steps_per_exchange),
            int(self.n_exchanges),
            int(self.seed),
        )
        self.s = final_spins
        self.T_ladder = T_ladder
        self.pt_energy = ene_2d
        self.pt_magn = magn_2d
        self.ene = ene_2d
        self.magn = magn_2d
        # Rebuild exchange list for API compatibility
        n_rep = len(T_ladder)
        self.pt_exchanges = []
        for ex_round in range(exchanges.shape[0]):
            start = ex_round % 2
            for i in range(start, n_rep - 1, 2):
                self.pt_exchanges.append(
                    (ex_round, i, i + 1, bool(exchanges[ex_round, i]))
                )

    # =========================================================================
    # Pybind11 cluster backends
    # =========================================================================

    def _run_pybind_wolff(self) -> None:
        """Run Wolff cluster algorithm via pybind11 native kernels."""
        mod = self._load_native_module()
        ni, nw, nptr = self._build_graph_csr()
        s_out, ene, magn, clusters = mod.wolff_sampling(
            self.s.astype(np.int8),
            ni, nw, nptr,
            self.field.astype(np.float64),
            float(self.T),
            int(self.steps),
            int(self.seed),
        )
        self.s = s_out
        self.ene = ene.tolist()
        self.magn = magn.tolist()
        self.cluster_sizes = clusters.tolist()

    def _run_pybind_sw(self) -> None:
        """Run Swendsen-Wang via pybind11 native kernels."""
        mod = self._load_native_module()
        ni, nw, nptr = self._build_graph_csr()
        s_out, ene, magn, clusters = mod.sw_sampling(
            self.s.astype(np.int8),
            ni, nw, nptr,
            self.field.astype(np.float64),
            float(self.T),
            int(self.steps),
            int(self.seed),
        )
        self.s = s_out
        self.ene = ene.tolist()
        self.magn = magn.tolist()
        self.cluster_sizes = clusters.tolist()

    # =========================================================================
    # CuPy cluster backends
    # =========================================================================

    def _run_cupy_wolff(self) -> None:
        """Run Wolff cluster algorithm via CuPy GPU kernels."""
        from ._cupy_ising import cupy_wolff

        ni, nw, nptr = self._build_graph_csr()
        final_spins, ene_trace, magn_trace, clusters = cupy_wolff(
            self.s.astype(np.int8),
            ni, nw, nptr,
            float(self.T),
            int(self.steps),
            int(self.seed),
        )
        self.s = final_spins
        self.ene = ene_trace
        self.magn = magn_trace
        self.cluster_sizes = clusters

    def _run_cupy_sw(self) -> None:
        """Run Swendsen-Wang via CuPy GPU kernels."""
        from ._cupy_ising import cupy_sw

        ni, nw, nptr = self._build_graph_csr()
        final_spins, ene_trace, magn_trace, clusters = cupy_sw(
            self.s.astype(np.int8),
            ni, nw, nptr,
            float(self.T),
            int(self.steps),
            int(self.seed),
        )
        self.s = final_spins
        self.ene = ene_trace
        self.magn = magn_trace
        self.cluster_sizes = clusters

    # =========================================================================
    # Topological / spectral subspace helpers
    # =========================================================================

    def _ensure_spectral_subspace(self, n_modes: int) -> None:
        """Compute partial Laplacian spectrum if not already cached.

        Calls ``self.sg.compute_k_eigvV(k=n_modes)`` which is a no-op
        when eigenvectors are already available.
        """
        # If full spectrum already computed, nothing to do
        if getattr(self.sg, "eigV", None) is not None:
            return
        self.sg.compute_k_eigvV(k=n_modes)

    def _build_subspace_matrix(self, n_modes: int) -> np.ndarray:
        """Return ``(n_modes, N)`` matrix of continuous eigenvectors.

        Row ``k`` is the *k*-th eigenvector (not binarized).  Uses
        ``self.sg.get_eigV(which=k)`` which handles NX/GT layout.
        """
        vecs = [
            self.sg.get_eigV(which=k, binarize=False)
            for k in range(n_modes)
        ]
        return np.vstack(vecs).astype(np.float64)   # (M, N)

    def _compute_rbim_energies(self, n_modes: int) -> np.ndarray:
        """Return ``(n_modes,)`` array of RBIM energies for binarized eigvecs."""
        self.sg.compute_rbim_energy_eigV_all()
        return np.array([
            self.sg.energy_eigV_RBIM[k]
            for k in range(n_modes)
        ], dtype=np.float64)

    def _best_eigenvector_seed(
        self, n_modes: int,
    ) -> tuple[np.ndarray, int]:
        """Return binarized spins and mode index of lowest-RBIM-energy eigvec."""
        rbim_E = self._compute_rbim_energies(n_modes)
        best_k = int(np.argmin(rbim_E))
        spins = self.sg.get_eigV(which=best_k, binarize=True).copy()
        # Ensure {-1, +1}  (sign(0) → 0 must become +1)
        spins[spins == 0] = 1
        return spins.astype(np.int8), best_k

    @staticmethod
    def _compute_spectral_field(
        subspace_vecs: np.ndarray,
        rbim_energies: np.ndarray,
        tau: float,
        field_strength: float,
    ) -> np.ndarray:
        """Construct TFCA external field from softmax-weighted eigvecs.

        Parameters
        ----------
        subspace_vecs : ndarray, shape ``(M, N)``
            Row-major eigenvectors.
        rbim_energies : ndarray, shape ``(M,)``
            RBIM energy per binarized eigenvector.
        tau : float
            Softmax temperature (lower → sharper weighting).
        field_strength : float
            Overall field magnitude scale.

        Returns
        -------
        ndarray, shape ``(N,)``
            External field ``h_i``.
        """
        # Normalize energies by |E_min| so gaps are O(0.01–0.1)
        E_min_abs = np.abs(rbim_energies.min())
        normed = rbim_energies / E_min_abs if E_min_abs > 0 else rbim_energies
        # softmax(-E_norm / tau): lower energy → higher weight
        logits = -normed / tau
        logits -= logits.max()                       # numerical stability
        weights = np.exp(logits)
        weights /= weights.sum()
        # h_i = field_strength * sum_k c_k * v_k[i]
        return field_strength * (weights @ subspace_vecs)  # (M,)@(M,N)→(N,)

    def _greedy_polish(
        self,
        spins: np.ndarray,
        max_sweeps: int = 50,
    ) -> np.ndarray:
        """T=0 Metropolis sweeps until no spin flips.

        Uses the graph edge structure to evaluate local energy changes.
        Returns polished spin configuration (does not modify ``self.s``).
        """
        s = spins.astype(np.float64).copy()
        edges = self.sg.get_edges_with_weights()
        # Precompute adjacency as dict-of-lists for fast node lookup
        adj: dict[int, list[tuple[int, float]]] = {
            i: [] for i in range(self.N)
        }
        for u, v, w in edges:
            adj[u].append((v, w))
            adj[v].append((u, w))

        for _ in range(max_sweeps):
            any_flip = False
            for node in range(self.N):
                h_eff = sum(w * s[nn] for nn, w in adj[node])
                dE = 2.0 * s[node] * h_eff
                if dE < 0:           # flipping reduces energy
                    s[node] = -s[node]
                    any_flip = True
            if not any_flip:
                break
        return s.astype(np.int8)

    # =========================================================================
    # Topological Metropolis  (py_topo_met)
    # =========================================================================

    def topological_metropolis_sampling(
        self, tqdm_on: bool = True,
    ) -> None:
        """MC in spectral coefficient space with mode-weighted selection.

        One MC step = *N* coefficient proposals (analogous to one lattice
        sweep), where *N* is the number of sites.  Energy and magnetization
        are recorded once per step.

        See :ref:`topological-algorithms` for algorithm details.
        """
        n_modes = min(self.topo_n_modes, self.N)
        self._ensure_spectral_subspace(n_modes)
        V = self._build_subspace_matrix(n_modes)       # (M, N)
        rbim_E = self._compute_rbim_energies(n_modes)   # (M,)

        # --- seed from best eigenvector ---
        best_spins, best_idx = self._best_eigenvector_seed(n_modes)
        coeffs = np.zeros(n_modes, dtype=np.float64)
        coeffs[best_idx] = 1.0
        field = V[best_idx].copy()                      # (N,) continuous
        spins = best_spins.astype(np.float64).copy()

        # --- edge arrays for fast energy ---
        edges = self.sg.get_edges_with_weights()
        E_current = float(compute_ising_pairwise_energy(spins, edges))

        # --- adaptive sigma per mode ---
        sigma = np.full(n_modes, self.topo_sigma_init, dtype=np.float64)
        accept_cnt = np.zeros(n_modes, dtype=np.int64)
        total_cnt = np.zeros(n_modes, dtype=np.int64)

        # --- mode selection weights ---
        mode_w = np.abs(rbim_E).copy()
        if mode_w.sum() == 0:
            mode_w[:] = 1.0
        mode_w /= mode_w.sum()

        self.ene: list[float] = []
        self.magn: list[float] = []
        best_E = E_current
        best_s = spins.copy()

        rng = np.random.default_rng(self.seed)
        iterator = (
            tqdm.tqdm(range(self.steps), desc="TopoMet")
            if tqdm_on else range(self.steps)
        )

        n_proposals = self.N  # 1 MC step = N proposals

        for step in iterator:
            for _ in range(n_proposals):
                # 1. select mode
                m = int(rng.choice(n_modes, p=mode_w))

                # 2. propose perturbation
                delta = float(rng.normal(0.0, sigma[m]))
                field_new = field + delta * V[m]

                # 3. binarize (keep old spin where field=0)
                s_new = np.sign(field_new)
                zeros = s_new == 0
                s_new[zeros] = spins[zeros]

                # 4. energy
                E_new = float(
                    compute_ising_pairwise_energy(s_new, edges)
                )
                dE = E_new - E_current

                # 5. Metropolis
                if dE <= 0 or (
                    self.T > 0
                    and rng.random() < np.exp(-dE / self.T)
                ):
                    coeffs[m] += delta
                    field = field_new
                    spins = s_new
                    E_current = E_new
                    accept_cnt[m] += 1

                total_cnt[m] += 1

            # --- record once per MC step ---
            self.ene.append(E_current)
            self.magn.append(float(np.sum(spins)))

            if E_current < best_E:
                best_E = E_current
                best_s = spins.copy()

            # adaptive sigma every chunk steps
            if (step + 1) % self.topo_chunk_size == 0:
                for mm in range(n_modes):
                    if total_cnt[mm] > 0:
                        rate = accept_cnt[mm] / total_cnt[mm]
                        if rate > 0.35:
                            sigma[mm] *= 1.1
                        elif rate < 0.25:
                            sigma[mm] *= 0.9
                accept_cnt[:] = 0
                total_cnt[:] = 0

                # optional greedy polish
                if self.topo_polish:
                    polished = self._greedy_polish(
                        spins.astype(np.int8),
                        max_sweeps=self.topo_polish_sweeps,
                    )
                    E_pol = float(
                        compute_ising_pairwise_energy(polished, edges)
                    )
                    if E_pol < best_E:
                        best_E = E_pol
                        best_s = polished.astype(np.float64).copy()

        self.s = spins.astype(np.int8)
        self.topo_met_coeffs = coeffs
        self.topo_met_best_spins = best_s.astype(np.int8)
        self.topo_met_best_energy = best_E

    # =========================================================================
    # Topological Field Cooling Annealing  (py_topo_fca)
    # =========================================================================

    def topological_fca_sampling(self, tqdm_on: bool = True) -> None:
        """Field-constrained quench: spectral field + T~0 Metropolis.

        Constructs ``h_i = field_strength * Σ_k c_k v_k[i]`` where
        ``c_k = softmax(-E_k / tau)`` from RBIM energies, then performs
        a sudden quench from a paramagnetic configuration: standard
        Metropolis at near-zero temperature lets the system relax under
        the topological field, reaching a structure-informed minimum.

        The total number of MC sweeps equals
        ``n_temperatures * steps_per_T`` (matching the SA step budget).
        """
        n_modes = min(self.topo_n_modes, self.N)
        self._ensure_spectral_subspace(n_modes)
        V = self._build_subspace_matrix(n_modes)
        rbim_E = self._compute_rbim_energies(n_modes)

        spectral_field = self._compute_spectral_field(
            V, rbim_E, self.topo_tau, self.topo_field_strength,
        )
        self.topo_fca_field = spectral_field.copy()

        # Sudden quench: Metropolis at T~0 with spectral field
        self.field = spectral_field
        saved_T, saved_steps = self.T, self.steps
        self.T = 1e-3
        self.steps = self.n_temperatures * self.steps_per_T
        self.metropolis_sampling(tqdm_on)
        self.T, self.steps = saved_T, saved_steps

    # =========================================================================
    # Pybind11 Topological FCA  (pb_topo_fca)
    # =========================================================================

    def _run_pybind_topo_fca(self) -> None:
        """Field-constrained quench: spectral field + T~0 Metropolis.

        Builds an eigenvector-weighted external field, then performs a
        sudden quench from a paramagnetic configuration: standard
        Metropolis at near-zero temperature lets the system relax under
        the topological field, reaching a structure-informed minimum.

        The total number of MC sweeps equals
        ``n_temperatures * steps_per_T`` (matching the SA step budget).
        """
        n_modes = min(self.topo_n_modes, self.N)
        self._ensure_spectral_subspace(n_modes)
        V = self._build_subspace_matrix(n_modes)
        rbim_E = self._compute_rbim_energies(n_modes)

        spectral_field = self._compute_spectral_field(
            V, rbim_E, self.topo_tau, self.topo_field_strength,
        )
        self.topo_fca_field = spectral_field.copy()
        self.field = spectral_field

        # Sudden quench: Metropolis at T~0 with spectral field
        saved_T, saved_steps = self.T, self.steps
        self.T = 1e-3
        self.steps = self.n_temperatures * self.steps_per_T
        self._run_pybind_met()
        self.T, self.steps = saved_T, saved_steps

    # =========================================================================
    # Pybind11 Topological Metropolis  (pb_topo_met)
    # =========================================================================

    def _run_pybind_topo_met(self) -> None:
        """Run Topological Metropolis via pybind11 native kernel."""
        mod = self._load_native_module()
        ni, nw, nptr = self._build_graph_csr()

        n_modes = min(self.topo_n_modes, self.N)
        self._ensure_spectral_subspace(n_modes)
        V = self._build_subspace_matrix(n_modes)
        rbim_E = self._compute_rbim_energies(n_modes)

        best_spins, best_idx = self._best_eigenvector_seed(n_modes)
        coeffs = np.zeros(n_modes, dtype=np.float64)
        coeffs[best_idx] = 1.0
        init_field = V[best_idx].copy()

        mode_w = np.abs(rbim_E).copy()
        if mode_w.sum() == 0:
            mode_w[:] = 1.0
        mode_w /= mode_w.sum()

        s_out, ene, magn, best_s, best_E, final_coeffs = (
            mod.topo_met_sampling(
                best_spins.astype(np.int8),
                ni, nw, nptr,
                V.astype(np.float64),
                coeffs.astype(np.float64),
                init_field.astype(np.float64),
                mode_w.astype(np.float64),
                float(self.T),
                int(self.steps),
                float(self.topo_sigma_init),
                int(self.topo_chunk_size),
                bool(self.topo_polish),
                int(self.topo_polish_sweeps),
                int(self.seed),
            )
        )
        self.s = s_out
        self.ene = ene.tolist()
        self.magn = magn.tolist()
        self.topo_met_coeffs = final_coeffs
        self.topo_met_best_spins = best_s
        self.topo_met_best_energy = float(best_E)

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

        # Resolve naming convention
        rl = _resolve_runlang(self.runlang)

        # Determine mode
        use_sa = sa_mode or self.sa_enabled
        use_pt = pt_mode or self.pt_enabled

        # ---- Pybind11 backend ----
        if rl.startswith("pb_"):
            if "topo_met" in rl:
                self._run_pybind_topo_met()
            elif "topo_fca" in rl:
                self._run_pybind_topo_fca()
            elif "wolff" in rl:
                self._run_pybind_wolff()
            elif "sw" in rl:
                self._run_pybind_sw()
            elif use_pt or rl.endswith("_pt") or "pt" in rl:
                self._run_pybind_pt()
            elif use_sa or rl.endswith("_sa") or "sa" in rl:
                self._run_pybind_sa()
            else:
                self._run_pybind_met()

        # ---- CuPy (GPU) backend ----
        elif rl.startswith("cu_"):
            if "topo_met" in rl:
                # CuPy topological not yet implemented — fall back to Python
                self.topological_metropolis_sampling(tqdm_on)
            elif "topo_fca" in rl:
                self.topological_fca_sampling(tqdm_on)
            elif "wolff" in rl:
                self._run_cupy_wolff()
            elif "sw" in rl:
                self._run_cupy_sw()
            elif use_pt or rl.endswith("_pt") or "pt" in rl:
                self._run_cupy_pt()
            elif use_sa or rl.endswith("_sa") or "sa" in rl:
                self._run_cupy_sa()
            else:
                self._run_cupy_met()

        # ---- C subprocess backend ----
        elif rl.startswith("c_") or self.runlang.upper().startswith("C"):
            self.build_cprogram_command()
            self.run_cprogram(verbose)
            if clean_export:
                self.remove_run_c_files()
                self.sg.remove_exported_files()

        # ---- Python backend ----
        else:
            if "topo_met" in rl:
                self.topological_metropolis_sampling(tqdm_on)
            elif "topo_fca" in rl:
                self.topological_fca_sampling(tqdm_on)
            elif "wolff" in rl:
                self.wolff_sampling(tqdm_on)
            elif "sw" in rl:
                self.swendsen_wang_sampling(tqdm_on)
            elif use_pt:
                self.parallel_tempering_sampling(tqdm_on)
            elif use_sa:
                self.simulated_annealing_sampling(tqdm_on)
            else:
                self.metropolis_sampling(tqdm_on)

    #
    def find_ising_clusters(self, import_cl: bool = False):
        if import_cl:
            path_ising = getattr(self.sg, 'path_ising', self.dynpath)
            std_fname = getattr(self.sg, 'std_fname', 'sg')
            for i in range(self.NoClust):
                self.Ising_clusters.append(
                    np.fromfile(
                        f"{path_ising}cl{i}_{std_fname}.bin",
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
            lnodes = list(range(self.sg.N))
            lnodes_tmp = lnodes[:]

            def recursive_search(seed, magn_i, clustertmp):
                neighs = np.array(self.sg.get_neighbor_indices(seed))
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
            if hasattr(self.sg, 'export_ising_clust'):
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

        # Lattice-specific mapping (requires side1/side2)
        from ...graphs.protocols import is_lattice_graph
        if not is_lattice_graph(self.sg):
            raise TypeError(
                "mapping_nodes_to_clusters requires a lattice graph "
                "with side1 and side2 attributes."
            )
        result_array = np.empty((self.sg.side1, self.sg.side2), dtype=object)
        self.result_array = result_array

        for i, sublist in enumerate(sorted_list):
            row, col = divmod(i, self.sg.side1)
            result_array[row, col] = sublist[1]
        self.mapping = result_array


from ...const import (
    MSG_P1 as DEFAULT_MC_P1,
    MSG_P2 as DEFAULT_MC_P2,
    MSG_P3 as DEFAULT_MC_P3,
    MSG_P4 as DEFAULT_MC_P4,
    MSG_ITERATIONS as DEFAULT_MC_ITERATIONS,
    MSG_FRACTION as DEFAULT_MC_FRACTION,
)

DEFAULT_MC_VARIANT = 'exp_clocks'  # exp_clocks is 2-16x faster (see benchmark_mc_variants.py)
DEFAULT_MC_STOCHASTIC = False
DEFAULT_MC_PERIODIC = False

# Serializer defaults (for parameter sweeps)
DEFAULT_MCG_SLAPLSPECT_SRUN_MATRIX_LIST = [
    (0.8, 0.6, 0.6, 0.8),  # Default symmetric matrix
    (0.9, 0.5, 0.5, 0.9),  # Stronger diagonal
    (0.7, 0.7, 0.7, 0.7),  # Uniform
]
DEFAULT_MCG_SLAPLSPECT_SRUN_FRACTION_LIST = [0.3, 0.4, 0.5]
DEFAULT_MCG_SLAPLSPECT_SRUN_ITERATIONS_LIST = [5, 6, 7]
DEFAULT_MCG_SLAPLSPECT_SRUN_P_LIST = [0.0, 0.1, 0.2]  # pflip values

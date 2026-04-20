"""Default values for the SignedRW walker CLI.

These mirror the library-level defaults in :mod:`lrgsglib.config.const`
but live here so the argparse layer can resolve them without importing
the statsys module.
"""

from ...const import (
    DEFAULT_SRW_COVERAGE_FRAC,
    DEFAULT_SRW_MAX_N_CROSS,
    DEFAULT_SRW_N_WALKERS,
    DEFAULT_SRW_RULE,
    DEFAULT_SRW_SEED,
    DEFAULT_SRW_START,
    DEFAULT_SRW_X_NODE,
    SRW_RULES,
    SRW_START_PROTOCOLS,
    SRW_X_NODE_BEHAVIORS,
)

DEFAULT_SRW_MODE = "single"
DEFAULT_SRW_SEED_B = DEFAULT_SRW_SEED + 1
DEFAULT_SRW_START_B = DEFAULT_SRW_START

SRW_MODES = ("single", "pair")

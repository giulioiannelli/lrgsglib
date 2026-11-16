from ...const import KREG_METHOD_NATIVE as RRG_SAMPLER_NATIVE
from ...const import KREG_METHOD_PAIRING as RRG_SAMPLER_PAIRING
from ...const import KREG_PAIRING_MAX_ATTEMPTS
from ...const import KREG_PHTABB as DEFAULT_RRG_WORKDIR  # noqa: F401

# RRG_TransCluster defaults (kernels/RRG.py, src/RRG_TransCluster.py)
RRG_K_MIN = 3  # k = 2 is a union of cycles
RRG_MAX_ATTEMPTS = 64  # redraws of a disconnected topology, then raise
RRG_PAIRING_K_MAX = 4  # --sampler auto: pairing up to this k, then native
# An explicit --sampler pairing is refused when the expected draws per graph
# exceed this (10x below the attempt cap): k >= 7. Above RRG_PAIRING_K_MAX
# it only warns (k = 5: ~400 draws, about 2 min per graph at N = 1e6).
RRG_PAIRING_MAX_EXPECTED_DRAWS = KREG_PAIRING_MAX_ATTEMPTS // 10
RRG_SAMPLER_AUTO = "auto"
# RRG_SAMPLER_PAIRING: exact uniform (graphs/_shared/_regular_pairing.py)
# RRG_SAMPLER_NATIVE: the engine's own generator. NX's is not uniform:
# measured mean triangle count 1.61 +- 0.04 vs 4/3 at N = 1000, k = 3.
RRG_SAMPLERS = (RRG_SAMPLER_AUTO, RRG_SAMPLER_PAIRING, RRG_SAMPLER_NATIVE)
DEFAULT_RRG_SAMPLER = RRG_SAMPLER_AUTO
DEFAULT_RRG_ENGINE = "gt"
DEFAULT_RRG_SEED = 20260930  # root seed: the s= token of the run dir
DEFAULT_RRG_NAVG = 1000  # -na: reps (disorder realizations) per run

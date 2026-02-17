from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments
#
# Merge serializer sweep args + L3D topology args + Ising dynamics args
# (for explicit passthrough to L3D_IsingDynamics.py, no unknown pattern)
optionalaction_args_dict = {
    **L3D_ISDYN_srun_opt_args,     # side1_list, pflip_linsp, Temp_linsp, slanzarv
    **L3D_opt_args,                 # geometry, pdil, mu, sigma, edge_weight, graph_engine
    **L3D_ISDYN_opt_args,           # SA/PT args, n_thermal, runlang, freq, etc.
    **L3D_ISDYN_srun_action_args,   # exec, print, nomail, short
}
#
parser = argparse.ArgumentParser(
    description=L3D_ISDYN_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
)
#
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)

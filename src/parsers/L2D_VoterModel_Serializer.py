from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments

# Merge serializer sweep args + geometry + voter program passthrough + actions.
# (No --save-results: the voter program persists per observable gate, so the
# save gates in L2D_VOTER_opt_args are forwarded directly.)
optionalaction_args_dict = {
    **L2D_VOTER_srun_opt_args,      # side1_list, pflip_linsp, voter axes, slanzarv, array
    **L2D_opt_args,                  # geometry, cell_type, number_of_averages, graph_engine
    **L2D_VOTER_opt_args,            # voter dynamics + run + save-gate passthrough
    **L2D_VOTER_srun_action_args,    # exec, print, nomail, short
}

parser = argparse.ArgumentParser(
    description=L2D_VOTER_srun_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    allow_abbrev=False,
)

for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)

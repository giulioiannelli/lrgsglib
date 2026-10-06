from lrgsglib.config.progargs import *
from parsers.shared import GraphArgumentParser

optionalaction_args_dict = {
    **SCSGeneralized_optional_args_dict,
    **SCS_TransCluster_optional_args_dict,
    **SCS_TransCluster_action_args_dict,
}

parser = GraphArgumentParser(
    description=SCS_TransCluster_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in SCS_TransCluster_args.items():
    parser.add_argument(k, **v)
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
# graph-tool path not ported yet
parser.refuse_engine("gt", "the SCS kernels build NetworkX graphs only")

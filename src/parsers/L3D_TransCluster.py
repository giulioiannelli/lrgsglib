from lrgsglib.config.progargs import *
from parsers.shared import GraphArgumentParser

optionalaction_args_dict = {
    **L3D_TransCluster_optional_args_dict,
    **L3D_TransCluster_action_args_dict,
}

parser = GraphArgumentParser(
    description=L3D_TransCluster_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)

for k, v in L3D_TransCluster_args.items():
    parser.add_argument(k, **v)
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
# graph-tool path not ported yet
parser.refuse_engine("gt", "the kernel builds NetworkX lattices only")

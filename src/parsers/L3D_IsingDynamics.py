from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import GraphArgumentParser, refuse_gt_ignored_options
from parsers.shared import parse_arguments
#
optionalaction_args_dict = {
    **L3D_opt_args,
    **L3D_ISDYN_opt_args,
    **L3D_ISDYN_action_args,
}
#
parser = GraphArgumentParser(
    description=L3D_ISDYN_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
#
for k, v in L3D_ISDYN_args.items():
    parser.add_argument(k, **v)
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
# graph-tool: refuse the options its prepare kernels would drop
refuse_gt_ignored_options(parser)
parser.refuse_engine(
    "gt",
    "--save-results reads the graph through NetworkX methods",
    when=lambda a: a.save_results,
)
parser.refuse_engine(
    "gt",
    "graph-tool 3D lattices ignore --pdil",
    when=lambda a: a.pdil != parser.get_default("pdil"),
)

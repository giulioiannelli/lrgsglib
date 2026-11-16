from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import GraphArgumentParser, refuse_gt_ignored_options

# Combine ER and IsingDynamics arguments
optionalaction_args_dict = {
    **ER_opt_args,
    **ER_ISDYN_opt_args,
    **ER_ISDYN_action_args
}

parser = GraphArgumentParser(
    description=ER_ISDYN_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Mandatory arguments
for k, v in ER_ISDYN_args.items():
    parser.add_argument(k, **v)

# Optional and action arguments
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
# graph-tool: refuse the options its prepare kernels would drop
refuse_gt_ignored_options(parser)
parser.refuse_engine(
    "gt",
    "the C backends fail on graph-tool ER graphs whose giant component is below N",
    when=lambda a: a.runlang.upper().startswith("C"),
)


def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    """Parse command line arguments and return namespace."""
    return parser.parse_args()

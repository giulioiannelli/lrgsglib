from lrgsglib.core import *
from lrgsglib.config.progargs import *

# Combine ER and IsingDynamics arguments
optionalaction_args_dict = {
    **ER_opt_args,
    **ER_ISDYN_opt_args,
    **ER_ISDYN_action_args
}

parser = argparse.ArgumentParser(
    description=ER_ISDYN_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Mandatory arguments
for k, v in ER_ISDYN_args.items():
    parser.add_argument(k, **v)

# Optional and action arguments
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)


def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    """Parse command line arguments and return namespace."""
    return parser.parse_args()

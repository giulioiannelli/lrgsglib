from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments
#
optionalaction_args_dict = {
    **L3D_opt_args,
    **L3D_ISDYN_opt_args,
    **L3D_ISDYN_action_args,
}
#
parser = argparse.ArgumentParser(
    description=L3D_ISDYN_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
#
for k, v in L3D_ISDYN_args.items():
    parser.add_argument(k, **v)
for k, v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)

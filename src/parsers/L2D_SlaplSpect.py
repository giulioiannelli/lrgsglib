# from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import parse_arguments as _parse_arguments
#
optionalaction_args_dict = {
    **L2D_opt_args,
    **L2D_SlaplSpect_optional_args_dict, 
    **action_args_dict,
    **L2D_SlaplSpect_action_args_dict,
}
#
parser = argparse.ArgumentParser(
    description=L2D_SlaplSpect_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#
for k, v in L2D_SlaplSpect_args.items():
    parser.add_argument(k, **v)
for k,v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)


def parse_arguments(parser: argparse.ArgumentParser) -> argparse.Namespace:
    args = _parse_arguments(parser)
    # Backward-compatible aliases expected by legacy kernels
    args.geo = getattr(args, "geo", args.geometry)
    args.workDir = getattr(args, "workDir", args.workdir)
    return args

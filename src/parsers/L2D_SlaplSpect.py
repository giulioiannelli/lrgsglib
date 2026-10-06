# from lrgsglib.core import *
from lrgsglib.config.progargs import *
from parsers.shared import GraphArgumentParser
#
optionalaction_args_dict = {
    **L2D_opt_args,
    **L2D_SlaplSpect_optional_args_dict,
    **L2D_SlaplSpect_action_args_dict,
    **action_args_dict,
}
#
parser = GraphArgumentParser(
    description=L2D_SlaplSpect_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#
for k, v in L2D_SlaplSpect_args.items():
    parser.add_argument(k, **v)
for k,v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
# graph-tool eigenvectors are stored node-major (N, k); eigvec_dist reads rows
parser.refuse_engine("gt", "-m eigvec_dist reads NetworkX-layout eigenvectors",
                     when=lambda a: a.mode == "eigvec_dist")

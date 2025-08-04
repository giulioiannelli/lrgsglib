from lrgsglib.proglib import *
#
optional_args_dict_tmp = L2D_ISDYN_opt_args.copy()
optional_args_dict_tmp[tuple(['-rl', '--runlang'])]['default'] = 'C3'

optionalaction_args_dict = {
    **L2D_opt_args,
    **L2D_Recon_optional_args_dict,
    **L2D_Recon_action_args_dict,
    **optional_args_dict_tmp, 
    **L2D_ISDYN_action_args
}

parser = argparse.ArgumentParser(
    description=L2D_Recon_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#
for k, v in L2D_Recon_args.items():
    parser.add_argument(k, **v)
for k,v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
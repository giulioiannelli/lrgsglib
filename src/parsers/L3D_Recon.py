from lrgsglib.proglib import *
#
optional_args_dict_tmp = L3D_ISDYN_opt_args.copy()
optional_args_dict_tmp[tuple(['-rl', '--runlang'])]['default'] = 'C3'

optionalaction_args_dict = {
    **L3D_opt_args,
    **L3D_Recon_optional_args_dict,
    **L3D_Recon_action_args_dict,
    **optional_args_dict_tmp, 
    **L3D_ISDYN_action_args
}

parser = argparse.ArgumentParser(
    description=L3D_Recon_description,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#
for k, v in L3D_Recon_args.items():
    parser.add_argument(k, **v)
for k,v in optionalaction_args_dict.items():
    parser.add_argument(*k, **v)
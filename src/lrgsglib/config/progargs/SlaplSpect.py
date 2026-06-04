from .common import *

from .phelp.SlaplSpect import *
from .defs.SlaplSpect import *

from .Lattice2D import *
from .MultiplicativeCascade import *

#
# Common SlaplSpect parameters (apply to all graph types)
#
SlaplSpect_progName = 'SlaplSpect'

# Generic spectral analysis parameters (shared across L2D, MCG, etc.)
SlaplSpect_common_optional_args_dict = {
    **common_opt_args,  # Inherit common optional args (out_suffix, workdir)
    tuple(['-bc', '--bins_count']): {
        'help': phelp_binsc,
        'type': int,
        'default': DEFAULT_BINSC
    },
    tuple(['-em', '--eigen_mode']): {
        'help': phelp_eigMode,
        'type': str,
        'default': DEFAULT_EIGMODE
    },
    tuple(['-hm', '--howmany']): {
        'help': phelp_howmany_eigs,
        'type': int,
        'default': DEFAULT_HOWMANY_EIGS
    },
    tuple(['-sf', '--save_frequency']): {
        'help': phelp_data_save_freq,
        'type': int,
        'default': DEFAULT_DATA_SAVE_FREQ
    },
    tuple(['--backend']): {
        'help': phelp_slaplspect_backend,
        'type': str,
        'default': DEFAULT_SLAPLSPECT_BACKEND,
        'choices': ['numpy', 'scipy', 'cupy']
    },
    tuple(['--keep-sparse']): {
        'help': phelp_keep_sparse,
        'type': lambda x: None if x.lower() == 'none' else x.lower() == 'true',
        'default': DEFAULT_KEEP_SPARSE,
        'metavar': 'BOOL|none'
    }
}

#
# L2D_SlaplSpect (Lattice2D spectral analysis)
#
L2D_SlaplSpect_progName = 'L2D_SlaplSpect'
L2D_SlaplSpect_progNameShrt = 'L2DSS'
L2D_SlaplSpect_description = f"""
    Computational resourses regarding the Signed Laplacian spectrum of 2D
    lattices: {L2D_SlaplSpect_progName}.py
"""

L2D_SlaplSpect_args = {**L2D_args}
L2D_SlaplSpect_optional_args_dict = {
    **SlaplSpect_common_optional_args_dict,  # Inherit common spectral parameters
    tuple(['-m', '--mode']): {
        'help': phelp_l2dsspect_mode,
        'type': str,
        'default': DEFAULT_L2DSSPECT_MODE
    }
}
L2D_SlaplSpect_action_args_dict = {}

# Serializer-specific configuration (L3D-style: only sweep + slurm flags here;
# everything else passes through to L2D_SlaplSpect.py via `parse_known_args`)
L2D_SlaplSpect_srun_description = (
    f"Serializer for {L2D_SlaplSpect_progName}.py"
)

DEFAULT_L2D_SLAPLSPECT_SRUN_L_LIST = [16, 32, 64, 128]
DEFAULT_L2D_SLAPLSPECT_SRUN_P_LIST = list(np.linspace(0.0, 0.3, 100))

L2D_SlaplSpect_srun_optional_args_dict = {
    tuple(['-m', '--mode']): {
        'help': (
            "Execution mode (e.g. 'eigvals', 'eigval_dist', 'eigvec_dist'); "
            "prefix with 'slanzarv_' to submit each child via slanzarv."
        ),
        'type': str,
        'default': f"slanzarv_{DEFAULT_L2DSSPECT_MODE}",
    },
    tuple(['--L-list']): {
        'help': (
            "Explicit lattice sizes "
            f"(default per-mode; falls back to {DEFAULT_L2D_SLAPLSPECT_SRUN_L_LIST})"
        ),
        'type': int,
        'nargs': '+',
        'default': None,
    },
    tuple(['--L-linsp']): {
        'help': "Semicolon-separated linspace tuples for lattice sizes",
        'type': parse_multiple_linspace,
        'default': None,
    },
    tuple(['--p-list']): {
        'help': "Explicit probabilities for flipped edges",
        'type': float,
        'nargs': '+',
        'default': None,
    },
    tuple(['--p-linsp']): {
        'help': "Semicolon-separated linspace tuples for flipped-edge probabilities",
        'type': parse_multiple_linspace,
        'default': None,
    },
    **srun_opt_args,
}

#
# MCG_SlaplSpect (MultiplicativeCascade spectral analysis)
#
MCG_SlaplSpect_progName = 'MCG_SlaplSpect'
MCG_SlaplSpect_progNameShrt = 'MCGSS'
MCG_SlaplSpect_description = f"""
    Computational resources regarding the Signed Laplacian spectrum of
    MultiplicativeCascadeGraph: {MCG_SlaplSpect_progName}.py
"""

MCG_SlaplSpect_args = {**MC_args}  # Inherits MC graph parameters
MCG_SlaplSpect_optional_args_dict = {
    **SlaplSpect_common_optional_args_dict,  # Inherit common spectral parameters
    tuple(['-m', '--mode']): {
        'help': phelp_mcgsspect_mode,
        'type': str,
        'default': DEFAULT_MCGSSPECT_MODE
    },
    # Entropy-specific parameters
    tuple(['--entropy-steps']): {
        'help': phelp_entropy_steps,
        'type': int,
        'default': DEFAULT_ENTROPY_STEPS
    },
    tuple(['--entropy-t1']): {
        'help': phelp_entropy_t1,
        'type': int,
        'default': DEFAULT_ENTROPY_T1
    },
    tuple(['--entropy-t2']): {
        'help': phelp_entropy_t2,
        'type': int,
        'default': DEFAULT_ENTROPY_T2
    },
    tuple(['--num-samples']): {
        'help': phelp_num_samples,
        'type': int,
        'default': DEFAULT_NUM_SAMPLES
    },
    tuple(['--entropy-seed']): {
        'help': phelp_entropy_seed,
        'type': int,
        'default': DEFAULT_ENTROPY_SEED
    },
    tuple(['--entropy-norm']): {
        'help': phelp_entropy_norm,
        'type': str,
        'default': DEFAULT_ENTROPY_NORM,
        'choices': ['standard', 'complement']
    },
    tuple(['--specific-heat-scale']): {
        'help': phelp_specific_heat_scale,
        'type': str,
        'default': DEFAULT_SPECIFIC_HEAT_SCALE,
        'choices': ['logN', 'none']
    }
}
MCG_SlaplSpect_action_args_dict = {}

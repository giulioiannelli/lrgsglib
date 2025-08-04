# -*- coding: utf-8 -*-
# ##############################################################################
#
#   L2D_Recon.py
#   (2D Lattice Reconstruction Module)
#
#   Part of the lrgsglib package for network topology analysis in lattice systems.
#   Copyright (C) 2025 Giulio Iannelli.
#   All rights reserved.
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the MIT License.
#   You should have received a copy of the MIT License along with this
#   program. If not, see <https://opensource.org/licenses/MIT>.
#
# ##############################################################################
#
#   Author:         Giulio Iannelli
#   Email:          giulioiannelli.w@gmail.com
#   Organization:   Centro di Ricerche Enrico Fermi (CREF), Rome, Italy
#   Version:        0.1.0 (Initial Release)
#   Created:        14/05/2025
#   Last Modified:  14/06/2025
#
# ##############################################################################
#
#   Description:
#   ------------
#   This script is dedicated to performing computational analyses for the
#   reconstruction of Laplacian eigenmodes in 2D lattice structures.
#   It implements robust simulation and reconstruction methodologies tailored
#   for various lattice geometries. The primary goal is to facilitate the
#   study of spin interactions and emergent network topology by utilizing
#   progressive, batched simulation runs, suitable for handling extensive
#   computational tasks.
#
#   Key Features:
#   -------------
#     - Advanced reconstruction of Laplacian eigenmodes.
#     - Flexible support for diverse 2D lattice geometries.
#     - Efficient progressive batched simulations for large-scale studies.
#     - Integrated command-line argument parsing for configurable runs.
#     - Custom logging for detailed execution tracking.
#
#   Usage:
#   ------
#   python L2D_Recon.py [options]
#   For detailed options, run: python L2D_Recon.py --help
#
# ##############################################################################
#
#   TO DO / Future Enhancements:
#   ----------------------------
#     - [x] Implement Chronometer functionality to track execution time.
#           (Assuming this was the next immediate task from previous state)
#     - [ ] Performance: Optimize memory usage, especially for very large
#           lattice configurations.
#     - [ ] Extensibility: Add support for additional or custom lattice
#           geometries and boundary conditions.
#     - [ ] Parallelism: Implement parallel processing capabilities
#           (e.g., using multiprocessing or Dask) to leverage multi-core
#           architectures for faster computations.
#     - [ ] Error Handling: Enhance error handling with more specific
#           exceptions and robust recovery mechanisms for failed simulations.
#     - [ ] Testing: Develop a comprehensive test suite (unit, integration)
#           to ensure code reliability and correctness.
#     - [ ] GPU Acceleration: Explore and integrate GPU acceleration
#           (e.g., via CuPy or Numba) for computationally intensive kernels.
#     - [ ] Checkpointing: Add functionality to save and resume
#           long-running simulations from intermediate checkpoints.
#     - [ ] Documentation: Auto-generate API documentation using tools
#           like Sphinx.
#     - [ ] Configuration: Allow simulation parameters to be loaded from a
#           configuration file (e.g., YAML, JSON).
#
# ##############################################################################

# Local application/library specific imports
from parsers.L2D_Recon import *
from kernels.L2D_Recon import *

# ##############################################################################
#   Main Program Logic
# ##############################################################################

def main():
    """
    Main execution function for the L2D_Recon module.

    Orchestrates the 2D lattice reconstruction process. This involves:
    1.  Parsing command-line arguments to obtain simulation parameters
        and configurations.
    2.  Initializing a dedicated logger for recording detailed operational
        progress and potential issues.
    3.  Invoking the core computation function (`compute_recon_prog_lattice_incr`)
        which performs the reconstruction progressively in batches, using the
        provided `prepare_lattice` function for lattice setup.
    """
    # Parse arguments from command line
    args = parse_arguments(parser)
    args_dict = vars(args)
    
    # Initialize logger with custom configuration
    locallog = initialize_custom_logger(args_dict, L2D_Recon_args)
    
    # Run reconstruction with progressive lattice increments
    compute_recon_prog_lattice_incr(args, locallog, prepare_lattice)

if __name__ == "__main__":
    main()
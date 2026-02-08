# ============================================================================
# LEGACY BUILD SYSTEM
# ============================================================================
# The primary build system is now CMake + scikit-build-core (see CMakeLists.txt).
# Use: pip install -e . --no-build-isolation
# Or:  pixi run build
#
# This Makefile is kept as a supported alternative for users who prefer it.
# ============================================================================
# use bash for shell
SHELL := /bin/bash
# python includes and libraries
PYTHON3 = $(shell which python3)
PYTHON_INC = $(shell python3 -m pybind11 --includes)
PYTHON_LIB = $(shell python3-config --ldflags)
# helper variables
DEBRIS = a.out *~ *.o *.so *.pyc *.pyo
# helper functions
print-%: ; @echo '$* = $($*)'
# conda environment name
CONDA_ENV_NAME ?= lrgsgenv
# *.mk files
include build/lrgsg-paths.mk
include build/conda-config.mk
include build/cprogn.mk
include build/cconfig.mk
# .env file
ENV_FILE := .env
ENV_PY := lrgsg_env.py
GEN_ENV_PY := generate-env.py
#
path-config: rootp-file echo-paths
env-config: dotenv-file py-env-file
#
basic-config: path-config env-config create-dirs chmod-scripts
conda-config: print-conda-prefix configure-conda-environment
full-config: basic-config conda-config
c-make: $(PROGS)
cpp-make: sub_make
all: full-config c-make cpp-make 

# Generic IsingSimulator pattern rule (for legacy and Metropolis variants)
$(LRGSG_STATSYS_ISING_BIN)/IsingSimulator%: $(LRGSG_STATSYS_ISING)/IsingSimulator%.c \
									$(PATH_SRCC_FILES) \
									$(PATH_SRCC_RBIM) \
									$(PATH_SFMT_FILES) \
									$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling IsingSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# Simulated Annealing variants (need LRGSG_sa)
$(LRGSG_STATSYS_ISING_BIN)/IsingSimulator3b: $(LRGSG_STATSYS_ISING)/IsingSimulator3b.c \
									$(PATH_SRCC_FILES) \
									$(PATH_SRCC_RBIM) \
									$(PATH_SRCC_SA) \
									$(PATH_SFMT_FILES) \
									$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling IsingSimulator3b (SA)...\n"
	$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# Parallel Tempering variants (need LRGSG_pt)
$(LRGSG_STATSYS_ISING_BIN)/IsingSimulator4b: $(LRGSG_STATSYS_ISING)/IsingSimulator4b.c \
									$(PATH_SRCC_FILES) \
									$(PATH_SRCC_RBIM) \
									$(PATH_SRCC_PT) \
									$(PATH_SFMT_FILES) \
									$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling IsingSimulator4b (PT)...\n"
	$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# special rule for VoterSimulator0
$(LRGSG_STATSYS_VM_BIN)/VoterSimulator0: $(LRGSG_STATSYS_VM)/VoterSimulator0.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_VM) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
		@printf "Compiling VoterSimulator0...\n"
		$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# pattern rule for VoterSimulator variants
$(LRGSG_STATSYS_VM_BIN)/VoterSimulator%: $(LRGSG_STATSYS_VM)/VoterSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_VM) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
		@printf "Compiling VoterSimulator%s...\n" "$*"
		$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# special rule for unified ContactSimulator (no suffix)
$(LRGSG_STATSYS_CP_BIN)/ContactSimulator: $(LRGSG_STATSYS_CP)/ContactSimulator.c \
                                                       $(PATH_SRCC_FILES) \
                                                       $(PATH_SRCC_CP) \
                                                       $(PATH_SFMT_FILES) \
                                                       $(PATH_SRCC_BINDYNSYS)
	@printf "Compiling ContactSimulator (unified)...\n"
	$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# pattern rule for ContactSimulator variants
$(LRGSG_STATSYS_CP_BIN)/ContactSimulator%: $(LRGSG_STATSYS_CP)/ContactSimulator%.c \
                                                       $(PATH_SRCC_FILES) \
                                                       $(PATH_SRCC_CP) \
                                                       $(PATH_SFMT_FILES) \
                                                       $(PATH_SRCC_BINDYNSYS)
	@printf "Compiling ContactSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# Kuramoto Model
$(LRGSG_STATSYS_KUR_BIN)/KuramotoSimulator%: $(LRGSG_STATSYS_KUR)/KuramotoSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_CONTDYNSYS) \
							$(PATH_SRCC_KUR) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling KuramotoSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_KUR) -o $@ $^ $(LMFLAG)

# Reaction-Diffusion
$(LRGSG_STATSYS_RD_BIN)/ReactionDiffusionSimulator%: $(LRGSG_STATSYS_RD)/ReactionDiffusionSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_CONTDYNSYS) \
							$(PATH_SRCC_RD) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling ReactionDiffusionSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_RD) -o $@ $^ $(LMFLAG)

# Coupled ODE
$(LRGSG_STATSYS_CODE_BIN)/CoupledODESimulator%: $(LRGSG_STATSYS_CODE)/CoupledODESimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_CONTDYNSYS) \
							$(PATH_SRCC_CODE) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling CoupledODESimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_CODE) -o $@ $^ $(LMFLAG)

# Potts Model
$(LRGSG_STATSYS_POTTS_BIN)/PottsSimulator%: $(LRGSG_STATSYS_POTTS)/PottsSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_VECDYNSYS) \
							$(PATH_SRCC_POTTS) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling PottsSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_POTTS) -o $@ $^ $(LMFLAG)

# XY Model (needs contdynsys + vecdynsys)
$(LRGSG_STATSYS_XY_BIN)/XYSimulator%: $(LRGSG_STATSYS_XY)/XYSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_CONTDYNSYS) \
							$(PATH_SRCC_VECDYNSYS) \
							$(PATH_SRCC_XY) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling XYSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_XY) -o $@ $^ $(LMFLAG)

# Heisenberg Model (needs contdynsys + vecdynsys)
$(LRGSG_STATSYS_HBERG_BIN)/HeisenbergSimulator%: $(LRGSG_STATSYS_HBERG)/HeisenbergSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_CONTDYNSYS) \
							$(PATH_SRCC_VECDYNSYS) \
							$(PATH_SRCC_HBERG) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling HeisenbergSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_HBERG) -o $@ $^ $(LMFLAG)

# MultiSpecies Model
$(LRGSG_STATSYS_MSPEC_BIN)/MultiSpeciesSimulator%: $(LRGSG_STATSYS_MSPEC)/MultiSpeciesSimulator%.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_VECDYNSYS) \
							$(PATH_SRCC_MSPEC) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling MultiSpeciesSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -I$(LRGSG_STATSYS_MSPEC) -o $@ $^ $(LMFLAG)

rootp-file:
	@echo "Creating .isrootf file..."
	@echo "YES" > .isrootf

dotenv-file: 
	@echo "Generating .env file..."
	@printf "%s\n" \
	  $(foreach V,$(LIST),LRGSG_$(V)=$(LRGSG_$(V))) \
	  LRGSG_LLIB=$(LRGSG_LLIB) \
	> $(ENV_FILE)

py-env-file:
	@echo "Generating $(ENV_PY) file..."
	@$(PYTHON3) $(LRGSG_TOOLS_PY)/$(GEN_ENV_PY)

create-dirs:
	@mkdir -p $(DIRS_TO_MAKE)

sub_make:
	$(foreach d,$(LRGSG_OBJ_DIRS),$(MAKE) -C $(d);)

clean-subdirs:
	$(foreach d,$(LRGSG_OBJ_DIRS),$(MAKE) -C $(d) clean;)

clean-progs:
	@printf "Removing binary files of C programs...\n"
	@echo "  $(notdir $(PROGS))"
	@rm -f $(PROGS)

clean-files:
	@rm -f $(RW_TARGET)
	@rm -f *.o main

clean-activations:
	@rm -rf $(ACTIVATE_D) $(DEACTIVATE_D)

clean: clean-progs clean-subdirs clean-files clean-activations

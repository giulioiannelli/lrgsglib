# ============================================================================
# LRGSGLIB BUILD SYSTEM
# ============================================================================
# The primary build system is now CMake + scikit-build-core (see CMakeLists.txt).
# Use: pip install -e . --no-build-isolation
# Or:  pixi run build
#
# This Makefile is kept as a supported alternative for users who prefer it.
# Run `make help` to see all available targets.
# ============================================================================
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
# formatting
_BOLD := $(shell tput bold 2>/dev/null)
_RST  := $(shell tput sgr0 2>/dev/null)
_CYAN := $(shell tput setaf 6 2>/dev/null)
_GRN  := $(shell tput setaf 2 2>/dev/null)
_YEL  := $(shell tput setaf 3 2>/dev/null)
_RED  := $(shell tput setaf 1 2>/dev/null)
_SEP  := $(_CYAN)────────────────────────────────────────────────────────$(_RST)

# ============================================================================
# TOP-LEVEL TARGETS
# ============================================================================
.PHONY: all full-config basic-config conda-config c-make cpp-make
.PHONY: path-config env-config help clean

all:  ## Full build: configure + compile C + compile C++
	@$(MAKE) --no-print-directory full-config
	@$(MAKE) --no-print-directory c-make-section
	@$(MAKE) --no-print-directory cpp-make
	@echo ""
	@echo "$(_SEP)"
	@echo "$(_BOLD)$(_GRN)  Build complete.$(_RST)"
	@echo "$(_SEP)"

c-make-section:
	@echo ""
	@echo "$(_SEP)"
	@echo "$(_BOLD)  C simulators$(_RST)"
	@echo "$(_SEP)"
	@$(MAKE) --no-print-directory c-make

path-config: rootp-file echo-paths  ## Set up path configuration
env-config: dotenv-file py-env-file  ## Generate .env and lrgsg_env.py

basic-config:  ## Paths + env + dirs
	@echo ""
	@echo "$(_SEP)"
	@echo "$(_BOLD)  Environment configuration$(_RST)"
	@echo "$(_SEP)"
	@$(MAKE) --no-print-directory path-config env-config create-dirs chmod-scripts

conda-config:  ## Configure conda env
	@echo ""
	@echo "$(_SEP)"
	@echo "$(_BOLD)  Conda environment$(_RST)"
	@echo "$(_SEP)"
	@$(MAKE) --no-print-directory print-conda-prefix configure-conda-environment

full-config: basic-config conda-config  ## Full environment configuration

c-make: $(PROGS)  ## Compile all C simulator binaries
	@echo "  $(_GRN)All C simulators compiled.$(_RST)"

cpp-make: sub_make  ## Compile all C++/pybind11 extensions

# ============================================================================
# C SIMULATOR BUILD RULES (link against pre-compiled object files)
# ============================================================================

# --- Ising Model -----------------------------------------------------------

# Generic IsingSimulator pattern rule (legacy and Metropolis variants)
$(LRGSG_STATSYS_ISING_BIN)/IsingSimulator%: $(LRGSG_STATSYS_ISING)/IsingSimulator%.c $(OBJS_ISING)
	@printf "  $(_YEL)[CC]$(_RST) IsingSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_ISING) -o $@ $< $(OBJS_ISING) $(LMFLAG)

# Simulated Annealing variants (need LRGSG_sa)
$(LRGSG_STATSYS_ISING_BIN)/IsingSimulator3b: $(LRGSG_STATSYS_ISING)/IsingSimulator3b.c $(OBJS_ISING_SA)
	@printf "  $(_YEL)[CC]$(_RST) IsingSimulator3b (SA)\n"
	@$(GCC) $(ALLFLAGS) $(INC_ISING) -o $@ $< $(OBJS_ISING_SA) $(LMFLAG)

# Parallel Tempering variants (need LRGSG_pt)
$(LRGSG_STATSYS_ISING_BIN)/IsingSimulator4b: $(LRGSG_STATSYS_ISING)/IsingSimulator4b.c $(OBJS_ISING_PT)
	@printf "  $(_YEL)[CC]$(_RST) IsingSimulator4b (PT)\n"
	@$(GCC) $(ALLFLAGS) $(INC_ISING) -o $@ $< $(OBJS_ISING_PT) $(LMFLAG)

# --- Voter Model -----------------------------------------------------------

$(LRGSG_STATSYS_VM_BIN)/VoterSimulator0: $(LRGSG_STATSYS_VM)/VoterSimulator0.c $(OBJS_VM)
	@printf "  $(_YEL)[CC]$(_RST) VoterSimulator0\n"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -o $@ $< $(OBJS_VM) $(LMFLAG)

$(LRGSG_STATSYS_VM_BIN)/VoterSimulator%: $(LRGSG_STATSYS_VM)/VoterSimulator%.c $(OBJS_VM)
	@printf "  $(_YEL)[CC]$(_RST) VoterSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -o $@ $< $(OBJS_VM) $(LMFLAG)

# --- Contact Process -------------------------------------------------------

$(LRGSG_STATSYS_CP_BIN)/ContactSimulator: $(LRGSG_STATSYS_CP)/ContactSimulator.c $(OBJS_CP)
	@printf "  $(_YEL)[CC]$(_RST) ContactSimulator (unified)\n"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -o $@ $< $(OBJS_CP) $(LMFLAG)

$(LRGSG_STATSYS_CP_BIN)/ContactSimulator%: $(LRGSG_STATSYS_CP)/ContactSimulator%.c $(OBJS_CP)
	@printf "  $(_YEL)[CC]$(_RST) ContactSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -o $@ $< $(OBJS_CP) $(LMFLAG)

# --- Kuramoto Model --------------------------------------------------------

$(LRGSG_STATSYS_KUR_BIN)/KuramotoSimulator%: $(LRGSG_STATSYS_KUR)/KuramotoSimulator%.c $(OBJS_KUR)
	@printf "  $(_YEL)[CC]$(_RST) KuramotoSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_KUR) -o $@ $< $(OBJS_KUR) $(LMFLAG)

# --- Reaction-Diffusion ----------------------------------------------------

$(LRGSG_STATSYS_RD_BIN)/ReactionDiffusionSimulator%: $(LRGSG_STATSYS_RD)/ReactionDiffusionSimulator%.c $(OBJS_RD)
	@printf "  $(_YEL)[CC]$(_RST) ReactionDiffusionSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_RD) -o $@ $< $(OBJS_RD) $(LMFLAG)

# --- Coupled ODE -----------------------------------------------------------

$(LRGSG_STATSYS_CODE_BIN)/CoupledODESimulator%: $(LRGSG_STATSYS_CODE)/CoupledODESimulator%.c $(OBJS_CODE)
	@printf "  $(_YEL)[CC]$(_RST) CoupledODESimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_CODE) -o $@ $< $(OBJS_CODE) $(LMFLAG)

# --- Potts Model -----------------------------------------------------------

$(LRGSG_STATSYS_POTTS_BIN)/PottsSimulator%: $(LRGSG_STATSYS_POTTS)/PottsSimulator%.c $(OBJS_POTTS)
	@printf "  $(_YEL)[CC]$(_RST) PottsSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_POTTS) -o $@ $< $(OBJS_POTTS) $(LMFLAG)

# --- XY Model (needs contdynsys + vecdynsys) --------------------------------

$(LRGSG_STATSYS_XY_BIN)/XYSimulator%: $(LRGSG_STATSYS_XY)/XYSimulator%.c $(OBJS_XY)
	@printf "  $(_YEL)[CC]$(_RST) XYSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_XY) -o $@ $< $(OBJS_XY) $(LMFLAG)

# --- Heisenberg Model (needs contdynsys + vecdynsys) ------------------------

$(LRGSG_STATSYS_HBERG_BIN)/HeisenbergSimulator%: $(LRGSG_STATSYS_HBERG)/HeisenbergSimulator%.c $(OBJS_HBERG)
	@printf "  $(_YEL)[CC]$(_RST) HeisenbergSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_HBERG) -o $@ $< $(OBJS_HBERG) $(LMFLAG)

# --- MultiSpecies Model ----------------------------------------------------

$(LRGSG_STATSYS_MSPEC_BIN)/MultiSpeciesSimulator%: $(LRGSG_STATSYS_MSPEC)/MultiSpeciesSimulator%.c $(OBJS_MSPEC)
	@printf "  $(_YEL)[CC]$(_RST) MultiSpeciesSimulator%s\n" "$*"
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_MSPEC) -o $@ $< $(OBJS_MSPEC) $(LMFLAG)

# ============================================================================
# CONFIGURATION TARGETS
# ============================================================================

rootp-file:
	@echo "YES" > .isrootf

dotenv-file:
	@printf "%s\n" \
	  $(foreach V,$(LIST),LRGSG_$(V)=$(LRGSG_$(V))) \
	  LRGSG_LLIB=$(LRGSG_LLIB) \
	> $(ENV_FILE)

py-env-file:
	@$(PYTHON3) $(LRGSG_TOOLS_PY)/$(GEN_ENV_PY)

create-dirs:
	@mkdir -p $(DIRS_TO_MAKE) $(OBJ_DIR)

# ============================================================================
# C++ / PYBIND11 EXTENSION TARGETS
# ============================================================================

sub_make:
	@echo ""
	@echo "$(_SEP)"
	@echo "$(_BOLD)  C++ / pybind11 extensions$(_RST)"
	@echo "$(_SEP)"
	@$(foreach d,$(LRGSG_OBJ_DIRS),\
		printf "  $(_YEL)[C++]$(_RST) %s\n" "$(notdir $(patsubst %/,%,$(d)))"; \
		$(MAKE) --no-print-directory -C $(d) 2>&1 | sed 's/^/       /'; \
	)

# ============================================================================
# CLEAN TARGETS
# ============================================================================

clean-subdirs:  ## Clean C++ extension build artifacts
	$(foreach d,$(LRGSG_OBJ_DIRS),$(MAKE) --no-print-directory -C $(d) clean;)

clean-progs:  ## Remove compiled C simulator binaries
	@printf "  Removing C simulator binaries...\n"
	@rm -f $(PROGS)

clean-obj:  ## Remove cached object files
	@printf "  Removing cached object files...\n"
	@rm -rf $(OBJ_DIR)

clean-files:
	@rm -f $(RW_TARGET)
	@rm -f *.o main

clean-activations:  ## Remove conda activation scripts
	@rm -rf $(ACTIVATE_D) $(DEACTIVATE_D)

clean: clean-progs clean-obj clean-subdirs clean-files clean-activations  ## Full clean

# ============================================================================
# HELP
# ============================================================================

help:  ## Show this help message
	@echo ""
	@echo "$(_BOLD)LRGSGLIB Build System$(_RST)"
	@echo "$(_SEP)"
	@echo ""
	@echo "$(_BOLD)Usage:$(_RST)  make [TARGET] [CONDA_ENV_NAME=name] [LRGSG_LLIB=path]"
	@echo ""
	@echo "$(_BOLD)Build targets:$(_RST)"
	@echo "  $(_GRN)all$(_RST)                 Full build: configure + compile C + compile C++"
	@echo "  $(_GRN)c-make$(_RST)              Compile C simulators only"
	@echo "  $(_GRN)cpp-make$(_RST)            Compile C++/pybind11 extensions only"
	@echo ""
	@echo "$(_BOLD)Configuration targets:$(_RST)"
	@echo "  $(_GRN)full-config$(_RST)         Full environment configuration"
	@echo "  $(_GRN)basic-config$(_RST)        Paths + env files + directories"
	@echo "  $(_GRN)conda-config$(_RST)        Configure conda environment"
	@echo "  $(_GRN)path-config$(_RST)         Generate path files only"
	@echo "  $(_GRN)env-config$(_RST)          Generate .env and lrgsg_env.py"
	@echo "  $(_GRN)echo-paths$(_RST)          Print all configured paths"
	@echo ""
	@echo "$(_BOLD)Clean targets:$(_RST)"
	@echo "  $(_GRN)clean$(_RST)               Full clean (all artifacts)"
	@echo "  $(_GRN)clean-progs$(_RST)         Remove C simulator binaries"
	@echo "  $(_GRN)clean-obj$(_RST)           Remove cached object files"
	@echo "  $(_GRN)clean-subdirs$(_RST)       Remove C++ extension build artifacts"
	@echo "  $(_GRN)clean-activations$(_RST)   Remove conda activation scripts"
	@echo ""
	@echo "$(_BOLD)Examples:$(_RST)"
	@echo "  make all                                    # Full build"
	@echo "  make all LRGSG_LLIB=\$$(pwd)/.. CONDA_ENV_NAME=lrgsgnb"
	@echo "  make c-make                                 # C simulators only"
	@echo "  make cpp-make                               # C++ extensions only"
	@echo "  make clean                                  # Remove all build artifacts"
	@echo ""

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
all: full-config c-make 

$(LRGSG_CCORE_BIN)/IsingSimulator%: $(LRGSG_RBIM_SIMC)/IsingSimulator%.c \
									$(PATH_SRCC_FILES) \
									$(PATH_SRCC_RBIM) \
									$(PATH_SFMT_FILES) \
									$(PATH_SRCC_BINDYNSYS)
	@printf "Compiling IsingSimulator%s...\n" "$*"
	$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# special rule for voter_model
$(LRGSG_CCORE_BIN)/voter_model: $(LRGSG_STATSYS_VM)/voter_model.c \
							$(PATH_SRCC_FILES) \
							$(PATH_SRCC_VM) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
		@printf "Compiling voter_model...\n"
		$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

# special rule for contact_process
$(LRGSG_CCORE_BIN)/contact_process: $(LRGSG_CCORE_STATSYS)/contactP/contact_process.c\
							$(PATH_SRCC_FILES) \
							$(PATH_SFMT_FILES) \
							$(PATH_SRCC_BINDYNSYS)
		@printf "Compiling contact_process...\n"
		$(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

rootp-file:
	@echo "Creating .isrootf file..."
	@echo "YES" > .isrootf

dotenv-file: 
	@echo "Generating .env file..."
	@printf "%s\n" \
	  $(foreach V,$(LIST),LRGSG_$(V)=$(LRGSG_$(V))) \
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
	@printf "Removing binary files of C programs in %s:\n" "$(LRGSG_CCORE_BIN)"
	@echo "  $(notdir $(PROGS))"
	@rm -f $(PROGS)

clean-files:
	@rm -f $(RW_TARGET)
	@rm -f *.o main

clean-activations:
	@rm -rf $(ACTIVATE_D) $(DEACTIVATE_D)

clean: clean-progs clean-subdirs clean-files clean-activations


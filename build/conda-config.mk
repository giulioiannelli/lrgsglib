# conda paths
CONDA_PREFIX	= $(shell conda info --base)/envs/$(CONDA_ENV_NAME)
CONDA_BIN		= $(CONDA_PREFIX)/bin
ACTIVATE_D		= $(CONDA_PREFIX)/etc/conda/activate.d
DEACTIVATE_D	= $(CONDA_PREFIX)/etc/conda/deactivate.d
#
CUSTOM_ACTIVATE_SCRIPT		= $(ACTIVATE_D)/custom_env_setup.sh
CUSTOM_DEACTIVATE_SCRIPT	= $(DEACTIVATE_D)/custom_env_cleanup.sh
#
export PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig
#
print-conda-prefix:
	@printf "  %-24s %s\n" "Conda prefix" "$(CONDA_PREFIX)"
	@printf "  %-24s %s\n" "Conda bin" "$(CONDA_BIN)"

generate_config_script:
	@chmod +x $(CONFIG_SCRIPT_GEN)
	@bash $(CONFIG_SCRIPT_GEN)

setup_conda_activate:
	@mkdir -p $(ACTIVATE_D)
	@echo 'source $(CONFIG_SCRIPT_PATH)' > $(CUSTOM_ACTIVATE_SCRIPT)
	@chmod +x $(CUSTOM_ACTIVATE_SCRIPT)

setup_conda_deactivate:
	@mkdir -p $(DEACTIVATE_D)
	@echo 'source $(UNCONFIG_SCRIPT_PATH)' > $(CUSTOM_DEACTIVATE_SCRIPT)
	@chmod +x $(CUSTOM_DEACTIVATE_SCRIPT)

chmod-scripts:
	@find $(LRGSG_TOOLS_SCRPT) -type f -name '*.sh' -exec chmod +x {} \;

configure-conda-environment: generate_config_script setup_conda_activate setup_conda_deactivate
	@# Find the gcc compiler binary
	$(eval GCC_BIN := $(shell find $(CONDA_BIN) -name 'x86_64-conda-linux-gnu-gcc' | head -n 1))
	@# Find the g++ compiler binary (c++, NOT cpp which is the preprocessor)
	$(eval GPP_BIN := $(shell find $(CONDA_BIN) -name 'x86_64-conda-linux-gnu-c++' | head -n 1))

	@# Remove existing gcc symlink if it exists and create a new one
	@if [ -L $(CONDA_BIN)/gcc ]; then \
	    rm $(CONDA_BIN)/gcc; \
	fi
	@ln -s $(GCC_BIN) $(CONDA_BIN)/gcc;

	@# Remove existing g++ symlink if it exists and create a new one
	@if [ -L $(CONDA_BIN)/g++ ]; then \
	    rm $(CONDA_BIN)/g++; \
	fi
	@ln -s $(GPP_BIN) $(CONDA_BIN)/g++;

	@printf "  %-24s %s\n" "GCC" "$(GCC_BIN)"
	@printf "  %-24s %s\n" "G++" "$(GPP_BIN)"
	@if [ -z "$(GCC_BIN)" ] || [ -z "$(GPP_BIN)" ]; then \
		echo "Error: Unable to find required compilers in $(CONDA_BIN)"; exit 1; \
	fi

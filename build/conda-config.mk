# conda paths
CONDA_PREFIX	= $(shell conda info --root)/envs/$(CONDA_ENV_NAME)
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
	@echo "Conda prefix: $(CONDA_PREFIX)"
	@echo "Conda bin: $(CONDA_BIN)"

generate_config_script:
	@echo "Generating config script..."
	@chmod +x $(CONFIG_SCRIPT_GEN)
	@bash $(CONFIG_SCRIPT_GEN)

setup_conda_activate:
	@echo "Creating activate.d directory..."
	@mkdir -p $(ACTIVATE_D)
	@echo "Creating custom activation script..."
	@echo 'source $(CONFIG_SCRIPT_PATH)' > $(CUSTOM_ACTIVATE_SCRIPT)
	@echo "Making custom activation script executable..."
	@chmod +x $(CUSTOM_ACTIVATE_SCRIPT)
	@echo "Setup conda environment activate.d complete."

setup_conda_deactivate:
	@echo "Creating deactivate.d directory..."
	@mkdir -p $(DEACTIVATE_D)
	@echo "Creating custom deactivation script..."
	@echo 'source $(UNCONFIG_SCRIPT_PATH)' > $(CUSTOM_DEACTIVATE_SCRIPT)
	@echo "Making custom deactivation script executable..."
	@chmod +x $(CUSTOM_DEACTIVATE_SCRIPT)
	@echo "Setup conda environment deactivate.d complete."

chmod-scripts:
	find $(LRGSG_TOOLS_SCRPT) -type f -name '*.sh' -exec chmod +x {} \;

configure-conda-environment: generate_config_script setup_conda_activate setup_conda_deactivate 
	@# Find the gcc compiler binary
	$(eval GCC_BIN := $(shell find $(CONDA_BIN) -name 'x86_64-conda-linux-gnu-gcc' | head -n 1))
	@# Find the g++ compiler binary
	$(eval GPP_BIN := $(shell find $(CONDA_BIN) -name 'x86_64-conda-linux-gnu-cpp' | head -n 1))

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

	@echo "Using GCC at $(GCC_BIN)"
	@echo "Using G++ at $(GPP_BIN)"
	@if [ -z "$(GCC_BIN)" ] || [ -z "$(GPP_BIN)" ]; then \
		echo "Error: Unable to find required compilers in $(CONDA_BIN)"; exit 1; \
	fi
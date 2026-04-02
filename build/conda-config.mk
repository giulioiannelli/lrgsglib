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

# --- Platform-aware compiler detection ---
# Detect platform and find the correct conda compiler prefix
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

ifeq ($(UNAME_S),Linux)
    COMPILER_PREFIX := x86_64-conda-linux-gnu
else ifeq ($(UNAME_S),Darwin)
    ifeq ($(UNAME_M),arm64)
        COMPILER_PREFIX := arm64-apple-darwin20.0.0
    else
        COMPILER_PREFIX := x86_64-apple-darwin13.4.0
    endif
else
    COMPILER_PREFIX :=
endif

configure-conda-environment: generate_config_script setup_conda_activate setup_conda_deactivate
	@# Find gcc: try conda cross-compiler first, then fall back to system gcc
	$(eval GCC_BIN := $(shell \
		if [ -n "$(COMPILER_PREFIX)" ]; then \
			find $(CONDA_BIN) -name '$(COMPILER_PREFIX)-gcc' 2>/dev/null | head -n 1; \
		fi))
	$(eval GCC_BIN := $(if $(GCC_BIN),$(GCC_BIN),$(shell which gcc 2>/dev/null)))

	@# Find g++: try conda cross-compiler first, then fall back to system g++
	$(eval GPP_BIN := $(shell \
		if [ -n "$(COMPILER_PREFIX)" ]; then \
			find $(CONDA_BIN) -name '$(COMPILER_PREFIX)-c++' 2>/dev/null | head -n 1; \
		fi))
	$(eval GPP_BIN := $(if $(GPP_BIN),$(GPP_BIN),$(shell which g++ 2>/dev/null)))

	@# Create gcc/g++ symlinks only if a conda compiler was found
	@if [ -n "$(COMPILER_PREFIX)" ] && [ -n "$(GCC_BIN)" ]; then \
		if [ -L $(CONDA_BIN)/gcc ]; then rm $(CONDA_BIN)/gcc; fi; \
		ln -s $(GCC_BIN) $(CONDA_BIN)/gcc; \
		if [ -L $(CONDA_BIN)/g++ ]; then rm $(CONDA_BIN)/g++; fi; \
		ln -s $(GPP_BIN) $(CONDA_BIN)/g++; \
	fi

	@printf "  %-24s %s\n" "Platform" "$(UNAME_S) $(UNAME_M)"
	@printf "  %-24s %s\n" "GCC" "$(GCC_BIN)"
	@printf "  %-24s %s\n" "G++" "$(GPP_BIN)"
	@if [ -z "$(GCC_BIN)" ] || [ -z "$(GPP_BIN)" ]; then \
		echo "Error: Unable to find required compilers (tried conda prefix '$(COMPILER_PREFIX)' and system)"; exit 1; \
	fi

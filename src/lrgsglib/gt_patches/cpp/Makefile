# gt_patches C++ extensions Makefile
# ===================================
# Builds graph-tool C++ extensions using Boost.Python bindings.
# Requires: graph-tool, boost-python, conda environment with CXX set.

# Use the conda environment compiler (x86_64-conda-linux-gnu-c++)
# This is critical because 'g++' in conda environments may be misconfigured.
CXX ?= $(shell which x86_64-conda-linux-gnu-c++ 2>/dev/null || echo g++)

# Python version for graph-tool pkg-config
PY3DOTVERSION := $(shell python -c "import sys; print(sys.version_info[0],sys.version_info[1],sep='.')")

# PKG_CONFIG_PATH for finding graph-tool
PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig

# Get graph-tool includes and libs separately (cleaner)
GT_CFLAGS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags graph-tool-py$(PY3DOTVERSION) 2>/dev/null)
GT_LIBS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --libs graph-tool-py$(PY3DOTVERSION) 2>/dev/null)

# Path to graph-tool core library (needed for linking)
GT_SITE_PACKAGES := $(CONDA_PREFIX)/lib/python$(PY3DOTVERSION)/site-packages/graph_tool
GT_CORE_LIB := $(GT_SITE_PACKAGES)/libgraph_tool_core.so

# Compiler flags (separate compile and link flags)
CXXFLAGS := -O3 -fopenmp -std=gnu++17 -Wall -fPIC $(GT_CFLAGS) -I$(CONDA_PREFIX)/include
# Link against boost_python and graph_tool_core for runtime symbols
# Note: We link directly against libgraph_tool_core.so for type_info symbols
LDFLAGS := -shared $(GT_LIBS) -L$(GT_SITE_PACKAGES) -lgraph_tool_core -Wl,-rpath,$(GT_SITE_PACKAGES)

# Output directory for compiled extensions (relative to this Makefile)
BUILD_DIR := .

# All targets
# Note: kcore is commented out - graph-tool has builtin gt.topology.kcore_decomposition()
# To re-enable, the kcore.cc needs to be updated to match current graph-tool dispatch API
TARGETS := $(BUILD_DIR)/triangular_lattice.so
# TARGETS += $(BUILD_DIR)/libkcore.so  # Disabled: needs API update

.PHONY: all clean info

all: $(TARGETS)
	@echo "Build complete. Extensions ready:"
	@ls -la $(BUILD_DIR)/*.so 2>/dev/null || echo "  (no .so files found)"

info:
	@echo "Build configuration:"
	@echo "  CXX           = $(CXX)"
	@echo "  PY3DOTVERSION = $(PY3DOTVERSION)"
	@echo "  CONDA_PREFIX  = $(CONDA_PREFIX)"
	@echo "  GT_CFLAGS     = $(GT_CFLAGS)"
	@echo "  GT_LIBS       = $(GT_LIBS)"
	@echo "  CXXFLAGS      = $(CXXFLAGS)"
	@echo "  LDFLAGS       = $(LDFLAGS)"

# Triangular lattice extension
$(BUILD_DIR)/triangular_lattice.so: triangular_lattice.cc triangular_lattice.hh
	@echo "Compiling triangular_lattice.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

# K-core decomposition extension
$(BUILD_DIR)/libkcore.so: kcore.cc kcore.hh
	@echo "Compiling libkcore.so..."
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $< -o $@

clean:
	@echo "Cleaning compiled extensions..."
	rm -f $(BUILD_DIR)/*.so $(BUILD_DIR)/*.o

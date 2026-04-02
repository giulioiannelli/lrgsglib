# gt-cpp.mk - Shared build configuration for graph-tool C++ extensions
# Include this from each graph's cpp/Makefile:
#   include ../../../../build/gt-cpp.mk
# or similar path based on depth

# Platform detection
UNAME_S := $(shell uname -s)

# Use the conda environment compiler (platform-aware)
ifeq ($(UNAME_S),Darwin)
    CXX ?= $(shell which clang++ 2>/dev/null || echo g++)
    SHARED_EXT := .so
    RPATH_FLAG = -Wl,-rpath,$(GT_SITE_PACKAGES)
else
    CXX ?= $(shell which x86_64-conda-linux-gnu-c++ 2>/dev/null || echo g++)
    SHARED_EXT := .so
    RPATH_FLAG = -Wl,-rpath,$(GT_SITE_PACKAGES)
endif

# Python version for graph-tool pkg-config
PY3DOTVERSION := $(shell python -c "import sys; print(sys.version_info[0],sys.version_info[1],sep='.')")

# PKG_CONFIG_PATH for finding graph-tool
PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig

# Get graph-tool includes and libs
GT_CFLAGS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --cflags graph-tool-py$(PY3DOTVERSION) 2>/dev/null)
GT_LIBS := $(shell PKG_CONFIG_PATH=$(PKG_CONFIG_PATH) pkg-config --libs graph-tool-py$(PY3DOTVERSION) 2>/dev/null)

# Path to graph-tool core library (needed for linking)
GT_SITE_PACKAGES := $(CONDA_PREFIX)/lib/python$(PY3DOTVERSION)/site-packages/graph_tool
GT_CORE_LIB := $(GT_SITE_PACKAGES)/libgraph_tool_core$(SHARED_EXT)

# Compiler flags
# Workaround: GCC 15 + old conda sysroot lacks timespec_get (see gcc15_compat.h)
GCC15_COMPAT := $(dir $(lastword $(MAKEFILE_LIST)))gcc15_compat.h
ifeq ($(UNAME_S),Darwin)
    CXXFLAGS := -O3 -std=gnu++17 -Wall -fPIC $(GT_CFLAGS) -I$(CONDA_PREFIX)/include
else
    CXXFLAGS := -O3 -fopenmp -std=gnu++17 -Wall -fPIC -include $(GCC15_COMPAT) $(GT_CFLAGS) -I$(CONDA_PREFIX)/include
endif

# Link flags - link against boost_python and graph_tool_core
LDFLAGS := -shared $(GT_LIBS) -L$(GT_SITE_PACKAGES) -lgraph_tool_core $(RPATH_FLAG)

# Output directory (current directory by default)
BUILD_DIR ?= .

# Info target for debugging
.PHONY: info
info:
	@echo "Build configuration:"
	@echo "  Platform      = $(UNAME_S)"
	@echo "  CXX           = $(CXX)"
	@echo "  PY3DOTVERSION = $(PY3DOTVERSION)"
	@echo "  CONDA_PREFIX  = $(CONDA_PREFIX)"
	@echo "  SHARED_EXT    = $(SHARED_EXT)"
	@echo "  GT_CFLAGS     = $(GT_CFLAGS)"
	@echo "  GT_LIBS       = $(GT_LIBS)"
	@echo "  CXXFLAGS      = $(CXXFLAGS)"
	@echo "  LDFLAGS       = $(LDFLAGS)"

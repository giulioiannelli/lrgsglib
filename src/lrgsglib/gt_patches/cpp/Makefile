
PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig


PY3DOTVERSION=$(shell python -c "import sys; print(sys.version_info[0],sys.version_info[1],sep='.')")
CXX := g++
CXXFLAGS := -O3 -fopenmp -std=gnu++17 -Wall -fPIC $(shell PKG_CONFIG_PATH=$(CONDA_PREFIX)/lib/pkgconfig pkg-config --cflags --libs graph-tool-py$(PY3DOTVERSION)) -I$(CONDA_PREFIX)/include -shared

all: triangular_lattice.so libkcore.so
	@echo "CONDA_PREFIX is: $(CONDA_PREFIX)"
	@echo "PKG_CONFIG_PATH is: $(PKG_CONFIG_PATH)"
	@echo "PY3DOTVERSION is: $(PY3DOTVERSION)"


triangular_lattice.so: triangular_lattice.cc
	$(CXX) $(CXXFLAGS) $< -o $@ 

libkcore.so: kcore.hh kcore.cc
	${CXX} ${CXXFLAGS} kcore.cc -o libkcore.so

clean:
	rm -f *.so

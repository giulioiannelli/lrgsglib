GCC := $(CONDA_PREFIX)/bin/gcc
CPP := $(CONDA_PREFIX)/bin/g++
CXXFLAGS      = -O3 -Wall -shared -std=c++14 -fPIC
GFLAGS        = -g
OFLAGS        = -O3
DSFMTFLAG     = -DSFMT_MEXP=19937
LMFLAG        = -lm
WFLAG1        = -Wall
WFLAG2        = -Wextra
OPENMPFLAG    = -fopenmp
WFLAGS        = ${WFLAG1} ${WFLAG2}
INC_PATHS    := -I$(LRGSG_LIB_CCORE) -I$(LRGSG_RBIM_SIMC) -I$(LRGSG_CCORE_SFMT)
ALLFLAGS      = ${GFLAGS} ${OFLAGS} ${WFLAGS} ${DSFMTFLAG} ${INC_PATHS}
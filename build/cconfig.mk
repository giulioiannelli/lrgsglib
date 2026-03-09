GCC := $(CONDA_PREFIX)/bin/gcc
CPP := $(CONDA_PREFIX)/bin/g++
# Workaround: GCC 15 + old conda sysroot lacks timespec_get (see gcc15_compat.h)
GCC15_COMPAT := $(LRGSG_BUILD)/gcc15_compat.h
CXXFLAGS      = -O3 -Wall -shared -std=c++14 -fPIC -include $(GCC15_COMPAT)
GFLAGS        = -g
OFLAGS        = -O3
DSFMTFLAG     = -DSFMT_MEXP=19937
LMFLAG        = -lm
WFLAG1        = -Wall
WFLAG2        = -Wextra
OPENMPFLAG    = -fopenmp
WFLAGS        = ${WFLAG1} ${WFLAG2}
INC_PATHS    := -I$(LRGSG_STATSYS_CCORE) -I$(LRGSG_CCORE_SFMT)
INC_ISING    := $(INC_PATHS) -I$(LRGSG_STATSYS_ISING)
ALLFLAGS      = ${GFLAGS} ${OFLAGS} ${WFLAGS} ${DSFMTFLAG}
# ============================================================================
# Object file cache — shared sources compiled once, linked many times
# ============================================================================
OBJ_DIR := $(LRGSG_BUILD)/obj
# Shared core objects (used by ALL simulators)
OBJ_UTILS     := $(OBJ_DIR)/LRGSG_utils.o
OBJ_SFMTRNG   := $(OBJ_DIR)/sfmtrng.o
OBJ_SFMT      := $(OBJ_DIR)/SFMT.o
OBJ_BINDYNSYS := $(OBJ_DIR)/LRGSG_bindynsys.o
OBJS_CORE     := $(OBJ_UTILS) $(OBJ_SFMTRNG) $(OBJ_SFMT) $(OBJ_BINDYNSYS)
# Per-model shared objects
OBJ_RBIM       := $(OBJ_DIR)/LRGSG_rbim.o
OBJ_SA         := $(OBJ_DIR)/LRGSG_sa.o
OBJ_PT         := $(OBJ_DIR)/LRGSG_pt.o
OBJ_VM         := $(OBJ_DIR)/LRGSG_vm.o
OBJ_CP         := $(OBJ_DIR)/LRGSG_cp.o
OBJ_CONTDYNSYS := $(OBJ_DIR)/LRGSG_contdynsys.o
OBJ_VECDYNSYS  := $(OBJ_DIR)/LRGSG_vecdynsys.o
OBJ_KURAMOTO   := $(OBJ_DIR)/LRGSG_kuramoto.o
OBJ_RD         := $(OBJ_DIR)/LRGSG_rd.o
OBJ_CODE       := $(OBJ_DIR)/LRGSG_code.o
OBJ_POTTS      := $(OBJ_DIR)/LRGSG_potts.o
OBJ_XY         := $(OBJ_DIR)/LRGSG_xy.o
OBJ_HBERG      := $(OBJ_DIR)/LRGSG_heisenberg.o
OBJ_MSPEC      := $(OBJ_DIR)/LRGSG_multispec.o
# Composite sets for convenience
OBJS_ISING     := $(OBJS_CORE) $(OBJ_RBIM)
OBJS_ISING_SA  := $(OBJS_ISING) $(OBJ_SA)
OBJS_ISING_PT  := $(OBJS_ISING) $(OBJ_PT)
OBJS_VM        := $(OBJS_CORE) $(OBJ_VM)
OBJS_CP        := $(OBJS_CORE) $(OBJ_CP)
OBJS_KUR       := $(OBJS_CORE) $(OBJ_CONTDYNSYS) $(OBJ_KURAMOTO)
OBJS_RD        := $(OBJS_CORE) $(OBJ_CONTDYNSYS) $(OBJ_RD)
OBJS_CODE      := $(OBJS_CORE) $(OBJ_CONTDYNSYS) $(OBJ_CODE)
OBJS_POTTS     := $(OBJS_CORE) $(OBJ_VECDYNSYS) $(OBJ_POTTS)
OBJS_XY        := $(OBJS_CORE) $(OBJ_CONTDYNSYS) $(OBJ_VECDYNSYS) $(OBJ_XY)
OBJS_HBERG     := $(OBJS_CORE) $(OBJ_CONTDYNSYS) $(OBJ_VECDYNSYS) $(OBJ_HBERG)
OBJS_MSPEC     := $(OBJS_CORE) $(OBJ_VECDYNSYS) $(OBJ_MSPEC)
# ============================================================================
# Object file compilation rules
# ============================================================================
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
# Core shared objects
$(OBJ_UTILS): $(LRGSG_STATSYS_CCORE)/LRGSG_utils.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
$(OBJ_SFMTRNG): $(LRGSG_STATSYS_CCORE)/sfmtrng.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
$(OBJ_SFMT): $(LRGSG_CCORE_SFMT)/SFMT.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
$(OBJ_BINDYNSYS): $(LRGSG_STATSYS_CCORE)/LRGSG_bindynsys.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
# Ising model objects
$(OBJ_RBIM): $(LRGSG_STATSYS_ISING)/LRGSG_rbim.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_ISING) -c $< -o $@
$(OBJ_SA): $(LRGSG_STATSYS_ISING)/LRGSG_sa.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_ISING) -c $< -o $@
$(OBJ_PT): $(LRGSG_STATSYS_ISING)/LRGSG_pt.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_ISING) -c $< -o $@
# Voter model
$(OBJ_VM): $(LRGSG_STATSYS_VM)/LRGSG_vm.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
# Contact process
$(OBJ_CP): $(LRGSG_STATSYS_CP)/LRGSG_cp.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
# Continuous dynamics
$(OBJ_CONTDYNSYS): $(LRGSG_STATSYS_CCORE)/LRGSG_contdynsys.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
$(OBJ_VECDYNSYS): $(LRGSG_STATSYS_CCORE)/LRGSG_vecdynsys.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -c $< -o $@
# Per-model objects
$(OBJ_KURAMOTO): $(LRGSG_STATSYS_KUR)/LRGSG_kuramoto.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_KUR) -c $< -o $@
$(OBJ_RD): $(LRGSG_STATSYS_RD)/LRGSG_rd.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_RD) -c $< -o $@
$(OBJ_CODE): $(LRGSG_STATSYS_CODE)/LRGSG_code.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_CODE) -c $< -o $@
$(OBJ_POTTS): $(LRGSG_STATSYS_POTTS)/LRGSG_potts.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_POTTS) -c $< -o $@
$(OBJ_XY): $(LRGSG_STATSYS_XY)/LRGSG_xy.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_XY) -c $< -o $@
$(OBJ_HBERG): $(LRGSG_STATSYS_HBERG)/LRGSG_heisenberg.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_HBERG) -c $< -o $@
$(OBJ_MSPEC): $(LRGSG_STATSYS_MSPEC)/LRGSG_multispec.c | $(OBJ_DIR)
	@$(GCC) $(ALLFLAGS) $(INC_PATHS) -I$(LRGSG_STATSYS_MSPEC) -c $< -o $@

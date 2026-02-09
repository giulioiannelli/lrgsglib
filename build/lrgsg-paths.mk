# project root (the lrgsglib submodule/package directory)
LRGSG_ROOT := $(shell pwd)
# library base directory (overridable via `make LRGSG_LLIB=/custom/path`)
# This is the outer project root when using lrgsglib as a submodule
# Default: same as LRGSG_ROOT for standalone usage
LRGSG_LLIB ?= $(LRGSG_ROOT)
# define every LRGSG_<NAME> := $(LRGSG_<PARENT>)/<dirname> 
# Note: DATA, IPYNB, and LOG are relative to LLIB (outer project)
# while all other paths are relative to ROOT (lrgsglib submodule)
LRGSG_PATHS := \
	BUILD:build:ROOT \
	DATA:data:LLIB \
	IPYNB:ipynb:LLIB \
	LOG:.log:LLIB \
	SRC:src:ROOT \
	TEST:test:ROOT \
	TOOLS:tools:ROOT \
	LIB:lrgsglib:SRC \
	TOOLS_SCRPT:bash:TOOLS \
	TOOLS_PY:py:TOOLS \
	LIB_GT_PATCHES:gt_patches:LIB \
	LIB_NX_PATCHES:nx_patches:LIB \
	LIB_STOCPROC:stocproc:LIB \
	GT_PATCHES_CPP:cpp:LIB_GT_PATCHES \
	LIB_STATSYS:statsys:LIB \
	STATSYS_CCORE:_ccore:LIB_STATSYS \
	CCORE_SFMT:SFMT:STATSYS_CCORE \
	STATSYS_ISING:IsingDynamics/ccore:LIB_STATSYS \
	STATSYS_ISING_BIN:bin:STATSYS_ISING \
	STATSYS_CP:ContactProcess/ccore:LIB_STATSYS \
	STATSYS_CP_BIN:bin:STATSYS_CP \
	STATSYS_VM:VoterModel/ccore:LIB_STATSYS \
	STATSYS_VM_BIN:bin:STATSYS_VM \
	STATSYS_SRW:SignedRW/ccore:LIB_STATSYS \
	STATSYS_KUR:KuramotoModel/ccore:LIB_STATSYS \
	STATSYS_KUR_BIN:bin:STATSYS_KUR \
	STATSYS_RD:ReactionDiffusionModel/ccore:LIB_STATSYS \
	STATSYS_RD_BIN:bin:STATSYS_RD \
	STATSYS_CODE:CoupledODEModel/ccore:LIB_STATSYS \
	STATSYS_CODE_BIN:bin:STATSYS_CODE \
	STATSYS_POTTS:PottsModel/ccore:LIB_STATSYS \
	STATSYS_POTTS_BIN:bin:STATSYS_POTTS \
	STATSYS_XY:XYModel/ccore:LIB_STATSYS \
	STATSYS_XY_BIN:bin:STATSYS_XY \
	STATSYS_HBERG:HeisenbergModel/ccore:LIB_STATSYS \
	STATSYS_HBERG_BIN:bin:STATSYS_HBERG \
	STATSYS_MSPEC:MultiSpeciesModel/ccore:LIB_STATSYS \
	STATSYS_MSPEC_BIN:bin:STATSYS_MSPEC \
	RBIM_BASE:base:STATSYS_ISING \
	RBIM_STORE:storer:STATSYS_ISING \
	SRW_LATT:Lattices:STATSYS_SRW \
	LIB_GRAPHS:graphs:LIB \
	GRAPHS_GT:gt:LIB_GRAPHS \
	GRAPHS_GT_BA:BarabasiAlbertGT:GRAPHS_GT \
	GRAPHS_GT_BA_CPP:cpp:GRAPHS_GT_BA \
	GRAPHS_GT_EBA:ExtendedBarabasiAlbertGT:GRAPHS_GT \
	GRAPHS_GT_EBA_CPP:cpp:GRAPHS_GT_EBA \
	GRAPHS_GT_DBA:DualBarabasiAlbertGT:GRAPHS_GT \
	GRAPHS_GT_DBA_CPP:cpp:GRAPHS_GT_DBA \
	GRAPHS_GT_FC:FullyConnectedGT:GRAPHS_GT \
	GRAPHS_GT_FC_CPP:cpp:GRAPHS_GT_FC \
	GRAPHS_GT_HK:HolmeKimGT:GRAPHS_GT \
	GRAPHS_GT_HK_CPP:cpp:GRAPHS_GT_HK \
	GRAPHS_GT_MC:MultiplicativeCascadeGT:GRAPHS_GT \
	GRAPHS_GT_MC_CPP:cpp:GRAPHS_GT_MC
# make paths
define mk_path
  $(eval LRGSG_$(word 1,$(subst :, ,$(1))) := \
         $(LRGSG_$(word 3,$(subst :, ,$(1))))/$(word 2,$(subst :, ,$(1))))
endef
$(eval $(foreach P,$(LRGSG_PATHS),$(call mk_path,$(P))))
# collect suffixes
LIST := $(foreach P,$(LRGSG_PATHS),$(word 1,$(subst :, ,$(P))))
# ALL_PATHS := $(sort $(foreach V,$(LIST),$(LRGSG_$(V))))
# print paths
.PHONY: echo-paths
echo-paths: $(ALL_PATHS)
	@echo "LRGSG_ROOT = $(LRGSG_ROOT)"
	@$(foreach V,$(LIST),echo LRGSG_$(V) = $(LRGSG_$(V));)
	@echo "LRGSG_LLIB = $(LRGSG_LLIB)"
#
LRGSG_OBJ_DIRS := $(LRGSG_GT_PATCHES_CPP) \
           $(LRGSG_SRW_LATT) \
           $(LRGSG_RBIM_BASE) \
           $(LRGSG_RBIM_STORE) \
           $(LRGSG_GRAPHS_GT_BA_CPP) \
           $(LRGSG_GRAPHS_GT_EBA_CPP) \
           $(LRGSG_GRAPHS_GT_DBA_CPP) \
           $(LRGSG_GRAPHS_GT_FC_CPP) \
           $(LRGSG_GRAPHS_GT_HK_CPP) \
           $(LRGSG_GRAPHS_GT_MC_CPP)
# Legacy aliases for backward compatibility
LRGSG_LIB_CCORE    := $(LRGSG_STATSYS_CCORE)
LRGSG_CCORE_BIN    := $(LRGSG_STATSYS_ISING_BIN)
LRGSG_CCORE_STATSYS := $(LRGSG_STATSYS_CCORE)
LRGSG_RBIM_SIMC    := $(LRGSG_STATSYS_ISING)
#
CONFIG_SCRIPT_GEN			= $(LRGSG_TOOLS_SCRPT)/generate_config.sh
CONFIG_SCRIPT_PATH			= $(LRGSG_TOOLS_SCRPT)/config_env.sh
UNCONFIG_SCRIPT_PATH		= $(LRGSG_TOOLS_SCRPT)/unconfig_env.sh
# 
#
DIRS_TO_MAKE = $(LRGSG_DATA) \
	$(LRGSG_STATSYS_ISING_BIN) \
	$(LRGSG_STATSYS_CP_BIN) \
	$(LRGSG_STATSYS_VM_BIN) \
	$(LRGSG_STATSYS_KUR_BIN) \
	$(LRGSG_STATSYS_RD_BIN) \
	$(LRGSG_STATSYS_CODE_BIN) \
	$(LRGSG_STATSYS_POTTS_BIN) \
	$(LRGSG_STATSYS_XY_BIN) \
	$(LRGSG_STATSYS_HBERG_BIN) \
	$(LRGSG_STATSYS_MSPEC_BIN)
#

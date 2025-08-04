# project root
LRGSG_ROOT := $(shell pwd)
# define every LRGSG_<NAME> := $(LRGSG_<PARENT>)/<dirname> 
LRGSG_PATHS := \
	BUILD:build:ROOT \
	DATA:data:ROOT \
	IPYNB:ipynb:ROOT \
	LOG:.log:ROOT \
	SRC:src:ROOT \
	TEST:test:ROOT \
	TOOLS:tools:ROOT \
	LIB:lrgsglib:SRC \
	TOOLS_SCRPT:bash:TOOLS \
	TOOLS_PY:py:TOOLS \
	LIB_CCORE:Ccore:LIB \
	LIB_GT_PATCHES:gt_patches:LIB \
	LIB_NX_PATCHES:nx_patches:LIB \
	LIB_STOCPROC:stocproc:LIB \
	CCORE_BIN:bin:LIB_CCORE \
	CCORE_SFMT:SFMT:LIB_CCORE \
	CCORE_STATSYS:statsys:LIB_CCORE \
	GT_PATCHES_CPP:cpp:LIB_GT_PATCHES \
        STATSYS_RBIM:RBIsingM:CCORE_STATSYS \
        STATSYS_SRW:signedRw:CCORE_STATSYS \
        STATSYS_VM:voterM:CCORE_STATSYS \
        STATSYS_CP:contactP:CCORE_STATSYS \
        RBIM_BASE:base:STATSYS_RBIM \
	RBIM_SIMC:simulatorC:STATSYS_RBIM \
	RBIM_STORE:storer:STATSYS_RBIM \
	SRW_LATT:Lattices:STATSYS_SRW
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
#
LRGSG_OBJ_DIRS := $(LRGSG_GT_PATCHES_CPP) \
           $(LRGSG_SRW_LATT) \
           $(LRGSG_RBIM_BASE) \
           $(LRGSG_RBIM_STORE)
#
CONFIG_SCRIPT_GEN			= $(LRGSG_TOOLS_SCRPT)/generate_config.sh
CONFIG_SCRIPT_PATH			= $(LRGSG_TOOLS_SCRPT)/config_env.sh
UNCONFIG_SCRIPT_PATH		= $(LRGSG_TOOLS_SCRPT)/unconfig_env.sh
# 
#
DIRS_TO_MAKE = $(LRGSG_DATA) $(LRGSG_CCORE_BIN)
#
